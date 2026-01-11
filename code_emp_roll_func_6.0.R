library(glmnet)
library(pls)

# ---------------------------
# metrics
# ---------------------------

parse_sasdate <- function(x, fmt = "%m/%d/%Y") {
  if (inherits(x, "Date")) return(x)
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  as.Date(as.character(x), format = fmt)
}

msfe <- function(y, yhat) mean((y - yhat)^2, na.rm = TRUE)

corr2 <- function(a, b) suppressWarnings(cor(a, b, use = "complete.obs"))

standardize_train_apply <- function(X_train, x_apply) {
  mu <- colMeans(X_train, na.rm = TRUE)
  sdv <- apply(X_train, 2, sd, na.rm = TRUE)
  sdv[sdv == 0 | is.na(sdv)] <- 1
  Xs <- sweep(X_train, 2, mu, "-")
  Xs <- sweep(Xs, 2, sdv, "/")
  xa <- (x_apply - mu) / sdv
  list(Xs = Xs, xa = xa, mu = mu, sd = sdv)
}

make_taus <- function(t0, t1, h, K = 5, min_gap = 12) {
  hi <- t1 - h
  if (hi < t0) return(integer(0))
  cand <- seq(t0, hi, by = min_gap)
  if (length(cand) == 0) return(integer(0))
  if (length(cand) <= K) return(cand)
  idx <- unique(round(seq(1, length(cand), length.out = K)))
  cand[idx]
}

prepare_predictors_S <- function(df_processed, date_col, drop_cols = character(0), date_fmt = "%m/%d/%Y") {
  df <- df_processed
  df[[date_col]] <- parse_sasdate(df[[date_col]], fmt = date_fmt)
  df <- df[!is.na(df[[date_col]]), , drop = FALSE]
  is_num <- sapply(df, is.numeric)
  x_cols <- names(df)[is_num]
  x_cols <- setdiff(x_cols, c(date_col, drop_cols))
  list(df = df, x_cols = x_cols)
}

build_targets_from_raw <- function(df_raw, date_col, ip_level_col, cpi_level_col, date_fmt = "%m/%d/%Y", h = 12) {
  df <- df_raw
  df[[date_col]] <- parse_sasdate(df[[date_col]], fmt = date_fmt)
  df <- df[!is.na(df[[date_col]]), , drop = FALSE]  # drop tcode row etc.
  df[[ip_level_col]]  <- as.numeric(df[[ip_level_col]])
  df[[cpi_level_col]] <- as.numeric(df[[cpi_level_col]])

  ip_t <- 100 * log(df[[ip_level_col]])
  pi_t <- 100 * log(df[[cpi_level_col]] / dplyr::lag(df[[cpi_level_col]], 12))

  z_ip  <- ip_t - dplyr::lag(ip_t,  h)
  z_cpi <- pi_t - dplyr::lag(pi_t, h)

  data.frame(
    date = df[[date_col]],
    ip_level = ip_t,
    pi_level = pi_t,
    w_ip_h = z_ip,
    w_pi_h = z_cpi,
    stringsAsFactors = FALSE
  )
}

align_S_W <- function(df_S, date_col_S, df_W, date_col_W = "date") {
  a <- df_S
  b <- df_W
  names(a)[names(a) == date_col_S] <- "date"
  names(b)[names(b) == date_col_W] <- "date"
  m <- merge(a, b, by = "date", all = FALSE)
  m <- m[order(m$date), , drop = FALSE]
  m
}

# ---------------------------
# Model fit/predict helpers
# ---------------------------

fit_predict_pcr <- function(Xs, y_train, x_now_std, r) {
  Xs0 <- Xs
  cm <- colMeans(Xs0, na.rm = TRUE)
  for (j in seq_len(ncol(Xs0))) {
    idx <- which(is.na(Xs0[, j]))
    if (length(idx) > 0) Xs0[idx, j] <- cm[j]
  }
  x0 <- x_now_std
  idxn <- which(is.na(x0))
  if (length(idxn) > 0) x0[idxn] <- cm[idxn]

  pc <- prcomp(Xs0, center = FALSE, scale. = FALSE)
  r <- min(r, ncol(pc$x))
  Z  <- pc$x[, 1:r, drop = FALSE]
  z_now <- as.numeric(x0 %*% pc$rotation[, 1:r, drop = FALSE])
  fit <- lm(y_train ~ Z)
  coef_fit <- coef(fit)
  as.numeric(coef_fit[1] + sum(coef_fit[-1] * z_now))
}

fit_predict_pls <- function(Xs, y_train, x_now_std, ncomp) {
  ncomp <- min(ncomp, ncol(Xs), nrow(Xs) - 1)
  if (ncomp < 1) return(NA_real_)
  dat <- data.frame(y = y_train, Xs)
  fit <- plsr(y ~ ., data = dat, ncomp = ncomp, validation = "none", scale = FALSE)
  pred <- predict(fit, newdata = as.data.frame(as.list(x_now_std)), ncomp = ncomp)
  as.numeric(pred)
}

fit_predict_glmnet <- function(Xs, y_train, x_now_std, alpha, lambda) {
  fit <- glmnet(x = Xs, y = y_train, alpha = alpha, lambda = lambda,
                standardize = FALSE, intercept = TRUE)
  as.numeric(predict(fit, newx = matrix(x_now_std, nrow = 1), s = lambda))
}

ridge_df <- function(Xs, lambda) {
  s <- svd(Xs, nu = 0, nv = 0)$d
  mu <- s^2
  sum(mu / (mu + lambda))
}

glmnet_k <- function(fit, lambda) {
  b <- as.matrix(coef(fit, s = lambda))[-1, 1]
  sum(abs(b) > 0)
}

# ---------------------------
# Time-series CV selector
# ---------------------------

cv_select_theta_ts <- function(method, X, w, h, taus,
                               lambda_grid = NULL, r_grid = NULL, alpha_grid = NULL,
                               verbose = FALSE) {
  
  K <- length(taus)
  if (K == 0) stop("No valid validation anchors (taus). Increase window or adjust K/min_gap.")
  
  val_loss_one <- function(predict_one) {
    errs <- rep(NA_real_, K)
    
    for (k in seq_len(K)) {
      tau <- taus[k]
      
      # ---- KEY FIX: training ends at (tau - h), not tau ----
      tr_end <- tau - h
      if (tr_end < 1) next  # too early, skip
      
      X_train <- X[1:tr_end, , drop = FALSE]
      y_train <- w[(1:tr_end) + h]       # = w[(h+1):tau]
      x_now   <- X[tau, ]                # info set at time tau
      
      std <- standardize_train_apply(X_train, x_now)
      Xs <- std$Xs
      xa <- std$xa
      
      yhat   <- predict_one(Xs, y_train, xa)
      y_true <- w[tau + h]               # validation target
      errs[k] <- (y_true - yhat)^2
    }
    
    mean(errs, na.rm = TRUE)
  }
  
  # ---- below: keep your original method branches unchanged ----
  if (method == "pcr") {
    stopifnot(!is.null(r_grid))
    losses <- rep(NA_real_, length(r_grid))
    for (i in seq_along(r_grid)) {
      r <- r_grid[i]
      losses[i] <- val_loss_one(function(Xs, y_train, xa) fit_predict_pcr(Xs, y_train, xa, r))
    }
    best_i <- which.min(losses)
    list(r = r_grid[best_i], cv_loss = losses[best_i], path = data.frame(r = r_grid, loss = losses))
  } else if (method == "pls") {
    stopifnot(!is.null(r_grid))
    losses <- rep(NA_real_, length(r_grid))
    for (i in seq_along(r_grid)) {
      r <- r_grid[i]
      losses[i] <- val_loss_one(function(Xs, y_train, xa) fit_predict_pls(Xs, y_train, xa, r))
    }
    best_i <- which.min(losses)
    list(ncomp = r_grid[best_i], cv_loss = losses[best_i], path = data.frame(ncomp = r_grid, loss = losses))
  } else if (method %in% c("ridge","lasso")) {
    stopifnot(!is.null(lambda_grid))
    alpha <- if (method == "ridge") 0 else 1
    losses <- rep(NA_real_, length(lambda_grid))
    for (i in seq_along(lambda_grid)) {
      lam <- lambda_grid[i]
      losses[i] <- val_loss_one(function(Xs, y_train, xa) fit_predict_glmnet(Xs, y_train, xa, alpha, lam))
    }
    best_i <- which.min(losses)
    list(lambda = lambda_grid[best_i], cv_loss = losses[best_i],
         path = data.frame(lambda = lambda_grid, loss = losses))
  } else if (method == "enet") {
    stopifnot(!is.null(lambda_grid), !is.null(alpha_grid), length(alpha_grid) > 0)
    best <- list(alpha = NA_real_, lambda = NA_real_, cv_loss = Inf)
    path_list <- vector("list", length(alpha_grid))
    for (a_i in seq_along(alpha_grid)) {
      a <- alpha_grid[a_i]
      losses <- rep(NA_real_, length(lambda_grid))
      for (i in seq_along(lambda_grid)) {
        lam <- lambda_grid[i]
        losses[i] <- val_loss_one(function(Xs, y_train, xa) fit_predict_glmnet(Xs, y_train, xa, a, lam))
      }
      j <- which.min(losses)
      path_list[[a_i]] <- data.frame(alpha = a, lambda = lambda_grid, loss = losses)
      if (losses[j] < best$cv_loss) {
        best$alpha <- a
        best$lambda <- lambda_grid[j]
        best$cv_loss <- losses[j]
      }
    }
    best$path <- do.call(rbind, path_list)
    best
  } else {
    stop("Unknown method: ", method)
  }
}

# ---------------------------
# Main rolling + CV driver
# ---------------------------

run_rolling_cv_forecast_v6 <- function(panel_SW, date_col = "date", x_cols, target = c("ip","cpi"),
                                       h = NULL, W = NULL, K = NULL, min_gap = NULL, step = NULL,
                                       r_grid = NULL, lambda_grid = NULL, alpha_grid = NULL, 
                                       models = NULL, verbose = TRUE) {
  
  # ---- require main to set core design params ----
  if (is.null(h)) stop("h must be provided from main.")
  if (is.null(W)) stop("W must be provided from main.")
  if (is.null(K)) stop("K must be provided from main.")
  if (is.null(min_gap)) stop("min_gap must be provided from main.")
  if (is.null(step)) stop("step must be provided from main.")
  
  # ---- require grids/models (or set minimal defaults if you prefer) ----
  if (is.null(r_grid)) stop("r_grid must be provided from main.")
  if (is.null(lambda_grid)) stop("lambda_grid must be provided from main.")
  if (is.null(alpha_grid)) stop("alpha_grid must be provided from main.")
  if (is.null(models)) stop("models must be provided from main.")
  
  target <- match.arg(target)

  dates <- panel_SW[[date_col]]
  X_full <- as.matrix(panel_SW[, x_cols, drop = FALSE])

  if (target == "ip") {
    w <- as.numeric(panel_SW$w_ip_h)
    level <- as.numeric(panel_SW$ip_level)
  } else {
    w <- as.numeric(panel_SW$w_pi_h)
    level <- as.numeric(panel_SW$pi_level)
  }
  
  if (target == "ip" && !("w_ip_h" %in% names(panel_SW))) stop("Missing w_ip_h in panel_SW.")
  if (target == "cpi" && !("w_pi_h" %in% names(panel_SW))) stop("Missing w_pi_h in panel_SW.")
  
  
  nT <- nrow(panel_SW)
  T_start <- max(W, 1)
  T_end <- nT - h
  if (T_end <= T_start) stop("Not enough sample after alignment for given W and h.")

  outer_Ts <- seq(T_start, T_end, by = step)
  n_out <- length(outer_Ts)

  out <- data.frame(
    date_pred   = dates[outer_Ts],
    date_target = dates[outer_Ts + h],
    z_true      = w[outer_Ts + h],
    level_t     = level[outer_Ts],
    y_true      = level[outer_Ts + h],
    stringsAsFactors = FALSE
  )

  for (m in models) {
    out[[paste0("zhat_", m)]] <- NA_real_
    out[[paste0("yhat_", m)]] <- NA_real_
  }
  out$yhat_rw <- out$level_t

  if ("pcr" %in% models) out$pcr_r <- NA_integer_
  if ("pls" %in% models) out$pls_ncomp <- NA_integer_
  if ("ridge" %in% models) { out$ridge_lambda <- NA_real_; out$ridge_df <- NA_real_ }
  if ("lasso" %in% models) { out$lasso_lambda <- NA_real_; out$lasso_k <- NA_integer_ }
  if ("enet" %in% models)  { out$enet_alpha <- NA_real_; out$enet_lambda <- NA_real_; out$enet_k <- NA_integer_ }

  for (ii in seq_along(outer_Ts)) {
    T0 <- outer_Ts[ii]
    t0 <- T0 - W + 1
    t1 <- T0
    if (verbose && (ii == 1 || ii %% max(1, floor(n_out/10)) == 0)) {
      message("Outer T=", T0, " (", as.character(dates[T0]), "), window=[", t0, ",", t1, "]")
    }

    X_win <- X_full[t0:t1, , drop = FALSE]
    w_win <- w[t0:t1]
    taus <- make_taus(t0 = 1, t1 = W, h = h, K = K, min_gap = min_gap)
    if (length(taus) == 0) next

    idx_tr <- 1:(W - h)
    X_tr_raw <- X_win[idx_tr, , drop = FALSE]
    y_tr <- w_win[idx_tr + h]
    x_now_raw <- X_win[W, ]
    std_full <- standardize_train_apply(X_tr_raw, x_now_raw)
    X_tr <- std_full$Xs
    x_now <- std_full$xa

    if ("pcr" %in% models) {
      sel <- cv_select_theta_ts("pcr", X_win, w_win, h = h, taus = taus, r_grid = r_grid)
      r <- sel$r
      out$pcr_r[ii] <- r
      zhat <- fit_predict_pcr(X_tr, y_tr, x_now, r)
      out$zhat_pcr[ii] <- zhat
      out$yhat_pcr[ii] <- out$level_t[ii] + zhat
    }

    if ("pls" %in% models) {
      sel <- cv_select_theta_ts("pls", X_win, w_win, h = h, taus = taus, r_grid = r_grid)
      nc <- sel$ncomp
      out$pls_ncomp[ii] <- nc
      zhat <- fit_predict_pls(X_tr, y_tr, x_now, nc)
      out$zhat_pls[ii] <- zhat
      out$yhat_pls[ii] <- out$level_t[ii] + zhat
    }

    if ("ridge" %in% models) {
      sel <- cv_select_theta_ts("ridge", X_win, w_win, h = h, taus = taus, lambda_grid = lambda_grid)
      lam <- sel$lambda
      out$ridge_lambda[ii] <- lam
      out$ridge_df[ii] <- ridge_df(X_tr, lam)
      zhat <- fit_predict_glmnet(X_tr, y_tr, x_now, alpha = 0, lambda = lam)
      out$zhat_ridge[ii] <- zhat
      out$yhat_ridge[ii] <- out$level_t[ii] + zhat
    }

    if ("lasso" %in% models) {
      sel <- cv_select_theta_ts("lasso", X_win, w_win, h = h, taus = taus, lambda_grid = lambda_grid)
      lam <- sel$lambda
      out$lasso_lambda[ii] <- lam
      fit <- glmnet(x = X_tr, y = y_tr, alpha = 1, lambda = lam, standardize = FALSE, intercept = TRUE)
      out$lasso_k[ii] <- glmnet_k(fit, lam)
      zhat <- as.numeric(predict(fit, newx = matrix(x_now, nrow = 1), s = lam))
      out$zhat_lasso[ii] <- zhat
      out$yhat_lasso[ii] <- out$level_t[ii] + zhat
    }

    if ("enet" %in% models) {
      sel <- cv_select_theta_ts("enet", X_win, w_win, h = h, taus = taus,
                                lambda_grid = lambda_grid, alpha_grid = alpha_grid)
      a <- sel$alpha; lam <- sel$lambda
      out$enet_alpha[ii] <- a
      out$enet_lambda[ii] <- lam
      fit <- glmnet(x = X_tr, y = y_tr, alpha = a, lambda = lam, standardize = FALSE, intercept = TRUE)
      out$enet_k[ii] <- glmnet_k(fit, lam)
      zhat <- as.numeric(predict(fit, newx = matrix(x_now, nrow = 1), s = lam))
      out$zhat_enet[ii] <- zhat
      out$yhat_enet[ii] <- out$level_t[ii] + zhat
    }
  }

  out$target <- target
  out
}
