## ================= 1) DGP (AFM: strict & weak) ==================
toeplitz_cor <- function(rho, n){
  idx <- 0:(n-1)
  outer(idx, idx, function(i,j) rho^abs(i-j))
}

dgp_gen_afm <- function(n, T, r_signal=1, r_noise=0, n_noise_affected=1, rho_cs=0.2, sigma_min=0.6, sigma_max=1.4){
  # Approximate Factor Model (AFM): X = F*Lambda' + common noise + idiosyncratic correlated noise
  F_signal <- matrix(rnorm(T*r_signal), nrow=T, ncol=r_signal)
  Lambda_signal <- matrix(rnorm(n*r_signal), nrow=n, ncol=r_signal)
  Lambda_signal <- Lambda_signal / sqrt(mean(Lambda_signal^2))
  C_signal <- F_signal %*% t(Lambda_signal)
  
  if(r_noise>0){
    G_noise <- matrix(rnorm(T*r_noise), nrow=T, ncol=r_noise)
    Lambda_noise_dense <- matrix(rnorm(n*r_noise), nrow=n, ncol=r_noise)
    if(n_noise_affected==n){
      Lambda_noise <- Lambda_noise_dense
    } else {
      Lambda_noise <- matrix(0, nrow=n, ncol=r_noise)
      affected_indices <- sample.int(n, n_noise_affected)
      Lambda_noise[affected_indices,] <- Lambda_noise_dense[affected_indices,]
    }
    E_common_noise <- G_noise %*% t(Lambda_noise)
  } else {
    E_common_noise <- matrix(0, nrow=T, ncol=n)
  }
  
  sig_i <- runif(n, sigma_min, sigma_max)
  D <- diag(sig_i)
  R <- toeplitz_cor(rho_cs, n)
  Psi <- D %*% R %*% D
  E_idio_noise <- MASS::mvrnorm(n=T, mu=rep(0, n), Sigma=Psi)
  
  X <- C_signal + E_common_noise + E_idio_noise
  list(X=X, F=F_signal, Lambda=Lambda_signal)
}

dgp_make_y_lead_from_F <- function(F, h, gamma=1, sigma_v=1){
  # Predictive regression DGP: y_{t+h} = gamma * F_t + v_{t+h}
  T <- length(F)
  v <- rnorm(T, sd=sigma_v)
  y <- rep(NA_real_, T)
  y[(h+1):T] <- gamma * F[1:(T-h)] + v[(h+1):T]
  list(y=y, idx=1:(T-h))
}

## ===================== 2) Metrics =====================
mn_ci <- function(z){
  # Robust CI: if all missing -> return NA triplet (prevents downstream errors)
  z <- z[is.finite(z)]
  k <- length(z)
  if(k==0) return(c(mean=NA_real_, l=NA_real_, u=NA_real_))
  m <- mean(z)
  s <- sd(z)
  se <- s/sqrt(k)
  q <- qt(0.975, df=max(1, k-1))
  c(mean=m, l=m-q*se, u=m+q*se)
}

proj_rowspace <- function(X){
  A <- X %*% t(X)
  P <- t(X) %*% MASS::ginv(A) %*% X
  (P+t(P))/2
}

angle_deg_subspace <- function(X, b, v){
  P <- proj_rowspace(X)
  vs <- as.numeric(P %*% v)
  nb <- sqrt(sum(b^2))
  nv <- sqrt(sum(vs^2))
  if(nb==0||nv==0) return(NA_real_)
  acos(min(1, max(0, abs(sum(b*vs)/(nb*nv))))) * 180/pi
}

frobenius_rel_subspace <- function(X, b, v){
  P <- proj_rowspace(X)
  vs <- as.numeric(P %*% v)
  if(sum(b^2)==0 || sum(vs^2)==0) return(NA_real_)
  Pb <- tcrossprod(b)/as.numeric(crossprod(b))
  Pv <- tcrossprod(vs)/as.numeric(crossprod(vs))
  A <- X %*% Pb
  B <- X %*% Pv
  num <- sqrt(sum((A-B)^2))
  den <- sqrt(sum(B^2))
  if(den==0) NA_real_ else num/den
}

subspace_distance <- function(A, B){
  qa <- qr.Q(qr(A))
  qb <- qr.Q(qr(B))
  cos_angles <- svd(t(qa) %*% qb)$d
  max_angle_rad <- acos(min(cos_angles[length(cos_angles)], 1.0))
  max_angle_rad * 180/pi
}

metrics_eigs <- function(Xtr_centered, k=5){
  cov_X <- cov(Xtr_centered)
  eig_vals <- eigen(cov_X, symmetric=TRUE, only.values=TRUE)$values
  k_use <- min(k, length(eig_vals))
  top <- eig_vals[seq_len(k_use)]
  if(k_use<k) top <- c(top, rep(NA_real_, k-k_use))
  PV <- eig_vals / sum(eig_vals)
  PV5 <- PV[1:5]
  if(length(PV5) < 5) PV5 <- c(PV5, rep(NA_real_, 5-length(PV5)))
  list(eig_1=top[1], eig_2=top[2], eig_3=top[3], eig_4=top[4], eig_5=top[5],
       PV1=PV5[1], PV2=PV5[2], PV3=PV5[3], PV4=PV5[4], PV5=PV5[5])
}

safe_cor <- function(a, b){
  a <- as.numeric(a)
  b <- as.numeric(b)
  ok <- is.finite(a) & is.finite(b)
  if(sum(ok) < 2) return(NA_real_)
  if(sd(a[ok])==0 || sd(b[ok])==0) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok]))
}

## ===================== 3) Models =====================
preprocess_train_test <- function(Xtr_raw, Xte_raw, scale_data=TRUE){
  # Center/scale using training moments only
  mu <- colMeans(Xtr_raw)
  if(scale_data){
    sdv <- apply(Xtr_raw, 2, sd)
    sdv[sdv==0] <- 1
  } else {
    sdv <- FALSE
  }
  list(Xtr=scale(Xtr_raw, center=mu, scale=sdv), Xte=scale(Xte_raw, center=mu, scale=sdv))
}

model_oracle <- function(Xtr_raw, Xte_raw, ytr, yte, truth, cfg){
  # Oracle benchmark: regress y on true factor F
  F_use <- truth$F_use
  Ttr <- length(ytr)
  Tuse <- length(F_use)
  Ftr <- F_use[1:Ttr]
  Fte <- F_use[(Ttr+1):Tuse]
  gamma_or <- as.numeric(sum(Ftr*ytr) / sum(Ftr^2))
  yhat_te <- gamma_or * Fte
  list(yhat_te=yhat_te, metrics=list(mse_or=mean((yte - yhat_te)^2)), extras=list(gamma_or=gamma_or))
}

model_pcr_bai <- function(Xtr_raw, Xte_raw, ytr, yte, truth, cfg){
  # PCR (Bai-Ng IC for r) + factor regression on estimated PCs
  estimate_r_bai_ng <- function(X, r_max=10){
    Xc <- scale(X, center=TRUE, scale=FALSE)
    sv <- svd(Xc)
    T <- nrow(Xc)
    n <- ncol(Xc)
    r_max_use <- min(r_max, max(0L, min(n, T)-1L))
    V <- sv$d^2 / (T*n)
    IC <- numeric(r_max_use+1L)
    for(r in 0:r_max_use){
      tailV <- if(r==0) V else V[(r+1L):length(V)]
      sigma2_r <- if(length(tailV)==0) mean(V) else mean(tailV)
      penalty <- r * (log(min(n, T)) * (n + T) / (n*T))
      IC[r+1L] <- log(sigma2_r) + penalty
    }
    which.min(IC) - 1L
  }
  
  bai_factors_and_b <- function(Xtr, Xte, ytr, r=1){
    sv <- svd(Xtr)
    U <- sv$u
    D <- sv$d
    V <- sv$v
    r <- min(r, ncol(V), nrow(U))
    if(r<1) return(NULL)
    Ur <- U[,1:r,drop=FALSE]
    Dr <- D[1:r]
    Vr <- V[,1:r,drop=FALSE]
    Ttr <- nrow(Xtr)
    Ftr_hat <- sqrt(Ttr) * Ur
    Fte_hat <- sqrt(Ttr) * Xte %*% sweep(Vr, 2, Dr, `/`)
    gamma_hat <- as.numeric(solve(crossprod(Ftr_hat), crossprod(Ftr_hat, ytr)))
    a <- as.numeric(gamma_hat / (sqrt(Ttr)*Dr))
    b_hat <- as.numeric(Vr %*% a)
    Lambda_hat <- Vr %*% diag(Dr[1:r], nrow=r) / sqrt(Ttr)
    list(Ftr=Ftr_hat, Fte=Fte_hat, gamma=gamma_hat, b=b_hat, Lambda_hat=Lambda_hat)
  }
  
  st <- preprocess_train_test(Xtr_raw, Xte_raw, scale_data=FALSE)
  Xtr <- st$Xtr
  Xte <- st$Xte
  r_hat <- estimate_r_bai_ng(Xtr, r_max=cfg$pcr_rmax)
  bf <- bai_factors_and_b(Xtr, Xte, ytr, r=r_hat)
  
  if(is.null(bf)){
    return(list(yhat_te=rep(NA_real_, length(yte)),
                metrics=list(mse_pcr=NA_real_, corr_pcr_or=NA_real_, frob_pcr=NA_real_, ang_pcr=NA_real_, r_hat=as.numeric(r_hat)),
                extras=list()))
  }
  
  yhat_te <- as.numeric(bf$Fte %*% bf$gamma)
  
  fr <- frobenius_rel_subspace(Xtr, bf$b, truth$Lambda)
  ang <- angle_deg_subspace(Xtr, bf$b, truth$Lambda)
  if(truth$r_noise>0){
    fr <- NA_real_
    ang <- subspace_distance(bf$Lambda_hat, truth$Lambda) # subspace-only comparison
  }
  
  list(yhat_te=yhat_te,
       metrics=list(mse_pcr=mean((yte - yhat_te)^2), corr_pcr_or=suppressWarnings(cor(yhat_te, truth$yhat_or)),
                    frob_pcr=fr, ang_pcr=ang, r_hat=as.numeric(r_hat)),
       extras=list())
}

model_pls <- function(Xtr_raw, Xte_raw, ytr, yte, truth, cfg){
  # PLS: choose ncomp by CV-MSE on training set; evaluate by test MSE like ridge/lasso
  st <- preprocess_train_test(Xtr_raw, Xte_raw, scale_data=TRUE)
  Xtr <- st$Xtr
  Xte <- st$Xte
  Ttr <- length(ytr)
  max_ncomp <- min(cfg$pls_max_ncomp, ncol(Xtr), max(1L, Ttr-1L))
  if(max_ncomp < 1L){
    return(list(yhat_te=rep(NA_real_, length(yte)), metrics=list(mse_pls=NA_real_, corr_pls_or=NA_real_, ncomp_best_pls=NA_real_), extras=list()))
  }
  
  df_tr <- data.frame(y=ytr, Xtr)
  fit <- pls::plsr(y ~ ., data=df_tr, ncomp=max_ncomp, validation="CV", method=cfg$pls_method)
  
  press <- fit$validation$PRESS
  if(is.null(press)){
    return(list(yhat_te=rep(NA_real_, length(yte)), metrics=list(mse_pls=NA_real_, corr_pls_or=NA_real_, ncomp_best_pls=NA_real_), extras=list()))
  }
  
  press_vec <- as.numeric(press)[seq_len(max_ncomp + 1L)]  # includes comp=0
  mse_cv_vec <- press_vec / Ttr                            # CV-MSE per ncomp (incl 0)
  mse_cv_vec <- mse_cv_vec[-1]                             # drop comp=0
  ncomp_best <- which.min(mse_cv_vec)
  
  yhat_te <- as.numeric(predict(fit, newdata=data.frame(Xte), ncomp=ncomp_best))
  list(yhat_te=yhat_te,
       metrics=list(mse_pls=mean((yte - yhat_te)^2), corr_pls_or=suppressWarnings(cor(yhat_te, truth$yhat_or)), ncomp_best_pls=as.numeric(ncomp_best)),
       extras=list())
}

model_ridge <- function(Xtr_raw, Xte_raw, ytr, yte, truth, cfg){
  # Ridge regression with CV lambda on standardized X
  st <- preprocess_train_test(Xtr_raw, Xte_raw, scale_data=TRUE)
  Xtr <- st$Xtr
  Xte <- st$Xte
  cv_fit <- glmnet::cv.glmnet(Xtr, ytr, alpha=0, intercept=TRUE, nfolds=cfg$cv_nfolds, type.measure=cfg$cv_type_measure)
  yhat_te <- as.numeric(predict(cv_fit, s="lambda.min", newx=Xte, type="response"))
  list(yhat_te=yhat_te,
       metrics=list(mse_ridge=mean((yte - yhat_te)^2), corr_ridge_or=suppressWarnings(cor(yhat_te, truth$yhat_or)), lambda_best_ridge=as.numeric(cv_fit$lambda.min)),
       extras=list())
}

model_lasso <- function(Xtr_raw, Xte_raw, ytr, yte, truth, cfg){
  # Lasso regression with CV lambda on standardized X (records sparsity)
  st <- preprocess_train_test(Xtr_raw, Xte_raw, scale_data=TRUE)
  Xtr <- st$Xtr
  Xte <- st$Xte
  cv_fit <- glmnet::cv.glmnet(Xtr, ytr, alpha=1, intercept=TRUE, nfolds=cfg$cv_nfolds, type.measure=cfg$cv_type_measure)
  yhat_te <- as.numeric(predict(cv_fit, s="lambda.min", newx=Xte, type="response"))
  beta <- coef(cv_fit, s="lambda.min")
  n_predictors <- sum(beta != 0) - 1
  list(yhat_te=yhat_te,
       metrics=list(mse_lasso=mean((yte - yhat_te)^2), corr_lasso_or=suppressWarnings(cor(yhat_te, truth$yhat_or)),
                    lambda_best_lasso=as.numeric(cv_fit$lambda.min), n_predictors_lasso=as.numeric(n_predictors)),
       extras=list())
}

model_enet <- function(Xtr_raw, Xte_raw, ytr, yte, truth, cfg){
  # Elastic Net: outer loop over alpha, inner CV over lambda 
  st <- preprocess_train_test(Xtr_raw, Xte_raw, scale_data=TRUE)
  Xtr <- st$Xtr
  Xte <- st$Xte
  
  best <- list(alpha=NA_real_, lambda=NA_real_, cvm=Inf, fit=NULL)
  for(a in cfg$enet_alpha_grid){
    cv_fit <- glmnet::cv.glmnet(Xtr, ytr, alpha=a, intercept=TRUE, nfolds=cfg$cv_nfolds, type.measure=cfg$cv_type_measure)
    cvm_min <- min(cv_fit$cvm, na.rm=TRUE)
    if(is.finite(cvm_min) && cvm_min < best$cvm){
      best$alpha <- a
      best$lambda <- as.numeric(cv_fit$lambda.min)
      best$cvm <- cvm_min
      best$fit <- cv_fit
    }
  }
  
  if(is.null(best$fit)){
    return(list(yhat_te=rep(NA_real_, length(yte)),
                metrics=list(mse_enet=NA_real_, corr_enet_or=NA_real_, lambda_best_enet=NA_real_, alpha_best_enet=NA_real_, n_predictors_enet=NA_real_),
                extras=list()))
  }
  
  yhat_te <- as.numeric(predict(best$fit, s="lambda.min", newx=Xte, type="response"))
  beta <- coef(best$fit, s="lambda.min")
  n_predictors <- sum(beta != 0) - 1
  list(yhat_te=yhat_te,
       metrics=list(mse_enet=mean((yte - yhat_te)^2), corr_enet_or=suppressWarnings(cor(yhat_te, truth$yhat_or)),
                    lambda_best_enet=best$lambda, alpha_best_enet=as.numeric(best$alpha), n_predictors_enet=as.numeric(n_predictors)),
       extras=list())
}

models_registry <- function(){ list(oracle=model_oracle, pcr=model_pcr_bai, pls=model_pls, ridge=model_ridge, lasso=model_lasso, enet=model_enet) }

## ===================== 4) Engine =====================
engine_split <- function(X0, y, test_frac, min_train){
  # Train/test split on aligned sample
  Tuse <- nrow(X0)
  Ttr <- max(min_train, floor((1-test_frac)*Tuse))
  if(Ttr>=Tuse) Ttr <- Tuse-1L
  list(Xtr_raw=X0[1:Ttr,,drop=FALSE], Xte_raw=X0[(Ttr+1):Tuse,,drop=FALSE], ytr=y[1:Ttr], yte=y[(Ttr+1):Tuse], Ttr=Ttr, Tuse=Tuse)
}

run_once <- function(n, T, rep_id, r_signal, r_noise, n_noise_affected, cfg, dgp_fun=dgp_gen_afm, model_funs=models_registry()){
  # One Monte Carlo replication for a given (n,T) and DGP params
  set.seed(cfg$seed_base + 11*n + 37*rep_id)
  
  data <- dgp_fun(n, T, r_signal=r_signal, r_noise=r_noise, n_noise_affected=n_noise_affected)
  X0 <- data$X
  F <- as.numeric(data$F)
  Lambda <- data$Lambda
  
  yy <- dgp_make_y_lead_from_F(F, h=cfg$h)
  idx <- yy$idx
  X0 <- X0[idx,,drop=FALSE]
  y <- yy$y[(cfg$h+1):T]
  
  sp <- engine_split(X0, y, test_frac=cfg$test_frac, min_train=cfg$min_train)
  Xtr_raw <- sp$Xtr_raw
  Xte_raw <- sp$Xte_raw
  ytr <- sp$ytr
  yte <- sp$yte
  
  truth <- list(F_use=F[idx], Lambda=Lambda, r_noise=r_noise)
  
  or_res <- model_funs$oracle(Xtr_raw, Xte_raw, ytr, yte, truth, cfg)
  truth$yhat_or <- or_res$yhat_te
  
  pcr_res <- model_funs$pcr(Xtr_raw, Xte_raw, ytr, yte, truth, cfg)
  pls_res <- model_funs$pls(Xtr_raw, Xte_raw, ytr, yte, truth, cfg)
  ridge_res <- model_funs$ridge(Xtr_raw, Xte_raw, ytr, yte, truth, cfg)
  lasso_res <- model_funs$lasso(Xtr_raw, Xte_raw, ytr, yte, truth, cfg)
  enet_res <- model_funs$enet(Xtr_raw, Xte_raw, ytr, yte, truth, cfg)
  
  corr_or_pcr <- safe_cor(or_res$yhat_te, pcr_res$yhat_te)
  corr_pls_pcr <- safe_cor(pls_res$yhat_te, pcr_res$yhat_te)
  corr_ridge_pcr <- safe_cor(ridge_res$yhat_te, pcr_res$yhat_te)
  corr_lasso_pcr <- safe_cor(lasso_res$yhat_te, pcr_res$yhat_te)
  corr_enet_pcr <- safe_cor(enet_res$yhat_te, pcr_res$yhat_te)
  
  st_pcr <- preprocess_train_test(Xtr_raw, Xte_raw, scale_data=FALSE)
  eigs <- metrics_eigs(st_pcr$Xtr, k=5) # structure of X only (same for all models)
  
  pred_obj <- list(n=n, T=T, rep_id=rep_id, yte=as.numeric(yte),
                   yhat_or=as.numeric(or_res$yhat_te),
                   yhat_pcr=as.numeric(pcr_res$yhat_te),
                   yhat_pls=as.numeric(pls_res$yhat_te),
                   yhat_ridge=as.numeric(ridge_res$yhat_te),
                   yhat_lasso=as.numeric(lasso_res$yhat_te),
                   yhat_enet=as.numeric(enet_res$yhat_te))
  
  out_vec <- c(n=n, T=T,
               mse_pcr=pcr_res$metrics$mse_pcr, mse_pls=pls_res$metrics$mse_pls, mse_ridge=ridge_res$metrics$mse_ridge, mse_lasso=lasso_res$metrics$mse_lasso, mse_enet=enet_res$metrics$mse_enet, mse_or=or_res$metrics$mse_or,
               corr_pcr_or=pcr_res$metrics$corr_pcr_or, corr_pls_or=pls_res$metrics$corr_pls_or, corr_ridge_or=ridge_res$metrics$corr_ridge_or, corr_lasso_or=lasso_res$metrics$corr_lasso_or, corr_enet_or=enet_res$metrics$corr_enet_or,
               corr_or_pcr=corr_or_pcr, corr_pls_pcr=corr_pls_pcr, corr_ridge_pcr=corr_ridge_pcr, corr_lasso_pcr=corr_lasso_pcr, corr_enet_pcr=corr_enet_pcr,
               frob_pcr=pcr_res$metrics$frob_pcr, ang_pcr=pcr_res$metrics$ang_pcr, r_hat=pcr_res$metrics$r_hat,
               ncomp_best_pls=pls_res$metrics$ncomp_best_pls,
               lambda_best_ridge=ridge_res$metrics$lambda_best_ridge,
               lambda_best_lasso=lasso_res$metrics$lambda_best_lasso, n_predictors_lasso=lasso_res$metrics$n_predictors_lasso,
               lambda_best_enet=enet_res$metrics$lambda_best_enet, alpha_best_enet=enet_res$metrics$alpha_best_enet, n_predictors_enet=enet_res$metrics$n_predictors_enet,
               eig_1=eigs$eig_1, eig_2=eigs$eig_2, eig_3=eigs$eig_3, eig_4=eigs$eig_4, eig_5=eigs$eig_5,
               PV1=eigs$PV1, PV2=eigs$PV2, PV3=eigs$PV3, PV4=eigs$PV4, PV5=eigs$PV5)
  
  attr(out_vec, "pred") <- pred_obj
  out_vec
} 

simulate_grid <- function(n_vals, k_ratio, R, r_signal, r_noise, noise_prop=NULL, noise_power=0, cfg, dgp_fun=dgp_gen_afm, model_funs=models_registry()){
  # Simulation loop over grid of n, with T = floor(k_ratio*n) (and a minimum floor)
  params_grid <- expand.grid(n=n_vals, rep_id=1:R)
  params_grid$T <- pmax(cfg$min_T, floor(k_ratio*params_grid$n))
  
  cat(sprintf(">> Running %d total simulations on %d workers...\n", nrow(params_grid), future::nbrOfWorkers()))
  
  results_list <- future.apply::future_lapply(1:nrow(params_grid), function(i){
    n_run <- params_grid$n[i]
    T_run <- params_grid$T[i]
    rep_id_run <- params_grid$rep_id[i]
    n_noise_calc <- round(n_run^noise_power)
    n_noise_run <- min(n_run, n_noise_calc)
    run_once(n=n_run, T=T_run, rep_id=rep_id_run, r_signal=r_signal, r_noise=r_noise, n_noise_affected=n_noise_run, cfg=cfg, dgp_fun=dgp_fun, model_funs=model_funs)
  }, future.seed=TRUE)
  
  cat(">> Simulations completed. Binding results...\n")
  out_df <- dplyr::bind_rows(results_list)
  for(cl in names(out_df)) out_df[[cl]] <- as.numeric(out_df[[cl]])
  attr(out_df, "pred") <- lapply(results_list, function(x) attr(x, "pred"))
  out_df
}

## ===================== 5) Aggregation =====================
aggregate_results <- function(out_df){
  # Aggregate by n: mean and 95% CI across replications
  agg <- do.call(rbind, lapply(split(out_df, out_df$n), function(d){
    c(n=d$n[1], T_m=mean(d$T, na.rm=TRUE),
      mn_ci(d$mse_pcr), mn_ci(d$mse_pls), mn_ci(d$mse_ridge), mn_ci(d$mse_lasso), mn_ci(d$mse_enet), mn_ci(d$mse_or),
      mn_ci(d$corr_pcr_or), mn_ci(d$corr_pls_or), mn_ci(d$corr_ridge_or), mn_ci(d$corr_lasso_or), mn_ci(d$corr_enet_or),
      mn_ci(d$corr_or_pcr), mn_ci(d$corr_pls_pcr), mn_ci(d$corr_ridge_pcr), mn_ci(d$corr_lasso_pcr), mn_ci(d$corr_enet_pcr),
      mn_ci(d$frob_pcr), mn_ci(d$ang_pcr), mn_ci(d$r_hat), mn_ci(d$ncomp_best_pls),
      mn_ci(d$lambda_best_ridge),
      mn_ci(d$lambda_best_lasso), mn_ci(d$n_predictors_lasso),
      mn_ci(d$lambda_best_enet), mn_ci(d$alpha_best_enet), mn_ci(d$n_predictors_enet),
      mn_ci(d$eig_1), mn_ci(d$eig_2), mn_ci(d$eig_3), mn_ci(d$eig_4), mn_ci(d$eig_5),
      mn_ci(d$PV1), mn_ci(d$PV2), mn_ci(d$PV3), mn_ci(d$PV4), mn_ci(d$PV5))
  }))
  agg <- as.data.frame(agg)
  names(agg) <- c("n","T_mean",
                  "mse_pcr_m","mse_pcr_l","mse_pcr_u",
                  "mse_pls_m","mse_pls_l","mse_pls_u",
                  "mse_ridge_m","mse_ridge_l","mse_ridge_u",
                  "mse_lasso_m","mse_lasso_l","mse_lasso_u",
                  "mse_enet_m","mse_enet_l","mse_enet_u",
                  "mse_or_m","mse_or_l","mse_or_u",
                  "corr_pcr_or_m","corr_pcr_or_l","corr_pcr_or_u",
                  "corr_pls_or_m","corr_pls_or_l","corr_pls_or_u",
                  "corr_ridge_or_m","corr_ridge_or_l","corr_ridge_or_u",
                  "corr_lasso_or_m","corr_lasso_or_l","corr_lasso_or_u",
                  "corr_enet_or_m","corr_enet_or_l","corr_enet_or_u",
                  "corr_or_pcr_m","corr_or_pcr_l","corr_or_pcr_u",
                  "corr_pls_pcr_m","corr_pls_pcr_l","corr_pls_pcr_u",
                  "corr_ridge_pcr_m","corr_ridge_pcr_l","corr_ridge_pcr_u",
                  "corr_lasso_pcr_m","corr_lasso_pcr_l","corr_lasso_pcr_u",
                  "corr_enet_pcr_m","corr_enet_pcr_l","corr_enet_pcr_u",
                  "frob_pcr_m","frob_pcr_l","frob_pcr_u",
                  "ang_pcr_m","ang_pcr_l","ang_pcr_u",
                  "r_hat_m","r_hat_l","r_hat_u",
                  "ncomp_pls_m","ncomp_pls_l","ncomp_pls_u",
                  "lambda_best_ridge_m","lambda_best_ridge_l","lambda_best_ridge_u",
                  "lambda_best_lasso_m","lambda_best_lasso_l","lambda_best_lasso_u",
                  "n_pred_lasso_m","n_pred_lasso_l","n_pred_lasso_u",
                  "lambda_best_enet_m","lambda_best_enet_l","lambda_best_enet_u",
                  "alpha_best_enet_m","alpha_best_enet_l","alpha_best_enet_u",
                  "n_pred_enet_m","n_pred_enet_l","n_pred_enet_u",
                  "eig_1_m","eig_1_l","eig_1_u",
                  "eig_2_m","eig_2_l","eig_2_u",
                  "eig_3_m","eig_3_l","eig_3_u",
                  "eig_4_m","eig_4_l","eig_4_u",
                  "eig_5_m","eig_5_l","eig_5_u",
                  "PV1_m","PV1_l","PV1_u",
                  "PV2_m","PV2_l","PV2_u",
                  "PV3_m","PV3_l","PV3_u",
                  "PV4_m","PV4_l","PV4_u",
                  "PV5_m","PV5_l","PV5_u")
  for(cl in names(agg)) agg[[cl]] <- as.numeric(agg[[cl]])
  agg
}

estimate_convergence_rates <- function(agg_df){
  agg_df$log_n <- log(agg_df$n)
  agg_df$log_T <- log(agg_df$T_mean)
  out <- list()
  
  add_rate <- function(metric_name, y){
    yy <- y
    ok <- is.finite(yy) & yy>0 & is.finite(agg_df$log_n) & is.finite(agg_df$log_T)
    if(sum(ok)>=3){
      fit_n <- lm(log(yy[ok]) ~ agg_df$log_n[ok])
      fit_T <- lm(log(yy[ok]) ~ agg_df$log_T[ok])
      out[[length(out)+1]] <<- data.frame(metric=metric_name, rate_vs_n=coef(fit_n)[2], rate_vs_T=coef(fit_T)[2], stringsAsFactors=FALSE)
    } else {
      out[[length(out)+1]] <<- data.frame(metric=metric_name, rate_vs_n=NA_real_, rate_vs_T=NA_real_, stringsAsFactors=FALSE)
    }
  }
  
  for(p in c("pcr","pls","ridge","lasso","enet")){
    col_m <- paste0("mse_", p, "_m")
    if(col_m %in% names(agg_df)){
      gap <- pmax(agg_df[[col_m]] - agg_df$mse_or_m, 1e-10)
      add_rate(paste0("Gap_", toupper(p), " (MSE)"), gap)
    }
  }
  
  if("ang_pcr_m" %in% names(agg_df)) add_rate("Angle_PCR", pmax(agg_df$ang_pcr_m, 1e-10))
  if("frob_pcr_m" %in% names(agg_df) && any(is.finite(agg_df$frob_pcr_m))) add_rate("Frobenius_PCR", pmax(agg_df$frob_pcr_m, 1e-10))
  
  if("lambda_best_ridge_m" %in% names(agg_df)) add_rate("Lambda_Ridge", pmax(agg_df$lambda_best_ridge_m, 1e-10))
  if("lambda_best_lasso_m" %in% names(agg_df)) add_rate("Lambda_Lasso", pmax(agg_df$lambda_best_lasso_m, 1e-10))
  if("n_pred_lasso_m" %in% names(agg_df)) add_rate("N_Predictors_Lasso", pmax(agg_df$n_pred_lasso_m, 1e-10))
  
  if("lambda_best_enet_m" %in% names(agg_df)) add_rate("Lambda_ENet", pmax(agg_df$lambda_best_enet_m, 1e-10))
  if("alpha_best_enet_m" %in% names(agg_df)) add_rate("Alpha_ENet", pmax(agg_df$alpha_best_enet_m, 1e-10))
  if("n_pred_enet_m" %in% names(agg_df)) add_rate("N_Predictors_ENet", pmax(agg_df$n_pred_enet_m, 1e-10))
  
  if("ncomp_pls_m" %in% names(agg_df)) add_rate("NComp_PLS", pmax(agg_df$ncomp_pls_m, 1e-10))
  if("r_hat_m" %in% names(agg_df)) add_rate("Rhat_PCR", pmax(agg_df$r_hat_m, 1e-10))
  
  for(k in 1:5){
    ck <- paste0("eig_", k, "_m")
    if(ck %in% names(agg_df)) add_rate(paste0("Eigen", k), pmax(agg_df[[ck]], 1e-10))
  }
  
  for(k in 1:5){
    pk <- paste0("PV", k, "_m")
    if(pk %in% names(agg_df)) add_rate(paste0("PV", k), pmax(agg_df[[pk]], 1e-10))
  }
  
  do.call(rbind, out)
}

## ===================== 6) Plots =====================
plot_eigenvalue_results <- function(agg){
  op <- par(mfrow=c(1,2), mar=c(4,4,3,1))
  on.exit(par(op))
  
  ylim_eigs <- range(c(agg$eig_1_l, agg$eig_1_u, agg$eig_2_l, agg$eig_2_u, agg$eig_5_l, agg$eig_5_u), na.rm=TRUE)
  plot(agg$n, agg$eig_1_m, log="xy", type="b", pch=19, lwd=2, ylim=ylim_eigs, xlab="n", ylab="Eigenvalue (log)", main="Eigenvalues Cov(X)")
  lines(agg$n, agg$eig_2_m, type="b", pch=2, lwd=1, col=2)
  lines(agg$n, agg$eig_3_m, type="b", pch=3, lwd=1, col=3)
  lines(agg$n, agg$eig_4_m, type="b", pch=4, lwd=1, col=4)
  lines(agg$n, agg$eig_5_m, type="b", pch=5, lwd=1, col=5)
  legend("topleft", legend=c("Eig 1","Eig 2","Eig 3","Eig 4","Eig 5"), col=1:5, pch=c(19,2:5), lty=1)
  
  if(all(paste0("PV", 1:5, "_m") %in% names(agg))){
    plot(agg$n, agg$PV1_m, log="x", type="b", pch=19, lwd=2,
         ylim=range(c(agg$PV1_l, agg$PV1_u, agg$PV2_l, agg$PV2_u, agg$PV3_l, agg$PV3_u, agg$PV4_l, agg$PV4_u, agg$PV5_l, agg$PV5_u), na.rm=TRUE),
         xlab="n", ylab="PV_k", main="Proportion of Variance (PV1–PV5)")
    segments(agg$n, agg$PV1_l, agg$n, agg$PV1_u)
    lines(agg$n, agg$PV2_m, type="b", pch=2, lwd=1, col=2); segments(agg$n, agg$PV2_l, agg$n, agg$PV2_u, col=2)
    lines(agg$n, agg$PV3_m, type="b", pch=3, lwd=1, col=3); segments(agg$n, agg$PV3_l, agg$n, agg$PV3_u, col=3)
    lines(agg$n, agg$PV4_m, type="b", pch=4, lwd=1, col=4); segments(agg$n, agg$PV4_l, agg$n, agg$PV4_u, col=4)
    lines(agg$n, agg$PV5_m, type="b", pch=5, lwd=1, col=5); segments(agg$n, agg$PV5_l, agg$n, agg$PV5_u, col=5)
    legend("topright", legend=c("PV1","PV2","PV3","PV4","PV5"), col=1:5, pch=c(19,2:5), lty=1)
  }
}

plot_model_results <- function(agg, prefix, hyper_key=NULL, hyper_label=NULL, has_sparsity=FALSE, sparsity_key=NULL, main_tag=NULL){
  # Generic 2x3 panel: MSE, Gap, Corr, Hyper, (optional) Sparsity, Corr vs PCR
  mse_m <- paste0("mse_", prefix, "_m")
  mse_l <- paste0("mse_", prefix, "_l")
  mse_u <- paste0("mse_", prefix, "_u")
  corr_m <- paste0("corr_", prefix, "_or_m")
  corr_l <- paste0("corr_", prefix, "_or_l")
  corr_u <- paste0("corr_", prefix, "_or_u")
  corr2_m <- paste0("corr_", prefix, "_pcr_m")
  corr2_l <- paste0("corr_", prefix, "_pcr_l")
  corr2_u <- paste0("corr_", prefix, "_pcr_u")
  if(!(mse_m %in% names(agg))) stop("Missing columns for prefix: ", prefix)
  ttl <- if(is.null(main_tag)) toupper(prefix) else main_tag
  
  op <- par(mfrow=c(2,3), mar=c(4,4,2,1))
  on.exit(par(op))
  
  plot(agg$n, agg[[mse_m]], log="x", pch=19,
       ylim=range(c(agg[[mse_l]], agg[[mse_u]], agg$mse_or_l, agg$mse_or_u), na.rm=TRUE),
       xlab="n", ylab="MSE (test)", main=paste0("MSE — ", ttl, " vs Oracle"))
  segments(agg$n, agg[[mse_l]], agg$n, agg[[mse_u]])
  lines(agg$n, agg[[mse_m]], lwd=2)
  lines(agg$n, agg$mse_or_m, lty=2)
  segments(agg$n, agg$mse_or_l, agg$n, agg$mse_or_u, lty=2)
  abline(h=1, lty=3)
  
  gap <- pmax(agg[[mse_m]] - agg$mse_or_m, 1e-9)
  plot(agg$n, gap, log="xy", pch=19, xlab="n", ylab="Gap MSE (log)", main=paste0("Gap: ", ttl, " − Oracle"))
  lines(agg$n, gap, lwd=2)
  
  if(corr_m %in% names(agg)){
    plot(agg$n, agg[[corr_m]], log="x", pch=19,
         ylim=range(c(0,1,agg[[corr_l]], agg[[corr_u]]), na.rm=TRUE),
         xlab="n", ylab=paste0("corr(ŷ_", ttl, ", ŷ_oracle)"), main=paste0("Forecast alignment (", ttl, ")"))
    segments(agg$n, agg[[corr_l]], agg$n, agg[[corr_u]])
    lines(agg$n, agg[[corr_m]], lwd=2)
  } else {
    plot.new(); title("Correlation (missing)")
  }
  
  if(!is.null(hyper_key)){
    hy_m <- paste0(hyper_key, "_m")
    hy_l <- paste0(hyper_key, "_l")
    hy_u <- paste0(hyper_key, "_u")
    if(hy_m %in% names(agg)){
      plot(agg$n, agg[[hy_m]], log="xy", pch=19, ylim=range(c(agg[[hy_l]], agg[[hy_u]]), na.rm=TRUE),
           xlab="n", ylab=if(is.null(hyper_label)) hyper_key else hyper_label, main=paste0("Hyperparameter (CV): ", ttl))
      segments(agg$n, agg[[hy_l]], agg$n, agg[[hy_u]])
      lines(agg$n, agg[[hy_m]], lwd=2)
    } else {
      plot.new(); title("Hyper (missing)")
    }
  } else {
    plot.new(); title("Hyper (none)")
  }
  
  if(has_sparsity && !is.null(sparsity_key)){
    sp_m <- paste0(sparsity_key, "_m")
    sp_l <- paste0(sparsity_key, "_l")
    sp_u <- paste0(sparsity_key, "_u")
    if(sp_m %in% names(agg)){
      plot(agg$n, agg[[sp_m]], log="xy", pch=19,
           ylim=range(c(1, agg[[sp_l]], agg[[sp_u]]), na.rm=TRUE),
           xlab="n", ylab="N. predictors (log)", main=paste0("Sparsity: ", ttl))
      segments(agg$n, agg[[sp_l]], agg$n, agg[[sp_u]])
      lines(agg$n, agg[[sp_m]], lwd=2)
    } else {
      plot.new(); title("Sparsity (missing)")
    }
  } else {
    plot.new()
  }
  
  if(corr2_m %in% names(agg)){
    plot(agg$n, agg[[corr2_m]], log="x", pch=19,
         ylim=range(c(0,1,agg[[corr2_l]], agg[[corr2_u]]), na.rm=TRUE),
         xlab="n", ylab=paste0("corr(ŷ_", ttl, ", ŷ_PCR)"), main=paste0("Forecast alignment vs PCR (", ttl, ")"))
    segments(agg$n, agg[[corr2_l]], agg$n, agg[[corr2_u]])
    lines(agg$n, agg[[corr2_m]], lwd=2)
  } else {
    plot.new(); title("Correlation vs PCR (missing)")
  }
}

# comparison plots
plot_mse_all <- function(agg, main_title){
  x <- agg$n
  ylimv <- range(c(agg$mse_pcr_m,
                   agg$mse_pls_m,
                   agg$mse_ridge_m,
                   agg$mse_lasso_m,
                   agg$mse_enet_m,
                   agg$mse_or_m), na.rm=TRUE)
  plot(x, agg$mse_pcr_m, log="x", type="b", pch=19, lwd=2, ylim=ylimv,
       xlab="n", ylab="MSE (test)", main=main_title)
  lines(x, agg$mse_pls_m,   type="b", pch=2, lwd=2, col=2)
  lines(x, agg$mse_ridge_m, type="b", pch=3, lwd=2, col=3)
  lines(x, agg$mse_lasso_m, type="b", pch=4, lwd=2, col=4)
  lines(x, agg$mse_enet_m,  type="b", pch=5, lwd=2, col=5)
  lines(x, agg$mse_or_m, lwd=2, lty=2, col=1)
  legend("topright",legend=c("PCR","PLS","Ridge","Lasso","ENet","Oracle"),
         col=c(1,2,3,4,5,1), pch=c(19,2,3,4,5,NA), lty=c(1,1,1,1,1,2),bty="n")
}

plot_gap_all <- function(agg, main_title){
  x <- agg$n
  g_pcr <- pmax(agg$mse_pcr_m - agg$mse_or_m, 1e-10)
  g_pls <- pmax(agg$mse_pls_m - agg$mse_or_m, 1e-10)
  g_rid <- pmax(agg$mse_ridge_m - agg$mse_or_m, 1e-10)
  g_las <- pmax(agg$mse_lasso_m - agg$mse_or_m, 1e-10)
  g_ene <- pmax(agg$mse_enet_m - agg$mse_or_m, 1e-10)
  ylimv <- range(c(g_pcr,g_pls,g_rid,g_las,g_ene), na.rm=TRUE)
  plot(x, g_pcr, log="xy", type="b", pch=19, lwd=2, ylim=ylimv,
       xlab="n", ylab="MSE gap (model − oracle, log)", main=main_title)
  lines(x, g_pls, type="b", pch=2, lwd=2, col=2)
  lines(x, g_rid, type="b", pch=3, lwd=2, col=3)
  lines(x, g_las, type="b", pch=4, lwd=2, col=4)
  lines(x, g_ene, type="b", pch=5, lwd=2, col=5)
  legend("topright",
         legend=c("PCR","PLS","Ridge","Lasso","ENet"),
         col=1:5, pch=c(19,2,3,4,5), lty=1, bty="n")
}

plot_corr_vs_pcr_all <- function(agg, main_title){
  x <- agg$n
  y_all <- c(
    agg$corr_or_pcr_m,
    agg$corr_pls_pcr_m,
    agg$corr_ridge_pcr_m,
    agg$corr_lasso_pcr_m,
    agg$corr_enet_pcr_m)
  y_min <- max(0, min(y_all, na.rm=TRUE) - 0.02)
  y_max <- min(1, max(y_all, na.rm=TRUE) + 0.02)
  plot(x, agg$corr_pls_pcr_m, log="x", type="b",
       pch=2, lwd=2, col=2,
       ylim=c(y_min, y_max),
       xlab="n", ylab="corr(ŷ_model, ŷ_PCR)",
       main=main_title)
  lines(x, agg$corr_ridge_pcr_m, type="b", pch=3, lwd=2, col=3)
  lines(x, agg$corr_lasso_pcr_m, type="b", pch=4, lwd=2, col=4)
  lines(x, agg$corr_enet_pcr_m,  type="b", pch=5, lwd=2, col=5)
  lines(x, agg$corr_or_pcr_m, lwd=2, lty=2, col=1)
  legend("bottomright",
         legend=c("PLS","Ridge","Lasso","ENet","Oracle"),
         col=c(2,3,4,5,1), pch=c(2,3,4,5,NA), lty=c(1,1,1,1,2), bty="n")
}

plot_convergence_all <- function(agg, main_title){
  x <- agg$n
  g_pcr <- pmax(agg$mse_pcr_m   - agg$mse_or_m, 1e-10)
  g_pls <- pmax(agg$mse_pls_m   - agg$mse_or_m, 1e-10)
  g_rid <- pmax(agg$mse_ridge_m - agg$mse_or_m, 1e-10)
  g_las <- pmax(agg$mse_lasso_m - agg$mse_or_m, 1e-10)
  g_ene <- pmax(agg$mse_enet_m  - agg$mse_or_m, 1e-10)
  plot(x, g_pcr, log="xy", type="b", pch=19, lwd=2, col=1,
       xlab="n", ylab="MSE gap to oracle (log)",
       main=main_title,
       ylim=range(c(g_pcr, g_pls, g_rid, g_las, g_ene), na.rm=TRUE))
  
  lines(x, g_pls, type="b", pch=2, lwd=2, col=2)
  lines(x, g_rid, type="b", pch=3, lwd=2, col=3)
  lines(x, g_las, type="b", pch=4, lwd=2, col=4)
  lines(x, g_ene, type="b", pch=5, lwd=2, col=5)
  lines(x, rep(1e-10, length(x)), lty=2, lwd=2, col=1)
  legend("topright",
         legend=c("PCR","PLS","Ridge","Lasso","ENet","Oracle baseline"),
         col=c(1,2,3,4,5,1),
         pch=c(19,2,3,4,5,NA),
         lty=c(1,1,1,1,1,2),
         bty="n")
}

# comparison PLS plots
plot_pls_proxy <- function(agg_strict, agg_weak,
                           main_left="PLS dimensionality selection",
                           main_right="Forecast alignment with oracle (PLS)"){
  
  stopifnot(exists("agg_strict"), exists("agg_weak"))
  merge_key <- intersect(agg_strict$n, agg_weak$n)
  S <- agg_strict[match(merge_key, agg_strict$n), ]
  W <- agg_weak  [match(merge_key, agg_weak$n), ]
  
  op <- par(mfrow=c(1,2), mar=c(4,4,3,1))
  on.exit(par(op), add=TRUE)
  
  ylim_nc <- range(c(1, S$ncomp_pls_l, S$ncomp_pls_u, W$ncomp_pls_l, W$ncomp_pls_u), na.rm=TRUE)
  plot(S$n, S$ncomp_pls_m, log="x", type="b", pch=19, lwd=2,
       ylim=ylim_nc, xlab="n",ylab="n_comp_best_pls",
       main="PLS components selection")
  segments(S$n, S$ncomp_pls_l, S$n, S$ncomp_pls_u)
  lines(W$n, W$ncomp_pls_m, type="b", pch=2, lwd=2, col=2)
  segments(W$n, W$ncomp_pls_l, W$n, W$ncomp_pls_u, col=2)
  abline(h=1, lty=2, lwd=1.5)
  legend("topright",legend=c("STRICT","WEAK"), col=c(1,2),pch=c(19,2,NA), lty=c(1,1), bty="n")

  ylim_corr <- range(c(S$corr_pls_or_l, S$corr_pls_or_u, W$corr_pls_or_l, W$corr_pls_or_u), na.rm=TRUE)
  plot(S$n, S$corr_pls_or_m, log="x", type="b", pch=19, lwd=2,
       ylim=ylim_corr, xlab="n", ylab="corr(ŷ_PLS, ŷ_oracle)",
       main="Forecast alignment (PLS)")
  segments(S$n, S$corr_pls_or_l, S$n, S$corr_pls_or_u)
  lines(W$n, W$corr_pls_or_m, type="b", pch=2, lwd=2, col=2)
  segments(W$n, W$corr_pls_or_l, W$n, W$corr_pls_or_u, col=2)
  abline(h=1, lty=2)
  legend("bottomright",legend=c("STRICT","WEAK"),col=c(1,2),pch=c(19,2), lty=1, bty="n")
}

