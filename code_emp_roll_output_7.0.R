rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(scales)


IN_DIR  <- "emp_outputs_roll"
OUT_DIR <- file.path(IN_DIR, "output_7.0")

RDS_FILE <- file.path(IN_DIR, "all_results_v6.rds")
IP_FILE  <- file.path(IN_DIR, "results_IP_v6.csv")
CPI_FILE <- file.path(IN_DIR, "results_CPI_v6.csv")

# Corr path settings
CORR_WINDOW  <- 120
CORR_MIN_OBS <- 24

# Subsample periods (user provided)
PERIODS <- list(
  full = c(as.Date("1900-01-01"), as.Date("2100-12-31")),
  s1   = c(as.Date("1971-01-01"), as.Date("1984-12-31")),
  s2   = c(as.Date("1985-01-01"), as.Date("2007-12-31")),
  s3   = c(as.Date("2008-01-01"), as.Date("2019-12-31"))
)

# Plot sizing 
PLOT_W_CORR  <- 8.0
PLOT_H_CORR  <- 4.8
PLOT_W_FCST  <- 8.0
PLOT_H_FCST  <- 5.5
PLOT_W_PARAM <- 8.0
PLOT_H_PARAM <- 4.5


# Line widths
LW_CORR  <- 0.55
LW_FCST  <- 0.40
LW_PARAM <- 0.55

# ---------------------------
# Helpers
# ---------------------------

dir_create <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

msfe <- function(y, yhat) mean((y - yhat)^2, na.rm = TRUE)

corr_safe <- function(a, b, min_obs = 10) {
  ok <- is.finite(a) & is.finite(b)
  if (sum(ok) < min_obs) return(NA_real_)
  suppressWarnings(cor(a[ok], b[ok]))
}

rolling_corr <- function(x, y, window, min_obs = 10) {
  n <- length(x)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    lo <- max(1, i - window + 1)
    hi <- i
    out[i] <- corr_safe(x[lo:hi], y[lo:hi], min_obs = min_obs)
  }
  out
}

expanding_corr <- function(x, y, min_obs = 10) {
  n <- length(x)
  out <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    out[i] <- corr_safe(x[1:i], y[1:i], min_obs = min_obs)
  }
  out
}

as_date_robust <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  suppressWarnings(as.Date(as.character(x)))
}

read_results_csv <- function(path) {
  if (!file.exists(path)) stop("Missing input file: ", path)
  df <- read_csv(path, show_col_types = FALSE) %>% as.data.frame()
  if ("date_pred" %in% names(df)) df$date_pred <- as_date_robust(df$date_pred)
  if ("date_target" %in% names(df)) df$date_target <- as_date_robust(df$date_target)
  df
}

read_results_rds <- function(path) {
  if (!file.exists(path)) return(NULL)
  obj <- readRDS(path)
  if (!is.list(obj)) stop("RDS file is not a list: ", path)

  ip_key  <- if ("IP"  %in% names(obj)) "IP"  else if ("ip"  %in% names(obj)) "ip"  else NULL
  cpi_key <- if ("CPI" %in% names(obj)) "CPI" else if ("cpi" %in% names(obj)) "cpi" else NULL
  if (is.null(ip_key) || is.null(cpi_key)) stop("RDS does not contain IP/CPI entries: ", path)

  ip_df  <- as.data.frame(obj[[ip_key]])
  cpi_df <- as.data.frame(obj[[cpi_key]])

  if ("date_pred" %in% names(ip_df))   ip_df$date_pred   <- as_date_robust(ip_df$date_pred)
  if ("date_target" %in% names(ip_df)) ip_df$date_target <- as_date_robust(ip_df$date_target)
  if ("date_pred" %in% names(cpi_df))   cpi_df$date_pred   <- as_date_robust(cpi_df$date_pred)
  if ("date_target" %in% names(cpi_df)) cpi_df$date_target <- as_date_robust(cpi_df$date_target)

  list(IP = ip_df, CPI = cpi_df)
}

require_cols <- function(df, cols, tag = "") {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0) stop("Missing columns ", tag, ": ", paste(miss, collapse = ", "))
}

# Save PNG only
save_plot_png <- function(p, file_base, w = 9, h = 4, dpi = 300) {
  ggsave(paste0(file_base, ".png"), p, width = w, height = h, dpi = dpi)
}

period_label <- function(name, rng) {
  paste0(name, " [", format(rng[1], "%Y-%m"), " to ", format(rng[2], "%Y-%m"), "]")
}

# ---------------------------
# Core outputs per target
# ---------------------------

make_outputs_for_target <- function(df, target_tag = c("IP","INF")) {
  target_tag <- match.arg(target_tag)
  target_label <- ifelse(target_tag == "IP", "IP (log level)", "Inflation (YoY)")

  fig_dir <- file.path(OUT_DIR, paste0("figures_", target_tag))
  tab_dir <- file.path(OUT_DIR, paste0("tables_", target_tag))
  dir_create(fig_dir); dir_create(tab_dir)

  base_cols <- c("date_pred", "date_target", "y_true",
                 "yhat_rw", "yhat_pcr", "yhat_pls", "yhat_ridge", "yhat_lasso", "yhat_enet")
  require_cols(df, base_cols, tag = paste0("(", target_tag, ")"))

  df <- df %>%
    mutate(date_pred = as_date_robust(date_pred),
           date_target = as_date_robust(date_target)) %>%
    arrange(date_pred) %>%
    mutate(target_tag = target_tag)

  # 1) Corr path to PCR (static rolling) using yhat (EXCLUDE RW)
  corr_df <- df %>%
    transmute(date_pred = date_pred,
              pcr = yhat_pcr, pls = yhat_pls, ridge = yhat_ridge, lasso = yhat_lasso, enet = yhat_enet)

  corr_static <- data.frame(date_pred = corr_df$date_pred)
  corr_static$corr_pls   <- rolling_corr(corr_df$pls,   corr_df$pcr, window = CORR_WINDOW, min_obs = CORR_MIN_OBS)
  corr_static$corr_ridge <- rolling_corr(corr_df$ridge, corr_df$pcr, window = CORR_WINDOW, min_obs = CORR_MIN_OBS)
  corr_static$corr_lasso <- rolling_corr(corr_df$lasso, corr_df$pcr, window = CORR_WINDOW, min_obs = CORR_MIN_OBS)
  corr_static$corr_enet  <- rolling_corr(corr_df$enet,  corr_df$pcr, window = CORR_WINDOW, min_obs = CORR_MIN_OBS)

  corr_static_long <- corr_static %>%
    pivot_longer(-date_pred, names_to = "model", values_to = "corr") %>%
    mutate(model = gsub("^corr_", "", model),
           model = factor(model, levels = c("pls","ridge","lasso","enet"),
                          labels = c("PLS","Ridge","Lasso","Elastic Net")))

  p1 <- ggplot(corr_static_long, aes(x = date_pred, y = corr, color = model)) +
    geom_line(linewidth = LW_CORR, alpha = 0.9, na.rm = TRUE) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    labs(title = paste0(target_label, ": Corr path to PCR"),
         subtitle = paste0("Corr of yhat with PCR; rolling window = ", CORR_WINDOW),
         x = "Prediction date (t)", y = "Corr(yhat_model, yhat_PCR)", color = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = c(0.98, 0.05),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.key.height = unit(0.35, "cm"),
          legend.text = element_text(size = 8),
          plot.title = element_text(face = "bold"))
  save_plot_png(p1, file.path(fig_dir, paste0("plot1_corr_to_pcr_static_", target_tag)),
                w = PLOT_W_CORR, h = PLOT_H_CORR)

  # 2) Corr path to PCR (cumulative) using yhat (EXCLUDE RW)
  corr_cum <- data.frame(date_pred = corr_df$date_pred)
  corr_cum$corr_pls   <- expanding_corr(corr_df$pls,   corr_df$pcr, min_obs = CORR_MIN_OBS)
  corr_cum$corr_ridge <- expanding_corr(corr_df$ridge, corr_df$pcr, min_obs = CORR_MIN_OBS)
  corr_cum$corr_lasso <- expanding_corr(corr_df$lasso, corr_df$pcr, min_obs = CORR_MIN_OBS)
  corr_cum$corr_enet  <- expanding_corr(corr_df$enet,  corr_df$pcr, min_obs = CORR_MIN_OBS)

  corr_cum_long <- corr_cum %>%
    pivot_longer(-date_pred, names_to = "model", values_to = "corr") %>%
    mutate(model = gsub("^corr_", "", model),
           model = factor(model, levels = c("pls","ridge","lasso","enet"),
                          labels = c("PLS","Ridge","Lasso","Elastic Net")))

  p2 <- ggplot(corr_cum_long, aes(x = date_pred, y = corr, color = model)) +
    geom_line(linewidth = LW_CORR, alpha = 0.9, na.rm = TRUE) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    labs(title = paste0(target_label, ": Corr path to PCR (cumulative)"),
         subtitle = "Expanding-window corr of yhat with PCR",
         x = "Prediction date (t)", y = "Cum Corr(yhat_model, yhat_PCR)", color = NULL) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank(),
          legend.position = c(0.98, 0.05),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
          legend.key.height = unit(0.35, "cm"),
          legend.text = element_text(size = 8),
          plot.title = element_text(face = "bold"))
  save_plot_png(p2, file.path(fig_dir, paste0("plot2_corr_to_pcr_cum_", target_tag)),
                w = PLOT_W_CORR, h = PLOT_H_CORR)

  # 3) Forecast paths (4 figures): thinner + dotted for forecasts
  make_forecast_path_plot <- function(model_key, model_label) {
    dfp <- df %>%
      transmute(date = date_target,
                True = y_true,
                RW = yhat_rw,
                PCR = yhat_pcr,
                Model = .data[[paste0("yhat_", model_key)]]) %>%
      pivot_longer(-date, names_to = "series", values_to = "value") %>%
      mutate(series = factor(series, levels = c("True","RW","PCR","Model")))

    color_map <- c(True = "black", RW = "blue", PCR = "green4", Model = "red")
    lty_map   <- c(True = "solid", RW = "dotted", PCR = "dotted", Model = "dotted")

    ggplot(dfp, aes(x = date, y = value, color = series, linetype = series)) +
      geom_line(linewidth = LW_FCST, na.rm = TRUE) +
      scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
      scale_color_manual(values = color_map) +
      scale_linetype_manual(values = lty_map) +
      labs(title = paste0(target_label, ": Forecast paths (Model = ", model_label, ")"),
           subtitle = "True (black solid) vs RW (blue dotted) vs PCR (green dotted) vs Model (red dotted)",
           x = "Target date (t+h)", y = "Level", color = NULL, linetype = NULL) +
      theme_minimal(base_size = 12) +
      theme(panel.grid.minor = element_blank(),
            plot.title = element_text(face = "bold"),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.box = "horizontal",
            legend.text = element_text(size = 9),
            plot.margin = margin(5.5, 5.5, 16, 5.5)) +
      guides(color = guide_legend(nrow = 1, byrow = TRUE),
             linetype = guide_legend(nrow = 1, byrow = TRUE))
  }

  save_plot_png(make_forecast_path_plot("pls","PLS"), file.path(fig_dir, paste0("plot3_forecast_path_PLS_", target_tag)),
                w = PLOT_W_FCST, h = PLOT_H_FCST)
  save_plot_png(make_forecast_path_plot("ridge","Ridge"), file.path(fig_dir, paste0("plot3_forecast_path_Ridge_", target_tag)),
                w = PLOT_W_FCST, h = PLOT_H_FCST)
  save_plot_png(make_forecast_path_plot("lasso","Lasso"), file.path(fig_dir, paste0("plot3_forecast_path_Lasso_", target_tag)),
                w = PLOT_W_FCST, h = PLOT_H_FCST)
  save_plot_png(make_forecast_path_plot("enet","Elastic Net"), file.path(fig_dir, paste0("plot3_forecast_path_Enet_", target_tag)),
                w = PLOT_W_FCST, h = PLOT_H_FCST)

  # 4) Hyperparameters paths: one parameter per plot (PNG)
  param_cols <- c("pcr_r","pls_ncomp",
                  "ridge_lambda","ridge_df",
                  "lasso_lambda","lasso_k",
                  "enet_alpha","enet_lambda","enet_k")
  param_cols_exist <- intersect(param_cols, names(df))
  if (length(param_cols_exist) > 0) {
    param_df <- df %>% select(date_pred, all_of(param_cols_exist)) %>% arrange(date_pred)
    
    for (pc in param_cols_exist) {
      dplot <- param_df %>% transmute(date = date_pred, value = .data[[pc]])
      use_log10 <- grepl("lambda", pc, ignore.case = TRUE)
      
      p <- ggplot(dplot, aes(x = date, y = value)) +
        geom_point(size = 0.35, alpha = 0.6, shape = 16, na.rm = TRUE)+
        scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
        labs(title = paste0(target_label, ": Hyperparameter path - ", pc),
             subtitle = "Scatter plot of window-specific selected hyperparameters",
             x = "Prediction date (t)", y = pc) +
        theme_minimal(base_size = 12) +
        theme(panel.grid.minor = element_blank(),
              plot.title = element_text(face = "bold"))
      if (use_log10) p <- p + scale_y_continuous(trans = "log10", labels = label_number())
      
      save_plot_png(p, file.path(fig_dir, paste0("plot4_hyperparam_", pc, "_", target_tag)),
                    w = PLOT_W_PARAM, h = PLOT_H_PARAM)
    }
    
    write_csv(param_df, file.path(tab_dir, paste0("table8_param_path_", target_tag, ".csv")))
  }
  

  # 5) Table: corr_path_to_pcr (all times) using yhat (EXCLUDE RW)
  corr_all <- tibble(
    model = c("PLS","Ridge","Lasso","Elastic Net"),
    corr_to_pcr = c(
      corr_safe(df$yhat_pls, df$yhat_pcr, min_obs = CORR_MIN_OBS),
      corr_safe(df$yhat_ridge, df$yhat_pcr, min_obs = CORR_MIN_OBS),
      corr_safe(df$yhat_lasso, df$yhat_pcr, min_obs = CORR_MIN_OBS),
      corr_safe(df$yhat_enet, df$yhat_pcr, min_obs = CORR_MIN_OBS)
    )
  )
  write_csv(corr_all, file.path(tab_dir, paste0("table5_corr_path_to_pcr_alltimes_", target_tag, ".csv")))

  # 6) Table: corr_subsamples (EXCLUDE RW)
  corr_sub <- lapply(names(PERIODS), function(pn) {
    rng <- PERIODS[[pn]]
    sub <- df %>% filter(date_pred >= rng[1], date_pred <= rng[2])
    tibble(period = pn,
           period_desc = period_label(pn, rng),
           corr_pls = corr_safe(sub$yhat_pls, sub$yhat_pcr, min_obs = CORR_MIN_OBS),
           corr_ridge = corr_safe(sub$yhat_ridge, sub$yhat_pcr, min_obs = CORR_MIN_OBS),
           corr_lasso = corr_safe(sub$yhat_lasso, sub$yhat_pcr, min_obs = CORR_MIN_OBS),
           corr_enet = corr_safe(sub$yhat_enet, sub$yhat_pcr, min_obs = CORR_MIN_OBS),
           n = nrow(sub))
  }) %>% bind_rows()
  write_csv(corr_sub, file.path(tab_dir, paste0("table6_corr_subsamples_", target_tag, ".csv")))

  # 7) Table: MSFE relative to RW (RW baseline retained)
  msfe_one_period <- function(sub_df, pn) {
    ms_rw <- msfe(sub_df$y_true, sub_df$yhat_rw)
    ms_pcr <- msfe(sub_df$y_true, sub_df$yhat_pcr)
    ms_pls <- msfe(sub_df$y_true, sub_df$yhat_pls)
    ms_ridge <- msfe(sub_df$y_true, sub_df$yhat_ridge)
    ms_lasso <- msfe(sub_df$y_true, sub_df$yhat_lasso)
    ms_enet <- msfe(sub_df$y_true, sub_df$yhat_enet)

    tibble(period = pn,
           msfe_rw = ms_rw,
           msfe_pcr = ms_pcr, rel_pcr = ms_pcr/ms_rw,
           msfe_pls = ms_pls, rel_pls = ms_pls/ms_rw,
           msfe_ridge = ms_ridge, rel_ridge = ms_ridge/ms_rw,
           msfe_lasso = ms_lasso, rel_lasso = ms_lasso/ms_rw,
           msfe_enet = ms_enet, rel_enet = ms_enet/ms_rw,
           n = nrow(sub_df))
  }

  msfe_tab <- lapply(names(PERIODS), function(pn) {
    rng <- PERIODS[[pn]]
    sub <- df %>% filter(date_pred >= rng[1], date_pred <= rng[2])
    msfe_one_period(sub, pn)
  }) %>% bind_rows() %>%
    left_join(tibble(period = names(PERIODS),
                     period_desc = sapply(names(PERIODS), function(pn) period_label(pn, PERIODS[[pn]]))),
              by = "period") %>%
    select(period, period_desc, everything())
  write_csv(msfe_tab, file.path(tab_dir, paste0("table7_msfe_relative_to_rw_", target_tag, ".csv")))

  # 8) Table: parameter summary (all models; subsample stats)
  #    Report mean + quartiles (p25/p75) + dispersion for key hyperparameters.
  summarize_param <- function(x) {
    x <- x[is.finite(x)]
    if (length(x) == 0) {
      return(tibble(n = 0L, mean = NA_real_, median = NA_real_, sd = NA_real_,
                    p25 = NA_real_, p75 = NA_real_, min = NA_real_, max = NA_real_))
    }
    tibble(
      n = as.integer(length(x)),
      mean = mean(x),
      median = median(x),
      sd = if (length(x) > 1) sd(x) else NA_real_,
      p25 = as.numeric(quantile(x, 0.25, type = 7)),
      p75 = as.numeric(quantile(x, 0.75, type = 7)),
      min = min(x),
      max = max(x)
    )
  }

  # Define which (model, parameter) pairs to summarize (only if the columns exist)
  param_map <- list(
    `PCR` = c("pcr_r"),
    `PLS` = c("pls_ncomp"),
    `Ridge` = c("ridge_lambda", "ridge_df"),
    `Lasso` = c("lasso_lambda", "lasso_k"),
    `Elastic Net` = c("enet_alpha", "enet_lambda", "enet_k")
  )

  # Collect stats by period, model, parameter
  stats_list <- list()
  for (pn in names(PERIODS)) {
    rng <- PERIODS[[pn]]
    sub <- df %>% filter(date_pred >= rng[1], date_pred <= rng[2])

    for (model_name in names(param_map)) {
      cols <- param_map[[model_name]]
      cols <- intersect(cols, names(sub))
      if (length(cols) == 0) next

      for (pc in cols) {
        s <- summarize_param(sub[[pc]])
        stats_list[[length(stats_list) + 1]] <- tibble(
          period = pn,
          model = model_name,
          param = pc
        ) %>% bind_cols(s)
      }
    }
  }

  if (length(stats_list) > 0) {
    param_stats <- bind_rows(stats_list) %>%
      left_join(
        tibble(
          period = names(PERIODS),
          period_desc = sapply(names(PERIODS), function(pn) period_label(pn, PERIODS[[pn]]))
        ),
        by = "period"
      ) %>%
      select(period, period_desc, model, param, everything())

    write_csv(param_stats, file.path(tab_dir, paste0("table9_stats_", target_tag, ".csv")))
  }

  invisible(TRUE)
}

# ---------------------------
# Run
# ---------------------------

dir_create(OUT_DIR)

rds_obj <- read_results_rds(RDS_FILE)
if (!is.null(rds_obj)) {
  message("Loaded RDS: ", RDS_FILE)
  ip_df  <- rds_obj$IP
  cpi_df <- rds_obj$CPI
} else {
  message("RDS not found. Falling back to CSV inputs.")
  ip_df  <- read_results_csv(IP_FILE)
  cpi_df <- read_results_csv(CPI_FILE)
}

message("Generating outputs for IP...")
make_outputs_for_target(ip_df, "IP")

message("Generating outputs for Inflation (YoY)...")
make_outputs_for_target(cpi_df, "INF")

message("All done. Outputs saved under: ", OUT_DIR)
