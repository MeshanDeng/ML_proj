rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# ------------------------------------------------------------
# Parameter distribution diagnostics (violin + box + jitter)
# Reads rolling results (preferred: all_results_v6.rds)
# Outputs PNG figures by target and parameter family
#
# File name: code_emp_roll_param_7.0.R
# ------------------------------------------------------------

# ---------------------------
# User configuration
# ---------------------------

IN_DIR  <- "emp_outputs_roll"
OUT_DIR <- file.path(IN_DIR, "param_distributions")

RDS_FILE <- file.path(IN_DIR, "all_results_v6.rds")

# Subsample periods (user provided)
PERIODS <- list(
  full = c(as.Date("1900-01-01"), as.Date("2100-12-31")),
  s1   = c(as.Date("1971-01-01"), as.Date("1984-12-31")),
  s2   = c(as.Date("1985-01-01"), as.Date("2007-12-31")),
  s3   = c(as.Date("2008-01-01"), as.Date("2019-12-31"))
)

PERIOD_LABELS <- c(
  full = "1971--2019",
  s1   = "1971--1984",
  s2   = "1985--2007",
  s3   = "2008--2019"
)

# Plot sizing
PLOT_W <- 9.0
PLOT_H <- 7.0
DPI    <- 300

# Jitter controls (reduce overplotting)
JITTER_WIDTH  <- 0.12
JITTER_HEIGHT <- 0.0
POINT_SIZE    <- 0.35
POINT_ALPHA   <- 0.10

# Violin/box controls
VIOLIN_ALPHA <- 0.35
BOX_WIDTH    <- 0.18
BOX_ALPHA    <- 0.55
BOX_SIZE     <- 0.25

# ---------------------------
# Helpers
# ---------------------------

dir_create <- function(p) if (!dir.exists(p)) dir.create(p, recursive = TRUE, showWarnings = FALSE)

as_date_robust <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  suppressWarnings(as.Date(as.character(x)))
}

read_results_rds <- function(path) {
  if (!file.exists(path)) stop("Missing RDS file: ", path)
  obj <- readRDS(path)
  if (!is.list(obj)) stop("RDS file is not a list: ", path)

  ip_key  <- if ("IP"  %in% names(obj)) "IP"  else if ("ip"  %in% names(obj)) "ip"  else NULL
  cpi_key <- if ("CPI" %in% names(obj)) "CPI" else if ("cpi" %in% names(obj)) "cpi" else NULL
  if (is.null(ip_key) || is.null(cpi_key)) stop("RDS does not contain IP/CPI entries: ", path)

  ip_df  <- as.data.frame(obj[[ip_key]])
  cpi_df <- as.data.frame(obj[[cpi_key]])

  if ("date_pred" %in% names(ip_df))   ip_df$date_pred   <- as_date_robust(ip_df$date_pred)
  if ("date_pred" %in% names(cpi_df))  cpi_df$date_pred  <- as_date_robust(cpi_df$date_pred)

  list(IP = ip_df, INF = cpi_df)
}

save_png <- function(p, file_base, w = 9, h = 6, dpi = 300) {
  ggsave(filename = paste0(file_base, ".png"), plot = p, width = w, height = h, dpi = dpi)
}

add_period <- function(df) {
  d <- df %>% mutate(date_pred = as_date_robust(date_pred))
  out <- lapply(names(PERIODS), function(pn) {
    rng <- PERIODS[[pn]]
    d %>%
      filter(date_pred >= rng[1], date_pred <= rng[2]) %>%
      mutate(period = pn, period_lbl = PERIOD_LABELS[[pn]])
  }) %>% bind_rows()
  out
}

theme_pub <- function() {
  theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
}

# ---------------------------
# Long-format builders
# ---------------------------

build_long_complexity <- function(df, target_label) {
  keep <- intersect(c("date_pred", "pcr_r", "pls_ncomp", "lasso_k", "enet_k"), names(df))
  if (!("date_pred" %in% keep)) stop("Missing date_pred for target: ", target_label)

  d <- df %>% select(all_of(keep))
  d <- add_period(d)
  d <- d %>%
    mutate(target = target_label) %>%
    pivot_longer(cols = -c(date_pred, period, period_lbl, target),
                 names_to = "param", values_to = "value") %>%
    mutate(
      method = case_when(
        param == "pcr_r"     ~ "PCR",
        param == "pls_ncomp" ~ "PLS",
        param == "lasso_k"   ~ "Lasso",
        param == "enet_k"    ~ "Elastic Net",
        TRUE ~ param
      ),
      param_family = case_when(
        param %in% c("pcr_r","pls_ncomp") ~ "Number of components",
        param %in% c("lasso_k","enet_k")  ~ "Sparsity (k)",
        TRUE ~ "Complexity"
      ),
      value = as.numeric(value)
    ) %>%
    filter(is.finite(value))

  if (nrow(d) == 0) return(NULL)
  d
}

build_long_shrinkage <- function(df, target_label) {
  keep <- intersect(c("date_pred", "ridge_lambda", "lasso_lambda", "enet_lambda", "enet_alpha"), names(df))
  if (!("date_pred" %in% keep)) stop("Missing date_pred for target: ", target_label)

  d <- df %>% select(all_of(keep))
  d <- add_period(d)
  d <- d %>%
    mutate(target = target_label) %>%
    pivot_longer(cols = -c(date_pred, period, period_lbl, target),
                 names_to = "param", values_to = "value") %>%
    mutate(
      method = case_when(
        param == "ridge_lambda" ~ "Ridge",
        param == "lasso_lambda" ~ "Lasso",
        param == "enet_lambda"  ~ "Elastic Net",
        param == "enet_alpha"   ~ "Elastic Net",
        TRUE ~ param
      ),
      param_family = case_when(
        grepl("lambda", param) ~ "Penalty (lambda)",
        param == "enet_alpha"  ~ "Mixing (alpha)",
        TRUE ~ "Shrinkage"
      ),
      value = as.numeric(value)
    ) %>%
    filter(is.finite(value))

  if (nrow(d) == 0) return(NULL)
  d
}

# ---------------------------
# Plotters
# ---------------------------

plot_violin_box_jitter <- function(d, title, ylab, log10_y = FALSE) {
  d <- d %>% mutate(
    target = factor(target, levels = c("INF (YoY)", "IP (log)")),
    period_lbl = factor(period_lbl, levels = c("1971--2019","1971--1984","1985--2007","2008--2019"))
  )

  p <- ggplot(d, aes(x = method, y = value)) +
    geom_violin(trim = TRUE, linewidth = 0.2, alpha = VIOLIN_ALPHA, na.rm = TRUE) +
    geom_boxplot(width = BOX_WIDTH, outlier.shape = NA, alpha = BOX_ALPHA,
                 linewidth = BOX_SIZE, na.rm = TRUE) +
    geom_jitter(width = JITTER_WIDTH, height = JITTER_HEIGHT,
                size = POINT_SIZE, alpha = POINT_ALPHA, na.rm = TRUE) +
    facet_grid(period_lbl ~ target, scales = "free_y") +
    labs(title = title, x = NULL, y = ylab) +
    theme_pub()

  if (log10_y) {
    p <- p + scale_y_continuous(trans = "log10", labels = label_number())
  }

  p
}

# ---------------------------
# Run
# ---------------------------

dir_create(OUT_DIR)

res <- read_results_rds(RDS_FILE)

ip_df  <- res$IP  %>% mutate(date_pred = as_date_robust(date_pred))
inf_df <- res$INF %>% mutate(date_pred = as_date_robust(date_pred))

# Target labels (edit here only if you want different text)
TARGET_IP  <- "IP (log)"
TARGET_INF <- "INF (YoY)"

# 1) Complexity distributions
d_comp_ip  <- build_long_complexity(ip_df,  TARGET_IP)
d_comp_inf <- build_long_complexity(inf_df, TARGET_INF)
d_comp <- bind_rows(d_comp_inf, d_comp_ip)

if (!is.null(d_comp) && nrow(d_comp) > 0) {
  d_components <- d_comp %>% filter(param_family == "Number of components")
  d_sparsity   <- d_comp %>% filter(param_family == "Sparsity (k)")

  if (nrow(d_components) > 0) {
    pC <- plot_violin_box_jitter(
      d_components,
      title = "Distribution of selected number of components",
      ylab  = "Selected components",
      log10_y = FALSE
    )
    save_png(pC, file.path(OUT_DIR, "fig_param_dist_components"), w = PLOT_W, h = PLOT_H, dpi = DPI)
  }

  if (nrow(d_sparsity) > 0) {
    pK <- plot_violin_box_jitter(
      d_sparsity,
      title = "Distribution of selected sparsity (k)",
      ylab  = "Number of selected predictors (k)",
      log10_y = FALSE
    )
    save_png(pK, file.path(OUT_DIR, "fig_param_dist_sparsity_k"), w = PLOT_W, h = PLOT_H, dpi = DPI)
  }
}

# 2) Shrinkage distributions (lambda and alpha)
d_sh_ip  <- build_long_shrinkage(ip_df,  TARGET_IP)
d_sh_inf <- build_long_shrinkage(inf_df, TARGET_INF)
d_sh <- bind_rows(d_sh_inf, d_sh_ip)

if (!is.null(d_sh) && nrow(d_sh) > 0) {
  d_lambda <- d_sh %>% filter(param_family == "Penalty (lambda)")
  d_alpha  <- d_sh %>% filter(param_family == "Mixing (alpha)")

  if (nrow(d_lambda) > 0) {
    pL <- plot_violin_box_jitter(
      d_lambda,
      title = "Distribution of selected penalty strength (lambda)",
      ylab  = "Lambda (log scale)",
      log10_y = TRUE
    )
    save_png(pL, file.path(OUT_DIR, "fig_param_dist_lambda"), w = PLOT_W, h = PLOT_H, dpi = DPI)
  }

  if (nrow(d_alpha) > 0) {
    pA <- plot_violin_box_jitter(
      d_alpha,
      title = "Distribution of selected elastic-net mixing parameter (alpha)",
      ylab  = "Alpha",
      log10_y = FALSE
    )
    save_png(pA, file.path(OUT_DIR, "fig_param_dist_alpha"), w = PLOT_W, h = PLOT_H, dpi = DPI)
  }
}

message("Done. Parameter distribution figures saved to: ", OUT_DIR)
