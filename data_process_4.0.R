library(readxl)
library(openxlsx)
library(dplyr)
library(tseries)
library(tidyr)
library(ggplot2)

# working directory
setwd("------")

file_name  <- "data.xlsx"
sheet_name <- "data_raw"

START_DATE  <- as.Date("1961-01-01")
DATE_FMT <- "%m/%d/%Y"

RAW_IP_COL  <- "INDPRO"
RAW_CPI_COL <- "CPIAUCSL"

H <- 12
L <- 12
small <- 1e-6

# Variables to drop
drop_vars  <- c("ACOGNO","ANDENOx","TWEXAFEGSMTHx","UMCSENTx","VIXCLSx")
drop_after <- c("NONBORRES")

# helper: normalize monthly dates to month-start
to_month_start <- function(d) as.Date(format(as.Date(d), "%Y-%m-01"))


# 1) Read raw
raw <- read_excel(file_name, sheet = sheet_name, col_names = FALSE)

# parse header rows
var_names <- as.character(unlist(raw[1, ]))
transcode <- suppressWarnings(as.numeric(unlist(raw[2, ])))
vars_all   <- var_names[-1]
tcodes_all <- transcode[-1]

# drop vars
keep_idx <- !(vars_all %in% drop_vars)
vars     <- vars_all[keep_idx]
tcodes   <- tcodes_all[keep_idx]

# extract sasdate + X
sasdate <- as.vector(unlist(raw[-c(1, 2), 1]))

X <- raw[-c(1, 2), -1]
X <- as.data.frame(lapply(X, function(z) suppressWarnings(as.numeric(z))))
colnames(X) <- vars_all
X <- X[, vars, drop = FALSE]

# tcode map
tcode_vec <- setNames(tcodes, vars)


# 2) Transform predictors by TCODE
X_trans <- X %>%
  mutate(across(all_of(vars),
                ~{
                  tc <- tcode_vec[[cur_column()]]
                  x  <- as.numeric(.x)
                  
                  # guard for log transforms with nonpositive values
                  if (tc %in% c(4, 5, 6) && any(x <= small, na.rm = TRUE)) {
                    return(rep(NA_real_, length(x)))
                  }
                  
                  d1  <- x - lag(x, L)
                  d2  <- d1 - lag(d1, L)
                  lx  <- log(x)
                  dl  <- lx - lag(lx, L)
                  ddl <- dl - lag(dl, L)
                  g   <- x / lag(x, L) - 1
                  dg  <- g - lag(g, L)
                  
                  switch(as.character(tc),
                         "1" = x, "2" = d1, "3" = d2, "4" = lx, "5" = dl, "6" = ddl, "7" = dg, NA_real_)
                }))

# drop first 24 and last 3
n_tot <- nrow(X_trans)
keep_idx <- 25:(n_tot - 3)

X_trans <- X_trans[keep_idx, , drop = FALSE]

# IMPORTANT: make sasdate a Date (month-start) to ensure perfect alignment later
sasdate_t_chr <- sasdate[keep_idx]
sasdate_t_date <- suppressWarnings(as.Date(sasdate_t_chr, format = DATE_FMT))
if (all(is.na(sasdate_t_date))) sasdate_t_date <- as.Date(sasdate_t_chr)
sasdate_t_date <- to_month_start(sasdate_t_date)

# write data_trans (sasdate as Date)
out_trans <- data.frame(sasdate = sasdate_t_date, X_trans, check.names = FALSE)

wb <- loadWorkbook(file_name)
if ("data_trans" %in% names(wb)) removeWorksheet(wb, "data_trans")
addWorksheet(wb, "data_trans")
writeData(wb, "data_trans", out_trans, startRow = 1, colNames = TRUE)


# 3) ADF tests (diagnostic)
adf_stationary <- function(x, alpha = 0.1, k = 1) {
  x <- x[is.finite(x)]
  if (length(x) < 30) return(list(stat = FALSE, p = NA_real_))
  if (sd(x) == 0)     return(list(stat = FALSE, p = NA_real_))
  res <- tryCatch(adf.test(x, k = k), error = function(e) NULL)
  if (is.null(res))   return(list(stat = FALSE, p = NA_real_))
  list(stat = (res$p.value < alpha), p = res$p.value)
}

adf_list <- lapply(vars, function(v) adf_stationary(out_trans[[v]], alpha = 0.1, k = 1))
adf_test <- data.frame(
  variable   = vars,
  stationary = sapply(adf_list, `[[`, "stat"),
  p_value    = sapply(adf_list, `[[`, "p"),
  check.names = FALSE
)

if ("adf_test" %in% names(wb)) removeWorksheet(wb, "adf_test")
addWorksheet(wb, "adf_test")
writeData(wb, "adf_test", adf_test)

saveWorkbook(wb, file_name, overwrite = TRUE)


# 4) IQR anomaly diagnostic + plots
iqr_anomaly <- function(x, k = 10) {
  med <- median(x, na.rm = TRUE)
  iqr <- IQR(x, na.rm = TRUE)
  if (iqr == 0 || !is.finite(iqr)) return(rep(FALSE, length(x)))
  abs(x - med) > k * iqr
}

anomaly_flag <- as.data.frame(lapply(X_trans, function(x) iqr_anomaly(x, k = 10)))

anomaly_summary <- data.frame(
  variable = names(anomaly_flag),
  has_anomaly = sapply(anomaly_flag, any, na.rm = TRUE),
  n_anomaly   = sapply(anomaly_flag, function(z) sum(z, na.rm = TRUE))
)

idx <- which(as.matrix(anomaly_flag), arr.ind = TRUE)
if (nrow(idx) == 0) {
  cat("No anomalies flagged (k = 10).\n")
} else {
  flagged_df <- data.frame(
    variable = colnames(X_trans)[idx[, 2]],
    date     = sasdate_t_date[idx[, 1]],
    value    = as.matrix(X_trans)[idx],
    row.names = NULL
  )
  print(flagged_df)
}

vars_k10 <- anomaly_summary %>%
  dplyr::filter(has_anomaly) %>%
  dplyr::pull(variable)
print(vars_k10)

out_dir <- "emp_data_plots"
if (!dir.exists(out_dir)) dir.create(out_dir)

date_vec <- sasdate_t_date

if (length(vars_k10) == 0) {
  cat("No flagged variables to plot (k = 10).\n")
} else {
  for (v in vars_k10) {
    df <- data.frame(
      date    = date_vec,
      value   = as.numeric(X_trans[[v]]),
      is_anom = as.logical(anomaly_flag[[v]])
    )
    
    p <- ggplot(df, aes(x = date, y = value)) +
      geom_line(linewidth = 0.5) +
      geom_point(
        data = df[df$is_anom, , drop = FALSE],
        aes(x = date, y = value),
        color = "red", size = 2
      ) +
      scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
      labs(
        title = paste0(v, " (transformed)"),
        subtitle = "Full sample time series with IQR anomalies (k = 10)",
        x = "Date",
        y = "Value"
      ) +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"),
            panel.grid.minor = element_blank())
    
    print(p)
    ggsave(
      filename = file.path(out_dir, paste0("ts_", v, "_k20.png")),
      plot = p, width = 10, height = 4, dpi = 300
    )
  }
  
  cat("Saved plots to folder:", out_dir, "\n")
}


# 5) Save data_predictors
data_predictors <- out_trans

# drop NONBORRES if present
drop_cols_exist <- intersect(drop_after, colnames(data_predictors))
if (length(drop_cols_exist) > 0) {
  data_predictors <- data_predictors %>% dplyr::select(-all_of(drop_cols_exist))
}

# drop raw target columns from predictors if they exist
drop_from_S <- intersect(c(RAW_IP_COL, RAW_CPI_COL), colnames(data_predictors))
if (length(drop_from_S) > 0) {
  data_predictors <- data_predictors %>% dplyr::select(-all_of(drop_from_S))
}

wb <- loadWorkbook(file_name)
if ("data_predictors" %in% names(wb)) removeWorksheet(wb, "data_predictors")
addWorksheet(wb, "data_predictors")
writeData(wb, "data_predictors", data_predictors, startRow = 1, colNames = TRUE)
saveWorkbook(wb, file_name, overwrite = TRUE)


# 6) Build data_target from raw levels
X_full_level <- raw[-c(1, 2), -1]
X_full_level <- as.data.frame(lapply(X_full_level, function(z) suppressWarnings(as.numeric(z))))
colnames(X_full_level) <- vars_all

if (!(RAW_IP_COL %in% colnames(X_full_level)))  stop("RAW_IP_COL not found in raw sheet: ", RAW_IP_COL)
if (!(RAW_CPI_COL %in% colnames(X_full_level))) stop("RAW_CPI_COL not found in raw sheet: ", RAW_CPI_COL)

ip_raw  <- as.numeric(X_full_level[[RAW_IP_COL]])
cpi_raw <- as.numeric(X_full_level[[RAW_CPI_COL]])

# raw dates -> Date -> month-start (CRITICAL)
date_raw <- suppressWarnings(as.Date(sasdate, format = DATE_FMT))
if (all(is.na(date_raw))) date_raw <- as.Date(sasdate)
date_raw <- to_month_start(date_raw)

# guard (rare) nonpositive values for logs
if (any(!is.na(ip_raw) & ip_raw <= 0)) warning("INDPRO has nonpositive values; ip_level will contain NA for those points.")
if (any(!is.na(cpi_raw) & cpi_raw <= 0)) warning("CPIAUCSL has nonpositive values; pi_level will contain NA for those points.")

ip_level <- 100 * log(ip_raw)
pi_level <- 100 * log(cpi_raw / dplyr::lag(cpi_raw, 12))

w_ip_h <- ip_level - dplyr::lag(ip_level, H)
w_pi_h <- pi_level - dplyr::lag(pi_level, H)

data_target_raw <- data.frame(
  date     = date_raw,
  ip_level = ip_level,
  pi_level = pi_level,
  w_ip_h   = w_ip_h,
  w_pi_h   = w_pi_h,
  check.names = FALSE
)

# predictors dates (already Date month-start)
data_predictors2 <- data_predictors %>%
  dplyr::mutate(date = to_month_start(as.Date(sasdate))) %>%
  dplyr::select(date, everything(), -sasdate)

# align targets to predictors time index
data_target <- dplyr::inner_join(
  data_target_raw %>% dplyr::select(date, ip_level, pi_level, w_ip_h, w_pi_h),
  data_predictors2 %>% dplyr::select(date),
  by = "date"
) %>%
  dplyr::filter(date >= START_DATE)

# strong requirement: no NA/Inf in targets
data_target <- data_target %>% dplyr::filter(if_all(everything(), ~ is.finite(.x)))

# extra hard check: common dates must be > 0
common_n <- length(intersect(data_predictors2$date, data_target$date))
if (common_n == 0) stop("No common dates between predictors and targets after normalization. Check DATE_FMT and month anchor.")

# write data_target
wb <- loadWorkbook(file_name)
if ("data_target" %in% names(wb)) removeWorksheet(wb, "data_target")
addWorksheet(wb, "data_target")
writeData(wb, "data_target", data_target, startRow = 1, colNames = TRUE)
saveWorkbook(wb, file_name, overwrite = TRUE)


# 7) Write final panel (predictors + targets)
panel_final <- dplyr::inner_join(data_predictors2, data_target, by = "date") %>%
  dplyr::filter(if_all(everything(), ~ is.finite(.x)))

if (nrow(panel_final) == 0) {
  stop("panel_final is empty after alignment/filters. Check START_DATE and raw availability.")
}

message("Final panel rows: ", nrow(panel_final), " | Date range: ",
        format(min(panel_final$date), "%Y-%m"), " to ", format(max(panel_final$date), "%Y-%m"))

wb <- loadWorkbook(file_name)
if ("data_processed_final" %in% names(wb)) removeWorksheet(wb, "data_processed_final")
addWorksheet(wb, "data_processed_final")
writeData(wb, "data_processed_final", panel_final, startRow = 1, colNames = TRUE)
saveWorkbook(wb, file_name, overwrite = TRUE)

message("Saved sheets: data_predictors, data_target, data_processed_final to workbook: ", file_name)

