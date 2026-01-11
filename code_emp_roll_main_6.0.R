rm(list = ls())

library(readxl)
library(dplyr)
library(glmnet)
library(pls)

source("code_emp_roll_func_6.0.R")

XLSX_FILE   <- "data.xlsx"
SHEET_PRED  <- "data_predictors"
SHEET_TARG  <- "data_target"

# Parameters
h       <- 12
W       <- 120
step    <- 1
K       <- 5
min_gap <- 12

r_grid      <- c(1, 3, 5, 10, 25, 50, 75)
lambda_grid <- exp(seq(-6, 6, length.out = 60))
alpha_grid  <- c(0.05, 0.25, 0.5, 0.75, 0.95)

model_list  <- c("pcr","pls","ridge","lasso","enet")

END_DATE <- as.Date("2019-12-01")

OUT_DIR <- "emp_outputs_roll"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
RDS_FILE <- file.path(OUT_DIR, "all_results_v6.rds")

verbose <- TRUE


# Robust date parsing + month anchor
to_month_start <- function(d) as.Date(format(as.Date(d), "%Y-%m-01"))

parse_date_robust <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  if (is.numeric(x)) return(as.Date(x, origin = "1899-12-30"))
  x_chr <- trimws(as.character(x))
  d <- suppressWarnings(as.Date(x_chr, format = "%m/%d/%Y"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x_chr, format = "%Y-%m-%d"))
  if (all(is.na(d))) d <- suppressWarnings(as.Date(x_chr))
  d
}


# Load predictors and targets
P_df <- as.data.frame(read_excel(path = XLSX_FILE, sheet = SHEET_PRED))
T_df <- as.data.frame(read_excel(path = XLSX_FILE, sheet = SHEET_TARG))

if (nrow(P_df) == 0) stop("Sheet ", SHEET_PRED, " is empty.")
if (nrow(T_df) == 0) stop("Sheet ", SHEET_TARG, " is empty.")

if (!("sasdate" %in% names(P_df))) stop("Predictor sheet missing column: sasdate")
if (!("date" %in% names(T_df)))     stop("Target sheet missing column: date")

need_targ <- c("ip_level","pi_level","w_ip_h","w_pi_h")
miss_targ <- setdiff(need_targ, names(T_df))
if (length(miss_targ) > 0) stop("Target sheet missing columns: ", paste(miss_targ, collapse = ", "))

P_df$sasdate <- to_month_start(parse_date_robust(P_df$sasdate))
T_df$date    <- to_month_start(parse_date_robust(T_df$date))

P_df <- P_df[!is.na(P_df$sasdate), , drop = FALSE]
T_df <- T_df[!is.na(T_df$date), , drop = FALSE]

if (nrow(P_df) == 0) stop("All predictor dates became NA after parsing. Check data_predictors$sasdate.")
if (nrow(T_df) == 0) stop("All target dates became NA after parsing. Check data_target$date.")

# truncate sample (paper sample truncation)
P_df <- subset(P_df, sasdate <= END_DATE)
T_df <- subset(T_df, date    <= END_DATE)


# Merge (bulletproof): string key YYYY-mm-dd
P_df$date_key <- format(P_df$sasdate, "%Y-%m-%d")
T_df$date_key <- format(T_df$date,    "%Y-%m-%d")

panel <- merge(P_df, T_df, by = "date_key", all = FALSE)

if (nrow(panel) == 0) {
  stop(
    "Merged panel is empty.\n",
    "Predictor date head: ", paste(head(sort(unique(P_df$date_key)), 5), collapse = ", "), "\n",
    "Target date head:    ", paste(head(sort(unique(T_df$date_key)), 5), collapse = ", "), "\n",
    "Intersection size:   ", length(intersect(unique(P_df$date_key), unique(T_df$date_key))), "\n"
  )
}


panel$date <- as.Date(panel$date_key)
panel <- panel %>%
  dplyr::select(
    -date_key,
    -dplyr::any_of(c("sasdate", "date_P", "date_T"))
  )
panel <- panel[order(panel$date), , drop = FALSE]


message("Merged panel rows: ", nrow(panel),
        " | Date range: ", format(min(panel$date), "%Y-%m"), " to ", format(max(panel$date), "%Y-%m"))



# Define predictors columns (prevent leakage)
target_cols <- c("ip_level","pi_level","w_ip_h","w_pi_h")

# predictors = all numeric columns except targets
is_num <- sapply(panel, is.numeric)
x_cols <- names(panel)[is_num]
x_cols <- setdiff(x_cols, target_cols)

# extra safety
if ("date" %in% x_cols) x_cols <- setdiff(x_cols, "date")
if (length(x_cols) == 0) stop("No predictor columns found after filtering. Check data_predictors sheet.")
if (length(intersect(x_cols, target_cols)) > 0) stop("Leakage: x_cols contains target columns.")

message("Predictors used: ", length(x_cols))


# Run rolling CV forecasts
all_results <- list()

message("\n=== Target: IP (log level) ===")
res_ip <- run_rolling_cv_forecast_v6(
  panel_SW = panel,
  date_col = "date",
  x_cols   = x_cols,
  target   = "ip",
  h = h, W = W, K = K, min_gap = min_gap, step = step,
  r_grid = r_grid,
  lambda_grid = lambda_grid,
  alpha_grid = alpha_grid,
  models = model_list,
  verbose = verbose
)
all_results[["IP"]] <- res_ip
write.csv(res_ip, file.path(OUT_DIR, "results_IP_v6.csv"), row.names = FALSE)

message("\n=== Target: CPI inflation level ===")
res_cpi <- run_rolling_cv_forecast_v6(
  panel_SW = panel,
  date_col = "date",
  x_cols   = x_cols,
  target   = "cpi",
  h = h, W = W, K = K, min_gap = min_gap, step = step,
  r_grid = r_grid,
  lambda_grid = lambda_grid,
  alpha_grid = alpha_grid,
  models = model_list,
  verbose = verbose
)
all_results[["CPI"]] <- res_cpi
write.csv(res_cpi, file.path(OUT_DIR, "results_CPI_v6.csv"), row.names = FALSE)

saveRDS(all_results, RDS_FILE)
message("\nDone. RDS saved: ", RDS_FILE)
message("Outputs in: ", OUT_DIR)
