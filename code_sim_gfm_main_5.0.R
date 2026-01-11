## ==================== 7) Main output (GFM mild & hard) =====================
set.seed(123)
library(MASS)
library(glmnet)
library(pls)
library(dplyr)
library(future)
library(future.apply)

plan(multisession, workers=max(1, parallel::detectCores()-2))

N_VALS    <- c(5,10,20,30,40,50,80,100,200,300,500,1000)
K_RATIO   <- 1
R_REP     <- 200
h         <- 1
test_frac <- 0.2

r_signal <- 1
r_noise <- 10
noise_power <- 0.7

cfg <- list(h=h, test_frac=test_frac, min_train=40, min_T=80, seed_base=700,
            cv_nfolds=10, cv_type_measure="mse",
            enet_alpha_grid=seq(0, 1, by=0.1),
            pls_max_ncomp=20, pls_method="simpls",
            pcr_rmax=20)

cfg_gfm_mild <- cfg; cfg_gfm_mild$alpha_vec <- seq(0.25, 0.45, length.out=r_noise)
cfg_gfm_hard <- cfg; cfg_gfm_hard$alpha_vec <- seq(0.75, 0.95, length.out=r_noise)

cat(">>> Starting simulation 'GFM-MILD'\n")
OUT_gfm_mild <- simulate_grid(n_vals=N_VALS, k_ratio=K_RATIO, R=R_REP, r_signal=r_signal, r_noise=r_noise, noise_power=noise_power, cfg=cfg_gfm_mild, dgp_fun=dgp_gen_gfm)
cat(">>> Simulation 'GFM-MILD' completed.\n")

cat(">>> Starting simulation 'GFM-HARD'\n")
OUT_gfm_hard <- simulate_grid(n_vals=N_VALS, k_ratio=K_RATIO, R=R_REP, r_signal=r_signal, r_noise=r_noise, noise_power=noise_power, cfg=cfg_gfm_hard, dgp_fun=dgp_gen_gfm)
cat(">>> Simulation 'GFM-HARD' completed.\n")

cat(">>> Analyzing results for 'GFM-MILD':\n")
agg_gfm_mild <- aggregate_results(OUT_gfm_mild)
print(agg_gfm_mild)
rates_gfm_mild <- estimate_convergence_rates(agg_gfm_mild)
cat(">>> Convergence rates 'GFM-MILD':\n")
print(rates_gfm_mild)

cat(">>> Analyzing results for 'GFM-HARD':\n")
agg_gfm_hard <- aggregate_results(OUT_gfm_hard)
agg_gfm_hard <- agg_gfm_hard[, !(names(agg_gfm_hard) %in% c("frob_pcr_m","frob_pcr_l","frob_pcr_u")), drop=FALSE]
print(agg_gfm_hard)
rates_gfm_hard <- estimate_convergence_rates(agg_gfm_hard)
cat(">>> Convergence rates 'GFM-HARD':\n")
print(rates_gfm_hard)

# save plots (png) to ./plots
plots_dir <- file.path(getwd(), "gfm_plots")
dir.create(plots_dir, showWarnings = FALSE)
png(filename = file.path(plots_dir, "plot_%03d.png"), width = 1600, height = 1000, res = 200)
on.exit(if(dev.cur() > 1) dev.off(), add=TRUE)

# plot generic panels by prefix
plot_eigenvalue_results(agg_gfm_mild)
plot_model_results(agg_gfm_mild, prefix="pcr", hyper_key="r_hat", hyper_label="r_hat ", main_tag="PCR")
plot_model_results(agg_gfm_mild, prefix="pls", hyper_key="ncomp_pls", hyper_label="ncomp ", main_tag="PLS")
plot_model_results(agg_gfm_mild, prefix="ridge", hyper_key="lambda_best_ridge", hyper_label="lambda ")
plot_model_results(agg_gfm_mild, prefix="lasso", hyper_key="lambda_best_lasso", hyper_label="lambda ", has_sparsity=TRUE, sparsity_key="n_pred_lasso")
plot_model_results(agg_gfm_mild, prefix="enet", hyper_key="lambda_best_enet", hyper_label="lambda ", has_sparsity=TRUE, sparsity_key="n_pred_enet", main_tag="ENet")

plot_eigenvalue_results(agg_gfm_hard)
plot_model_results(agg_gfm_hard, prefix="pcr", hyper_key="r_hat", hyper_label="r_hat ", main_tag="PCR")
plot_model_results(agg_gfm_hard, prefix="pls", hyper_key="ncomp_pls", hyper_label="ncomp ", main_tag="PLS")
plot_model_results(agg_gfm_hard, prefix="ridge", hyper_key="lambda_best_ridge", hyper_label="lambda ")
plot_model_results(agg_gfm_hard, prefix="lasso", hyper_key="lambda_best_lasso", hyper_label="lambda ", has_sparsity=TRUE, sparsity_key="n_pred_lasso")
plot_model_results(agg_gfm_hard, prefix="enet", hyper_key="lambda_best_enet", hyper_label="lambda ", has_sparsity=TRUE, sparsity_key="n_pred_enet", main_tag="ENet")

# plot Enet alpha
plot(agg_gfm_mild$n, agg_gfm_mild$alpha_best_enet_m, log="x", pch=19,
     xlab="n", ylab="alpha_best_enet", main="ENet alpha (CV) vs n - GFM mild")
lines(agg_gfm_mild$n, agg_gfm_mild$alpha_best_enet_m, lwd=2)
plot(agg_gfm_hard$n, agg_gfm_hard$alpha_best_enet_m, log="x", pch=19,
     xlab="n", ylab="alpha_best_enet", main="ENet alpha (CV) vs n - GFM hard")
lines(agg_gfm_hard$n, agg_gfm_hard$alpha_best_enet_m, lwd=2)


# plot agg comparison
plot_mse_all(agg_gfm_mild, main_title="GFM-MILD: Test MSE (all methods)")
plot_corr_vs_pcr_all(agg_gfm_mild, main_title="GFM-MILD: corr(model, PCR)")
plot_convergence_all(agg_gfm_mild, main_title="GFM-MILD: Convergence rate")
plot_gap_all(agg_gfm_mild, main_title="GFM-MILD: Gap MSE (model − oracle)")

plot_mse_all(agg_gfm_hard, main_title="GFM-HARD: Test MSE (all methods)")
plot_corr_vs_pcr_all(agg_gfm_hard, main_title="GFM-HARD: corr(model, PCR)")
plot_convergence_all(agg_gfm_hard, main_title="GFM-HARD: Convergence rate")
plot_gap_all(agg_gfm_hard, main_title="GFM-HARD: Gap MSE (model − oracle)")

plot_pls_proxy(agg_mild=agg_gfm_mild, agg_hard=agg_gfm_hard)


if(dev.cur() > 1) dev.off()
plan(sequential)

# write tables (CSV) to ./tables
out_dir <- file.path(getwd(), "gfm_tables")
dir.create(out_dir, showWarnings = FALSE)

write.table(agg_gfm_mild, file = file.path(out_dir, "agg_gfm_mild.csv"), sep = ",", row.names = FALSE, quote = FALSE)
write.table(rates_gfm_mild, file = file.path(out_dir, "rates_gfm_mild.csv"), sep = ",", row.names = FALSE, quote = FALSE)
write.table(agg_gfm_hard, file = file.path(out_dir, "agg_gfm_hard.csv"), sep = ",", row.names = FALSE, quote = FALSE)
write.table(rates_gfm_hard, file = file.path(out_dir, "rates_gfm_hard.csv"), sep = ",", row.names = FALSE, quote = FALSE)


# Save cached simulation results
cache_dir <- file.path(getwd(), "gfm_cache")
dir.create(cache_dir, showWarnings = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

saveRDS(
  list(
    OUT_gfm_mild = OUT_gfm_mild,
    OUT_gfm_hard = OUT_gfm_hard,
    PRED_gfm_mild = attr(OUT_gfm_mild, "pred"),
    PRED_gfm_hard = attr(OUT_gfm_hard, "pred"),
    agg_gfm_mild = agg_gfm_mild,
    rates_gfm_mild = rates_gfm_mild,
    agg_gfm_hard = agg_gfm_hard,
    rates_gfm_hard = rates_gfm_hard,
    cfg_gfm_mild = cfg_gfm_mild,
    cfg_gfm_hard = cfg_gfm_hard,
    meta = list(N_VALS=N_VALS, K_RATIO=K_RATIO, R_REP=R_REP, h=h, test_frac=test_frac, r_noise=r_noise, noise_power=noise_power, saved_at=Sys.time())
  ),
  file = file.path(cache_dir, paste0("ALL_", run_id, ".rds"))
)

cat("Saved to: ", file.path(cache_dir, paste0("ALL_", run_id, ".rds")), "\n")

## Load cached GFM simulation results to skip the long run and only redo analysis/plots
# obj <- readRDS("gfm_cache/ALL_20260105_123456.rds")

# OUT_gfm_mild <- obj$OUT_gfm_mild
# OUT_gfm_hard <- obj$OUT_gfm_hard
# PRED_gfm_mild <- obj$PRED_gfm_mild
# PRED_gfm_hard <- obj$PRED_gfm_hard
# agg_gfm_mild <- obj$agg_gfm_mild
# rates_gfm_mild <- obj$rates_gfm_mild
# agg_gfm_hard <- obj$agg_gfm_hard
# rates_gfm_hard <- obj$rates_gfm_hard
# cfg_gfm_mild <- obj$cfg_gfm_mild
# cfg_gfm_hard <- obj$cfg_gfm_hard

# plot Enet aplha-single
plots_dir <- file.path(getwd(), "gfm_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

w <- 1600/3
h <- 1000/2
res <- 200

png(filename = file.path(plots_dir, "enet_alpha_mild_singlepanel.png"),
    width = w, height = h, res = res)
par(mar=c(4,4,2,1), cex=0.6)
plot(agg_gfm_mild$n, agg_gfm_mild$alpha_best_enet_m,
     log = "x", pch = 19, xlab = "n", ylab = "alpha_best_enet",
     main = "Hyperparameter (CV): ENet")
lines(agg_gfm_mild$n, agg_gfm_mild$alpha_best_enet_m, lwd = 2)
dev.off()


png(filename = file.path(plots_dir, "enet_alpha_hard_singlepanel.png"),
    width = w, height = h, res = res)
par(mar=c(4,4,2,1), cex=0.6)
plot(agg_gfm_hard$n, agg_gfm_hard$alpha_best_enet_m,
     log = "x", pch = 19, xlab = "n", ylab = "alpha_best_enet",
     main = "Hyperparameter (CV): ENet")
lines(agg_gfm_hard$n, agg_gfm_hard$alpha_best_enet_m, lwd = 2)
dev.off()

cat("Saved to:\n",
    file.path(plots_dir, "enet_alpha_mild_singlepanel.png"), "\n",
    file.path(plots_dir, "enet_alpha_hard_singlepanel.png"), "\n")


