## =============== 7) Main output (AFM: strict & weak) ================
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

cfg <- list(h=h, test_frac=test_frac, min_train=40, min_T=80, seed_base=700,
            cv_nfolds=10, cv_type_measure="mse",
            enet_alpha_grid=seq(0, 1, by=0.1),
            pls_max_ncomp=10, pls_method="simpls",
            pcr_rmax=10)

cat(">>> Starting simulation 'Approximate STRUCTURE'\n")
OUT_strict <- simulate_grid(n_vals=N_VALS, k_ratio=K_RATIO, R=R_REP, r_signal=1, r_noise=0, noise_power=0, cfg=cfg)
cat(">>> Simulation 'Approximate STRUCTURE' completed.\n")

cat(">>> Starting simulation 'WEAK STRUCTURE'\n")
OUT_weak <- simulate_grid(n_vals=N_VALS, k_ratio=K_RATIO, R=R_REP, r_signal=1, r_noise=2, noise_power=0.7, cfg=cfg)
cat(">>> Simulation 'WEAK STRUCTURE' completed.\n")

cat(">>> Analyzing results for 'Approximate STRUCTURE':\n")
agg_strict <- aggregate_results(OUT_strict)
print(agg_strict)
rates_strict <- estimate_convergence_rates(agg_strict)
cat(">>> Convergence rates 'Approximate STRUCTURE':\n")
print(rates_strict)

cat(">>> Analyzing results for 'WEAK STRUCTURE':\n")
agg_weak <- aggregate_results(OUT_weak)
agg_weak <- agg_weak[, !(names(agg_weak) %in% c("frob_pcr_m","frob_pcr_l","frob_pcr_u")), drop=FALSE]
print(agg_weak)
rates_weak <- estimate_convergence_rates(agg_weak)
cat(">>> Convergence rates 'WEAK STRUCTURE':\n")
print(rates_weak)

# save plots (png) to ./plots
plots_dir <- file.path(getwd(), "afm_plots")
dir.create(plots_dir, showWarnings = FALSE)
png(filename = file.path(plots_dir, "plot_%03d.png"), width = 1600, height = 1000, res = 200)
on.exit(if(dev.cur() > 1) dev.off(), add=TRUE)

# plot generic panels by prefix
plot_eigenvalue_results(agg_strict)
plot_model_results(agg_strict, prefix="pcr", hyper_key="r_hat", hyper_label="r_hat ", main_tag="PCR")
plot_model_results(agg_strict, prefix="pls", hyper_key="ncomp_pls", hyper_label="ncomp ", main_tag="PLS")
plot_model_results(agg_strict, prefix="ridge", hyper_key="lambda_best_ridge", hyper_label="lambda ")
plot_model_results(agg_strict, prefix="lasso", hyper_key="lambda_best_lasso", hyper_label="lambda ", has_sparsity=TRUE, sparsity_key="n_pred_lasso")
plot_model_results(agg_strict, prefix="enet", hyper_key="lambda_best_enet", hyper_label="lambda ", has_sparsity=TRUE, sparsity_key="n_pred_enet", main_tag="ENet")

plot_eigenvalue_results(agg_weak)
plot_model_results(agg_weak, prefix="pcr", hyper_key="r_hat", hyper_label="r_hat ", main_tag="PCR")
plot_model_results(agg_weak, prefix="pls", hyper_key="ncomp_pls", hyper_label="ncomp ", main_tag="PLS")
plot_model_results(agg_weak, prefix="ridge", hyper_key="lambda_best_ridge", hyper_label="lambda ")
plot_model_results(agg_weak, prefix="lasso", hyper_key="lambda_best_lasso", hyper_label="lambda ", has_sparsity=TRUE, sparsity_key="n_pred_lasso")
plot_model_results(agg_weak, prefix="enet", hyper_key="lambda_best_enet", hyper_label="lambda ", has_sparsity=TRUE, sparsity_key="n_pred_enet", main_tag="ENet")

# plot ENet alpha
plot(agg_strict$n, agg_strict$alpha_best_enet_m, log="x", pch=19, xlab="n ", ylab="alpha_best_enet", main="ENet alpha (CV) vs n - strict")
lines(agg_strict$n, agg_strict$alpha_best_enet_m, lwd=2)
plot(agg_weak$n, agg_weak$alpha_best_enet_m, log="x", pch=19, xlab="n ", ylab="alpha_best_enet", main="ENet alpha (CV) vs n — weak")
lines(agg_weak$n, agg_weak$alpha_best_enet_m, lwd=2)

# plot agg comparison
plot_mse_all(agg_strict, main_title="STRICT: Test MSE (all methods)")
plot_corr_vs_pcr_all(agg_strict, main_title="STRICT: corr(model, PCR)")
plot_convergence_all(agg_strict, main_title="STRICT: Convergence rate")
plot_gap_all(agg_strict, main_title="STRICT: Gap MSE (model − oracle)")

plot_mse_all(agg_weak, main_title="WEAK: Test MSE (all methods)")
plot_corr_vs_pcr_all(agg_weak, main_title="WEAK: corr(model, PCR)")
plot_convergence_all(agg_weak, main_title="WEAK: Convergence rate")
plot_gap_all(agg_weak, main_title="WEAK: Gap MSE (model − oracle)")

plot_pls_proxy(agg_strict, agg_weak)


if(dev.cur() > 1) dev.off()
plan(sequential)


# write tables (CSV) to ./tables
out_dir <- file.path(getwd(), "afm_tables")
dir.create(out_dir, showWarnings = FALSE)

write.table(agg_strict,  file = file.path(out_dir, "agg_strict.csv"), sep = ",", row.names = FALSE, quote = FALSE)
write.table(rates_strict, file = file.path(out_dir, "rates_strict.csv"), sep = ",", row.names = FALSE, quote = FALSE)
write.table(agg_weak,    file = file.path(out_dir, "agg_weak.csv"), sep = ",", row.names = FALSE, quote = FALSE)
write.table(rates_weak,  file = file.path(out_dir, "rates_weak.csv"), sep = ",", row.names = FALSE, quote = FALSE)

# Save cached simulation results
cache_dir <- file.path(getwd(), "afm_cache")
dir.create(cache_dir, showWarnings = FALSE)
run_id <- format(Sys.time(), "%Y%m%d_%H%M%S")

saveRDS(
  list(
    OUT_strict = OUT_strict,
    OUT_weak   = OUT_weak,
    PRED_strict = attr(OUT_strict, "pred"),
    PRED_weak   = attr(OUT_weak, "pred"),
    agg_strict = agg_strict,
    rates_strict = rates_strict,
    agg_weak   = agg_weak,
    rates_weak = rates_weak,
    cfg = cfg,
    meta = list(
      N_VALS = N_VALS, K_RATIO = K_RATIO, R_REP = R_REP,
      h = h, test_frac = test_frac,
      saved_at = Sys.time()
    )
  ),
  file = file.path(cache_dir, paste0("ALL_", run_id, ".rds"))
)

cat("Saved to: ", file.path(cache_dir, paste0("ALL_", run_id, ".rds")), "\n")

## Load cached simulation results to skip the long run and only redo analysis/plots
#obj <- readRDS("afm_cache/ALL_20260106_104519.rds")  

#OUT_strict <- obj$OUT_strict
#OUT_weak   <- obj$OUT_weak
#PRED_strict <- obj$PRED_strict
#PRED_weak   <- obj$PRED_weak
#agg_strict <- obj$agg_strict
#rates_strict <- obj$rates_strict
#agg_weak   <- obj$agg_weak
#rates_weak <- obj$rates_weak
#plot_files <- obj$plot_files
#cfg <- obj$cfg


# plot ENet alpha-single
plots_dir <- file.path(getwd(), "afm_plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

w <- 1600/3
h <- 1000/2
res <- 200

png(filename = file.path(plots_dir, "enet_alpha_strict_singlepanel.png"),
    width = w, height = h, res = res)

par(mar = c(4,4,2,1), cex = 0.6)
plot(agg_strict$n, agg_strict$alpha_best_enet_m,
     log = "x", pch = 19,xlab = "n", ylab = "alpha_best_enet",
     main = "Hyperparameter (CV): ENet")
lines(agg_strict$n, agg_strict$alpha_best_enet_m, lwd = 2)

dev.off()


png(filename = file.path(plots_dir, "enet_alpha_weak_singlepanel.png"),
    width = w, height = h, res = res)

par(mar = c(4,4,2,1), cex = 0.6)
plot(agg_weak$n, agg_weak$alpha_best_enet_m,
     log = "x", pch = 19, xlab = "n", ylab = "alpha_best_enet",
     main = "Hyperparameter (CV): ENet")
lines(agg_weak$n, agg_weak$alpha_best_enet_m, lwd = 2)

dev.off()

cat("Saved to:\n",
    file.path(plots_dir, "enet_alpha_strict_singlepanel.png"), "\n",
    file.path(plots_dir, "enet_alpha_weak_singlepanel.png"), "\n")


