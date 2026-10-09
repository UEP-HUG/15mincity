# ============================================================================
# Check 03: is the "30-minute threshold / reversal" a sum-to-zero artifact?
#
# Refits the Model 3 linear and RW2 non-linear Gamma models for non-occupational
# AM and leisure-time MVPA (positive part), exactly as in utils.R but with
# config = TRUE so that joint posterior samples of the RW2 curve can be drawn.
# Then:
#   - re-anchors the curve at a reference PT (15 min) instead of the
#     sum-to-zero centering,
#   - tests whether the slope changes sign (monotonicity, low vs high PT slope),
#   - compares linear vs RW2 fit (DIC, WAIC, mean log CPO).
#
# Reuses src/utils.R without modifying it. Outputs aggregate results only to
# results/checks/03_rw2/.
#
# Run from the repository root:
#   Rscript src/checks/03_rw2_threshold.R
# ============================================================================
set.seed(12345)
suppressPackageStartupMessages({library(INLA); library(sf); library(dplyr)})
source("./src/utils.R")

out_dir <- "./results/checks/03_rw2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
REF_MIN <- 15
N_SAMPLES <- 1000

# --- Data and model components (identical to "4. INLA Models.R") ------------
data_list <- load_data()
dfs <- prepare_datasets(data_list$df, data_list$df_urban)
meshes <- create_meshes(dfs$coords, dfs$coords_35)
spde_models <- create_spde_models(meshes$mesh1, meshes$mesh1_35)
inla_components <- create_indices_and_projections(
  spde_models, meshes, dfs$coords, dfs$coords_35, dfs$coords_mobility,
  dfs$coords_mobility_35, dfs$coords_mvpa, dfs$coords_mvpa_35)
covariates_list <- prepare_covariates(dfs)
inla_stacks <- create_inla_stacks(dfs, inla_components, covariates_list)

pt_mean <- mean(dfs$df$overall_15min_city_proximity_time)
pt_sd <- sd(dfs$df$overall_15min_city_proximity_time)

# Same settings as run_inla_model_gamma(), plus config = TRUE
fit_gamma <- function(formula, stack) {
  inla(formula, family = "gamma", control.family = list(link = "log"),
       data = inla.stack.data(stack),
       control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
       control.predictor = list(A = inla.stack.A(stack)),
       control.mode = list(restart = TRUE),
       control.inla = list(int.strategy = "eb", diagonal = 1e-6,
                           tolerance = 1e-4, h = 0.001),
       verbose = FALSE)
}

fit_stats <- function(m, label) {
  cpo <- m$cpo$cpo[m$cpo$cpo > 0]
  data.frame(model = label, dic = m$dic$dic, waic = m$waic$waic,
             mean_log_cpo = mean(log(cpo), na.rm = TRUE),
             n_cpo_failure = sum(m$cpo$failure > 0, na.rm = TRUE))
}

analyse_curve <- function(m, outcome, df_pos) {
  re_name <- "inla.group(pt_all, n = 20)"
  ids <- m$summary.random[[re_name]]$ID
  minutes <- ids * pt_sd + pt_mean
  k <- length(ids)

  # Observations per RW2 node (density behind the sum-to-zero centering)
  grp <- inla.group(df_pos$overall_15min_city_proximity_time_std, n = 20)
  n_per_node <- as.integer(table(factor(grp, levels = ids)))

  # Joint posterior samples of the RW2 node values
  smp <- inla.posterior.sample(N_SAMPLES, m)
  rn <- rownames(smp[[1]]$latent)
  idx <- grep(paste0("^", gsub("([().])", "\\\\\\1", re_name), ":"), rn)
  f <- sapply(smp, function(s) s$latent[idx, 1])   # k x N_SAMPLES

  ref <- which.min(abs(minutes - REF_MIN))
  f_anch <- sweep(f, 2, f[ref, ], "-")
  pct <- 100 * (exp(f_anch) - 1)

  curve <- data.frame(
    outcome = outcome, node = seq_len(k), pt_min = round(minutes, 2),
    n_obs = n_per_node,
    centered_mean_pct = round(100 * (exp(m$summary.random[[re_name]]$mean) - 1), 2),
    anchored_mean_pct = round(rowMeans(pct), 2),
    anchored_lo_pct = round(apply(pct, 1, quantile, 0.025), 2),
    anchored_hi_pct = round(apply(pct, 1, quantile, 0.975), 2))

  # Slope diagnostics
  d <- apply(f, 2, diff)                                   # (k-1) x N
  seg_mid <- (minutes[-1] + minutes[-k]) / 2
  lo <- seg_mid < 30; hi <- seg_mid >= 30
  slope_lo <- colMeans(d[lo, , drop = FALSE])
  slope_hi <- colMeans(d[hi, , drop = FALSE])
  centered_cross <- {
    cm <- m$summary.random[[re_name]]$mean
    j <- which(diff(sign(cm)) != 0)[1]
    if (is.na(j)) NA else minutes[j] + (0 - cm[j]) * (minutes[j + 1] - minutes[j]) / (cm[j + 1] - cm[j])
  }
  diag <- data.frame(
    outcome = outcome,
    pt_range_min = round(min(minutes), 2), pt_range_max = round(max(minutes), 2),
    pt_range_midpoint = round((min(minutes) + max(minutes)) / 2, 2),
    centered_curve_zero_crossing = round(centered_cross, 2),
    p_monotone_increasing = mean(apply(d, 2, function(x) all(x > 0))),
    p_monotone_decreasing = mean(apply(d, 2, function(x) all(x < 0))),
    mean_slope_below30_per_node = signif(mean(slope_lo), 3),
    mean_slope_above30_per_node = signif(mean(slope_hi), 3),
    p_slope_sign_differs_below_vs_above30 = mean(sign(slope_lo) != sign(slope_hi)),
    n_segments_slope_CI_excludes0_pos = sum(apply(d, 1, quantile, 0.025) > 0),
    n_segments_slope_CI_excludes0_neg = sum(apply(d, 1, quantile, 0.975) < 0))
  list(curve = curve, diag = diag)
}

outcomes <- list(
  mobility = list(stack = inla_stacks$stk_mobility,
                  spde = spde_models$spde_mobility_continuous,
                  df = dfs$df_mobility),
  mvpa = list(stack = inla_stacks$stk_mvpa,
              spde = spde_models$spde_mvpa_continuous,
              df = dfs$df_mvpa))

stats_all <- list(); curves <- list(); diags <- list()
for (o in names(outcomes)) {
  cat(sprintf("[%s] Fitting %s ...\n", format(Sys.time(), "%H:%M"), o))
  fm <- define_model_formulas(outcomes[[o]]$spde)
  m_lin <- fit_gamma(fm$f_m3_pt, outcomes[[o]]$stack)
  m_nl <- fit_gamma(fm$f_m3_pt_nonlinear, outcomes[[o]]$stack)
  stats_all[[o]] <- rbind(fit_stats(m_lin, paste0(o, "_M3_linear")),
                          fit_stats(m_nl, paste0(o, "_M4_rw2")))
  stats_all[[o]]$pt_coef_linear <- c(m_lin$summary.fixed["pt_all", "mean"], NA)
  res <- analyse_curve(m_nl, o, outcomes[[o]]$df)
  curves[[o]] <- res$curve; diags[[o]] <- res$diag
  print(stats_all[[o]]); print(res$diag)
}

write.csv(do.call(rbind, stats_all), file.path(out_dir, "fit_comparison.csv"), row.names = FALSE)
write.csv(do.call(rbind, curves), file.path(out_dir, "rw2_curves_anchored_15min.csv"), row.names = FALSE)
write.csv(do.call(rbind, diags), file.path(out_dir, "rw2_slope_diagnostics.csv"), row.names = FALSE)

# --- Figure: centered (as published) vs anchored at 15 min -------------------
png(file.path(out_dir, "rw2_centered_vs_anchored.png"), width = 2000, height = 900, res = 200)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))
for (o in names(curves)) {
  cv <- curves[[o]]
  yl <- range(c(cv$anchored_lo_pct, cv$anchored_hi_pct, cv$centered_mean_pct))
  plot(cv$pt_min, cv$anchored_mean_pct, type = "n", ylim = yl,
       xlab = "Proximity time (min)", ylab = "% change in outcome",
       main = ifelse(o == "mvpa", "Leisure-time MVPA", "Non-occupational AM"),
       adj = 0, cex.main = 1, font.main = 1)
  polygon(c(cv$pt_min, rev(cv$pt_min)), c(cv$anchored_lo_pct, rev(cv$anchored_hi_pct)),
          col = adjustcolor("#2a78d6", 0.2), border = NA)
  lines(cv$pt_min, cv$anchored_mean_pct, col = "#2a78d6", lwd = 2)
  lines(cv$pt_min, cv$centered_mean_pct, col = "#eb6834", lwd = 2, lty = 2)
  abline(h = 0, col = "grey60", lwd = 0.8); abline(v = REF_MIN, col = "grey60", lty = 3)
  legend("topleft", bty = "n", cex = 0.8, lwd = 2, lty = c(1, 2),
         col = c("#2a78d6", "#eb6834"),
         legend = c("Anchored at 15 min (95% CrI)", "Sum-to-zero centering (as published)"))
}
dev.off()
cat("Done.\n")
