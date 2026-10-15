# ============================================================================
# Check 05: corrected Model 4 (RW2) and convenience-frontier test
#
# For the four hurdle parts (AM and MVPA, binary and Gamma), fits:
#   - Model 3 (linear PT),
#   - Model 4 (RW2 on 20 PT groups, default prior, as in utils.R),
#   - Model 4 with a looser PC prior on the RW2 precision (prior sensitivity;
#     scale.model = TRUE, P(sigma > 1) = 0.01).
# For each RW2 fit: curve anchored at the node nearest 15 min (% change for
# Gamma, OR for binary), posterior probability of monotonicity, slope below vs
# above 30 min. Fit comparison: DIC, WAIC, mean log CPO.
#
# Same settings as utils.R (Gamma: int.strategy = "eb"; binomial: defaults),
# plus config = TRUE for posterior sampling. utils.R is sourced, not modified.
# Outputs aggregate results only to results/checks/05_model4/.
#
# Run from the repository root:
#   Rscript src/checks/05_model4_frontier.R
# ============================================================================
set.seed(12345)
suppressPackageStartupMessages({library(INLA); library(sf); library(dplyr)})
inla.setOption(fmesher.evolution.warn = FALSE)
source("./src/utils.R")

out_dir <- "./results/checks/05_model4"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
REF_MIN <- 15
N_SAMPLES <- 1000
RE_NAME <- "inla.group(pt_all, n = 20)"

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

# RW2 with PC prior (prior sensitivity); same covariates as f_m3_pt_nonlinear
pc_formula <- function(spde) {
  y ~ 0 + b0 +
    f(inla.group(pt_all, n = 20), model = "rw2", scale.model = TRUE,
      hyper = list(prec = list(prior = "pc.prec", param = c(1, 0.01)))) +
    Age + Sex + work_status + educ_lvl + smoking_status + bmi_cat + season +
    appt_year + f(spatial.field, model = spde)
}

fit_part <- function(formula, stack, family) {
  common <- list(formula = formula, data = inla.stack.data(stack),
                 control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE, config = TRUE),
                 control.predictor = list(A = inla.stack.A(stack)),
                 control.mode = list(restart = TRUE), verbose = FALSE)
  if (family == "gamma") {
    do.call(inla, c(common, list(family = "gamma", control.family = list(link = "log"),
                                 control.inla = list(int.strategy = "eb", diagonal = 1e-6,
                                                     tolerance = 1e-4, h = 0.001))))
  } else {
    do.call(inla, c(common, list(family = "binomial", control.family = list(link = "logit"))))
  }
}

fit_stats <- function(m, part, model) {
  cpo <- m$cpo$cpo[m$cpo$cpo > 0]
  data.frame(part = part, model = model, dic = m$dic$dic, waic = m$waic$waic,
             mean_log_cpo = mean(log(cpo), na.rm = TRUE),
             pt_linear_coef = if ("pt_all" %in% rownames(m$summary.fixed))
               m$summary.fixed["pt_all", "mean"] else NA)
}

analyse_curve <- function(m, part, prior, family) {
  ids <- m$summary.random[[RE_NAME]]$ID
  minutes <- ids * pt_sd + pt_mean
  k <- length(ids)
  smp <- inla.posterior.sample(N_SAMPLES, m)
  rn <- rownames(smp[[1]]$latent)
  idx <- grep(paste0("^", gsub("([().])", "\\\\\\1", RE_NAME), ":"), rn)
  f <- sapply(smp, function(s) s$latent[idx, 1])
  ref <- which.min(abs(minutes - REF_MIN))
  eff <- exp(sweep(f, 2, f[ref, ], "-"))                     # ratio vs reference
  scale_fn <- if (family == "gamma") function(x) 100 * (x - 1) else identity
  curve <- data.frame(part = part, prior = prior, node = seq_len(k),
                      pt_min = round(minutes, 2), ref_min = round(minutes[ref], 2),
                      effect_mean = round(scale_fn(rowMeans(eff)), 3),
                      effect_lo = round(scale_fn(apply(eff, 1, quantile, 0.025)), 3),
                      effect_hi = round(scale_fn(apply(eff, 1, quantile, 0.975)), 3),
                      effect_scale = ifelse(family == "gamma", "% change", "odds ratio"))
  d <- apply(f, 2, diff)
  seg_mid <- (minutes[-1] + minutes[-k]) / 2
  s_lo <- colMeans(d[seg_mid < 30, , drop = FALSE])
  s_hi <- colMeans(d[seg_mid >= 30, , drop = FALSE])
  diag <- data.frame(part = part, prior = prior,
                     rw2_sd_posterior_mean = round(1 / sqrt(m$summary.hyperpar[
                       grep("Precision for inla.group", rownames(m$summary.hyperpar)), "mean"]), 4),
                     p_monotone_increasing = mean(apply(d, 2, function(x) all(x > 0))),
                     p_monotone_decreasing = mean(apply(d, 2, function(x) all(x < 0))),
                     mean_slope_below30 = signif(mean(s_lo), 3),
                     mean_slope_above30 = signif(mean(s_hi), 3),
                     p_slope_sign_differs = mean(sign(s_lo) != sign(s_hi)),
                     p_slope_steeper_below30 = mean(abs(s_lo) > abs(s_hi)))
  list(curve = curve, diag = diag)
}

parts <- list(
  am_binary   = list(stack = inla_stacks$stk_mobility_binary, spde = spde_models$spde_mobility_binary, family = "binomial"),
  am_gamma    = list(stack = inla_stacks$stk_mobility, spde = spde_models$spde_mobility_continuous, family = "gamma"),
  mvpa_binary = list(stack = inla_stacks$stk_mvpa_binary, spde = spde_models$spde_mvpa_binary, family = "binomial"),
  mvpa_gamma  = list(stack = inla_stacks$stk_mvpa, spde = spde_models$spde_mvpa_continuous, family = "gamma"))

stats <- list(); curves <- list(); diags <- list()
for (p in names(parts)) {
  P <- parts[[p]]
  fm <- define_model_formulas(P$spde)
  cat(sprintf("[%s] %s: linear\n", format(Sys.time(), "%H:%M"), p))
  m_lin <- fit_part(fm$f_m3_pt, P$stack, P$family)
  cat(sprintf("[%s] %s: RW2 default prior\n", format(Sys.time(), "%H:%M"), p))
  m_def <- fit_part(fm$f_m3_pt_nonlinear, P$stack, P$family)
  cat(sprintf("[%s] %s: RW2 PC prior\n", format(Sys.time(), "%H:%M"), p))
  m_pc <- fit_part(pc_formula(P$spde), P$stack, P$family)
  stats[[p]] <- rbind(fit_stats(m_lin, p, "M3_linear"),
                      fit_stats(m_def, p, "M4_rw2_default"),
                      fit_stats(m_pc, p, "M4_rw2_pc"))
  a <- analyse_curve(m_def, p, "default", P$family)
  b <- analyse_curve(m_pc, p, "pc", P$family)
  curves[[p]] <- rbind(a$curve, b$curve)
  diags[[p]] <- rbind(a$diag, b$diag)
  print(stats[[p]]); print(diags[[p]])
  write.csv(do.call(rbind, stats), file.path(out_dir, "fit_comparison.csv"), row.names = FALSE)
  write.csv(do.call(rbind, curves), file.path(out_dir, "rw2_curves_anchored.csv"), row.names = FALSE)
  write.csv(do.call(rbind, diags), file.path(out_dir, "rw2_diagnostics.csv"), row.names = FALSE)
  rm(m_lin, m_def, m_pc); gc()
}

# --- Figure: anchored curves, 2 x 2, default vs PC prior --------------------
titles <- c(am_binary = "AM participation (OR)", am_gamma = "AM time (% change)",
            mvpa_binary = "MVPA participation (OR)", mvpa_gamma = "MVPA time (% change)")
png(file.path(out_dir, "model4_anchored_curves.png"), width = 2000, height = 1700, res = 200)
par(mfrow = c(2, 2), mar = c(4.5, 4.5, 2.5, 1))
for (p in names(curves)) {
  cv <- curves[[p]]; d0 <- cv[cv$prior == "default", ]; d1 <- cv[cv$prior == "pc", ]
  null_line <- ifelse(d0$effect_scale[1] == "odds ratio", 1, 0)
  yl <- range(c(cv$effect_lo, cv$effect_hi))
  plot(d0$pt_min, d0$effect_mean, type = "n", ylim = yl, xlab = "Proximity time (min)",
       ylab = d0$effect_scale[1], main = titles[[p]], adj = 0, cex.main = 1, font.main = 1)
  polygon(c(d0$pt_min, rev(d0$pt_min)), c(d0$effect_lo, rev(d0$effect_hi)),
          col = adjustcolor("#2a78d6", 0.2), border = NA)
  lines(d0$pt_min, d0$effect_mean, col = "#2a78d6", lwd = 2)
  lines(d1$pt_min, d1$effect_lo, col = "#eb6834", lwd = 1, lty = 3)
  lines(d1$pt_min, d1$effect_hi, col = "#eb6834", lwd = 1, lty = 3)
  lines(d1$pt_min, d1$effect_mean, col = "#eb6834", lwd = 2, lty = 2)
  abline(h = null_line, col = "grey60", lwd = 0.8); abline(v = d0$ref_min[1], col = "grey60", lty = 3)
  legend("topleft", bty = "n", cex = 0.75, lwd = 2, lty = c(1, 2), col = c("#2a78d6", "#eb6834"),
         legend = c("Default prior (95% CrI)", "PC prior (dotted: 95% CrI)"))
}
dev.off()
cat("Done.\n")
