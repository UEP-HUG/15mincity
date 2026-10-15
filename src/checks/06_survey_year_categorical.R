# ============================================================================
# Check 06: survey year as linear (current) vs categorical in Model 3
#
# Refits Model 3 for the four hurdle parts (AM/MVPA x binary/Gamma) with
# survey year entered (a) linearly, as in utils.R, and (b) as a factor
# (one level per survey year). Compares the PT coefficient and model fit.
# utils.R is sourced, not modified. Aggregate outputs only:
# results/checks/06_survey_year/.
#
# Run from the repository root:
#   Rscript src/checks/06_survey_year_categorical.R
# ============================================================================
set.seed(12345)
suppressPackageStartupMessages({library(INLA); library(sf); library(dplyr)})
inla.setOption(fmesher.evolution.warn = FALSE)
source("./src/utils.R")

out_dir <- "./results/checks/06_survey_year"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_list <- load_data()
dfs <- prepare_datasets(data_list$df, data_list$df_urban)
meshes <- create_meshes(dfs$coords, dfs$coords_35)
spde_models <- create_spde_models(meshes$mesh1, meshes$mesh1_35)
inla_components <- create_indices_and_projections(
  spde_models, meshes, dfs$coords, dfs$coords_35, dfs$coords_mobility,
  dfs$coords_mobility_35, dfs$coords_mvpa, dfs$coords_mvpa_35)
covariates_list <- prepare_covariates(dfs)
# Add a categorical version of survey year (levels shared across subsets)
yr_levels <- sort(unique(dfs$df$appointment_year))
for (nm in names(covariates_list)) {
  covariates_list[[nm]]$appt_year_f <- factor(covariates_list[[nm]]$appt_year, levels = yr_levels)
}
inla_stacks <- create_inla_stacks(dfs, inla_components, covariates_list)

f_year_cat <- function(spde) {
  y ~ 0 + b0 + pt_all + Age + Sex + work_status + educ_lvl + smoking_status +
    bmi_cat + season + appt_year_f + f(spatial.field, model = spde)
}

fit_part <- function(formula, stack, family) {
  common <- list(formula = formula, data = inla.stack.data(stack),
                 control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
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

summarise_fit <- function(m, part, family, year_spec) {
  b <- m$summary.fixed["pt_all", c("mean", "0.025quant", "0.975quant")]
  tr <- if (family == "gamma") function(x) 100 * (exp(x) - 1) else exp
  data.frame(part = part, year_spec = year_spec,
             pt_log_coef = round(b[[1]], 4),
             pt_effect = round(tr(b[[1]]), 3), pt_lo = round(tr(b[[2]]), 3), pt_hi = round(tr(b[[3]]), 3),
             effect_scale = ifelse(family == "gamma", "% change per SD", "OR per SD"),
             dic = round(m$dic$dic, 2), waic = round(m$waic$waic, 2))
}

parts <- list(
  am_binary   = list(stack = inla_stacks$stk_mobility_binary, spde = spde_models$spde_mobility_binary, family = "binomial"),
  am_gamma    = list(stack = inla_stacks$stk_mobility, spde = spde_models$spde_mobility_continuous, family = "gamma"),
  mvpa_binary = list(stack = inla_stacks$stk_mvpa_binary, spde = spde_models$spde_mvpa_binary, family = "binomial"),
  mvpa_gamma  = list(stack = inla_stacks$stk_mvpa, spde = spde_models$spde_mvpa_continuous, family = "gamma"))

res <- list()
for (p in names(parts)) {
  P <- parts[[p]]
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M"), p))
  m_lin <- fit_part(define_model_formulas(P$spde)$f_m3_pt, P$stack, P$family)
  m_cat <- fit_part(f_year_cat(P$spde), P$stack, P$family)
  res[[p]] <- rbind(summarise_fit(m_lin, p, P$family, "linear"),
                    summarise_fit(m_cat, p, P$family, "categorical"))
  print(res[[p]])
}
write.csv(do.call(rbind, res), file.path(out_dir, "survey_year_linear_vs_categorical.csv"), row.names = FALSE)
cat("Done.\n")
