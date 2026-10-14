# ============================================================================
# Check 07: Model 3 with current working status added (four hurdle parts)
# ============================================================================
# Occupational category records the last occupation, also for participants not
# currently in paid work. This check adds current working status (1/0) to
# Model 3 and compares the PT estimates with the main model.
# Run from the repository root:
#   Rscript -e 'pdf(NULL); source("src/checks/07_working_status.R")'
# ============================================================================
set.seed(12345)
source("./src/load_packages.R")
source("./src/utils.R")

out_dir <- "./results/checks/07_working_status"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_list <- load_data()
dfs <- prepare_datasets(data_list$df, data_list$df_urban)
meshes <- create_meshes(dfs$coords, dfs$coords_35)
spde <- create_spde_models(meshes$mesh1, meshes$mesh1_35)
comps <- create_indices_and_projections(spde, meshes, dfs$coords, dfs$coords_35,
                                        dfs$coords_mobility, dfs$coords_mobility_35,
                                        dfs$coords_mvpa, dfs$coords_mvpa_35)
stacks <- create_inla_stacks(dfs, comps, prepare_covariates(dfs))

parts <- list(
  am_binary   = list(stack = stacks$stk_mobility_binary, spde = spde$spde_mobility_binary, family = "binomial"),
  am_gamma    = list(stack = stacks$stk_mobility, spde = spde$spde_mobility_continuous, family = "gamma"),
  mvpa_binary = list(stack = stacks$stk_mvpa_binary, spde = spde$spde_mvpa_binary, family = "binomial"),
  mvpa_gamma  = list(stack = stacks$stk_mvpa, spde = spde$spde_mvpa_continuous, family = "gamma"))

m3_formula <- function(spde, extra = "") {
  as.formula(paste("y ~ 0 + b0 + pt_all + Age + Sex + work_status", extra,
                   "+ educ_lvl + smoking_status + bmi_cat + season + appt_year + f(spatial.field, model = spde)"),
             env = environment())
}

res <- list()
for (p in names(parts)) {
  P <- parts[[p]]
  fit <- function(f, label) if (P$family == "gamma") run_inla_model_gamma(f, P$stack, model_name = label) else run_inla_model_binomial(f, P$stack, model_name = label)
  tr <- if (P$family == "gamma") function(x) 100 * (exp(x) - 1) else exp
  for (v in c("main", "with working status")) {
    m <- fit(m3_formula(P$spde, if (v == "main") "" else "+ work_current"), paste(p, v))
    b <- m$summary.fixed["pt_all", c("mean", "0.025quant", "0.975quant")]
    w <- if (v == "main") c(NA, NA, NA) else unlist(m$summary.fixed["work_current1", c("mean", "0.025quant", "0.975quant")])
    res[[paste(p, v)]] <- data.frame(part = p, model = v,
      pt_est = tr(b[[1]]), pt_lo = tr(b[[2]]), pt_hi = tr(b[[3]]),
      work_est = tr(w[1]), work_lo = tr(w[2]), work_hi = tr(w[3]),
      dic = m$dic$dic, waic = m$waic$waic,
      scale = ifelse(P$family == "gamma", "% change per SD", "OR per SD"))
    rm(m); gc()
  }
  write.csv(do.call(rbind, res), file.path(out_dir, "working_status_check.csv"), row.names = FALSE)
}
print(do.call(rbind, res), digits = 3)
