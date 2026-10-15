# ============================================================================
# 15-Minute City and Physical Activity: sensitivity and supplementary models
# ============================================================================
# Model 3 (linear PT, fully adjusted, SPDE spatial field) for the four hurdle
# parts (AM and MVPA, participation and time), refitted under:
#   A. PT x period interaction (2015-2024 vs 2005-2014), full sample
#   B. Sample restricted to 2015-2024
#   C. Category-specific PT (nine categories, one at a time), full sample
# The 35-75 age-restricted analysis stays in "4. INLA Models.R".
#
# Run from the repository root, after "4. INLA Models.R" inputs exist:
#   Rscript "src/5. Sensitivity and supplementary models.R"
# ============================================================================
set.seed(12345)
source("./src/load_packages.R")
source("./src/utils.R")

out_dir <- "./results/all_pois/sensitivity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_list <- load_data()
df_full <- data_list$df

# Data, mesh, SPDE and stacks for a (sub)sample, standardized on the full sample
build_components <- function(df_in) {
  dfs <- prepare_datasets(df_in, data_list$df_urban, ref_df = df_full)
  meshes <- create_meshes(dfs$coords, dfs$coords_35)
  spde <- create_spde_models(meshes$mesh1, meshes$mesh1_35)
  comps <- create_indices_and_projections(spde, meshes, dfs$coords, dfs$coords_35,
                                          dfs$coords_mobility, dfs$coords_mobility_35,
                                          dfs$coords_mvpa, dfs$coords_mvpa_35)
  stacks <- create_inla_stacks(dfs, comps, prepare_covariates(dfs))
  list(
    am_binary   = list(stack = stacks$stk_mobility_binary, spde = spde$spde_mobility_binary, family = "binomial"),
    am_gamma    = list(stack = stacks$stk_mobility, spde = spde$spde_mobility_continuous, family = "gamma"),
    mvpa_binary = list(stack = stacks$stk_mvpa_binary, spde = spde$spde_mvpa_binary, family = "binomial"),
    mvpa_gamma  = list(stack = stacks$stk_mvpa, spde = spde$spde_mvpa_continuous, family = "gamma"))
}

fit_part <- function(formula, P, label) {
  if (P$family == "gamma") run_inla_model_gamma(formula, P$stack, model_name = label)
  else run_inla_model_binomial(formula, P$stack, model_name = label)
}

# Coefficient on the reporting scale: % change per SD (gamma) or OR per SD (binomial)
extract_effect <- function(m, coef, family, analysis, part, term = coef) {
  b <- m$summary.fixed[coef, c("mean", "0.025quant", "0.975quant")]
  tr <- if (family == "gamma") function(x) 100 * (exp(x) - 1) else exp
  data.frame(analysis = analysis, part = part, term = term,
             estimate = tr(b[[1]]), lower = tr(b[[2]]), upper = tr(b[[3]]),
             scale = ifelse(family == "gamma", "% change per SD", "OR per SD"))
}

# ----------------------------------------------------------------------------
# A. PT x period interaction, and C. category-specific PT (full sample)
# ----------------------------------------------------------------------------
parts_full <- build_components(df_full)
res_main <- list(); res_period <- list(); res_cat <- list()
categories <- c("physical", "culture", "education", "healthcare", "services",
                "transport", "outdoor", "supplies", "restaurant")

for (p in names(parts_full)) {
  P <- parts_full[[p]]
  fm <- define_model_formulas(P$spde)

  m <- fit_part(fm$f_m3_pt, P, paste(p, "M3"))
  res_main[[p]] <- extract_effect(m, "pt_all", P$family, "main (overall PT)", p, "overall")

  m <- fit_part(fm$f_m3_pt_period, P, paste(p, "M3 x period"))
  res_period[[p]] <- rbind(
    extract_effect(m, "pt_all", P$family, "PT x period", p, "PT, 2005-2014"),
    extract_effect(m, "pt_all:period", P$family, "PT x period", p, "ratio 2015-2024 vs 2005-2014"))
  write.csv(do.call(rbind, res_period), file.path(out_dir, "sens_pt_by_period.csv"), row.names = FALSE)

  for (k in categories) {
    m <- fit_part(define_category_formula(P$spde, paste0("pt_", k)), P, paste(p, "M3", k))
    res_cat[[paste(p, k)]] <- extract_effect(m, paste0("pt_", k), P$family, "category PT", p, k)
    rm(m); gc()
  }
  write.csv(rbind(do.call(rbind, res_main), do.call(rbind, res_cat)),
            file.path(out_dir, "supp_category_pt.csv"), row.names = FALSE)
}

# ----------------------------------------------------------------------------
# B. Sample restricted to 2015-2024 (appointment_year = year - 1997)
# ----------------------------------------------------------------------------
parts_recent <- build_components(df_full[df_full$appointment_year >= 18, ])
res_recent <- list()
for (p in names(parts_recent)) {
  P <- parts_recent[[p]]
  m <- fit_part(define_model_formulas(P$spde)$f_m3_pt, P, paste(p, "M3 2015-2024"))
  res_recent[[p]] <- extract_effect(m, "pt_all", P$family, "2015-2024 only", p, "overall")
  rm(m); gc()
}
write.csv(rbind(do.call(rbind, res_main), do.call(rbind, res_recent)),
          file.path(out_dir, "sens_2015_2024.csv"), row.names = FALSE)

# ----------------------------------------------------------------------------
# Supplementary figure: category-specific PT vs overall PT, four parts
# ----------------------------------------------------------------------------
cat_res <- rbind(do.call(rbind, res_main), do.call(rbind, res_cat))
part_titles <- c(am_binary = "A) Non-occupational AM, participation (OR per SD)",
                 am_gamma = "B) Non-occupational AM, time (% change per SD)",
                 mvpa_binary = "C) Leisure-time MVPA, participation (OR per SD)",
                 mvpa_gamma = "D) Leisure-time MVPA, time (% change per SD)")
cat_res$panel <- factor(part_titles[cat_res$part], levels = part_titles)
cat_labels <- c(overall = "Composite", physical = "Physical activity facilities", culture = "Culture",
                education = "Education", healthcare = "Healthcare", services = "Services",
                transport = "Transport", outdoor = "Outdoor spaces", supplies = "Daily supplies",
                restaurant = "Restaurants")
cat_res$composite <- cat_res$term == "overall"
cat_res$term <- factor(cat_labels[cat_res$term], levels = rev(cat_labels[c("overall", categories)]))
cat_res$null <- ifelse(cat_res$scale == "OR per SD", 1, 0)
supp_plot <- ggplot(cat_res, aes(x = estimate, y = term)) +
  geom_vline(aes(xintercept = null), linetype = "dashed", colour = "gray50") +
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, colour = "gray30") +
  geom_point(aes(shape = composite), size = 2) +
  scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 16), guide = "none") +
  facet_wrap(~ panel, scales = "free_x", ncol = 2) +
  labs(x = "Association per SD increase in proximity time", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(strip.text = element_text(hjust = 0), panel.grid.minor = element_blank())
ggsave(file.path(out_dir, "supp_category_pt_forest.png"), supp_plot, width = 10, height = 7, dpi = 300)
ggsave(file.path(out_dir, "supp_category_pt_forest.pdf"), supp_plot, width = 10, height = 7)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_R.txt"))
cat("Done.\n")
