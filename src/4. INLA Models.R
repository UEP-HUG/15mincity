# ============================================================================
# 15-Minute City and Physical Activity Analysis
# ============================================================================
# Set working directory and source utilities
setwd('/Users/david/Dropbox/PhD/GitHub/15min_city/')
set.seed(12345)

# Load required packages
source("./src/load_packages.R") 
source("./src/utils.R")
# Set result folder path
result_folder <- './results/all_pois'
# ============================================================================
# 1. Data Loading and Preparation
# ============================================================================
df <- read.csv("./data/confidential/gdf_final_1200.csv")
df
# df_urban <- read.csv("./data/gdf_final_1200_urban.csv")
ch_shp <- st_read('./data/input/canton_gecontour.GeoJSON', quiet = TRUE)
canton_ge <- st_read('./data/input/canton_ge.GeoJSON', quiet = TRUE)
communes_shp <- st_read("./data/input/communes_ge.GeoJSON", quiet = TRUE)


data_list <- load_data()
dfs <- prepare_datasets(data_list$df, data_list$df_urban) 

# ============================================================================
# 2. Mesh Construction and SPDE Models
# ============================================================================
meshes <- create_meshes(dfs$coords, dfs$coords_35)

spde_models <- create_spde_models(meshes$mesh1, meshes$mesh1_35)

inla_components <- create_indices_and_projections(spde_models, meshes, dfs$coords, dfs$coords_35, dfs$coords_mobility, dfs$coords_mobility_35, dfs$coords_mvpa, dfs$coords_mvpa_35)

# ============================================================================
# 3. INLA Stack Creation
# ============================================================================
covariates_list <- prepare_covariates(dfs)

inla_stacks <- create_inla_stacks(dfs, inla_components, covariates_list)

# ============================================================================
# 4. Model Formulas
# ============================================================================

mvpa_binary_formulas <- define_model_formulas(spde_models$spde_mvpa_binary)
mvpa_continuous_formulas <- define_model_formulas(spde_models$spde_mvpa_continuous)

mobility_binary_formulas <- define_model_formulas(spde_models$spde_mobility_binary)
mobility_continuous_formulas <- define_model_formulas(spde_models$spde_mobility_continuous)


mobility_binary_35_formulas <- define_model_formulas(spde_models$spde_mobility_binary_35)
mobility_continuous_35_formulas <- define_model_formulas(spde_models$spde_mobility_continuous_35)
mvpa_binary_35_formulas <- define_model_formulas(spde_models$spde_mvpa_binary_35)
mvpa_continuous_35_formulas <- define_model_formulas(spde_models$spde_mvpa_continuous_35)

# seden_urban_formulas <- define_model_formulas(spde_models$spde_seden_urban)
# energy_urban_formulas <- define_model_formulas(spde_models$spde_energy_urban)
# mvpa_urban_formulas <- define_model_formulas(spde_models$spde_mvpa_urban)
# mobility_urban_formulas <- define_model_formulas(spde_models$spde_mobility_urban)
# sensitivity_urban_formulas <- define_model_formulas(spde_models$spde_urban_35)

# ============================================================================
# 5. Model Fitting
# ============================================================================
#sedentary_models <- fit_sedentary_models(seden_formulas, inla_stacks)
#energy_models <- fit_energy_models(energy_formulas, inla_stacks)

mobility_binary_models <- fit_mobility_binary_models(mobility_binary_formulas, inla_stacks)
mobility_continuous_models <- fit_mobility_continuous_models(mobility_continuous_formulas, inla_stacks)
mvpa_binary_models <- fit_mvpa_binary_models(mvpa_binary_formulas, inla_stacks)
mvpa_continuous_models <- fit_mvpa_continuous_models(mvpa_continuous_formulas, inla_stacks)


## Sensitivity analyses
mobility_binary_35_models <- fit_mobility_35_binary_models(mobility_binary_35_formulas, inla_stacks)
mobility_continous_35_models <- fit_mobility_35_continuous_models(mobility_continuous_35_formulas, inla_stacks)
mvpa_binary_35_models <- fit_mvpa_35_binary_models(mvpa_binary_35_formulas, inla_stacks)
mvpa_continous_35_models <- fit_mvpa_35_continuous_models(mvpa_continuous_35_formulas, inla_stacks)


mvpa_continous_35_models$IM3_mvpa_pt_35$summary.fixed
# sedentary_urban_models <- fit_sedentary_models(seden_urban_formulas, inla_stacks)
# mobility_urban_models <- fit_mobility_models(mobility_urban_formulas, inla_stacks)
# energy_urban_models <- fit_energy_models(energy_urban_formulas, inla_stacks)
# mvpa_urban_models <- fit_mvpa_models(mvpa_urban_formulas, inla_stacks)
# sensitivity_urban_models <- fit_sensitivity_models(sensitivity_urban_formulas, inla_stacks)

pt_all_effect <- mobility_continuous_models$IM3_mobility_nonlinear$summary.random$`inla.group(pt_all, n = 20)`
pt_all_effect

plot(pt_all_effect$mean ~ pt_all_effect$ID, type = "l", 
     xlab = "Proximity Time (standardized)", ylab = "Effect on Outcome")
# Add confidence intervals
lines(pt_all_effect$ID, pt_all_effect$mean + 1.96*pt_all_effect$sd, lty = 2)
lines(pt_all_effect$ID, pt_all_effect$mean - 1.96*pt_all_effect$sd, lty = 2)



mean_pt <- mean(dfs$df_urban$overall_15min_city_proximity_time)
sd_pt <- sd(dfs$df_urban$overall_15min_city_proximity_time)

actual_time_z_minus1 <- mean_pt - 1*sd_pt
actual_time_z_0 <- mean_pt
actual_time_z_0
actual_time_z_1 <- mean_pt + 1*sd_pt
actual_time_z_1
actual_time_z_4 <- mean_pt + 4*sd_pt
actual_time_z_4


# ============================================================================
# 6. Results Processing and Visualization
# ============================================================================
# results_sedentary <- process_results_sedentary(sedentary_models)
# results_sedentary

results_mobility <- process_results_active_mobility(mobility_continuous_models)
results_mobility

results_mobility_binary <- process_results_binary_mobility(mobility_binary_models)
results_mobility_binary


results_mvpa <- process_results_active_mobility(mvpa_continuous_models)
results_mvpa
results_mvpa_binary <- process_results_binary_mobility(mvpa_binary_models)
results_mvpa_binary


results_mvpa_35 <- process_results_active_mobility(mvpa_continous_35_models)


results_mvpa_35 <- process_results_sensitivity(mvpa_continous_35_models)
results_mobility_35 <- process_results_sensitivity(mobility_continous_35_models)

results_mvpa_35
results_mobility_35

# Save results
save_results(results_mvpa,'results_mvpa_continuous.csv')
save_results(results_mvpa_binary,'results_mvpa_binary.csv')

save_results(results_mobility,'results_mobility_continuous.csv')
save_results(results_mobility_binary,'results_mobility_binary.csv')

save_results(results_mvpa_35,'results_mpva_continuous_35.csv')
save_results(results_mobility_35,'results_mobility_continuous_35.csv')



source("./src/utils.R")

create_effect_plots(mvpa_continuous_models, 'Medium & vigorous PA (standardized)', 'pt_all')
create_effect_plots(mobility_continuous_models,'Active mobility time (standardized)', 'pt_all')

figure_folder<- './results/all_pois/figures'


# Create and save the effect plot

mvpa_effect_plot <- create_effect_plots(
  mvpa_continuous_models, 
  "moderate & vigorous PA", 
  "pt_all", 
  model_names = c("Model 1", "Model 2", "Model 3",'Age interaction'),
  save_path = figure_folder
)
mvpa_effect_plot


mobility_effect_plot <- create_effect_plots(
  mobility_continuous_models, 
  "Active mobility", 
  "pt_all", 
  model_names = c("Model 1", "Model 2", "Model 3",'Age interaction'),
  save_path = figure_folder
)
mobility_effect_plot

mobility_binary_effect_plot <- create_effect_plots(
  mobility_binary_models, 
  "Active mobility", 
  "pt_all", 
  model_names = c("Model 1", "Model 2", "Model 3",'Age interaction'),
  save_path = figure_folder
)


 
mobility_continuous_models_all <- mobility_continuous_models[!grepl("urban", names(mobility_continuous_models))]
mobility_binary_models_all <- mobility_binary_models[!grepl("urban", names(mobility_binary_models))]
mvpa_continuous_models_all <- mvpa_continuous_models[!grepl("urban", names(mvpa_continuous_models))]
mvpa_binary_models_all <- mvpa_binary_models[!grepl("urban", names(mvpa_binary_models))]

hist(mobility_continuous_models_all$IM3_mobility_nonlinear$cpo$pit, breaks=20, main="PIT Histogram")
plot(mobility_continuous_models_all$IM3_mobility_nonlinear$cpo$cpo, ylab = "CPO values", main = "CPO Diagnostics")


source("./src/utils.R")

# Collect all model results in a list
all_models <- list(
  mobility = mobility_continuous_models_all,
  mvpa = mvpa_continuous_models_all
)



all_models_binary <- list(
  mobility_binary = mobility_binary_models_all,
  mvpa_binary = mvpa_binary_models_all
)

# Create and save the combined forest plot
combined_plot <- create_forest_plot(
  all_models, 
  parameter = "pt_all",
  save_path = "./results/all_pois/figures"
)

combined_plot

# Display the plot
print(combined_plot)

mean_pt <- mean(dfs$df$overall_15min_city_proximity_time)
sd_pt <- sd(dfs$df$overall_15min_city_proximity_time)


source('src/utils.R')
 # After fitting all models
plot_nonlinear_effects(list(
  energy = energy_models,
  seden = sedentary_models,
  mvpa = mvpa_models,
  mobility = mobility_continuous_models,
  mobility_binary = mobility_binary_models
), pt_mean = mean_pt,
   pt_sd = sd_pt)

plot_nonlinear_effects(list(
  mvpa = mvpa_continuous_models,
  mobility = mobility_continuous_models
), pt_mean = mean_pt,
pt_sd = sd_pt)

mobility_continuous_models_all$IM3_mobility_nonlinear
# Example usage
sweet_spot_viz <- create_sweet_spot_visualization(
  mobility_model = mobility_continuous_models_all$IM3_mobility_nonlinear,
  mvpa_model = mvpa_models_all$IM3_mvpa_nonlinear,
  pt_mean = mean(dfs$df$overall_15min_city_proximity_time, na.rm = TRUE),  # Your actual mean proximity time
  pt_sd = sd(dfs$df$overall_15min_city_proximity_time, na.rm = TRUE),     # Your actual SD of proximity time
  baseline_mobility = mean(dfs$df$Mobility..commute...personal.time...standardized..min.day., na.rm = TRUE), # Mean active mobility time
  output_path = "./results/all_pois/figures/sweet_spot_visualization.png"
)


sd(dfs$df$Mobility..commute...personal.time...standardized..min.day., na.rm = TRUE)



