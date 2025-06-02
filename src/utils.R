# ============================================================================
# 1. Data Loading and Preparation
# ============================================================================

load_data <- function() {
  # Load CSV data
  df <- read.csv("./data/gdf_final_1200.csv")
  df_urban <- read.csv("./data/gdf_final_1200_urban.csv")
  
  # Load spatial data
  ch_shp <- st_read('./data/canton_gecontour.gpkg', quiet = TRUE)
  canton_ge <- st_read('./data/canton_ge.gpkg', quiet = TRUE)
  communes_shp <- st_read("./data/communes_ge.gpkg", quiet = TRUE)
  
  # Create boundary segments for INLA
  bdry <- inla.sp2segment(canton_ge)
  bdry$loc <- inla.mesh.map(bdry$loc)
  
  return(list(
    df = df, 
    df_urban = df_urban, 
    ch_shp = ch_shp, 
    canton_ge = canton_ge, 
    communes_shp = communes_shp,
    bdry = bdry
  ))
}

df$smoking_status

prepare_datasets <- function(df, df_urban) {
  #Factorize age_group
  df$age_group <- as.factor(df$age_group)
  df_urban$age_group <- as.factor(df_urban$age_group)
  
  # Create age-filtered datasets
  df_35 <- df[df$age >= 34, ]
  df_urban_35 <- df_urban[df_urban$age >= 34, ]
  
  # Standardize continuous variables
  continuous_vars <- c("age", 'overall_15min_city_proximity_time', 'overall_15min_city_pa_proximity_time')
  continuous_vars_std <- c("age_std", 'overall_15min_city_proximity_time_std', 'overall_15min_city_pa_proximity_time_std')
  df[, continuous_vars_std] <- scale(df[, continuous_vars])
  df_35[, continuous_vars_std] <- scale(df_35[, continuous_vars])
  df_urban[, continuous_vars_std] <- scale(df_urban[, continuous_vars])
  df_urban_35[, continuous_vars_std] <- scale(df_urban_35[, continuous_vars])
  
  # Convert to spatial objects
  df <- st_as_sf(df, coords = c("E", "N"))
  df_35 <- st_as_sf(df_35, coords = c("E", "N"))
  df_urban <- st_as_sf(df_urban, coords = c("E", "N"))
  df_urban_35 <- st_as_sf(df_urban_35, coords = c("E", "N"))
  
  # Set coordinate reference system
  st_crs(df) <- "EPSG:2056"
  st_crs(df_35) <- "EPSG:2056"
  st_crs(df_urban) <- "EPSG:2056"
  st_crs(df_urban_35) <- "EPSG:2056"
  
  # Create subset of individuals with positive mobility time
  df_mobility <- df[df$'use_active_mobility_binary' ==1, ]
  df_mobility_35 <- df_35[df_35$'use_active_mobility_binary' ==1, ]
  
  df_mvpa <- df[df$'leisure_mvpa_binary' ==1, ]
  df_mvpa_35 <- df_35[df_35$'leisure_mvpa_binary' ==1, ]
  
  df_urban_mobility <- df_urban[df_urban$use_active_mobility_binary ==1, ]
  df_urban_mobility_35 <- df_urban_35[df_urban_35$use_active_mobility_binary ==1, ]
  
  # Extract coordinates
  coords <- st_coordinates(df)
  coords_35 <- st_coordinates(df_35)
  
  coords_mobility <- st_coordinates(df_mobility)
  coords_urban_mobility <- st_coordinates(df_urban_mobility)
  coords_mobility_35 <- st_coordinates(df_mobility_35)
  coords_urban_mobility_35 <- st_coordinates(df_urban_mobility_35)
  
  coords_urban <- st_coordinates(df_urban)
  coords_urban_35 <- st_coordinates(df_urban_35)
  
  coords_mvpa = st_coordinates(df_mvpa)
  coords_mvpa_35 = st_coordinates(df_mvpa_35)
  
  
  return(list(
    df = df,
    df_mobility = df_mobility,
    df_mvpa = df_mvpa,
    df_35 = df_35, 
    df_mobility_35 = df_mobility_35,
    df_mvpa_35 = df_mvpa_35,
    
    df_urban = df_urban,
    df_urban_mobility = df_urban_mobility,
    df_urban_35 = df_urban_35,
    df_urban_mobility_35 = df_urban_mobility_35,
    
    coords = coords,
    coords_35 = coords_35,
    coords_urban = coords_urban,
    coords_urban_35 = coords_urban_35,
    
    coords_mobility = coords_mobility,
    coords_mobility_35 = coords_mobility_35,
    coords_urban_mobility = coords_urban_mobility,
    coords_urban_mobility_35 = coords_urban_mobility_35,
    coords_mvpa = coords_mvpa,
    coords_mvpa_35 = coords_mvpa_35
  ))
}


# ============================================================================
# 2. Mesh Construction and SPDE Models
# ============================================================================


create_meshes <- function(coords, coords_35) {
  # Create non-convex hulls for each dataset
  non_convex_bdry <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
  non_convex_bdry_35 <- inla.nonconvex.hull(coords_35, -0.03, -0.05, resolution = c(100, 100))
  

  # Create primary meshes (fine resolution)
  mesh1 <- inla.mesh.2d(boundary = non_convex_bdry,
                        max.edge = c(400, 1200),
                        offset = c(600, 3000))
  
  mesh1_35 <- inla.mesh.2d(boundary = non_convex_bdry_35,
                           max.edge = c(400, 1200),
                           offset = c(600, 3000))

  # Return all meshes in a list for convenience
  return(list(
    mesh1 = mesh1,
    mesh1_35 = mesh1_35

  ))
}

create_spde_models <- function(mesh, mesh_35) {
  # Create SPDE models with parameter-specific priors

  
  spde_mvpa_continuous <- inla.spde2.pcmatern(mesh, alpha = 2,
                                              constr = TRUE,
                                              prior.range = c(1125, 0.1),
                                              prior.sigma = c(0.9, 0.1))
  
  spde_mvpa_binary <- inla.spde2.pcmatern(mesh, alpha = 2,
                                          constr = TRUE,
                                          prior.range = c(1125, 0.1),
                                          prior.sigma = c(1, 0.5))
  
  spde_mvpa_continuous_35 <- inla.spde2.pcmatern(mesh_35, alpha = 2,
                                                 constr = TRUE,
                                                 prior.range = c(1125, 0.1),
                                                 prior.sigma = c(0.9, 0.1))
  
  spde_mvpa_binary_35 <- inla.spde2.pcmatern(mesh_35, alpha = 2,
                                             constr = TRUE,
                                             prior.range = c(1125, 0.1),
                                             prior.sigma = c(1, 0.5))
  
  spde_mobility_binary <- inla.spde2.pcmatern(mesh, alpha = 2,
                                              constr = TRUE,
                                              prior.range = c(1125, 0.1),
                                              prior.sigma = c(1, 0.5))
  
  spde_mobility_binary_35 <- inla.spde2.pcmatern(mesh_35, alpha = 2,
                                                 constr = TRUE,
                                                 prior.range = c(1125, 0.1),
                                                 prior.sigma = c(1, 0.5))
  
  spde_mobility_continuous <- inla.spde2.pcmatern(mesh, alpha = 2,
                                                  constr = TRUE,
                                                  prior.range = c(1125, 0.1),
                                                  prior.sigma = c(0.7, 0.1))
  
  spde_mobility_continuous_35 <- inla.spde2.pcmatern(mesh_35, alpha = 2,
                                                     constr = TRUE,
                                                     prior.range = c(1125, 0.1),
                                                     prior.sigma = c(0.7, 0.1))
  
  
  # Return all SPDE models
  return(list(
    spde_mvpa_continuous = spde_mvpa_continuous,
    spde_mvpa_binary = spde_mvpa_binary,
    spde_mvpa_continuous_35 = spde_mvpa_continuous_35,
    spde_mvpa_binary_35 = spde_mvpa_binary_35,
    
    spde_mobility_binary = spde_mobility_binary,
    spde_mobility_binary_35 = spde_mobility_binary_35,
    spde_mobility_continuous = spde_mobility_continuous,
    spde_mobility_continuous_35 = spde_mobility_continuous_35
  ))
}

create_indices_and_projections <- function(spde_models, meshes, coords, coords_35, coords_mobility, coords_mobility_35, coords_mvpa, coords_mvpa_35, coop = NULL) {
  # Create index sets for each SPDE model
  iset_mvpa_binary <- inla.spde.make.index(name = "spatial.field", n.spde = spde_models$spde_mvpa_binary$n.spde)
  iset_mvpa_binary_35 <- inla.spde.make.index(name = "spatial.field", n.spde = spde_models$spde_mvpa_binary_35$n.spde)
  iset_mvpa_continuous <- inla.spde.make.index(name = "spatial.field", n.spde = spde_models$spde_mvpa_continuous$n.spde)
  iset_mvpa_continuous_35 <- inla.spde.make.index(name = "spatial.field", n.spde = spde_models$spde_mvpa_continuous_35$n.spde)

  
  iset_mobility_binary <- inla.spde.make.index(name = "spatial.field", n.spde = spde_models$spde_mobility_binary$n.spde)
  iset_mobility_binary_35 <- inla.spde.make.index(name = "spatial.field", n.spde = spde_models$spde_mobility_binary_35$n.spde)
  iset_mobility_continuous <- inla.spde.make.index(name = "spatial.field", n.spde = spde_models$spde_mobility_continuous$n.spde)
  iset_mobility_continuous_35 <- inla.spde.make.index(name = "spatial.field", n.spde = spde_models$spde_mobility_continuous_35$n.spde)
  

  
  # Create projection matrices (A matrices)
  A <- inla.spde.make.A(mesh = meshes$mesh1, loc = as.matrix(coords))
  A_35 <- inla.spde.make.A(mesh = meshes$mesh1_35, loc = as.matrix(coords_35))

  A_mobility <- inla.spde.make.A(mesh = meshes$mesh1, loc = as.matrix(coords_mobility))
  A_mobility_35 <- inla.spde.make.A(mesh = meshes$mesh1_35, loc = as.matrix(coords_mobility_35))
  
  A_mvpa <- inla.spde.make.A(mesh = meshes$mesh1, loc = as.matrix(coords_mvpa))
  A_mvpa_35 <- inla.spde.make.A(mesh = meshes$mesh1_35, loc = as.matrix(coords_mvpa_35))
  
  
  # Create projection matrix for prediction locations if provided
  Ap <- NULL
  if (!is.null(coop)) {
    Ap <- inla.spde.make.A(mesh = meshes$mesh1, loc = coop)
  }
  
  # Return indices and projections
  return(list(
    indices = list(
      iset_mvpa_continuous = iset_mvpa_continuous,
      iset_mvpa_continuous_35 = iset_mvpa_continuous_35,
      iset_mvpa_binary = iset_mvpa_binary,
      iset_mvpa_binary_35 = iset_mvpa_binary_35,
      iset_mobility_binary = iset_mobility_binary,
      iset_mobility_binary_35 = iset_mobility_binary_35,
      iset_mobility_continuous = iset_mobility_continuous,
      iset_mobility_continuous_35 = iset_mobility_continuous_35
    ),
    projections = list(
      A = A,
      A_35 = A_35,
      Ap = Ap,
      A_mobility = A_mobility, 
      A_mobility_35 = A_mobility_35,
      A_mvpa = A_mvpa,
      A_mvpa_35 = A_mvpa_35

    )
  ))
}

# ============================================================================
# 3. INLA Stack Creation
# ============================================================================


create_covariates_df <- function(df) {
  data.frame(
    b0 = 1,
    appt_year = df$appointment_year,
    pt_all = df$overall_15min_city_proximity_time_std,
    pt_pa_pois = df$overall_15min_city_pa_proximity_time_std,
    Age = df$age_std,
    Age_group = df$age_group,
    Sex = df$SEX_F,
    work_status = df$employment_coded,
    educ_lvl = df$education_coded,
    smoking_status = df$smoking_status,
    bmi_cat = df$bmi_category,
    season = df$season
  )
}

prepare_covariates <- function(dfs) {
  # Extract covariates for each dataset
  covariates <- create_covariates_df(dfs$df)
  covariates_35 <- create_covariates_df(dfs$df_35)
  # covariates_urban <- create_covariates_df(dfs$df_urban)
  # covariates_urban_35 <- create_covariates_df(dfs$df_urban_35)
  
  covariates_mobility <- create_covariates_df(dfs$df_mobility)
  covariates_mobility_35 <- create_covariates_df(dfs$df_mobility_35)
  covariates_mvpa <- create_covariates_df(dfs$df_mvpa)
  covariates_mvpa_35 <- create_covariates_df(dfs$df_mvpa_35)
  
  # covariates_urban_mobility <- create_covariates_df(dfs$df_urban_mobility)
  # covariates_urban_mobility_35 <- create_covariates_df(dfs$df_urban_mobility_35)
  # 
  return(list(
    covariates = covariates,
    covariates_35 = covariates_35,
    covariates_mvpa = covariates_mvpa,
    covariates_mvpa_35 = covariates_mvpa_35,
    covariates_mobility = covariates_mobility,
    covariates_mobility_35 = covariates_mobility_35

  ))
}


create_inla_stacks <- function(dfs, inla_components, covariates_list) {
  # Create stacks for full dataset

  stk_mvpa <- create_inla_stack("Leisure.time.MVPA..standardized..min.day.",
                                dfs$df_mvpa, inla_components$projections$A_mvpa,
                                inla_components$indices$iset_mvpa_continuous,
                                covariates_list$covariates_mvpa, "mvpa_continuous")
  stk_mvpa_binary <- create_inla_stack("leisure_mvpa_binary",
                                dfs$df, inla_components$projections$A,
                                inla_components$indices$iset_mvpa_binary,
                                covariates_list$covariates, "mvpa_binary")

  stk_mobility <- create_inla_stack("Mobility..commute...personal.time...standardized..min.day.",
                                    dfs$df_mobility, inla_components$projections$A_mobility,
                                    inla_components$indices$iset_mobility_continuous,
                                    covariates_list$covariates_mobility, "mobility_continuous",
                                    log_transform = FALSE)
  
  stk_mobility_binary <- create_inla_stack("use_active_mobility_binary",
                                           dfs$df, inla_components$projections$A,
                                           inla_components$indices$iset_mobility_binary,
                                           covariates_list$covariates, "mobility_binary",
                                           log_transform = FALSE)
  
  # Create stacks for age ≥ 35 dataset

  stk_mvpa_35 <- create_inla_stack("Leisure.time.MVPA..standardized..min.day.",
                                   dfs$df_mvpa_35, inla_components$projections$A_mvpa_35,
                                   inla_components$indices$iset_mvpa_continuous_35,
                                   covariates_list$covariates_mvpa_35, "mvpa_continuous_35") #Replace by MVPA to go back to previous results
  
  stk_mvpa_binary_35 <- create_inla_stack("Leisure.time.MVPA..standardized..min.day.",
                                   dfs$df_35, inla_components$projections$A_35,
                                   inla_components$indices$iset_mvpa_binary_35,
                                   covariates_list$covariates_35, "mvpa_binary_35")

  stk_mobility_35 <- create_inla_stack("Mobility..commute...personal.time...standardized..min.day.",
                                       dfs$df_mobility_35, inla_components$projections$A_mobility_35 ,
                                       inla_components$indices$iset_mobility_continuous_35,
                                       covariates_list$covariates_mobility_35, "mobility_continuous_35",
                                       log_transform = FALSE)
  
  stk_mobility_binary_35 <- create_inla_stack("use_active_mobility_binary",
                                              dfs$df_35, inla_components$projections$A_35 ,
                                              inla_components$indices$iset_mobility_binary_35,
                                              covariates_list$covariates_35, "mobility_binary_35",
                                              log_transform = FALSE)
  

  
  # Return all stacks
  return(list(
    stk_mvpa = stk_mvpa,
    stk_mvpa_binary = stk_mvpa_binary,
    stk_mvpa_35 = stk_mvpa_35,
    stk_mvpa_binary_35 = stk_mvpa_binary_35,
    
    stk_mobility = stk_mobility,
    stk_mobility_binary = stk_mobility_binary,
    stk_mobility_35 = stk_mobility_35,
    stk_mobility_binary_35 = stk_mobility_binary_35
    

  ))
}

# ============================================================================
# 4. Model Formulas
# ============================================================================

define_model_formulas <- function(spde) {
  # Basic spatial model
  f_m0_pt <- y ~ 0 + b0 + f(spatial.field, model = spde)
  # Models with proximity time
  f_m1_pt <- y ~ 0 + b0 + pt_all + Age + Sex + f(spatial.field, model = spde)
  f_m2_pt <- y ~ 0 + b0 + pt_all + Age + Sex + work_status + educ_lvl + smoking_status + bmi_cat + f(spatial.field, model = spde)
  f_m3_pt <- y ~ 0 + b0 + pt_all + Age + Sex + work_status + educ_lvl + smoking_status + bmi_cat + season + appt_year + f(spatial.field, model = spde)
  # Non-linear proximity time model
  f_m3_pt_nonlinear <- y ~ 0 + b0 + f(inla.group(pt_all, n=20), model="rw2") + Age + Sex + work_status + educ_lvl + smoking_status + bmi_cat + season + appt_year + f(spatial.field, model = spde)
  f_m3_pt_pa_nonlinear <- y ~ 0 + b0 + f(inla.group(pt_pa_pois, n=20), model="rw2") + Age + Sex + work_status + educ_lvl + smoking_status + bmi_cat + season + appt_year + f(spatial.field, model = spde)
  # Modified model with interaction
  # f_m3_pt_age_interaction <- y ~ 0 + b0 + f(inla.group(pt_all, n=20), model="rw2", Age_group) + Sex + work_status + educ_lvl + smoking_status + bmi_cat + season + appt_year + f(spatial.field, model = spde)
  formula_text <- "y ~ 0 + b0 + pt_all * Age_group + Sex + work_status + educ_lvl + smoking_status + bmi_cat + season + appt_year + f(spatial.field, model=spde)"
  f_m3_pt_age_interaction <- as.formula(formula_text)
  
  
  # Return all formulas
  return(list(
    f_m0_pt = f_m0_pt, 
    f_m1_pt = f_m1_pt,
    f_m2_pt = f_m2_pt,
    f_m3_pt = f_m3_pt,
    f_m3_pt_nonlinear = f_m3_pt_nonlinear,
    f_m3_pt_pa_nonlinear = f_m3_pt_pa_nonlinear,
    f_m3_pt_age_interaction = f_m3_pt_age_interaction
  ))
}


# ============================================================================
# 5. Model Fitting
# ============================================================================

run_inla_model_gamma <- function(formula, stack, model_name = NULL) {
  # Extract stack name if model_name is not provided
  if (is.null(model_name)) {
    # Try to get the name from the stack
    stack_name <- deparse(substitute(stack))
    # Extract relevant part of the formula to identify model type
    formula_str <- deparse(formula)[1]
    if (grepl("pois_all", formula_str)) {
      model_type <- "POI density"
    } else if (grepl("pt_all", formula_str)) {
      model_type <- "proximity time"
    } else if (grepl("pois_phys", formula_str)) {
      model_type <- "physical POI"
    } else if (grepl("pois_out", formula_str)) {
      model_type <- "outdoor POI"
    } else if (grepl("pois_trans", formula_str)) {
      model_type <- "transport POI"
    } else if (grepl("nonlinear", formula_str)) {
      model_type <- "nonlinear"
    } else {
      model_type <- "base"
    }
    model_name <- paste0(stack_name, " (", model_type, ")")
  }
  
  # Print start message with timestamp
  cat(sprintf("\n[%s] Starting Gamma model: %s\n", format(Sys.time(), "%H:%M:%S"), model_name))
  
  # Get the number of observations
  n_obs <- nrow(inla.stack.data(stack))
  cat(sprintf("  - Number of observations: %d\n", n_obs))
  
  # Record start time for measuring duration
  start_time <- Sys.time()
  
  # Run the model with error handling
  result <- tryCatch({
    # First attempt - standard gamma model with stability controls
    cat("  Attempting standard gamma model with numerical stability controls...\n")
    inla(formula,
         family = "gamma",  # Gamma distribution for right-skewed data
         control.family = list(link = "log"),  # Log link function
         data = inla.stack.data(stack),
         control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
         control.predictor = list(A = inla.stack.A(stack)),
         control.mode = list(restart = TRUE),
         # Add controls for numerical stability
         control.inla = list(
           int.strategy = "eb",       # Use empirical Bayes
           diagonal = 1e-6,           # Add small value to diagonal for stability
           tolerance = 1e-4,          # Increase tolerance
           h = 0.001                  # Step size in numerical derivatives
         ),
         verbose = FALSE)
  }, error = function(e) {
    # Log the error
    cat(sprintf("  ERROR in first attempt: %s\n", e$message))
    
    # Second attempt - try with lognormal family
    cat("  Trying with lognormal family instead...\n")
    
    tryCatch({
      inla(formula,
           family = "lognormal",      # Try lognormal instead of gamma
           data = inla.stack.data(stack),
           control.predictor = list(A = inla.stack.A(stack)),
           control.compute = list(dic = TRUE, waic = TRUE),
           control.mode = list(restart = TRUE),
           control.inla = list(int.strategy = "eb"),
           verbose = FALSE)
    }, error = function(e2) {
      cat(sprintf("  ERROR in second attempt: %s\n", e2$message))
      
      # Third attempt - simplify the model by removing some terms
      cat("  Trying with simplified formula (removing spatial component)...\n")
      
      # Extract fixed effects part of the formula
      formula_str <- deparse(formula)
      fixed_part <- sub("f\\(spatial\\.field.*$", "", formula_str)
      simplified_formula <- as.formula(paste0(fixed_part, "1)"))
      
      tryCatch({
        inla(simplified_formula,
             family = "gamma",
             control.family = list(link = "log"),
             data = as.data.frame(inla.stack.data(stack)),  # Convert to data frame
             control.compute = list(dic = TRUE),
             verbose = FALSE)
      }, error = function(e3) {
        cat(sprintf("  ERROR in third attempt: %s\n", e3$message))
        cat("  All model fitting attempts failed.\n")
        return(NULL)
      })
    })
  })
  
  # Check if model fitting succeeded
  if (is.null(result)) {
    cat("  Model fitting FAILED for all attempts. Returning NULL.\n")
    return(NULL)
  }
  
  # Calculate and print duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  # Print completion message with model fit statistics
  cat(sprintf("  [%s] Completed model: %s\n", format(Sys.time(), "%H:%M:%S"), model_name))
  cat(sprintf("  - Duration: %.2f minutes\n", as.numeric(duration)))
  
  # Check if DIC and WAIC are available (they might not be in fallback models)
  if (!is.null(result$dic) && !is.null(result$dic$dic)) {
    cat(sprintf("  - DIC: %.2f\n", result$dic$dic))
  } else {
    cat("  - DIC: Not available\n")
  }
  
  if (!is.null(result$waic) && !is.null(result$waic$waic)) {
    cat(sprintf("  - WAIC: %.2f\n", result$waic$waic))
  } else {
    cat("  - WAIC: Not available\n")
  }
  
  # Return model result
  return(result)
}


run_inla_model_binomial <- function(formula, stack, model_name = NULL) {
  # Extract stack name if model_name is not provided
  if (is.null(model_name)) {
    # Try to get the name from the stack
    stack_name <- deparse(substitute(stack))
    # Extract relevant part of the formula to identify model type
    formula_str <- deparse(formula)[1]
    if (grepl("pois_all", formula_str)) {
      model_type <- "POI density"
    } else if (grepl("pt_all", formula_str)) {
      model_type <- "proximity time"
    } else if (grepl("pois_phys", formula_str)) {
      model_type <- "physical POI"
    } else if (grepl("pois_out", formula_str)) {
      model_type <- "outdoor POI"
    } else if (grepl("pois_trans", formula_str)) {
      model_type <- "transport POI"
    } else if (grepl("nonlinear", formula_str)) {
      model_type <- "nonlinear"
    } else {
      model_type <- "base"
    }
    model_name <- paste0(stack_name, " (", model_type, ")")
  }
  
  # Print start message with timestamp
  cat(sprintf("\n[%s] Starting Binomial model: %s\n", format(Sys.time(), "%H:%M:%S"), model_name))
  
  # Get the number of observations
  n_obs <- nrow(inla.stack.data(stack))
  cat(sprintf("  - Number of observations: %d\n", n_obs))
  
  # Record start time for measuring duration
  start_time <- Sys.time()
  
  # Run the model with Gamma distribution and log link
  result <- inla(formula,
                 family = "binomial",  # Binomial
                 control.family = list(link = "logit"),  # Log link function
                 data = inla.stack.data(stack),
                 control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
                 control.predictor = list(A = inla.stack.A(stack)),
                 # Add initial values for hyperparameters to help convergence
                 control.mode = list(restart = TRUE))
  
  
  # Calculate and print duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  # Print completion message with model fit statistics
  cat(sprintf("  [%s] Completed Binomial model: %s\n", format(Sys.time(), "%H:%M:%S"), model_name))
  cat(sprintf("  - Duration: %.2f minutes\n", as.numeric(duration)))
  cat(sprintf("  - DIC: %.2f\n", result$dic$dic))
  cat(sprintf("  - WAIC: %.2f\n", result$waic$waic))
  
  # Return model result
  return(result)
}

run_inla_model <- function(formula, stack, model_name = NULL) {
  # Extract stack name if model_name is not provided
  if (is.null(model_name)) {
    # Try to get the name from the stack
    stack_name <- deparse(substitute(stack))
    # Extract relevant part of the formula to identify model type
    formula_str <- deparse(formula)[1]
    if (grepl("pois_all", formula_str)) {
      model_type <- "POI density"
    } else if (grepl("pt_all", formula_str)) {
      model_type <- "proximity time"
    } else if (grepl("pois_phys", formula_str)) {
      model_type <- "physical POI"
    } else if (grepl("pois_out", formula_str)) {
      model_type <- "outdoor POI"
    } else if (grepl("pois_trans", formula_str)) {
      model_type <- "transport POI"
    } else if (grepl("nonlinear", formula_str)) {
      model_type <- "nonlinear"
    } else {
      model_type <- "base"
    }
    model_name <- paste0(stack_name, " (", model_type, ")")
  }
  
  # Print start message with timestamp
  cat(sprintf("\n[%s] Starting model: %s\n", format(Sys.time(), "%H:%M:%S"), model_name))
  
  # Get the number of observations
  n_obs <- nrow(inla.stack.data(stack))
  cat(sprintf("  - Number of observations: %d\n", n_obs))
  
  # Record start time for measuring duration
  start_time <- Sys.time()
  
  # Run the model
  result <- inla(formula,
                 family = "gaussian",
                 data = inla.stack.data(stack),
                 control.compute = list(cpo = TRUE, dic = TRUE, waic = TRUE),
                 control.predictor = list(A = inla.stack.A(stack)),
                 control.mode = list(restart = TRUE))
  
  # Calculate and print duration
  end_time <- Sys.time()
  duration <- difftime(end_time, start_time, units = "mins")
  
  # Print completion message with model fit statistics
  cat(sprintf("  [%s] Completed model: %s\n", format(Sys.time(), "%H:%M:%S"), model_name))
  cat(sprintf("  - Duration: %.2f minutes\n", as.numeric(duration)))
  cat(sprintf("  - DIC: %.2f\n", result$dic$dic))
  cat(sprintf("  - WAIC: %.2f\n", result$waic$waic))
  
  # Return model result
  return(result)
}

fit_sedentary_models <- function(formulas, stacks) {
  # Initialize empty list for models
  models <- list()
  
  # Define model configurations to try
  model_configs <- list(
    IM0_seden_pt = list(formula = formulas$f_m0_pt, stack = stacks$stk_seden),
    IM1_seden_pt = list(formula = formulas$f_m1_pt, stack = stacks$stk_seden),
    IM2_seden_pt = list(formula = formulas$f_m2_pt, stack = stacks$stk_seden),
    IM3_seden_pt = list(formula = formulas$f_m3_pt, stack = stacks$stk_seden),
    IM3_seden_nonlinear = list(formula = formulas$f_m3_pt_nonlinear, stack = stacks$stk_seden),
    
    IM0_seden_urban_pt = list(formula = formulas$f_m0_urban_pt, stack = stacks$stk_urban_seden),
    IM1_seden_urban_pt = list(formula = formulas$f_m1_urban_pt, stack = stacks$stk_urban_seden),
    IM2_seden_urban_pt = list(formula = formulas$f_m2_urban_pt, stack = stacks$stk_urban_seden),
    IM3_seden_urban_pt = list(formula = formulas$f_m3_urban_pt, stack = stacks$stk_urban_seden),
    IM3_seden_urban_nonlinear = list(formula = formulas$f_m3_urban_pt_nonlinear, stack = stacks$stk_urban_seden)
  )
  
  # Try each model configuration
  for (model_name in names(model_configs)) {
    config <- model_configs[[model_name]]
    cat(sprintf("\nAttempting to fit model: %s\n", model_name))
    
    # Use tryCatch to handle errors
    result <- tryCatch({
      run_inla_model(config$formula, config$stack, model_name = model_name)
    }, error = function(e) {
      cat(sprintf("ERROR fitting model %s: %s\n", model_name, e$message))
      return(NULL)  # Return NULL for failed models
    })
    
    # Only add successful model fits to the results list
    if (!is.null(result)) {
      models[[model_name]] <- result
      # Save model state for reproducibility
      saveRDS(list(model = result, formula = config$formula),
              file = file.path("./results/all_pois/model_states", 
                               paste0(model_name, "_state.rds")))
      cat(sprintf("Successfully fitted model: %s\n", model_name))
    }
  }
  
  
  return(models)
}


fit_mvpa_models <- function(formulas, stacks) {
  # Fit models for MVPA
  models <- list()
  
  # Define a list of model configurations to try
  model_configs <- list(
    IM1_mvpa_pt = list(formula = formulas$f_m1_pt, stack = stacks$stk_mvpa),
    IM2_mvpa_pt = list(formula = formulas$f_m2_pt, stack = stacks$stk_mvpa),
    IM3_mvpa_pt = list(formula = formulas$f_m3_pt, stack = stacks$stk_mvpa),
    IM3_mvpa_nonlinear = list(formula = formulas$f_m3_pt_nonlinear, stack = stacks$stk_mvpa),
    IM3_mvpa_pa_nonlinear = list(formula = formulas$f_m3_pt_pa_nonlinear, stack = stacks$stk_mvpa),
    IM3_mvpa_age_interaction = list(formula = formulas$f_m3_pt_age_interaction, stack = stacks$stk_mvpa),
    
    IM1_mvpa_urban_pt = list(formula = formulas$f_m1_urban_pt, stack = stacks$stk_urban_mvpa),
    IM2_mvpa_urban_pt = list(formula = formulas$f_m2_urban_pt, stack = stacks$stk_urban_mvpa),
    IM3_mvpa_urban_pt = list(formula = formulas$f_m3_urban_pt, stack = stacks$stk_urban_mvpa),
    IM3_mvpa_urban_nonlinear = list(formula = formulas$f_m3_urban_pt_nonlinear, stack = stacks$stk_urban_mvpa),
    IM3_mvpa_urban_pa_nonlinear = list(formula = formulas$f_m3_urban_pt_pa_nonlinear, stack = stacks$stk_urban_mvpa),
    IM3_mvpa_urban_age_interaction = list(formula = formulas$f_m3_urban_pt_age_interaction, stack = stacks$stk_mvpa)
    
  )
  
  # Try each model configuration
  for (model_name in names(model_configs)) {
    config <- model_configs[[model_name]]
    cat(sprintf("\nAttempting to fit model: %s\n", model_name))
    
    # Use tryCatch to handle errors
    result <- tryCatch({
      run_inla_model(config$formula, config$stack, model_name = model_name)
    }, error = function(e) {
      cat(sprintf("ERROR fitting model %s: %s\n", model_name, e$message))
      return(NULL)  # Return NULL for failed models
    })
    
    # Only add successful model fits to the results list
    if (!is.null(result)) {
      models[[model_name]] <- result
      saveRDS(list(model = result, formula = config$formula),
              file = file.path("./results/all_pois/model_states", 
                               paste0(model_name, "_state.rds")))
      cat(sprintf("Successfully fitted model: %s\n", model_name))
    }
  }
  
  # Check if any models were successfully fit
  if (length(models) == 0) {
    warning("No models were successfully fit!")
  } else {
    cat(sprintf("\nSuccessfully fit %d out of %d models\n", 
                length(models), length(model_configs)))
  }
  
  return(models)
}



fit_mvpa_binary_models <- function(formulas, stacks) {
  # Fit models for energy expenditure
  models <- list(
    IM1_mvpa_pt = run_inla_model_binomial(formulas$f_m1_pt, stacks$stk_mvpa_binary),
    IM2_mvpa_pt = run_inla_model_binomial(formulas$f_m2_pt, stacks$stk_mvpa_binary),
    IM3_mvpa_pt = run_inla_model_binomial(formulas$f_m3_pt, stacks$stk_mvpa_binary),
    IM3_mvpa_nonlinear = run_inla_model_binomial(formulas$f_m3_pt_nonlinear, stacks$stk_mvpa_binary),
    IM3_mvpa_pa_nonlinear = run_inla_model_binomial(formulas$f_m3_pt_pa_nonlinear, stacks$stk_mvpa_binary),
    IM3_mvpa_age_interaction = run_inla_model_binomial(formulas$f_m3_pt_age_interaction, stacks$stk_mvpa_binary)
  )
  return(models)
}

fit_mvpa_continuous_models <- function(formulas, stacks) {
  # Fit models for energy expenditure
  models <- list(
    IM1_mvpa_pt = run_inla_model_gamma(formulas$f_m1_pt, stacks$stk_mvpa),
    IM2_mvpa_pt = run_inla_model_gamma(formulas$f_m2_pt, stacks$stk_mvpa),
    IM3_mvpa_pt = run_inla_model_gamma(formulas$f_m3_pt, stacks$stk_mvpa),
    IM3_mvpa_nonlinear = run_inla_model_gamma(formulas$f_m3_pt_nonlinear, stacks$stk_mvpa),
    IM3_mvpa_pa_nonlinear = run_inla_model_gamma(formulas$f_m3_pt_pa_nonlinear, stacks$stk_mvpa),
    IM3_mvpa_age_interaction = run_inla_model_gamma(formulas$f_m3_pt_age_interaction, stacks$stk_mvpa)
  )
  return(models)
}

fit_mobility_binary_models <- function(formulas, stacks) {
  # Fit models for energy expenditure
  models <- list(
    IM1_mobility_pt = run_inla_model_binomial(formulas$f_m1_pt, stacks$stk_mobility_binary),
    IM2_mobility_pt = run_inla_model_binomial(formulas$f_m2_pt, stacks$stk_mobility_binary),
    IM3_mobility_pt = run_inla_model_binomial(formulas$f_m3_pt, stacks$stk_mobility_binary),
    IM3_mobility_nonlinear = run_inla_model_binomial(formulas$f_m3_pt_nonlinear, stacks$stk_mobility_binary),
    IM3_mobility_pa_nonlinear = run_inla_model_binomial(formulas$f_m3_pt_pa_nonlinear, stacks$stk_mobility_binary),
    IM3_mobility_age_interaction = run_inla_model_binomial(formulas$f_m3_pt_age_interaction, stacks$stk_mobility_binary)

  )
  return(models)
}

fit_mobility_continuous_models <- function(formulas, stacks) {
  # Fit models for energy expenditure
  models <- list(
    IM1_mobility_pt = run_inla_model_gamma(formulas$f_m1_pt, stacks$stk_mobility),
    IM2_mobility_pt = run_inla_model_gamma(formulas$f_m2_pt, stacks$stk_mobility),
    IM3_mobility_pt = run_inla_model_gamma(formulas$f_m3_pt, stacks$stk_mobility),
    IM3_mobility_nonlinear = run_inla_model_gamma(formulas$f_m3_pt_nonlinear, stacks$stk_mobility),
    IM3_mobility_pa_nonlinear = run_inla_model_gamma(formulas$f_m3_pt_pa_nonlinear, stacks$stk_mobility),
    IM3_mobility_age_interaction = run_inla_model_gamma(formulas$f_m3_pt_age_interaction, stacks$stk_mobility)

  )
  
  return(models)
}

fit_sensitivity_models <- function(formulas, stacks) {
  models <- list(
    IM3_mvpa_pt_35 = run_inla_model_gamma(formulas$f_m3_pt, stacks$stk_mvpa_35),
    IM3_mvpa_nonlinear_35 = run_inla_model_gamma(formulas$f_m3_pt_nonlinear, stacks$stk_mvpa_35),
    
    IM3_mvpa_binary_pt_35 = run_inla_model_binomial(formulas$f_m3_pt, stacks$stk_mvpa_binary_35),

    IM3_mobility_pt_35 = run_inla_model_gamma(formulas$f_m3_pt, stacks$stk_mobility_35),
    IM3_mobility_binary_pt_35 = run_inla_model_binomial(formulas$f_m3_pt, stacks$stk_mobility_binary_35),
    IM3_mobility_nonlinear_35 = run_inla_model_gamma(formulas$f_m3_pt_nonlinear, stacks$stk_mobility_35),
    IM3_mobility_binary_nonlinear_35 = run_inla_model_binomial(formulas$f_m3_pt_nonlinear, stacks$stk_mobility_binary_35)

  )
  
  return(models)
}

fit_mvpa_35_continuous_models <- function(formulas, stacks) {
  models <- list(
    IM3_mvpa_pt_35 = run_inla_model_gamma(formulas$f_m3_pt, stacks$stk_mvpa_35),
    IM3_mvpa_nonlinear_35 = run_inla_model_gamma(formulas$f_m3_pt_nonlinear, stacks$stk_mvpa_35)
  )
  
  return(models)
}

fit_mvpa_35_binary_models <- function(formulas, stacks) {
  models <- list(
    IM3_mvpa_binary_pt_35 = run_inla_model_binomial(formulas$f_m3_pt, stacks$stk_mvpa_binary_35)
  )
  
  return(models)
}

fit_mobility_35_continuous_models <- function(formulas, stacks) {
  models <- list(
    IM3_mobility_pt_35 = run_inla_model_gamma(formulas$f_m3_pt, stacks$stk_mobility_35),
    IM3_mobility_nonlinear_35 = run_inla_model_gamma(formulas$f_m3_pt_nonlinear, stacks$stk_mobility_35)
    
  )
  
  return(models)
}

fit_mobility_35_binary_models <- function(formulas, stacks) {
  models <- list(
    IM3_mobility_binary_pt_35 = run_inla_model_binomial(formulas$f_m3_pt, stacks$stk_mobility_binary_35),
    IM3_mobility_binary_nonlinear_35 = run_inla_model_binomial(formulas$f_m3_pt_nonlinear, stacks$stk_mobility_binary_35)
  )
  
  return(models)
}


# ============================================================================
# 6. Results Processing and Visualization
# ============================================================================

process_results_sedentary <- function(models) {
  # Helper functions for safely extracting values
  safe_extract_fixed <- function(model, coef_name) {
    if (is.null(model) || is.null(model$summary.fixed) || 
        !(coef_name %in% rownames(model$summary.fixed))) {
      return(rep(NA, 4))  # Return NAs if coefficient is missing
    }
    return(model$summary.fixed[coef_name, c(1,2,3,5)])
  }
  
  safe_extract_dic <- function(model) {
    if (is.null(model) || is.null(model$dic) || is.null(model$dic$dic)) {
      return(NA)  # Return NA if DIC is missing
    }
    return(model$dic$dic)
  }
  
  safe_mean_log_cpo <- function(model) {
    if (is.null(model) || is.null(model$cpo) || !is.numeric(model$cpo$cpo) || length(model$cpo$cpo) == 0) {
      return(NA)  # Return NA if CPO is missing or invalid
    }
    valid_cpo <- model$cpo$cpo[model$cpo$cpo > 0]
    if (length(valid_cpo) == 0) {
      return(NA)  # Return NA if no valid CPO values
    }
    return(mean(log(valid_cpo)))
  }
  
  # Initialize results list
  results_list_seden <- list()
  
  # Only try to extract results for models that exist and have valid data
  if (!is.null(models$IM1_seden_pt)) {
    results_list_seden$Model1_seden_PT <- c(
      safe_extract_fixed(models$IM1_seden_pt, 'pt_all'),
      safe_extract_dic(models$IM1_seden_pt),
      safe_mean_log_cpo(models$IM1_seden_pt)
    )
  }
  
  if (!is.null(models$IM2_seden_pt)) {
    results_list_seden$Model2_seden_PT <- c(
      safe_extract_fixed(models$IM2_seden_pt, 'pt_all'),
      safe_extract_dic(models$IM2_seden_pt),
      safe_mean_log_cpo(models$IM2_seden_pt)
    )
  }
  
  if (!is.null(models$IM3_seden_pt)) {
    results_list_seden$Model3_seden_PT <- c(
      safe_extract_fixed(models$IM3_seden_pt, 'pt_all'),
      safe_extract_dic(models$IM3_seden_pt),
      safe_mean_log_cpo(models$IM3_seden_pt)
    )
  }
  
  if (!is.null(models$IM3_seden_nonlinear)) {
    # Note: Nonlinear models need special handling - they're included here for completeness
    # but will need further processing to extract meaningful coefficients
    cat("Note: Nonlinear model detected but not processed in summary table\n")
  }
  
  # Check if any models were processed
  if (length(results_list_seden) == 0) {
    warning("No valid sedentary time models to process!")
    return(NULL)
  }
  
  # Combine results and create data frame
  combined_results_seden <- do.call(rbind, results_list_seden)
  
  combined_results_seden_clean <- data.frame(
    Model = rownames(combined_results_seden),
    Estimate = round(as.numeric(combined_results_seden[,1]), 3),
    SD = round(as.numeric(combined_results_seden[,2]), 3),
    '95_CI' = sprintf("(%.3f - %.3f)", 
                      as.numeric(combined_results_seden[,3]), 
                      as.numeric(combined_results_seden[,4])),
    DIC = round(as.numeric(combined_results_seden[,5]), 3),
    Mean_log_CPO = round(as.numeric(combined_results_seden[,6]), 3)
  )
  
  # Add informative column names
  colnames(combined_results_seden_clean) <- c('Model', "Estimate", "SD", "95% CI", "DIC", "Mean log CPO")
  
  return(combined_results_seden_clean)
}

process_results_energy <- function(models) {
  # Extract results for energy expenditure models
  results_list_energy <- list(
    #Model1_energy_POI = c(models$IM1_energy$summary.fixed[c('pois_all'),c(1,2,3,5)], 
    #                      models$IM1_energy$dic$dic, 
    #                      mean(log(models$IM1_energy$cpo$cpo[models$IM1_energy$cpo$cpo > 0]))),
    #Model2_energy_POI = c(models$IM2_energy$summary.fixed[c('pois_all'),c(1,2,3,5)], 
    #                      models$IM2_energy$dic$dic, 
    #                      mean(log(models$IM2_energy$cpo$cpo[models$IM2_energy$cpo$cpo > 0]))),
    #Model3_energy_POI = c(models$IM3_energy$summary.fixed[c('pois_all'),c(1,2,3,5)], 
    #                      models$IM3_energy$dic$dic, 
    #                      mean(log(models$IM3_energy$cpo$cpo[models$IM3_energy$cpo$cpo > 0]))),
    Model1_energy_PT = c(models$IM1_energy_pt$summary.fixed[c('pt_all'),c(1,2,3,5)], 
                         models$IM1_energy_pt$dic$dic, 
                         mean(log(models$IM1_energy_pt$cpo$cpo[models$IM1_energy_pt$cpo$cpo > 0]))),
    Model2_energy_PT = c(models$IM2_energy_pt$summary.fixed[c('pt_all'),c(1,2,3,5)], 
                         models$IM2_energy_pt$dic$dic, 
                         mean(log(models$IM2_energy_pt$cpo$cpo[models$IM2_energy_pt$cpo$cpo > 0]))),
    Model3_energy_PT = c(models$IM3_energy_pt$summary.fixed[c('pt_all'),c(1,2,3,5)], 
                         models$IM3_energy_pt$dic$dic, 
                         mean(log(models$IM3_energy_pt$cpo$cpo[models$IM3_energy_pt$cpo$cpo > 0])))
    #Model3_energy_POI_phys = c(models$IM3_energy_phys$summary.fixed[c('pois_phys'),c(1,2,3,5)], 
    #                           models$IM3_energy_phys$dic$dic, 
    #                           mean(log(models$IM3_energy_phys$cpo$cpo[models$IM3_energy_phys$cpo$cpo > 0]))),
    #Model3_energy_POI_out = c(models$IM3_energy_out$summary.fixed[c('pois_out'),c(1,2,3,5)], 
    #                          models$IM3_energy_out$dic$dic, 
    #                          mean(log(models$IM3_energy_out$cpo$cpo[models$IM3_energy_out$cpo$cpo > 0]))),
    #Model3_energy_POI_trans = c(models$IM3_energy_trans$summary.fixed[c('pois_trans'),c(1,2,3,5)], 
    #                            models$IM3_energy_trans$dic$dic, 
    #                            mean(log(models$IM3_energy_trans$cpo$cpo[models$IM3_energy_trans$cpo$cpo > 0])))
  )
  
  combined_results_energy <- do.call(rbind, results_list_energy)
  
  combined_results_energy_clean <- data.frame(
    Model = rownames(combined_results_energy),
    Estimate = round(as.numeric(combined_results_energy[,1]), 3),
    Percent_Change = round(100 * (exp(as.numeric(combined_results_energy[,1])) - 1), 1),
    SD = round(as.numeric(combined_results_energy[,2]), 3),
    '95_CI' = sprintf("(%.1f%% - %.1f%%)", 
                      100 * (exp(as.numeric(combined_results_energy[,3])) - 1),
                      100 * (exp(as.numeric(combined_results_energy[,4])) - 1)),
    DIC = round(as.numeric(combined_results_energy[,5]), 3),
    Mean_log_CPO = round(as.numeric(combined_results_energy[,6]), 3)
  )
  
  # Add informative column names
  colnames(combined_results_energy_clean) <- c('Model', "Estimate", 'Percent change', 
                                               "SD", "95% CI", "DIC", "Mean log CPO")
  
  return(combined_results_energy_clean)
}

process_results_mvpa <- function(models) {
  # Create helper functions for safely extracting results
  safe_extract_fixed <- function(model, coef_name) {
    if (is.null(model) || is.null(model$summary.fixed) || 
        !(coef_name %in% rownames(model$summary.fixed))) {
      return(rep(NA, 4))  # Return NAs if coefficient is missing
    }
    return(model$summary.fixed[coef_name, c(1,2,3,5)])
  }
  
  safe_extract_dic <- function(model) {
    if (is.null(model) || is.null(model$dic) || is.null(model$dic$dic)) {
      return(NA)  # Return NA if DIC is missing
    }
    return(model$dic$dic)
  }
  
  safe_mean_log_cpo <- function(model) {
    if (is.null(model) || is.null(model$cpo) || !is.numeric(model$cpo$cpo) || length(model$cpo$cpo) == 0) {
      return(NA)  # Return NA if CPO is missing or invalid
    }
    valid_cpo <- model$cpo$cpo[model$cpo$cpo > 0]
    if (length(valid_cpo) == 0) {
      return(NA)  # Return NA if no valid CPO values
    }
    return(mean(log(valid_cpo)))
  }
  
  # Initialize results list
  results_list_mvpa <- list()
  
  # Only try to extract results for models that exist and have valid data
  if (!is.null(models$IM1_mvpa_pt)) {
    results_list_mvpa$Model1_mvpa_PT <- c(
      safe_extract_fixed(models$IM1_mvpa_pt, 'pt_all'),
      safe_extract_dic(models$IM1_mvpa_pt),
      safe_mean_log_cpo(models$IM1_mvpa_pt)
    )
  }
  
  if (!is.null(models$IM2_mvpa_pt)) {
    results_list_mvpa$Model2_mvpa_PT <- c(
      safe_extract_fixed(models$IM2_mvpa_pt, 'pt_all'),
      safe_extract_dic(models$IM2_mvpa_pt),
      safe_mean_log_cpo(models$IM2_mvpa_pt)
    )
  }
  
  if (!is.null(models$IM3_mvpa_pt)) {
    results_list_mvpa$Model3_mvpa_PT <- c(
      safe_extract_fixed(models$IM3_mvpa_pt, 'pt_all'),
      safe_extract_dic(models$IM3_mvpa_pt),
      safe_mean_log_cpo(models$IM3_mvpa_pt)
    )
  }
  
  # Add other models as needed using the same pattern
  
  # Check if any models were processed
  if (length(results_list_mvpa) == 0) {
    warning("No valid MVPA models to process!")
    return(NULL)
  }
  
  # Combine results and create data frame
  combined_results_mvpa <- do.call(rbind, results_list_mvpa)
  
  combined_results_mvpa_clean <- data.frame(
    Model = rownames(combined_results_mvpa),
    Estimate = round(as.numeric(combined_results_mvpa[,1]), 3),
    SD = round(as.numeric(combined_results_mvpa[,2]), 3),
    '95_CI' = sprintf("(%.3f - %.3f)", 
                      as.numeric(combined_results_mvpa[,3]), 
                      as.numeric(combined_results_mvpa[,4])),
    DIC = round(as.numeric(combined_results_mvpa[,5]), 3),
    Mean_log_CPO = round(as.numeric(combined_results_mvpa[,6]), 3)
  )
  
  # Add informative column names
  colnames(combined_results_mvpa_clean) <- c('Model', "Estimate", "SD", "95% CI", "DIC", "Mean log CPO")
  
  return(combined_results_mvpa_clean)
}

process_results_binary_mobility <- function(models) {
  # Extract results for binary mobility models with binomial distribution
  # For binomial models with logit link, the effects are on the log-odds scale
  results_list <- list()
  
  # Define helper functions for safely extracting values
  safe_extract_fixed <- function(model, coef_name) {
    if (is.null(model) || is.null(model$summary.fixed) || 
        !(coef_name %in% rownames(model$summary.fixed))) {
      return(rep(NA, 4))  # Return NAs if coefficient is missing
    }
    return(model$summary.fixed[coef_name, c(1,2,3,5)])
  }
  
  safe_extract_dic <- function(model) {
    if (is.null(model) || is.null(model$dic) || is.null(model$dic$dic)) {
      return(NA)  # Return NA if DIC is missing
    }
    return(model$dic$dic)
  }
  
  safe_mean_log_cpo <- function(model) {
    tryCatch({
      if (is.null(model) || is.null(model$cpo) || !is.numeric(model$cpo$cpo) || length(model$cpo$cpo) == 0) {
        return(NA)  # Return NA if CPO is missing or invalid
      }
      valid_cpo <- model$cpo$cpo[model$cpo$cpo > 0]
      if (length(valid_cpo) == 0) {
        return(NA)  # Return NA if no valid CPO values
      }
      return(mean(log(valid_cpo)))
    }, error = function(e) {
      cat(sprintf("Error in CPO calculation: %s\n", e$message))
      return(NA)
    })
  }
  
  # Extract information for each model
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    # Skip NULL models
    if (is.null(model)) {
      cat(sprintf("Skipping NULL model: %s\n", model_name))
      next
    }
    
    # Check if the model is a binomial model
    if (is.null(model$family) || model$family != "binomial" || 
        is.null(model$link) || model$link != "logit") {
      cat(sprintf("Warning: Model %s is not a binomial logit-link model\n", model_name))
      # You could skip or adjust processing here
    }
    
    # Determine which coefficient to extract based on model name (update patterns)
    if (grepl("_POI$|_poi$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pois_all'
    } else if (grepl("_PT$|_pt$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pt_all'
    } else if (grepl("_POI_phys$|_poi_phys$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pois_phys'
    } else if (grepl("_POI_out$|_poi_out$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pois_out'
    } else if (grepl("_POI_trans$|_poi_trans$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pois_trans'
    } else if (grepl("_nonlinear$", model_name, ignore.case = TRUE)) {
      cat(sprintf("Skipping non-linear model: %s (requires special handling)\n", model_name))
      next  # Skip non-linear models as they require special handling
    } else {
      cat(sprintf("Unknown model type: %s\n", model_name))
      next  # Skip if we can't determine coefficient
    }
    
    # Extract model results safely
    tryCatch({
      # For binomial models with logit link, the exponentiated effects represent odds ratios
      results_list[[model_name]] <- c(
        # Raw coefficient on log-odds scale (with safety checks)
        safe_extract_fixed(model, coef_name),
        # DIC and CPO (with safety checks)
        safe_extract_dic(model),
        safe_mean_log_cpo(model)
      )
      
      cat(sprintf("Successfully processed model: %s\n", model_name))
    }, error = function(e) {
      cat(sprintf("Error processing model %s: %s\n", model_name, e$message))
    })
  }
  
  # Combine results
  if (length(results_list) > 0) {
    combined_results <- do.call(rbind, results_list)
    
    # Create clean dataframe with results
    combined_results_clean <- data.frame(
      Model = rownames(combined_results),
      Estimate_logit = round(as.numeric(combined_results[,1]), 3),
      Odds_Ratio = round(exp(as.numeric(combined_results[,1])), 3),
      SD = round(as.numeric(combined_results[,2]), 3),
      '95_CI' = sprintf("(%.3f - %.3f)", 
                        exp(as.numeric(combined_results[,3])),
                        exp(as.numeric(combined_results[,4]))),
      DIC = round(as.numeric(combined_results[,5]), 3),
      Mean_log_CPO = round(as.numeric(combined_results[,6]), 3)
    )
    
    # Add informative column names
    colnames(combined_results_clean) <- c('Model', "Log-Odds Estimate", 'Odds Ratio', 
                                          "SD", "95% CI (OR)", "DIC", "Mean log CPO")
    
    cat(sprintf("Successfully processed %d models\n", nrow(combined_results_clean)))
    return(combined_results_clean)
  } else {
    warning("No valid binary mobility models to process!")
    return(NULL)
  }
}

process_results_active_mobility <- function(models) {
  # Extract results for active mobility models with gamma distribution
  # For gamma models with log link, the effects are multiplicative
  results_list <- list()
  
  # Define helper functions for safely extracting values
  safe_extract_fixed <- function(model, coef_name) {
    if (is.null(model) || is.null(model$summary.fixed) || 
        !(coef_name %in% rownames(model$summary.fixed))) {
      return(rep(NA, 4))  # Return NAs if coefficient is missing
    }
    return(model$summary.fixed[coef_name, c(1,2,3,5)])
  }
  
  safe_extract_dic <- function(model) {
    if (is.null(model) || is.null(model$dic) || is.null(model$dic$dic)) {
      return(NA)  # Return NA if DIC is missing
    }
    return(model$dic$dic)
  }
  
  safe_mean_log_cpo <- function(model) {
    tryCatch({
      if (is.null(model) || is.null(model$cpo) || !is.numeric(model$cpo$cpo) || length(model$cpo$cpo) == 0) {
        return(NA)  # Return NA if CPO is missing or invalid
      }
      valid_cpo <- model$cpo$cpo[model$cpo$cpo > 0]
      if (length(valid_cpo) == 0) {
        return(NA)  # Return NA if no valid CPO values
      }
      return(mean(log(valid_cpo)))
    }, error = function(e) {
      cat(sprintf("Error in CPO calculation: %s\n", e$message))
      return(NA)
    })
  }
  
  # Extract information for each model
  for (model_name in names(models)) {
    model <- models[[model_name]]
    if (is.null(model$family) || model$family != "gamma" || 
        is.null(model$link) || model$link != "log") {
      cat(sprintf("Warning: Model %s is not a gamma log-link model\n", model_name))
      # You could skip or adjust processing here
    }
    
    # Skip NULL models
    if (is.null(model)) {
      cat(sprintf("Skipping NULL model: %s\n", model_name))
      next
    }
    
    # Determine which coefficient to extract based on model name (update patterns)
    if (grepl("_POI$|_poi$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pois_all'
    } else if (grepl("_PT$|_pt$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pt_all'
    } else if (grepl("_POI_phys$|_poi_phys$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pois_phys'
    } else if (grepl("_POI_out$|_poi_out$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pois_out'
    } else if (grepl("_POI_trans$|_poi_trans$", model_name, ignore.case = TRUE)) {
      coef_name <- 'pois_trans'
    } else if (grepl("_nonlinear$", model_name, ignore.case = TRUE)) {
      cat(sprintf("Skipping non-linear model: %s (requires special handling)\n", model_name))
      next  # Skip non-linear models as they require special handling
    } else {
      cat(sprintf("Unknown model type: %s\n", model_name))
      next  # Skip if we can't determine coefficient
    }
    
    # Extract model results safely
    tryCatch({
      # For gamma models with log link, the exponentiated effects represent multiplicative effects
      results_list[[model_name]] <- c(
        # Raw coefficient on log scale (with safety checks)
        safe_extract_fixed(model, coef_name),
        # DIC and CPO (with safety checks)
        safe_extract_dic(model),
        safe_mean_log_cpo(model)
      )
      
      cat(sprintf("Successfully processed model: %s\n", model_name))
    }, error = function(e) {
      cat(sprintf("Error processing model %s: %s\n", model_name, e$message))
    })
  }
  
  # Combine results
  if (length(results_list) > 0) {
    combined_results <- do.call(rbind, results_list)
    
    # Create clean dataframe with results
    combined_results_clean <- data.frame(
      Model = rownames(combined_results),
      Estimate_log = round(as.numeric(combined_results[,1]), 3),
      Percent_Change = round(100 * (exp(as.numeric(combined_results[,1])) - 1), 1),
      SD = round(as.numeric(combined_results[,2]), 3),
      '95_CI' = sprintf("(%.1f%% - %.1f%%)", 
                        100 * (exp(as.numeric(combined_results[,3])) - 1),
                        100 * (exp(as.numeric(combined_results[,4])) - 1)),
      DIC = round(as.numeric(combined_results[,5]), 3),
      Mean_log_CPO = round(as.numeric(combined_results[,6]), 3)
    )
    
    # Add informative column names
    colnames(combined_results_clean) <- c('Model', "Log Estimate", 'Percent change', 
                                          "SD", "95% CI", "DIC", "Mean log CPO")
    
    cat(sprintf("Successfully processed %d models\n", nrow(combined_results_clean)))
    return(combined_results_clean)
  } else {
    warning("No valid active mobility models to process!")
    return(NULL)
  }
}

process_results_sensitivity <- function(models) {
  # Helper functions for safely extracting values
  safe_extract_fixed <- function(model, coef_name) {
    if (is.null(model) || is.null(model$summary.fixed) || 
        !(coef_name %in% rownames(model$summary.fixed))) {
      return(rep(NA, 4))  # Return NAs if coefficient is missing
    }
    return(model$summary.fixed[coef_name, c(1,2,3,5)])
  }
  
  safe_extract_dic <- function(model) {
    if (is.null(model) || is.null(model$dic) || is.null(model$dic$dic)) {
      return(NA)  # Return NA if DIC is missing
    }
    return(model$dic$dic)
  }
  
  safe_mean_log_cpo <- function(model) {
    tryCatch({
      if (is.null(model) || is.null(model$cpo) || !is.numeric(model$cpo$cpo) || length(model$cpo$cpo) == 0) {
        return(NA)  # Return NA if CPO is missing or invalid
      }
      valid_cpo <- model$cpo$cpo[model$cpo$cpo > 0]
      if (length(valid_cpo) == 0) {
        return(NA)  # Return NA if no valid CPO values
      }
      return(mean(log(valid_cpo)))
    }, error = function(e) {
      cat(sprintf("Error in CPO calculation: %s\n", e$message))
      return(NA)
    })
  }
  
  # Initialize results list
  results_list_sensitivity <- list()
  
  # Model configurations to try to process
  model_configs <- list(
    Model1_energy_PT_35 = list(model = "IM1_energy_pt_35", coef = "pt_all"),
    Model2_energy_PT_35 = list(model = "IM2_energy_pt_35", coef = "pt_all"),
    Model3_energy_PT_35 = list(model = "IM3_energy_pt_35", coef = "pt_all"),
    
    Model1_seden_PT_35 = list(model = "IM1_seden_pt_35", coef = "pt_all"),
    Model2_seden_PT_35 = list(model = "IM2_seden_pt_35", coef = "pt_all"),
    Model3_seden_PT_35 = list(model = "IM3_seden_pt_35", coef = "pt_all"),
    
    Model1_mvpa_PT_35 = list(model = "IM1_mvpa_pt_35", coef = "pt_all"),
    Model2_mvpa_PT_35 = list(model = "IM2_mvpa_pt_35", coef = "pt_all"),
    Model3_mvpa_PT_35 = list(model = "IM3_mvpa_pt_35", coef = "pt_all"),
    
    Model1_mobility_PT_35 = list(model = "IM1_mobility_pt_35", coef = "pt_all"),
    Model2_mobility_PT_35 = list(model = "IM2_mobility_pt_35", coef = "pt_all"),
    Model3_mobility_PT_35 = list(model = "IM3_mobility_pt_35", coef = "pt_all")
  )
  
  # Process each model configuration
  for (result_name in names(model_configs)) {
    config <- model_configs[[result_name]]
    model_name <- config$model
    coef_name <- config$coef
    
    # Check if model exists
    if (!model_name %in% names(models) || is.null(models[[model_name]])) {
      cat(sprintf("Model %s not found or is NULL, skipping\n", model_name))
      next
    }
    
    model <- models[[model_name]]
    
    # Try to extract results
    tryCatch({
      results_list_sensitivity[[result_name]] <- c(
        safe_extract_fixed(model, coef_name),
        safe_extract_dic(model),
        safe_mean_log_cpo(model)
      )
      cat(sprintf("Successfully processed %s\n", result_name))
    }, error = function(e) {
      cat(sprintf("Error processing %s: %s\n", result_name, e$message))
    })
  }
  
  # Check if we have any results to process
  if (length(results_list_sensitivity) == 0) {
    warning("No sensitivity models were successfully processed!")
    return(NULL)
  }
  
  # Combine results and create data frame
  combined_results <- do.call(rbind, results_list_sensitivity)
  
  # Handle potential NA values in percentage calculations
  calculate_percent_change <- function(value) {
    if (is.na(value)) return(NA)
    return(100 * (exp(value) - 1))
  }
  
  # Create clean data frame with results
  combined_results_clean <- data.frame(
    Model = rownames(combined_results),
    Estimate = round(as.numeric(combined_results[,1]), 3),
    Percent_Change = round(sapply(as.numeric(combined_results[,1]), calculate_percent_change), 1),
    SD = round(as.numeric(combined_results[,2]), 3),
    '95_CI' = sprintf("(%.1f%% - %.1f%%)", 
                      sapply(as.numeric(combined_results[,3]), calculate_percent_change),
                      sapply(as.numeric(combined_results[,4]), calculate_percent_change)),
    DIC = round(as.numeric(combined_results[,5]), 3),
    Mean_log_CPO = round(as.numeric(combined_results[,6]), 3)
  )
  
  # Add informative column names
  colnames(combined_results_clean) <- c('Model', "Estimate", 'Pct change', 
                                        "SD", "95% CI", "DIC", "Mean log CPO")
  
  cat(sprintf("Successfully processed %d sensitivity models\n", nrow(combined_results_clean)))
  return(combined_results_clean)
}


save_results <- function(results, filename) {
  # Create results directory if it doesn't exist
  result_folder <- './results/all_pois'
  if (!dir.exists(result_folder)) {
    dir.create(result_folder, recursive = TRUE)
  }
  
  # Save results to CSV
  write.csv(results, file.path(result_folder, filename), row.names = FALSE)
}

visualize_spatial_field <- function(model, mesh, canton_ge, resolution = c(200, 200)) {
  # Force garbage collection before starting
  gc()
  
  # Create projector with lower resolution
  rang <- apply(mesh$loc[, c(1, 2)], 2, range)
  proj <- inla.mesh.projector(mesh,
                              xlim = rang[, 1], 
                              ylim = rang[, 2],
                              dims = resolution)  # Reduced from 500x500
  
  # Project spatial field values
  mean_s <- inla.mesh.project(proj, model$summary.random$spatial.field$mean)
  sd_s <- inla.mesh.project(proj, model$summary.random$spatial.field$sd)
  
  # Create data frames more efficiently
  coords <- expand.grid(x = proj$x, y = proj$y)
  plot_data_mean <- cbind(coords, z = as.vector(mean_s))
  plot_data_sd <- cbind(coords, z = as.vector(sd_s))
  
  # Remove large objects no longer needed
  rm(mean_s, sd_s, coords)
  gc()
  
  # Create mean field plot
  mean_plot <- ggplot(plot_data_mean, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c(name = "Mean") +
    theme_minimal() +
    labs(title = "Posterior Mean of Spatial Random Field",
         x = "Longitude", y = "Latitude") +
    coord_equal()
  
  # Remove large object no longer needed
  rm(plot_data_mean)
  gc()
  
  # Create standard deviation plot
  sd_plot <- ggplot(plot_data_sd, aes(x = x, y = y, fill = z)) +
    geom_raster() +
    scale_fill_viridis_c(name = "SD", option = "magma") +
    theme_minimal() +
    labs(title = "Posterior Standard Deviation of Spatial Random Field",
         x = "Longitude", y = "Latitude") +
    coord_equal()
  
  # Remove large object no longer needed
  rm(plot_data_sd)
  gc()
  
  # Return plots combined as a patchwork object
  return(mean_plot + sd_plot)
}

plot_nonlinear_effects <- function(models, save_path = "./results/all_pois/figures/", 
                                   pt_mean = NULL, pt_sd = NULL, data = NULL) {
  # Create directory if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  # Define outcomes and colors
  outcomes <- list(
    energy = list(name = "Energy Expenditure", model = models$energy$IM3_energy_nonlinear, 
                  color = "#E64B35", ylab = "Effect on log(Energy Expenditure)"),
    seden = list(name = "Sedentary Time", model = models$seden$IM3_seden_nonlinear, 
                 color = "#4DBBD5", ylab = "Effect on Sedentary Time (min/day)"),
    mvpa = list(name = "MVPA", model = models$mvpa$IM3_mvpa_nonlinear, 
                color = "#feb24c", ylab = "% Change"),
    mobility = list(name = "Active Mobility", model = models$mobility$IM3_mobility_nonlinear, 
                    color = "#2c7fb8", ylab = "% Change"),
    mobility_binary = list(name = "Active Mobility", model = models$mobility_binary$IM3_mobility_nonlinear, 
                           color = "#2c7fb8", ylab = "Effect on Active Mobility"),
    energy_urban = list(name = "Energy Expenditure", model = models$energy$IM3_energy_urban_nonlinear, 
                        color = "#E64B35", ylab = "Effect on log(Energy Expenditure)"),
    seden_urban = list(name = "Sedentary Time", model = models$seden$IM3_seden_urban_nonlinear, 
                       color = "#4DBBD5", ylab = "Effect on Sedentary Time (min/day)"),
    mvpa_urban = list(name = "MVPA", model = models$mvpa$IM3_mvpa_urban_nonlinear, 
                      color = "#feb24c", ylab = "Effect on MVPA (min/day)"),
    mobility_urban = list(name = "Active Mobility", model = models$mobility$IM3_mobility_urban_nonlinear, 
                          color = "#2c7fb8", ylab = "Effect on Active Mobility"),
    mobility_binary_urban = list(name = "Active Mobility", model = models$mobility_binary$IM3_mobility_urban_nonlinear, 
                                 color = "#2c7fb8", ylab = "Effect on Active Mobility")
  )
  
  # For removing log scale, we need baseline/mean value of active mobility
  baseline_mobility <- if (!is.null(data)) mean(data$mobility, na.rm = TRUE) else 100
  
  # Plot for each outcome
  for (outcome_key in names(outcomes)) {
    outcome <- outcomes[[outcome_key]]
    if (is.null(outcome$model)) {
      cat(sprintf("No nonlinear model available for %s\n", outcome$name))
      next
    }
    
    # Extract effect estimates
    effect_data <- outcome$model$summary.random$`inla.group(pt_all, n = 20)`
    if (is.null(effect_data)) {
      cat(sprintf("No nonlinear effect data found for %s\n", outcome$name))
      next
    }
    
    # Create plot - using PNG for publication quality
    png(file.path(save_path, sprintf("nonlinear_effect_%s.png", outcome_key)), 
        width = 800, height = 600, res = 120)
    
    par(mar = c(5, 5, 2, 2), family = "Helvetica")
    
    # If mean and sd are provided, create transformed x-axis
    if (!is.null(pt_mean) && !is.null(pt_sd)) {
      # Transform standardized values back to original scale
      original_x <- effect_data$ID * pt_sd + pt_mean
      
      # Check if we need to transform from log-scale (for mobility outcomes)
      is_gamma_outcome <- grepl("mobility|mvpa", outcome_key, ignore.case = TRUE)
      
      if (is_gamma_outcome) {
        # Transform log-scale effects to actual scale
        # exp(log(baseline) + effect) - baseline
        effect_mean <- baseline_mobility * (exp(effect_data$mean) - 1)
        effect_lower <- baseline_mobility * (exp(effect_data$mean - 1.96*effect_data$sd) - 1)
        effect_upper <- baseline_mobility * (exp(effect_data$mean + 1.96*effect_data$sd) - 1)
        
        ylab <- "% Change"
      } else {
        effect_mean <- effect_data$mean
        effect_lower <- effect_data$mean - 1.96*effect_data$sd
        effect_upper <- effect_data$mean + 1.96*effect_data$sd
        
        ylab <- outcome$ylab
      }
      
      # Get data density for rug plot if data is available
      has_data <- !is.null(data) && "pt_all" %in% names(data)
      
      # Create the main plot
      plot(effect_mean ~ original_x, type = "n", 
           xlab = "Proximity Time (min)",
           ylab = ylab,
           cex.lab = 1.2, cex.axis = 1.1,
           las = 1, # Horizontal axis labels
           xaxs = "i", yaxs = "i") # Tight axes
      
      # Add grid
      grid(lty = "dotted", col = "gray90")
      
      # Add confidence intervals
      polygon(c(original_x, rev(original_x)), 
              c(effect_upper, rev(effect_lower)),
              col = adjustcolor(outcome$color, alpha.f = 0.2), border = NA)
      
      # Add line for mean effect
      lines(original_x, effect_mean, col = outcome$color, lwd = 2)
      
      # Add rug plot for data density if available
      if (has_data) {
        rug(data$pt_all, side = 1, col = adjustcolor("black", alpha.f = 0.3))
      }
      
      # Add reference line at y = 0
      abline(h = 0, lty = 2, col = "gray50")
      
      # Find where effect crosses zero (if it does)
      zero_cross <- NULL
      for (i in 1:(length(original_x) - 1)) {
        if ((effect_mean[i] >= 0 && effect_mean[i+1] <= 0) || 
            (effect_mean[i] <= 0 && effect_mean[i+1] >= 0)) {
          # Simple linear interpolation to find crossing point
          prop <- abs(effect_mean[i]) / (abs(effect_mean[i]) + abs(effect_mean[i+1]))
          zero_cross <- original_x[i] + prop * (original_x[i+1] - original_x[i])
          break
        }
      }
      
      # Add vertical line at zero crossing
      if (!is.null(zero_cross)) {
        abline(v = zero_cross, lty = 3, col = "gray50")
        text(zero_cross, min(effect_lower) * 0.9, 
             sprintf("%.1f", zero_cross), pos = 4, cex = 0.8)
      }
      
      # Add annotations for positive/negative regions
      if (any(effect_mean > 0) && any(effect_mean < 0)) {
        text(min(original_x) + (max(original_x) - min(original_x))*0.05, 
             max(effect_upper) * 0.8, 
             "Positive effect", cex = 0.9)
        text(max(original_x) - (max(original_x) - min(original_x))*0.05, 
             min(effect_lower) * 0.8, 
             "Negative effect", cex = 0.9)
      }
      
      
    } else {
      # Original plot without transformation - similar logic but with standardized values
      # (Code omitted for brevity, but would follow same pattern as above)
      plot(effect_data$mean ~ effect_data$ID, type = "l", 
           col = outcome$color, lwd = 2,
           xlab = "Proximity Time (min)", 
           ylab = outcome$ylab,
           cex.lab = 1.2, cex.axis = 1.1)
      
      # Add confidence intervals
      polygon(c(effect_data$ID, rev(effect_data$ID)), 
              c(effect_data$mean + 1.96*effect_data$sd, 
                rev(effect_data$mean - 1.96*effect_data$sd)),
              col = adjustcolor(outcome$color, alpha.f = 0.2), border = NA)
      
      # Add reference line at y = 0
      abline(h = 0, lty = 2, col = "gray50")
      
    }
    
    
    
    dev.off()
    cat(sprintf("Saved nonlinear effect plot for %s\n", outcome$name))
    
    # Also save data for reproducibility
    result_data <- effect_data
    if (!is.null(pt_mean) && !is.null(pt_sd)) {
      result_data$original_x <- original_x
      
      if (grepl("mobility", outcome_key, fixed = TRUE)) {
        result_data$non_log_effect <- baseline_mobility * (exp(effect_data$mean) - 1)
        result_data$non_log_lower <- baseline_mobility * (exp(effect_data$mean - 1.96*effect_data$sd) - 1)
        result_data$non_log_upper <- baseline_mobility * (exp(effect_data$mean + 1.96*effect_data$sd) - 1)
      }
    }
    
    write.csv(result_data, 
              file.path(save_path, sprintf("nonlinear_effect_%s_data.csv", outcome_key)),
              row.names = FALSE)
  }
}

plot_predicted_values <- function(model, df, canton_ge, variable_name) {
  # Add predicted values to the dataframe
  df[[paste0(variable_name, "_pred")]] <- model$summary.fitted.values[, "mean"]
  
  # Create plot
  plot <- ggplot() + 
    geom_sf(data = canton_ge) +
    geom_sf(data = df, aes(col = .data[[paste0(variable_name, "_pred")]])) +
    scale_color_viridis() +
    theme_minimal() +
    labs(title = paste0("Predicted ", variable_name),
         color = "Predicted Value")
  
  return(plot)
}

perform_model_diagnostics <- function(model, df, coords, outcome_var) {
  # Extract fitted values and residuals
  fitted_values <- model$summary.fitted.values$mean[1:nrow(df)]
  residuals <- df[[outcome_var]] - fitted_values
  standardized_residuals <- residuals / sqrt(model$summary.fitted.values$sd[1:nrow(df)]^2)
  
  # Create a new plotting device to avoid parameter conflicts
  dev.new(width = 10, height = 8, noRStudioGD = TRUE)
  
  # Set plot parameters safely
  par(mfrow = c(2, 2), mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  # Residuals vs Fitted
  plot(fitted_values, residuals,
       xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs Fitted")
  abline(h = 0, col = "red", lty = 2)
  
  # Q-Q plot
  qqnorm(standardized_residuals)
  qqline(standardized_residuals, col = "red")
  
  # Scale-Location plot
  plot(fitted_values, sqrt(abs(standardized_residuals)),
       xlab = "Fitted values", ylab = "√|Standardized residuals|",
       main = "Scale-Location")
  
  # Histogram of residuals for the fourth panel
  hist(standardized_residuals, main = "Histogram of Residuals", 
       xlab = "Standardized Residuals", breaks = 30)
  
  # Title for the overall plot
  mtext(paste("Diagnostics for", outcome_var), outer = TRUE, line = 0)
  
  # Return diagnostics summary without resetting parameters
  return(list(
    residuals = residuals,
    standardized_residuals = standardized_residuals,
    fitted_values = fitted_values,
    mean_residual = mean(residuals),
    sd_residual = sd(residuals)
  ))
}

diagnose_inla_gamma <- function(model, stack, df, outcome_var, plot_title = NULL) {
  
  # Extract fitted and linear predictor values correctly
  idx <- inla.stack.index(stack, "est")$data
  fitted_values <- model$summary.fitted.values$mean[idx]
  linear_pred <- model$summary.linear.predictor$mean[idx]
  
  # Precision parameter extraction
  phi_name <- grep("Gamma observations|gamma observations|precision",
                   rownames(model$summary.hyperpar), value = TRUE, ignore.case = TRUE)[1]
  phi <- model$summary.hyperpar[phi_name, "mean"]
  if (is.na(phi)) stop("Precision parameter extraction failed.")
  
  # Observed values
  observed <- df[[outcome_var]]
  
  # Deviance residuals for Gamma-log model
  deviance_residuals <- sign(observed - fitted_values) *
    sqrt(2 * (log(ifelse(observed == 0, 1, observed / fitted_values)) - 
                (observed - fitted_values) / fitted_values))
  
  # Pearson residuals
  pearson_residuals <- (observed - fitted_values) / sqrt(fitted_values^2 / phi)
  
  # Plots
  par(mfrow = c(2,2), mar=c(4,4,2,1))
  
  plot(linear_pred, deviance_residuals, pch=20,
       xlab="Linear predictor", ylab="Deviance residuals")
  abline(h=0, col="red")
  lines(lowess(linear_pred, deviance_residuals), col="blue")
  
  hist(deviance_residuals, breaks=30, main="Histogram of deviance residuals",
       xlab="Deviance residuals")
  
  plot(linear_pred, sqrt(abs(pearson_residuals)), pch=20,
       xlab="Linear predictor", ylab="√|Pearson residuals|")
  lines(lowess(linear_pred, sqrt(abs(pearson_residuals))), col="blue")
  
  # Mean-Variance check
  means <- tapply(fitted_values, cut(fitted_values, breaks=20), mean)
  vars <- tapply((observed - fitted_values)^2, cut(fitted_values, breaks=20), mean)
  plot(means, vars, pch=20, xlab="Mean fitted value", ylab="Variance")
  curve(x^2/phi, add=TRUE, col="red", lwd=2)
  
  if (!is.null(plot_title)) {
    mtext(plot_title, outer=TRUE, line=-2, cex=1.2)
  }
  
  # Return values
  return(list(
    deviance_residuals = deviance_residuals,
    pearson_residuals = pearson_residuals,
    fitted_values = fitted_values,
    linear_predictor = linear_pred,
    dispersion_estimate = 1/phi
  ))
}

# ============================================================================
# 7. Main Analysis Pipeline
# ============================================================================

run_analysis <- function() {
  # 1. Load data
  data_list <- load_data()
  
  # 2. Prepare datasets
  dfs <- prepare_datasets(data_list$df, data_list$df_urban)
  
  # 3. Create meshes
  meshes <- create_meshes(dfs$coords, dfs$coords_35, dfs$coords_urban, dfs$coords_urban_35)
  
  # 4. Create SPDE models
  spde_models <- create_spde_models(meshes$mesh1, meshes$mesh1_35, 
                                    meshes$mesh1_urban, meshes$mesh1_urban_35)
  
  # 5. Create indices and projections
  inla_components <- create_indices_and_projections(
    spde_models, meshes, 
    dfs$coords, dfs$coords_35, 
    dfs$coords_urban, dfs$coords_urban_35
  )
  
  # 6. Prepare covariates
  covariates_list <- prepare_covariates(dfs)
  
  # 7. Create INLA stacks
  stacks <- create_inla_stacks(dfs, inla_components, covariates_list)
  
  # 8. Define model formulas
  formulas <- define_model_formulas()
  
  # 9. Fit models
  seden_models <- fit_sedentary_models(formulas, stacks)
  mvpa_models <- fit_mvpa_models(formulas, stacks)
  energy_models <- fit_energy_models(formulas, stacks)
  sensitivity_models <- fit_sensitivity_models(formulas, stacks)
  
  # 10. Process results
  seden_results <- process_results_sedentary(seden_models)
  mvpa_results <- process_results_mvpa(mvpa_models)
  energy_results <- process_results_energy(energy_models)
  sensitivity_results <- process_results_sensitivity(sensitivity_models)
  
  # 11. Save results
  save_results(seden_results, "model_comparison_sedentary.csv")
  save_results(mvpa_results, "model_comparison_mvpa.csv")
  save_results(energy_results, "model_comparison_energy.csv")
  save_results(sensitivity_results, "model_comparison_energy_35.csv")
  
  # 12. Visualize results for selected models
  spatial_plots <- list(
    seden_spatial = visualize_spatial_field(seden_models$IM3_seden, meshes$mesh1, data_list$canton_ge),
    mvpa_spatial = visualize_spatial_field(mvpa_models$IM3_mvpa, meshes$mesh1, data_list$canton_ge),
    energy_spatial = visualize_spatial_field(energy_models$IM3_energy, meshes$mesh1, data_list$canton_ge)
  )
  
  prediction_plots <- list(
    seden_pred = plot_predicted_values(seden_models$IM3_seden, dfs$df, data_list$canton_ge, "Sedentary"),
    mvpa_pred = plot_predicted_values(mvpa_models$IM3_mvpa, dfs$df, data_list$canton_ge, "MVPA"),
    energy_pred = plot_predicted_values(energy_models$IM3_energy, dfs$df, data_list$canton_ge, "Energy")
  )
  
  # 13. Return all results and plots
  return(list(
    models = list(
      seden = seden_models,
      mvpa = mvpa_models,
      energy = energy_models,
      sensitivity = sensitivity_models
    ),
    results = list(
      seden = seden_results,
      mvpa = mvpa_results,
      energy = energy_results,
      sensitivity = sensitivity_results
    ),
    plots = list(
      spatial = spatial_plots,
      prediction = prediction_plots
    ),
    data = dfs,
    spatial_data = list(
      canton_ge = data_list$canton_ge,
      communes_shp = data_list$communes_shp
    ),
    meshes = meshes
  ))
}

# ============================================================================
# 8. Additional Utility Functions
# ============================================================================

compare_models <- function(model_list, model_names = NULL) {
  # Compare DIC values
  dic_values <- sapply(model_list, function(m) m$dic$dic)
  if (!is.null(model_names)) {
    names(dic_values) <- model_names
  }
  
  # Compare WAIC values
  waic_values <- sapply(model_list, function(m) m$waic$waic)
  if (!is.null(model_names)) {
    names(waic_values) <- model_names
  }
  
  # Compare CPO values
  cpo_values <- sapply(model_list, function(m) {
    mean(log(m$cpo$cpo[m$cpo$cpo > 0]))
  })
  if (!is.null(model_names)) {
    names(cpo_values) <- model_names
  }
  
  # Return comparison table
  comparison <- data.frame(
    Model = if (!is.null(model_names)) model_names else names(model_list),
    DIC = dic_values,
    WAIC = waic_values,
    Mean_log_CPO = cpo_values
  )
  
  return(comparison)
}

extract_fixed_effects <- function(model, terms = NULL) {
  # Extract all fixed effects
  fixed_effects <- model$summary.fixed
  
  # Filter specific terms if requested
  if (!is.null(terms)) {
    fixed_effects <- fixed_effects[rownames(fixed_effects) %in% terms, ]
  }
  
  return(fixed_effects)
}

extract_spatial_range <- function(model) {
  # Extract posterior for spatial range parameter
  if ("spatial.field" %in% names(model$summary.hyperpar)) {
    # For matern models
    range_index <- grep("Range for spatial.field", rownames(model$summary.hyperpar))
    if (length(range_index) > 0) {
      return(model$summary.hyperpar[range_index, ])
    }
  }
  return(NULL)
}

create_effect_plots <- function(model_list, outcome, parameter, model_names = NULL, 
                                save_path = NULL, width = 8, height = 6, dpi = 300) {
  # Filter out nonlinear models
  linear_models <- model_list[!grepl("nonlinear", names(model_list), ignore.case = TRUE)]
  
  if (length(linear_models) == 0) {
    warning("No linear models available after filtering out nonlinear models.")
    return(NULL)
  }
  
  # Extract parameter estimates from each model
  estimates <- sapply(linear_models, function(m) m$summary.fixed[parameter, "mean"])
  lower_ci <- sapply(linear_models, function(m) m$summary.fixed[parameter, "0.025quant"])
  upper_ci <- sapply(linear_models, function(m) m$summary.fixed[parameter, "0.975quant"])
  
  # Create dataframe for plotting
  if (is.null(model_names)) {
    model_names <- names(linear_models)
  } else if (length(model_names) > length(linear_models)) {
    # Adjust model_names if it includes nonlinear models
    model_names <- model_names[!grepl("nonlinear", names(model_list), ignore.case = TRUE)]
  }
  
  plot_data <- data.frame(
    Model = factor(model_names, levels = model_names),
    Estimate = estimates,
    Lower = lower_ci,
    Upper = upper_ci
  )
  
  # Add color coding based on CI crossing zero
  plot_data$color_code <- ifelse(plot_data$Lower <= 0 & plot_data$Upper >= 0, "grey",  # CI crosses 0
                                 ifelse(plot_data$Lower > 0, "red",                    # CI entirely above 0
                                        "blue"))                                        # CI entirely below 0
  
  # Create plot with color coding
  effect_plot <- ggplot(plot_data, aes(y = Model, x = Estimate, xmin = Lower, xmax = Upper, color = color_code)) +
    geom_point(size = 3) +
    geom_errorbar(width = 0.2) +
    scale_color_identity() +  # Use the color values directly
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = paste("Effect of", variable_labels[parameter],'on', outcome),
         x = "", y = "Estimate with 95% CI") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black")
  
  # Save the plot if a path is provided
  if (!is.null(save_path)) {
    file_name <- file.path(save_path, 
                           paste0("effect_", gsub(" ", "_", parameter), "_on_", 
                                  gsub(" ", "_", outcome), ".png"))
    
    ggsave(filename = file_name, plot = effect_plot, 
           width = width, height = height, dpi = dpi)
    
    cat(sprintf("Plot saved to: %s\n", file_name))
  }
  
  return(effect_plot)
}

create_forest_plot <- function(all_models, parameter = "pt_all", save_path = NULL) {
  library(ggplot2)
  library(grid)  # Required for gtable operations
  library(gtable)  # Required for manipulating the plot layout
  
  # Create directory if needed
  if (!is.null(save_path) && !dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  # Define outcome configurations
  outcome_configs <- list(
    energy = list(
      models = all_models$energy,
      title = "Total Energy Expenditure",
      exponentiate = TRUE,
      color = "darkblue", # Red
      transform_label = "% Change per SD increase"
    ),
    mobility = list(
      models = all_models$mobility,
      title = "A) Non-occupational AM",
      exponentiate = TRUE,
      color = "#2c7fb8", # Blue
      transform_label = "% Change per SD increase"
    ),
    mobility_binary = list(
      models = all_models$mobility_binary,
      title = "Active Mobility (Binomial)",
      exponentiate = TRUE,
      show_as_percent = FALSE,  # New flag
      color = "red", # Blue
      transform_label = "Odds Ratio"
    ),
    
    mvpa = list(
      models = all_models$mvpa,
      title = "B) Leisure-time MVPA",
      exponentiate = TRUE,
      color = "#feb24c", # Green
      transform_label = "% Change per SD increase"
    ),
    sedentary = list(
      models = all_models$sedentary,
      title = "Home Sedentary Time",
      exponentiate = FALSE,
      color = "darkred", # Purple
      transform_label = "Minutes per SD increase"
    )
  )
  
  
  
  
  # Initialize data frame
  all_data <- data.frame()
  
  # Create model descriptions - improved for clarity
  model_descriptions <- c(
    "1" = "Model 1: Unadjusted",
    "2" = "Model 2: Adjusted for demographics",
    "3" = "Model 3: Fully adjusted"
  )
  
  # Process each outcome
  for (outcome_name in names(outcome_configs)) {
    config <- outcome_configs[[outcome_name]]
    models <- config$models
    
    # Filter linear models
    linear_models <- models[!grepl("nonlinear", names(models), ignore.case = TRUE) & 
                              grepl("pt", names(models), ignore.case = TRUE)]
    if (length(linear_models) == 0) next
    
    # Process each model
    for (model_idx in 1:length(linear_models)) {
      model_name <- names(linear_models)[model_idx]
      model <- linear_models[[model_name]]
      
      if (is.null(model) || !parameter %in% rownames(model$summary.fixed)) next
      
      # Extract model number
      model_num <- gsub(".*([0-9])_.*", "\\1", model_name)
      
      # Extract estimates
      est <- model$summary.fixed[parameter, "mean"]
      lower <- model$summary.fixed[parameter, "0.025quant"]
      upper <- model$summary.fixed[parameter, "0.975quant"]
      
      # Transform if needed
      if (config$exponentiate) {
        if (!is.null(config$show_as_percent) && !config$show_as_percent) {
          # Just exponentiate for odds ratios
          est <- exp(est)
          lower <- exp(lower)
          upper <- exp(upper)
        } else {
          # Convert to percentage change
          est <- 100 * (exp(est) - 1)
          lower <- 100 * (exp(lower) - 1)
          upper <- 100 * (exp(upper) - 1)
        }
      }
      # Add to data frame
      all_data <- rbind(all_data, data.frame(
        Outcome = config$title,
        OutcomeOrder = match(outcome_name, names(outcome_configs)),
        Model = ifelse(model_num %in% names(model_descriptions),
                       model_descriptions[model_num],
                       paste("Model", model_num)),
        ModelNum = as.numeric(model_num),
        Estimate = est,
        Lower = lower,
        Upper = upper,
        Color = config$color,
        TransformLabel = config$transform_label,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  if (nrow(all_data) == 0) {
    warning("No data available for plotting")
    return(NULL)
  }
  
  # Organize data
  all_data <- all_data[order(all_data$OutcomeOrder, -all_data$ModelNum), ]
  
  # Calculate plot dimensions
  plot_height <- 2 + (length(unique(all_data$Outcome)) * 2)
  
  # Calculate appropriate x-axis limits
  padding_factor <- 0.2
  x_min <- min(all_data$Lower) * (1 + sign(min(all_data$Lower)) * padding_factor)
  x_max <- max(all_data$Upper) * (1 + sign(max(all_data$Upper)) * padding_factor)
  
  # Add model number as standalone column for y-axis labels
  all_data$ModelLabel <- paste("Model", all_data$ModelNum)
  
  # Create plot - using only ModelLabel for y-axis
  p <- ggplot(all_data, aes(x = Estimate, y = reorder(ModelLabel, desc(ModelNum)))) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper, color = Outcome), height = 0.3, size = 0.7) +
    geom_point(aes(color = Outcome), size = 3, shape = 18) +
    facet_wrap(~ Outcome, ncol = 1, scales = "free_y", strip.position = "top") +
    scale_color_manual(values = unique(all_data$Color)) +
    scale_y_discrete(labels = function(x) gsub(".*\\.", "", x)) +
    scale_x_continuous(
      limits = c(-8.3, 8.3),
      breaks = seq(from = -8, to = 8, by = 1)
    ) +
    labs(
      x = "Effect Estimate",
      y = NULL,
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.placement = "outside",
      strip.text = element_text(
        face = "bold", 
        size = 12,
        hjust=0
      ),
      strip.background = element_rect(fill = "white", color = NA),
      panel.spacing.y = unit(1.5, "lines"),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      plot.caption = element_text(hjust = 0, size = 10, face = "italic"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_line(color = "gray90"),
      axis.title.x = element_text(size = 12, margin = margin(t = 10)),
      axis.text = element_text(size = 11),
      plot.margin = margin(0.5, 1.5, 0.5, 0.5, "cm"),
      panel.border = element_rect(color = "black", fill = NA),
      panel.spacing = unit(0.5, "lines")
    )
  
  # Create a clean version without transform labels to extract layout information
  p_clean <- p
  
  # Get the gtable representation of the plot
  gt <- ggplot_gtable(ggplot_build(p_clean))
  
  # Add transform labels uniquely (one per outcome facet)
  unique_outcomes <- unique(all_data$Outcome)
  
  # Create a new grob list for the transform labels
  transform_grobs <- list()
  
  # Loop through outcomes and add transform labels
  # for (i in 1:length(unique_outcomes)) {
  #   outcome <- unique_outcomes[i]
  #   data_subset <- all_data[all_data$Outcome == outcome, ]
  #   transform_label <- unique(data_subset$TransformLabel)
  #   print(transform_label)
  #   
  #   # Find the panel for this outcome
  #   panel_row <- which(unique_outcomes == outcome)
  #   
  #   # Add transform label as a separate transparent layer
  #   p <- p + annotate(
  #     "text",
  #     x = 8,                 # Position at right side (adjust as needed)
  #     y = nrow(data_subset) / 2,  # Center vertically in facet
  #     label = transform_label,
  #     color = unique(data_subset$Color),
  #     hjust = 1,               # Right-align
  #     vjust = 0.5,             # Center vertically
  #     size = 3,                # Appropriate size
  #     fontface = "italic",     # Italicize for distinction
  #     alpha = 0.85             # Slightly transparent
  #   )
  # }
  
  # Save plot
  if (!is.null(save_path)) {
    filename <- file.path(save_path, "forest_plot_proximity_time.png")
    ggsave(filename, p, width = 10, height = plot_height, dpi = 300, units = "in")
    
    filename_pdf <- file.path(save_path, "forest_plot_proximity_time.pdf")
    ggsave(filename_pdf, p, width = 10, height = plot_height, units = "in")
  }
  
  return(p)
}




###################

create_inla_stack <- function(response_var, df, A_matrix, iset, covariates_df, tag_name, log_transform = FALSE) {
  # Transform response variable if log_transform is TRUE
  y_value <- if(log_transform) {
    log(df[[response_var]])
  } else {
    df[[response_var]]
  }
  
  inla.stack(
    data = list(y = y_value),
    A = list(1, A_matrix),
    effects = list(
      covariates_df,
      s = iset
    ),
    tag = tag_name
  )
}

create_sweet_spot_visualization <- function(mobility_model, mvpa_model, 
                                            pt_mean = NULL, pt_sd = NULL,
                                            baseline_mobility = NULL,
                                            output_path = NULL, 
                                            width = 10, height = 8) {
  library(ggplot2)
  library(dplyr)
  
  # Extract non-linear effects
  mobility_effects <- mobility_model$summary.random$`inla.group(pt_all, n = 20)`
  mvpa_effects <- mvpa_model$summary.random$`inla.group(pt_all, n = 20)`
  
  # Transform standardized values to original scale
  mobility_effects$original_x <- mobility_effects$ID * pt_sd + pt_mean
  mvpa_effects$original_x <- mvpa_effects$ID * pt_sd + pt_mean
  
  # Process mobility effects - convert to minutes per day
  mobility_df <- data.frame(
    x = mobility_effects$original_x,
    effect = baseline_mobility * (exp(mobility_effects$mean) - 1),
    lower = baseline_mobility * (exp(mobility_effects$mean - 1.96*mobility_effects$sd) - 1),
    upper = baseline_mobility * (exp(mobility_effects$mean + 1.96*mobility_effects$sd) - 1),
    type = "Active Mobility"
  )
  
  # Process MVPA effects
  mvpa_df <- data.frame(
    x = mvpa_effects$original_x,
    effect = mvpa_effects$mean,
    lower = mvpa_effects$mean - 1.96 * mvpa_effects$sd, 
    upper = mvpa_effects$mean + 1.96 * mvpa_effects$sd,
    type = "MVPA"
  )
  
  # Sort data frames by x
  mobility_df <- mobility_df %>% arrange(x)
  mvpa_df <- mvpa_df %>% arrange(x)
  
  # Create common x-grid with 100 points for smooth interpolation
  x_min <- min(c(mobility_df$x, mvpa_df$x), na.rm=TRUE)
  x_max <- max(c(mobility_df$x, mvpa_df$x), na.rm=TRUE)
  x_grid <- seq(x_min, x_max, length.out = 100)
  
  # Interpolate mobility values
  mobility_interp <- approx(mobility_df$x, mobility_df$effect, xout = x_grid)
  mobility_lower_interp <- approx(mobility_df$x, mobility_df$lower, xout = x_grid)
  mobility_upper_interp <- approx(mobility_df$x, mobility_df$upper, xout = x_grid)
  
  # Interpolate MVPA values
  mvpa_interp <- approx(mvpa_df$x, mvpa_df$effect, xout = x_grid)
  mvpa_lower_interp <- approx(mvpa_df$x, mvpa_df$lower, xout = x_grid)
  mvpa_upper_interp <- approx(mvpa_df$x, mvpa_df$upper, xout = x_grid)
  
  # Create smooth total effect dataframe
  total_df <- data.frame(
    x = x_grid,
    effect = mobility_interp$y + mvpa_interp$y,
    lower = mobility_lower_interp$y + mvpa_lower_interp$y,
    upper = mobility_upper_interp$y + mvpa_upper_interp$y,
    type = "Total Physical Activity"
  )
  
  # Find zero crossings using the interpolated data
  find_zero_crossing <- function(x, y) {
    for (i in 1:(length(y)-1)) {
      if (!is.na(y[i]) && !is.na(y[i+1])) {
        if ((y[i] >= 0 && y[i+1] <= 0) || (y[i] <= 0 && y[i+1] >= 0)) {
          # Linear interpolation
          prop <- abs(y[i]) / (abs(y[i]) + abs(y[i+1]))
          return(x[i] + prop * (x[i+1] - x[i]))
        }
      }
    }
    return(NULL)
  }
  
  mobility_zero_cross <- find_zero_crossing(mobility_interp$x, mobility_interp$y)
  mvpa_zero_cross <- find_zero_crossing(mvpa_interp$x, mvpa_interp$y)
  
  # Find sweet spot (max of total effect)
  if(any(!is.na(total_df$effect))) {
    sweet_spot_idx <- which.max(total_df$effect)
    sweet_spot_x <- total_df$x[sweet_spot_idx]
    max_total <- total_df$effect[sweet_spot_idx]
  } else {
    sweet_spot_x <- NULL
    max_total <- NULL
  }
  
  # Create plot for mobility
  p1 <- ggplot(mobility_df, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#2c7fb8", alpha = 0.2) +
    geom_line(aes(y = effect), color = "#2c7fb8", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Effect of Proximity Time on Active Mobility",
         x = "Proximity Time (min)",
         y = "Effect on Active Mobility (min/day)") +
    theme_light() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  # Add zero crossing line for mobility if found
  if(!is.null(mobility_zero_cross)) {
    p1 <- p1 + 
      geom_vline(xintercept = mobility_zero_cross, linetype = "dotted", color = "#2c7fb8") +
      annotate("text", x = mobility_zero_cross, y = min(mobility_df$lower, na.rm = TRUE),
               label = sprintf("%.1f min", mobility_zero_cross),
               color = "#2c7fb8", hjust = -0.1, size = 3.5)
  }
  
  # Create plot for MVPA
  p2 <- ggplot(mvpa_df, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#feb24c", alpha = 0.2) +
    geom_line(aes(y = effect), color = "#feb24c", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Effect of Proximity Time on MVPA",
         x = "Proximity Time (min)",
         y = "Effect on MVPA (min/day)") +
    theme_light() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  # Add zero crossing line for MVPA if found
  if(!is.null(mvpa_zero_cross)) {
    p2 <- p2 + 
      geom_vline(xintercept = mvpa_zero_cross, linetype = "dotted", color = "#feb24c") +
      annotate("text", x = mvpa_zero_cross, y = min(mvpa_df$lower, na.rm = TRUE),
               label = sprintf("%.1f min", mvpa_zero_cross),
               color = "#feb24c", hjust = -0.1, size = 3.5)
  }
  
  # Create plot for total effect - using interpolated data
  p3 <- ggplot(total_df, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), fill = "darkgreen", alpha = 0.15) +
    geom_line(aes(y = effect), color = "darkgreen", linewidth = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    labs(title = "Total Effect on Physical Activity (Active Mobility + MVPA)",
         x = "Proximity Time (min)",
         y = "Combined Effect (min/day)") +
    theme_light() +
    theme(
      plot.title = element_text(hjust = 0.5),
      axis.title = element_text(face = "bold")
    )
  
  # Add sweet spot if found
  if(!is.null(sweet_spot_x) && !is.null(max_total)) {
    p3 <- p3 + 
      geom_vline(xintercept = sweet_spot_x, linetype = "dotted", color = "darkgreen") +
      annotate("text", x = sweet_spot_x, y = max_total,
               label = sprintf("Sweet Spot: %.1f min\n(+%.1f min/day)", sweet_spot_x, max_total),
               color = "darkgreen", hjust = -0.1, size = 3.5)
  }
  
  # Combine plots if patchwork is available
  if(requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    combined_plot <- p1 / p2 / p3 +
      plot_annotation(
        title = "Effects of Proximity Time on Physical Activity (min/day)",
        subtitle = "Identifying optimal proximity time for combined physical activity",
        theme = theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12)
        )
      )
    
    # Save plot
    if(!is.null(output_path)) {
      ggsave(output_path, combined_plot, width = width, height = height, dpi = 300)
      if(!is.null(sweet_spot_x)) {
        cat(sprintf("Identified sweet spot at proximity time: %.2f minutes\n", sweet_spot_x))
        cat(sprintf("Maximum total physical activity effect: %.2f minutes/day\n", max_total))
      }
    }
    
    return(combined_plot)
  } else {
    # Save individual plots if patchwork is not available
    if(!is.null(output_path)) {
      mobility_path <- gsub("\\.([^.]*)$", "_mobility.\\1", output_path)
      mvpa_path <- gsub("\\.([^.]*)$", "_mvpa.\\1", output_path)
      total_path <- gsub("\\.([^.]*)$", "_total.\\1", output_path)
      
      ggsave(mobility_path, p1, width = width, height = height/3, dpi = 300)
      ggsave(mvpa_path, p2, width = width, height = height/3, dpi = 300)
      ggsave(total_path, p3, width = width, height = height/3, dpi = 300)
      
      if(!is.null(sweet_spot_x)) {
        cat(sprintf("Identified sweet spot at proximity time: %.2f minutes\n", sweet_spot_x))
        cat(sprintf("Maximum total physical activity effect: %.2f minutes/day\n", max_total))
      }
    }
    
    return(list(mobility = p1, mvpa = p2, total = p3))
  }
}

variable_labels <- c(
  b0 = "Intercept",
  pois_all = 'Number of POIs (All) per 1,000 inhabitants',
  'MVPA' = 'Moderate|Vigorous physical activity',
  'Age' = 'Age',
  'Sex' = 'Sex',
  'work_statusHousewife/Househusband' = 'Work status - Housewife/Househusband',
  'work_statusManual self-employed worker' = 'Work status - Manual self-employed worker',
  'work_statusManual worker' = 'Work status - Manual worker',
  'work_statusNon-manual manager' = 'Work status - Non-manual manager',
  'work_statusNon-manual worker' = 'Work status - Non-manual worker',
  'educ_lvlPrimary' = 'Education level - Primary',
  'educ_lvlSecondary' = 'Tertiary',
  'smoking_statusEx-smoker' = 'Smoking status - Ex-smoker',
  'smoking_statusNever-smoker' = 'Smoking status - Never-smoker',
  'bmi_catObesity' = 'BMI category - Obesity',
  'bmi_catOverweight' = 'BMI category - Overweight',
  'bmi_catUnderweight' = 'BMI category - Underweight',
  'seasonSpring' = 'Season - Spring',
  'seasonSummer' = 'Season - Summer',
  'seasonWinter' = 'Season - Winter',
  'appt_year' = 'Year of participation',
  'Age:Sex' = 'Age x Sex',
  'pt_all' = 'Proximity time'
)