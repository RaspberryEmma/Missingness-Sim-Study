# ****************************************
# Missingness Simulation Study
#
# Simulation procedure for confounder-handling
# with missingness considerations
#
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 01/02/2026
# ****************************************


# ----- Missingness mechanisms and missingness handling -----

# See: https://www.rdocumentation.org/packages/mice/versions/3.17.0/topics/mice
apply_stacked_MI <- function(data = NULL, num_datasets = NULL, repetitions = NULL, imp_method = NULL) {
  imp_data <- NULL
  
  # find all vars containing missingness to be imputed
  if (binary_Z) {
    vars_with_missingness <- colnames(data)[ apply(data, 2, anyNA) ]
    for (var in vars_with_missingness) {
      data[, var] <- as.factor(data[, var])
    }
  }
  
  capture.output(                      # suppress command line output
    imp <- mice(data,
                m      = num_datasets, # number of imputations
                maxit  = repetitions,  # number of iterations
                method = imp_method)
  )
  
  imp_data <- complete(imp,
                       action  = "long", # stacked
                       include = FALSE)  # do not include original data
  
  imp_data <- subset(imp_data, select=-c(.imp, .id))
  
  return (imp_data)
}

run_stacked_MI_simulation <- function(n_scenario = NULL,
                                      n_obs      = NULL,
                                      n_rep      = NULL,
                                      
                                      Z_correlation     = NULL,
                                      Z_subgroups       = NULL,
                                      target_r_sq_X     = NULL,
                                      target_r_sq_Y     = NULL,
                                      causal            = NULL,
                                      
                                      num_total_conf  = NULL,
                                      num_meas_conf   = NULL,
                                      num_unmeas_conf = NULL,
                                      
                                      vars_to_make_unmeasured = NULL,
                                      vars_to_censor          = NULL) {

  # ----- Recording results -----
  
  var_sel_methods <- c("fully_adjusted", "unadjusted", "two_step_lasso", "two_step_lasso_X", "two_step_lasso_union")
  
  results_methods <- c("causal_true_value", "causal_estimate", "causal_bias", "causal_bias_proportion", "causal_coverage",
                       "open_paths", "blocked_paths", "proportion_paths",
                       "empirical_SE", "model_SE")
  
  var_names                         <- c("Y", "X", paste('Z', c(1:num_total_conf), sep=''))
  var_names_except_Y                <- var_names[ !var_names == 'Y']
  var_names_except_Y_with_intercept <- c("(Intercept)", var_names_except_Y)
  
  n_variables       <- length(var_names)
  n_var_sel_methods <- length(var_sel_methods)
  n_results         <- length(results_methods)
  
  FULL_coefs         <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, (n_variables - 1 + 1), n_rep),
                              dimnames = list(var_sel_methods, var_names_except_Y_with_intercept, 1:n_rep) )
  MNAR_coefs         <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, (n_variables - 1 + 1), n_rep),
                              dimnames = list(var_sel_methods, var_names_except_Y_with_intercept, 1:n_rep) )
  MCAR_coefs         <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, (n_variables - 1 + 1), n_rep),
                              dimnames = list(var_sel_methods, var_names_except_Y_with_intercept, 1:n_rep) )
  
  FULL_cov_selection <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, (n_variables - 1 + 1), n_rep),
                              dimnames = list(var_sel_methods, var_names_except_Y_with_intercept, 1:n_rep) )
  MNAR_cov_selection <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, (n_variables - 1 + 1), n_rep),
                              dimnames = list(var_sel_methods, var_names_except_Y_with_intercept, 1:n_rep) )
  MCAR_cov_selection <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, (n_variables - 1 + 1), n_rep),
                              dimnames = list(var_sel_methods, var_names_except_Y_with_intercept, 1:n_rep) )
  
  FULL_results       <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, n_results, n_rep),
                              dimnames = list(var_sel_methods, results_methods, 1:n_rep))
  MNAR_results       <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, n_results, n_rep),
                              dimnames = list(var_sel_methods, results_methods, 1:n_rep))
  MCAR_results       <- array(data     = NaN,
                              dim      = c(n_var_sel_methods, n_results, n_rep),
                              dimnames = list(var_sel_methods, results_methods, 1:n_rep))
  
  sample_size_table <- array(data     = NaN,
                             dim      = c(3, 2),
                             dimnames = list(c("None", "MNAR", "MCAR"), c("complete_cases", "sample_size_after_handling")))
  
  
  # ----- Simulation procedure ------
  
  message("\n\n\n ***** Running simulation procedure, scenario ", n_scenario, " *****")
  
  # track time
  current_time  <- as.numeric(Sys.time())*1000
  remaining_rep <- n_rep + 1
  
  for (repetition in c(1:n_rep)) {
    # track time
    remaining_rep   <- remaining_rep - 1
    previous_time   <- current_time
    current_time    <- as.numeric(Sys.time())*1000
    remaining_time  <- remaining_rep * (current_time - previous_time)
    expected_finish <- signif(((current_time + remaining_time)/1000), digits = 13) %>% as.POSIXct()
    
    print(paste0("Running repetition ", repetition, "/", n_rep, ", Expected to finish running at: ", expected_finish))
  
    FULL_dataset <- generate_dataset(n_obs      = n_obs,
                                     n_rep      = n_rep,
                                     
                                     Z_correlation     = Z_correlation,
                                     Z_subgroups       = Z_subgroups,
                                     target_r_sq_X     = target_r_sq_X,
                                     target_r_sq_Y     = target_r_sq_Y,
                                     causal            = causal,
                                     
                                     num_total_conf  = num_total_conf,
                                     num_meas_conf   = num_meas_conf,
                                     num_unmeas_conf = num_unmeas_conf,
                                     
                                     vars_to_make_unmeasured = vars_to_make_unmeasured,
                                     vars_to_censor          = vars_to_censor,
                                     var_names               = var_names)
    
    MNAR_data       <- apply_MNAR_missingness(FULL_dataset, vars_to_censor = vars_to_censor)
    MNAR_dataset    <- MNAR_data[[1]]
    MNAR_psel       <- MNAR_data[[2]]
    MNAR_censorship <- MNAR_data[[3]]
    
    MCAR_data       <- apply_MCAR_missingness(FULL_dataset, vars_to_censor = vars_to_censor)
    MCAR_dataset    <- MCAR_data[[1]]
    MCAR_psel       <- MCAR_data[[2]]
    MCAR_censorship <- MCAR_data[[3]]
    
    message("\n\n Missingness Mechanisms")
    print("MNAR p(selection into sample)")
    print(summary(as.factor(MNAR_psel)))
    
    print("MCAR p(selection into sample)")
    print(summary(as.factor(MCAR_psel)))
    
    # apply stacked MI
    # imp_method = norm -> MI using Bayesian linear regression
    # see: section "Details" of https://www.rdocumentation.org/packages/mice/versions/3.17.0/topics/mice
    handled_MNAR_dataset <- apply_stacked_MI(data         = MNAR_dataset,
                                             num_datasets = 5,
                                             repetitions  = 20,
                                             imp_method   = "norm")
    handled_MCAR_dataset <- apply_stacked_MI(data         = MCAR_dataset,
                                             num_datasets = 5,
                                             repetitions  = 20,
                                             imp_method   = "norm")
    
    # record sample sizes before and after missingness handling is applied
    sample_size_table["None", "complete_cases"]             <- dim(FULL_dataset)[1]
    sample_size_table["None", "sample_size_after_handling"] <- dim(FULL_dataset)[1]
    sample_size_table["MNAR", "complete_cases"]             <- dim(MNAR_dataset[complete.cases(MNAR_dataset), ])[1]
    sample_size_table["MNAR", "sample_size_after_handling"] <- dim(handled_MNAR_dataset)[1]
    sample_size_table["MCAR", "complete_cases"]             <- dim(MCAR_dataset[complete.cases(MCAR_dataset), ])[1]
    sample_size_table["MCAR", "sample_size_after_handling"] <- dim(handled_MCAR_dataset)[1]
    
    # cut-up versions of the data as needed
    X_FULL_dataset <- subset(FULL_dataset, select=-c(Y))
    Z_FULL_dataset <- subset(FULL_dataset, select=-c(Y, X))
    Y_FULL_column  <- subset(FULL_dataset, select=c(Y))
    X_FULL_column  <- subset(FULL_dataset, select=c(X))
    
    X_handled_MNAR_dataset <- subset(handled_MNAR_dataset, select=-c(Y))
    Z_handled_MNAR_dataset <- subset(handled_MNAR_dataset, select=-c(Y, X))
    Y_handled_MNAR_column  <- subset(handled_MNAR_dataset, select=c(Y))
    X_handled_MNAR_column  <- subset(handled_MNAR_dataset, select=c(X))
    
    X_handled_MCAR_dataset <- subset(handled_MCAR_dataset, select=-c(Y))
    Z_handled_MCAR_dataset <- subset(handled_MCAR_dataset, select=-c(Y, X))
    Y_handled_MCAR_column  <- subset(handled_MCAR_dataset, select=c(Y))
    X_handled_MCAR_column  <- subset(handled_MCAR_dataset, select=c(X))
    
    
    # fully adjusted model
    fully_adjusted_FULL_model  <- lm("Y ~ .", data = FULL_dataset)
    fully_adjusted_MNAR_model  <- lm("Y ~ .", data = handled_MNAR_dataset)
    fully_adjusted_MCAR_model  <- lm("Y ~ .", data = handled_MCAR_dataset)
    
    # unadjusted model
    unadjusted_FULL_model  <- lm("Y ~ X", data = FULL_dataset)
    unadjusted_MNAR_model  <- lm("Y ~ X", data = handled_MNAR_dataset)
    unadjusted_MCAR_model  <- lm("Y ~ X", data = handled_MCAR_dataset)
    
    # LASSO (outcome Y)
    
    # full dataset
    
    cv_lasso_FULL_model <- cv.glmnet(x = data.matrix(X_FULL_dataset), y = data.matrix(Y_FULL_column), alpha=1)
    lambda              <- cv_lasso_FULL_model$lambda.min
    lasso_FULL_model    <- glmnet(x = data.matrix(X_FULL_dataset), y = data.matrix(Y_FULL_column), alpha=1, lambda=lambda)
    
    lasso_FULL_coefs          <- as.vector(lasso_FULL_model$beta)
    names(lasso_FULL_coefs)   <- rownames(lasso_FULL_model$beta)
    lasso_FULL_vars_selected  <- union(c('X'), names(lasso_FULL_coefs[lasso_FULL_coefs != 0.0])) # always select X
    lasso_FULL_vars_selected  <- lasso_FULL_vars_selected[lasso_FULL_vars_selected != "(Intercept)"] # exclude intercept
    two_step_lasso_FULL_model <- lm(make_model_formula(vars_selected = lasso_FULL_vars_selected), data = FULL_dataset)
    
    # MNAR
    
    cv_lasso_MNAR_model <- cv.glmnet(x = data.matrix(X_handled_MNAR_dataset), y = data.matrix(Y_handled_MNAR_column), alpha=1)
    lambda              <- cv_lasso_MNAR_model$lambda.min
    lasso_MNAR_model    <- glmnet(x = data.matrix(X_handled_MNAR_dataset), y = data.matrix(Y_handled_MNAR_column), alpha=1, lambda=lambda)
    
    lasso_MNAR_coefs          <- as.vector(lasso_MNAR_model$beta)
    names(lasso_MNAR_coefs)   <- rownames(lasso_MNAR_model$beta)
    lasso_MNAR_vars_selected  <- union(c('X'), names(lasso_MNAR_coefs[lasso_MNAR_coefs != 0.0])) # always select X
    lasso_MNAR_vars_selected  <- lasso_MNAR_vars_selected[lasso_MNAR_vars_selected != "(Intercept)"] # exclude intercept
    two_step_lasso_MNAR_model <- lm(make_model_formula(vars_selected = lasso_MNAR_vars_selected), data = handled_MNAR_dataset)
    
    # MCAR
    
    cv_lasso_MCAR_model <- cv.glmnet(x = data.matrix(X_handled_MCAR_dataset), y = data.matrix(Y_handled_MCAR_column), alpha=1)
    lambda              <- cv_lasso_MCAR_model$lambda.min
    lasso_MCAR_model    <- glmnet(x = data.matrix(X_handled_MCAR_dataset), y = data.matrix(Y_handled_MCAR_column), alpha=1, lambda=lambda)
    
    lasso_MCAR_coefs          <- as.vector(lasso_MCAR_model$beta)
    names(lasso_MCAR_coefs)   <- rownames(lasso_MCAR_model$beta)
    lasso_MCAR_vars_selected  <- union(c('X'), names(lasso_MCAR_coefs[lasso_MCAR_coefs != 0.0])) # always select X
    lasso_MCAR_vars_selected  <- lasso_MCAR_vars_selected[lasso_MCAR_vars_selected != "(Intercept)"] # exclude intercept
    two_step_lasso_MCAR_model <- lm(make_model_formula(vars_selected = lasso_MCAR_vars_selected), data = handled_MCAR_dataset)
    
    # LASSO_X (exposure X)
      
    # full dataset
    
    cv_lasso_X_FULL_model <- cv.glmnet(x = data.matrix(Z_FULL_dataset), y = data.matrix(X_FULL_column), alpha=1)
    lambda                <- cv_lasso_X_FULL_model$lambda.min
    lasso_X_FULL_model    <- glmnet(x = data.matrix(Z_FULL_dataset), y = data.matrix(X_FULL_column), alpha=1, lambda=lambda)
    
    lasso_X_FULL_coefs          <- as.vector(lasso_X_FULL_model$beta)
    names(lasso_X_FULL_coefs)   <- rownames(lasso_X_FULL_model$beta)
    lasso_X_FULL_vars_selected  <- union(c('X'), names(lasso_X_FULL_coefs[lasso_X_FULL_coefs != 0.0])) # always select X
    lasso_X_FULL_vars_selected  <- lasso_X_FULL_vars_selected[lasso_X_FULL_vars_selected != "(Intercept)"] # exclude intercept
    two_step_lasso_X_FULL_model <- lm(make_model_formula(vars_selected = lasso_X_FULL_vars_selected), data = FULL_dataset)
    
    # MNAR
    
    cv_lasso_X_MNAR_model <- cv.glmnet(x = data.matrix(Z_handled_MNAR_dataset), y = data.matrix(X_handled_MNAR_column), alpha=1)
    lambda                <- cv_lasso_X_MNAR_model$lambda.min
    lasso_X_MNAR_model    <- glmnet(x = data.matrix(Z_handled_MNAR_dataset), y = data.matrix(X_handled_MNAR_column), alpha=1, lambda=lambda)
    
    lasso_X_MNAR_coefs          <- as.vector(lasso_X_MNAR_model$beta)
    names(lasso_X_MNAR_coefs)   <- rownames(lasso_X_MNAR_model$beta)
    lasso_X_MNAR_vars_selected  <- union(c('X'), names(lasso_X_MNAR_coefs[lasso_X_MNAR_coefs != 0.0])) # always select X
    lasso_X_MNAR_vars_selected  <- lasso_X_MNAR_vars_selected[lasso_X_MNAR_vars_selected != "(Intercept)"] # exclude intercept
    two_step_lasso_X_MNAR_model <- lm(make_model_formula(vars_selected = lasso_X_MNAR_vars_selected), data = handled_MNAR_dataset)
    
    # MCAR
    
    cv_lasso_X_MCAR_model <- cv.glmnet(x = data.matrix(Z_handled_MCAR_dataset), y = data.matrix(X_handled_MCAR_column), alpha=1)
    lambda                <- cv_lasso_X_MCAR_model$lambda.min
    lasso_X_MCAR_model    <- glmnet(x = data.matrix(Z_handled_MCAR_dataset), y = data.matrix(X_handled_MCAR_column), alpha=1, lambda=lambda)
    
    lasso_X_MCAR_coefs          <- as.vector(lasso_X_MCAR_model$beta)
    names(lasso_X_MCAR_coefs)   <- rownames(lasso_X_MCAR_model$beta)
    lasso_X_MCAR_vars_selected  <- union(c('X'), names(lasso_X_MCAR_coefs[lasso_X_MCAR_coefs != 0.0])) # always select X
    lasso_X_MCAR_vars_selected  <- lasso_X_MCAR_vars_selected[lasso_X_MCAR_vars_selected != "(Intercept)"] # exclude intercept
    two_step_lasso_X_MCAR_model <- lm(make_model_formula(vars_selected = lasso_X_MCAR_vars_selected), data = handled_MCAR_dataset)
    
    # LASSO_union (exposure X and outcome Y)
    
    # full dataset
    lasso_union_FULL_vars_selected  <- union(lasso_FULL_vars_selected, lasso_X_FULL_vars_selected)
    two_step_lasso_union_FULL_model <- lm(make_model_formula(vars_selected = lasso_union_FULL_vars_selected), data = FULL_dataset)
    
    # MNAR
    lasso_union_MNAR_vars_selected  <- union(lasso_MNAR_vars_selected, lasso_X_MNAR_vars_selected)
    two_step_lasso_union_MNAR_model <- lm(make_model_formula(vars_selected = lasso_union_MNAR_vars_selected), data = handled_MNAR_dataset)
    
    # MCAR
    lasso_union_MCAR_vars_selected  <- union(lasso_MCAR_vars_selected, lasso_X_MCAR_vars_selected)
    two_step_lasso_union_MCAR_model <- lm(make_model_formula(vars_selected = lasso_union_MCAR_vars_selected), data = handled_MCAR_dataset)
    
    
    # ------ Record covariate selection ------
    
    FULL_cov_selection["fully_adjusted", , repetition]       <- fill_in_cov_selection(fully_adjusted_FULL_model$coefficients, var_names_except_Y_with_intercept)
    FULL_cov_selection["unadjusted", , repetition]           <- fill_in_cov_selection(unadjusted_FULL_model$coefficients, var_names_except_Y_with_intercept)
    FULL_cov_selection["two_step_lasso", , repetition]       <- fill_in_cov_selection(two_step_lasso_FULL_model$coefficients, var_names_except_Y_with_intercept)
    FULL_cov_selection["two_step_lasso_X", , repetition]     <- fill_in_cov_selection(two_step_lasso_X_FULL_model$coefficients, var_names_except_Y_with_intercept)
    FULL_cov_selection["two_step_lasso_union", , repetition] <- fill_in_cov_selection(two_step_lasso_union_FULL_model$coefficients, var_names_except_Y_with_intercept)
    
    MNAR_cov_selection["fully_adjusted", , repetition]       <- fill_in_cov_selection(fully_adjusted_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    MNAR_cov_selection["unadjusted", , repetition]           <- fill_in_cov_selection(unadjusted_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    MNAR_cov_selection["two_step_lasso", , repetition]       <- fill_in_cov_selection(two_step_lasso_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    MNAR_cov_selection["two_step_lasso_X", , repetition]     <- fill_in_cov_selection(two_step_lasso_X_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    MNAR_cov_selection["two_step_lasso_union", , repetition] <- fill_in_cov_selection(two_step_lasso_union_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    
    MCAR_cov_selection["fully_adjusted", , repetition]       <- fill_in_cov_selection(fully_adjusted_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    MCAR_cov_selection["unadjusted", , repetition]           <- fill_in_cov_selection(unadjusted_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    MCAR_cov_selection["two_step_lasso", , repetition]       <- fill_in_cov_selection(two_step_lasso_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    MCAR_cov_selection["two_step_lasso_X", , repetition]     <- fill_in_cov_selection(two_step_lasso_X_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    MCAR_cov_selection["two_step_lasso_union", , repetition] <- fill_in_cov_selection(two_step_lasso_union_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    
    
    
    # ------ Record coefficients ------
    
    FULL_coefs["fully_adjusted", , repetition]       <- fill_in_blanks(fully_adjusted_FULL_model$coefficients, var_names_except_Y_with_intercept)
    FULL_coefs["unadjusted", , repetition]           <- fill_in_blanks(unadjusted_FULL_model$coefficients, var_names_except_Y_with_intercept)
    FULL_coefs["two_step_lasso", , repetition]       <- fill_in_blanks(two_step_lasso_FULL_model$coefficients, var_names_except_Y_with_intercept)
    FULL_coefs["two_step_lasso_X", , repetition]     <- fill_in_blanks(two_step_lasso_X_FULL_model$coefficients, var_names_except_Y_with_intercept)
    FULL_coefs["two_step_lasso_union", , repetition] <- fill_in_blanks(two_step_lasso_union_FULL_model$coefficients, var_names_except_Y_with_intercept)
    
    MNAR_coefs["fully_adjusted", , repetition]       <- fill_in_blanks(fully_adjusted_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    MNAR_coefs["unadjusted", , repetition]           <- fill_in_blanks(unadjusted_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    MNAR_coefs["two_step_lasso", , repetition]       <- fill_in_blanks(two_step_lasso_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    MNAR_coefs["two_step_lasso_X", , repetition]     <- fill_in_blanks(two_step_lasso_X_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    MNAR_coefs["two_step_lasso_union", , repetition] <- fill_in_blanks(two_step_lasso_union_MNAR_model$coefficients, var_names_except_Y_with_intercept)
    
    MCAR_coefs["fully_adjusted", , repetition]       <- fill_in_blanks(fully_adjusted_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    MCAR_coefs["unadjusted", , repetition]           <- fill_in_blanks(unadjusted_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    MCAR_coefs["two_step_lasso", , repetition]       <- fill_in_blanks(two_step_lasso_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    MCAR_coefs["two_step_lasso_X", , repetition]     <- fill_in_blanks(two_step_lasso_X_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    MCAR_coefs["two_step_lasso_union", , repetition] <- fill_in_blanks(two_step_lasso_union_MCAR_model$coefficients, var_names_except_Y_with_intercept)
    
    
    
    # ------ Record results ------
    
    FULL_results["fully_adjusted", "causal_true_value", repetition]      <- causal
    FULL_results["fully_adjusted", "causal_estimate", repetition]        <- unname(fully_adjusted_FULL_model$coefficients['X'])
    FULL_results["fully_adjusted", "causal_bias", repetition]            <- (unname(fully_adjusted_FULL_model$coefficients['X']) - causal)
    FULL_results["fully_adjusted", "causal_bias_proportion", repetition] <- ((unname(fully_adjusted_FULL_model$coefficients['X']) - causal)/causal)
    FULL_results["fully_adjusted", "causal_coverage", repetition]        <- estimate_within_CI(model = fully_adjusted_FULL_model)
    FULL_results["fully_adjusted", "open_paths", repetition]             <- num_total_conf
    FULL_results["fully_adjusted", "blocked_paths", repetition]          <- num_meas_conf
    FULL_results["fully_adjusted", "proportion_paths", repetition]       <- num_meas_conf / num_total_conf
    FULL_results["fully_adjusted", "empirical_SE", repetition]           <- NaN
    FULL_results["fully_adjusted", "model_SE", repetition]               <- (coef(summary(fully_adjusted_FULL_model))[, "Std. Error"])['X']
    
    FULL_results["unadjusted", "causal_true_value", repetition]      <- causal
    FULL_results["unadjusted", "causal_estimate", repetition]        <- unname(unadjusted_FULL_model$coefficients['X'])
    FULL_results["unadjusted", "causal_bias", repetition]            <- (unname(unadjusted_FULL_model$coefficients['X']) - causal)
    FULL_results["unadjusted", "causal_bias_proportion", repetition] <- ((unname(unadjusted_FULL_model$coefficients['X']) - causal)/causal)
    FULL_results["unadjusted", "causal_coverage", repetition]        <- estimate_within_CI(model = unadjusted_FULL_model)
    FULL_results["unadjusted", "open_paths", repetition]             <- num_total_conf
    FULL_results["unadjusted", "blocked_paths", repetition]          <- 0.0
    FULL_results["unadjusted", "proportion_paths", repetition]       <- 0.0 / num_total_conf
    FULL_results["unadjusted", "empirical_SE", repetition]           <- NaN
    FULL_results["unadjusted", "model_SE", repetition]               <- (coef(summary(unadjusted_FULL_model))[, "Std. Error"])['X']
    
    FULL_results["two_step_lasso", "causal_true_value", repetition]      <- causal
    FULL_results["two_step_lasso", "causal_estimate", repetition]        <- unname(two_step_lasso_FULL_model$coefficients['X'])
    FULL_results["two_step_lasso", "causal_bias", repetition]            <- (unname(two_step_lasso_FULL_model$coefficients['X']) - causal)
    FULL_results["two_step_lasso", "causal_bias_proportion", repetition] <- ((unname(two_step_lasso_FULL_model$coefficients['X']) - causal)/causal)
    FULL_results["two_step_lasso", "causal_coverage", repetition]        <- estimate_within_CI(model = two_step_lasso_FULL_model)
    FULL_results["two_step_lasso", "open_paths", repetition]             <- num_total_conf
    FULL_results["two_step_lasso", "blocked_paths", repetition]          <- length(lasso_FULL_vars_selected[lasso_FULL_vars_selected != 'X'])
    FULL_results["two_step_lasso", "proportion_paths", repetition]       <- length(lasso_FULL_vars_selected[lasso_FULL_vars_selected != 'X']) / num_total_conf
    FULL_results["two_step_lasso", "empirical_SE", repetition]           <- NaN
    FULL_results["two_step_lasso", "model_SE", repetition]               <- (coef(summary(two_step_lasso_FULL_model))[, "Std. Error"])['X']
    
    FULL_results["two_step_lasso_X", "causal_true_value", repetition]      <- causal
    FULL_results["two_step_lasso_X", "causal_estimate", repetition]        <- unname(two_step_lasso_X_FULL_model$coefficients['X'])
    FULL_results["two_step_lasso_X", "causal_bias", repetition]            <- (unname(two_step_lasso_X_FULL_model$coefficients['X']) - causal)
    FULL_results["two_step_lasso_X", "causal_bias_proportion", repetition] <- ((unname(two_step_lasso_X_FULL_model$coefficients['X']) - causal)/causal)
    FULL_results["two_step_lasso_X", "causal_coverage", repetition]        <- estimate_within_CI(model = two_step_lasso_X_FULL_model)
    FULL_results["two_step_lasso_X", "open_paths", repetition]             <- num_total_conf
    FULL_results["two_step_lasso_X", "blocked_paths", repetition]          <- length(lasso_X_FULL_vars_selected[lasso_X_FULL_vars_selected != 'X'])
    FULL_results["two_step_lasso_X", "proportion_paths", repetition]       <- length(lasso_X_FULL_vars_selected[lasso_X_FULL_vars_selected != 'X']) / num_total_conf
    FULL_results["two_step_lasso_X", "empirical_SE", repetition]           <- NaN
    FULL_results["two_step_lasso_X", "model_SE", repetition]               <- (coef(summary(two_step_lasso_X_FULL_model))[, "Std. Error"])['X']
    
    FULL_results["two_step_lasso_union", "causal_true_value", repetition]      <- causal
    FULL_results["two_step_lasso_union", "causal_estimate", repetition]        <- unname(two_step_lasso_union_FULL_model$coefficients['X'])
    FULL_results["two_step_lasso_union", "causal_bias", repetition]            <- (unname(two_step_lasso_union_FULL_model$coefficients['X']) - causal)
    FULL_results["two_step_lasso_union", "causal_bias_proportion", repetition] <- ((unname(two_step_lasso_union_FULL_model$coefficients['X']) - causal)/causal)
    FULL_results["two_step_lasso_union", "causal_coverage", repetition]        <- estimate_within_CI(model = two_step_lasso_union_FULL_model)
    FULL_results["two_step_lasso_union", "open_paths", repetition]             <- num_total_conf
    FULL_results["two_step_lasso_union", "blocked_paths", repetition]          <- length(lasso_union_FULL_vars_selected[lasso_union_FULL_vars_selected != 'X'])
    FULL_results["two_step_lasso_union", "proportion_paths", repetition]       <- length(lasso_union_FULL_vars_selected[lasso_union_FULL_vars_selected != 'X']) / num_total_conf
    FULL_results["two_step_lasso_union", "empirical_SE", repetition]           <- NaN
    FULL_results["two_step_lasso_union", "model_SE", repetition]               <- (coef(summary(two_step_lasso_union_FULL_model))[, "Std. Error"])['X']
    
    MNAR_results["fully_adjusted", "causal_true_value", repetition]      <- causal
    MNAR_results["fully_adjusted", "causal_estimate", repetition]        <- unname(fully_adjusted_MNAR_model$coefficients['X'])
    MNAR_results["fully_adjusted", "causal_bias", repetition]            <- (unname(fully_adjusted_MNAR_model$coefficients['X']) - causal)
    MNAR_results["fully_adjusted", "causal_bias_proportion", repetition] <- ((unname(fully_adjusted_MNAR_model$coefficients['X']) - causal)/causal)
    MNAR_results["fully_adjusted", "causal_coverage", repetition]        <- estimate_within_CI(model = fully_adjusted_MNAR_model)
    MNAR_results["fully_adjusted", "open_paths", repetition]             <- num_total_conf
    MNAR_results["fully_adjusted", "blocked_paths", repetition]          <- num_meas_conf
    MNAR_results["fully_adjusted", "proportion_paths", repetition]       <- num_meas_conf / num_total_conf
    MNAR_results["fully_adjusted", "empirical_SE", repetition]           <- NaN
    MNAR_results["fully_adjusted", "model_SE", repetition]               <- (coef(summary(fully_adjusted_MNAR_model))[, "Std. Error"])['X']
    
    MNAR_results["unadjusted", "causal_true_value", repetition]      <- causal
    MNAR_results["unadjusted", "causal_estimate", repetition]        <- unname(unadjusted_MNAR_model$coefficients['X'])
    MNAR_results["unadjusted", "causal_bias", repetition]            <- (unname(unadjusted_MNAR_model$coefficients['X']) - causal)
    MNAR_results["unadjusted", "causal_bias_proportion", repetition] <- ((unname(unadjusted_MNAR_model$coefficients['X']) - causal)/causal)
    MNAR_results["unadjusted", "causal_coverage", repetition]        <- estimate_within_CI(model = unadjusted_MNAR_model)
    MNAR_results["unadjusted", "open_paths", repetition]             <- num_total_conf
    MNAR_results["unadjusted", "blocked_paths", repetition]          <- 0.0
    MNAR_results["unadjusted", "proportion_paths", repetition]       <- 0.0 / num_total_conf
    MNAR_results["unadjusted", "empirical_SE", repetition]           <- NaN
    MNAR_results["unadjusted", "model_SE", repetition]               <- (coef(summary(unadjusted_MNAR_model))[, "Std. Error"])['X']
    
    MNAR_results["two_step_lasso", "causal_true_value", repetition]      <- causal
    MNAR_results["two_step_lasso", "causal_estimate", repetition]        <- unname(two_step_lasso_MNAR_model$coefficients['X'])
    MNAR_results["two_step_lasso", "causal_bias", repetition]            <- (unname(two_step_lasso_MNAR_model$coefficients['X']) - causal)
    MNAR_results["two_step_lasso", "causal_bias_proportion", repetition] <- ((unname(two_step_lasso_MNAR_model$coefficients['X']) - causal)/causal)
    MNAR_results["two_step_lasso", "causal_coverage", repetition]        <- estimate_within_CI(model = two_step_lasso_MNAR_model)
    MNAR_results["two_step_lasso", "open_paths", repetition]             <- num_total_conf
    MNAR_results["two_step_lasso", "blocked_paths", repetition]          <- length(lasso_MNAR_vars_selected[lasso_MNAR_vars_selected != 'X'])
    MNAR_results["two_step_lasso", "proportion_paths", repetition]       <- length(lasso_MNAR_vars_selected[lasso_MNAR_vars_selected != 'X']) / num_total_conf
    MNAR_results["two_step_lasso", "empirical_SE", repetition]           <- NaN
    MNAR_results["two_step_lasso", "model_SE", repetition]               <- (coef(summary(two_step_lasso_MNAR_model))[, "Std. Error"])['X']
    
    MNAR_results["two_step_lasso_X", "causal_true_value", repetition]      <- causal
    MNAR_results["two_step_lasso_X", "causal_estimate", repetition]        <- unname(two_step_lasso_X_MNAR_model$coefficients['X'])
    MNAR_results["two_step_lasso_X", "causal_bias", repetition]            <- (unname(two_step_lasso_X_MNAR_model$coefficients['X']) - causal)
    MNAR_results["two_step_lasso_X", "causal_bias_proportion", repetition] <- ((unname(two_step_lasso_X_MNAR_model$coefficients['X']) - causal)/causal)
    MNAR_results["two_step_lasso_X", "causal_coverage", repetition]        <- estimate_within_CI(model = two_step_lasso_X_MNAR_model)
    MNAR_results["two_step_lasso_X", "open_paths", repetition]             <- num_total_conf
    MNAR_results["two_step_lasso_X", "blocked_paths", repetition]          <- length(lasso_X_MNAR_vars_selected[lasso_X_MNAR_vars_selected != 'X'])
    MNAR_results["two_step_lasso_X", "proportion_paths", repetition]       <- length(lasso_X_MNAR_vars_selected[lasso_X_MNAR_vars_selected != 'X']) / num_total_conf
    MNAR_results["two_step_lasso_X", "empirical_SE", repetition]           <- NaN
    MNAR_results["two_step_lasso_X", "model_SE", repetition]               <- (coef(summary(two_step_lasso_X_MNAR_model))[, "Std. Error"])['X']
    
    MNAR_results["two_step_lasso_union", "causal_true_value", repetition]      <- causal
    MNAR_results["two_step_lasso_union", "causal_estimate", repetition]        <- unname(two_step_lasso_union_MNAR_model$coefficients['X'])
    MNAR_results["two_step_lasso_union", "causal_bias", repetition]            <- (unname(two_step_lasso_union_MNAR_model$coefficients['X']) - causal)
    MNAR_results["two_step_lasso_union", "causal_bias_proportion", repetition] <- ((unname(two_step_lasso_union_MNAR_model$coefficients['X']) - causal)/causal)
    MNAR_results["two_step_lasso_union", "causal_coverage", repetition]        <- estimate_within_CI(model = two_step_lasso_union_MNAR_model)
    MNAR_results["two_step_lasso_union", "open_paths", repetition]             <- num_total_conf
    MNAR_results["two_step_lasso_union", "blocked_paths", repetition]          <- length(lasso_union_MNAR_vars_selected[lasso_union_MNAR_vars_selected != 'X'])
    MNAR_results["two_step_lasso_union", "proportion_paths", repetition]       <- length(lasso_union_MNAR_vars_selected[lasso_union_MNAR_vars_selected != 'X']) / num_total_conf
    MNAR_results["two_step_lasso_union", "empirical_SE", repetition]           <- NaN
    MNAR_results["two_step_lasso_union", "model_SE", repetition]               <- (coef(summary(two_step_lasso_union_MNAR_model))[, "Std. Error"])['X']
    
    MCAR_results["fully_adjusted", "causal_true_value", repetition]      <- causal
    MCAR_results["fully_adjusted", "causal_estimate", repetition]        <- unname(fully_adjusted_MCAR_model$coefficients['X'])
    MCAR_results["fully_adjusted", "causal_bias", repetition]            <- (unname(fully_adjusted_MCAR_model$coefficients['X']) - causal)
    MCAR_results["fully_adjusted", "causal_bias_proportion", repetition] <- ((unname(fully_adjusted_MCAR_model$coefficients['X']) - causal)/causal)
    MCAR_results["fully_adjusted", "causal_coverage", repetition]        <- estimate_within_CI(model = fully_adjusted_MCAR_model)
    MCAR_results["fully_adjusted", "open_paths", repetition]             <- num_total_conf
    MCAR_results["fully_adjusted", "blocked_paths", repetition]          <- num_meas_conf
    MCAR_results["fully_adjusted", "proportion_paths", repetition]       <- num_meas_conf / num_total_conf
    MCAR_results["fully_adjusted", "empirical_SE", repetition]           <- NaN
    MCAR_results["fully_adjusted", "model_SE", repetition]               <- (coef(summary(fully_adjusted_MCAR_model))[, "Std. Error"])['X']
    
    MCAR_results["unadjusted", "causal_true_value", repetition]      <- causal
    MCAR_results["unadjusted", "causal_estimate", repetition]        <- unname(unadjusted_MCAR_model$coefficients['X'])
    MCAR_results["unadjusted", "causal_bias", repetition]            <- (unname(unadjusted_MCAR_model$coefficients['X']) - causal)
    MCAR_results["unadjusted", "causal_bias_proportion", repetition] <- ((unname(unadjusted_MCAR_model$coefficients['X']) - causal)/causal)
    MCAR_results["unadjusted", "causal_coverage", repetition]        <- estimate_within_CI(model = unadjusted_MCAR_model)
    MCAR_results["unadjusted", "open_paths", repetition]             <- num_total_conf
    MCAR_results["unadjusted", "blocked_paths", repetition]          <- 0.0
    MCAR_results["unadjusted", "proportion_paths", repetition]       <- 0.0 / num_total_conf
    MCAR_results["unadjusted", "empirical_SE", repetition]           <- NaN
    MCAR_results["unadjusted", "model_SE", repetition]               <- (coef(summary(unadjusted_MCAR_model))[, "Std. Error"])['X']
    
    MCAR_results["two_step_lasso", "causal_true_value", repetition]      <- causal
    MCAR_results["two_step_lasso", "causal_estimate", repetition]        <- unname(two_step_lasso_MCAR_model$coefficients['X'])
    MCAR_results["two_step_lasso", "causal_bias", repetition]            <- (unname(two_step_lasso_MCAR_model$coefficients['X']) - causal)
    MCAR_results["two_step_lasso", "causal_bias_proportion", repetition] <- ((unname(two_step_lasso_MCAR_model$coefficients['X']) - causal)/causal)
    MCAR_results["two_step_lasso", "causal_coverage", repetition]        <- estimate_within_CI(model = two_step_lasso_MCAR_model)
    MCAR_results["two_step_lasso", "open_paths", repetition]             <- num_total_conf
    MCAR_results["two_step_lasso", "blocked_paths", repetition]          <- length(lasso_MCAR_vars_selected[lasso_MCAR_vars_selected != 'X'])
    MCAR_results["two_step_lasso", "proportion_paths", repetition]       <- length(lasso_MCAR_vars_selected[lasso_MCAR_vars_selected != 'X']) / num_total_conf
    MCAR_results["two_step_lasso", "empirical_SE", repetition]           <- NaN
    MCAR_results["two_step_lasso", "model_SE", repetition]               <- (coef(summary(two_step_lasso_MCAR_model))[, "Std. Error"])['X']
    
    MCAR_results["two_step_lasso_X", "causal_true_value", repetition]      <- causal
    MCAR_results["two_step_lasso_X", "causal_estimate", repetition]        <- unname(two_step_lasso_X_MCAR_model$coefficients['X'])
    MCAR_results["two_step_lasso_X", "causal_bias", repetition]            <- (unname(two_step_lasso_X_MCAR_model$coefficients['X']) - causal)
    MCAR_results["two_step_lasso_X", "causal_bias_proportion", repetition] <- ((unname(two_step_lasso_X_MCAR_model$coefficients['X']) - causal)/causal)
    MCAR_results["two_step_lasso_X", "causal_coverage", repetition]        <- estimate_within_CI(model = two_step_lasso_X_MCAR_model)
    MCAR_results["two_step_lasso_X", "open_paths", repetition]             <- num_total_conf
    MCAR_results["two_step_lasso_X", "blocked_paths", repetition]          <- length(lasso_X_MCAR_vars_selected[lasso_X_MCAR_vars_selected != 'X'])
    MCAR_results["two_step_lasso_X", "proportion_paths", repetition]       <- length(lasso_X_MCAR_vars_selected[lasso_X_MCAR_vars_selected != 'X']) / num_total_conf
    MCAR_results["two_step_lasso_X", "empirical_SE", repetition]           <- NaN
    MCAR_results["two_step_lasso_X", "model_SE", repetition]               <- (coef(summary(two_step_lasso_X_MCAR_model))[, "Std. Error"])['X']
    
    MCAR_results["two_step_lasso_union", "causal_true_value", repetition]      <- causal
    MCAR_results["two_step_lasso_union", "causal_estimate", repetition]        <- unname(two_step_lasso_union_MCAR_model$coefficients['X'])
    MCAR_results["two_step_lasso_union", "causal_bias", repetition]            <- (unname(two_step_lasso_union_MCAR_model$coefficients['X']) - causal)
    MCAR_results["two_step_lasso_union", "causal_bias_proportion", repetition] <- ((unname(two_step_lasso_union_MCAR_model$coefficients['X']) - causal)/causal)
    MCAR_results["two_step_lasso_union", "causal_coverage", repetition]        <- estimate_within_CI(model = two_step_lasso_union_MCAR_model)
    MCAR_results["two_step_lasso_union", "open_paths", repetition]             <- num_total_conf
    MCAR_results["two_step_lasso_union", "blocked_paths", repetition]          <- length(lasso_union_MCAR_vars_selected[lasso_union_MCAR_vars_selected != 'X'])
    MCAR_results["two_step_lasso_union", "proportion_paths", repetition]       <- length(lasso_union_MCAR_vars_selected[lasso_union_MCAR_vars_selected != 'X']) / num_total_conf
    MCAR_results["two_step_lasso_union", "empirical_SE", repetition]           <- NaN
    MCAR_results["two_step_lasso_union", "model_SE", repetition]               <- (coef(summary(two_step_lasso_union_MCAR_model))[, "Std. Error"])['X']
    
  } # repetitions loop
  
  # "level" coefficient for effect of Zs on X and Y before subgroups
  beta_X  <- beta_X_formula(num_total_conf = num_total_conf, # Z on X
                            target_r_sq_X  = target_r_sq_X,
                            Z_correlation  = Z_correlation)
  
  # variance of the error term for Y
  var_Y <- determine_subgroup_var_Y(num_total_conf = num_total_conf,
                                    beta_X         = beta_X,
                                    causal         = causal,
                                    Z_correlation  = Z_correlation,
                                    target_r_sq_Y  = target_r_sq_Y)
  
  TRUE_coefs        <- c(beta_X_subgroups_formula(beta_X = beta_X),
                         beta_Y_subgroups_formula(beta_X = beta_X),
                         causal,
                         determine_subgroup_var_error_Y(var_Y = var_Y, target_r_sq_Y  = target_r_sq_Y))
  
  names(TRUE_coefs) <- c("Z on X (Subgroup Low X, Low Y)",
                         "Z on X (Subgroup Low X, High Y)",
                         "Z on X (Subgroup High X, Low Y)",
                         "Z on X (Subgroup High X, High Y)",
                         
                         "Z on Y (Subgroup Low X, Low Y)",
                         "Z on Y (Subgroup Low X, High Y)",
                         "Z on Y (Subgroup High X, Low Y)",
                         "Z on Y (Subgroup High X, High Y)",
                         
                         "Causal effect X on Y",
                         "Var in error for Y")
  
  return (list(FULL_cov_selection,
               MNAR_cov_selection,
               MCAR_cov_selection,
               
               TRUE_coefs,
               FULL_coefs,
               MNAR_coefs,
               MCAR_coefs,
               
               FULL_results,
               MNAR_results,
               MCAR_results,
               
               sample_size_table,
               
               FULL_dataset,
               handled_MNAR_dataset,
               handled_MCAR_dataset
  ))
}


