# ****************************************
# Missingness Simulation Study
#
# Simulation procedure for confounder-handling
# with missingness considerations
#
# Emma Tarmey
#
# Started:          30/01/2026
# Most Recent Edit: 01/02/2026
# ****************************************


# ----- Complete data generation mechanism -----


# datasets representative of our DAG
generate_dataset <- function(n_obs      = NULL,
                             n_rep      = NULL,
                             
                             Z_correlation     = NULL,
                             Z_subgroups       = NULL,
                             target_r_sq_X     = NULL,
                             target_r_sq_Y     = NULL,
                             causal            = NULL,
                             
                             binary_X          = NULL,
                             binary_Y          = NULL,
                             binary_Z          = NULL,
                             
                             num_total_conf  = NULL,
                             num_meas_conf   = NULL,
                             num_unmeas_conf = NULL,
                             
                             vars_to_make_unmeasured = NULL,
                             vars_to_censor          = NULL,
                             var_names               = NULL) {
  
  # coefficient values for DAG
  alpha <- sqrt(Z_correlation)     # U on Z
  beta  <- sqrt(1 - Z_correlation) # error_Z on Z
  
  beta_X  <- beta_X_formula(num_total_conf = num_total_conf, # Z on X
                            target_r_sq_X  = target_r_sq_X,
                            Z_correlation  = Z_correlation)
  
  beta_Xs <- beta_X_subgroups_formula(beta_X = beta_X) # subgroup Z on X
  beta_Ys <- beta_Y_subgroups_formula(beta_X = beta_X) # subgroup Z on Y
  
  # variance of the error term for Y
  var_Y <- determine_subgroup_var_Y(num_total_conf = num_total_conf,
                                    beta_X         = beta_X,
                                    causal         = causal,
                                    Z_correlation  = Z_correlation,
                                    target_r_sq_Y  = target_r_sq_Y)
  
  var_error_Y <- determine_subgroup_var_error_Y(var_Y          = var_Y,
                                                target_r_sq_Y  = target_r_sq_Y)
  
  dataset <- data.frame(matrix(NaN, nrow = n_obs, ncol = length(var_names)))
  colnames(dataset) <- var_names
  
  # shared prior U for all Z_i
  prior_U <- rnorm(n = n_obs, mean = 0, sd = 1)
  
  # generate all confounders Z_i
  # low, low subgroup
  Z1  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z2  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z3  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z4  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z5  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z6  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z7  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z8  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  
  # low, high subgroup
  Z9  <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z10 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z11 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z12 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z13 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z14 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z15 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z16 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  
  # high, low subgroup
  Z17 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z18 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z19 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z20 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z21 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z22 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z23 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z24 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  
  # high, high subgroup
  Z25 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z26 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z27 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z28 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z29 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z30 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z31 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  Z32 <- (alpha * prior_U) + (beta * rnorm(n = n_obs, mean = 0, sd = 1))
  
  # generate exposure X
  X <- (beta_Xs[1] * (Z1  + Z2  + Z3  + Z4  + Z5  + Z6  + Z7  + Z8))  +
    (beta_Xs[2] * (Z9  + Z10 + Z11 + Z12 + Z13 + Z14 + Z15 + Z16)) +
    (beta_Xs[3] * (Z17 + Z18 + Z19 + Z20 + Z21 + Z22 + Z23 + Z24)) +
    (beta_Xs[4] * (Z25 + Z26 + Z27 + Z28 + Z29 + Z30 + Z31 + Z32))
  
  # generate outcome Y
  Y <- (beta_Ys[1] * (Z1  + Z2  + Z3  + Z4  + Z5  + Z6  + Z7  + Z8))  +
    (beta_Ys[2] * (Z9  + Z10 + Z11 + Z12 + Z13 + Z14 + Z15 + Z16)) +
    (beta_Ys[3] * (Z17 + Z18 + Z19 + Z20 + Z21 + Z22 + Z23 + Z24)) +
    (beta_Ys[4] * (Z25 + Z26 + Z27 + Z28 + Z29 + Z30 + Z31 + Z32))
  
  # add error term to X if X is not binary
  if (!binary_X) {
    error_X <- rnorm(n = n_obs, mean = 0, sd = 1)
    X <- X + error_X
  }
  
  # add error term to Y if Y is not binary
  if (!binary_Y) {
    error_Y <- rnorm(n = n_obs, mean = 0, sd = sqrt(var_error_Y))
    Y <- Y + error_Y
  }
  
  # binarize X if binary
  if (binary_X) {
    # NB: intercept term of logit expression controls prevalence (mean) of binary var
    logit_prob_X  <- X - 1.40                                   # interpret existing values as logit(probability)
    prob_X        <- inverse_logit(logit_prob_X)                # apply inverse to obtain prob values
    binary_vals_X <- rbinom(n = n_obs, size = 1, prob = prob_X) # re-sample to obtain X
    X             <- binary_vals_X                              # write binary values over previous continuous values
  }
  
  # add causal effect (X on Y)
  Y <- Y + (causal * X) 
  
  # binarize Y if binary
  if (binary_Y) {
    # NB: intercept term of logit expression controls prevalence (mean) of binary var
    # common Y: intercept = 1.45; rare Y: intercept = 3.71
    logit_prob_Y  <- Y - 1.45                                    # interpret existing values as logit(probability)
    prob_Y        <- inverse_logit(logit_prob_Y)                # apply inverse to obtain prob values
    binary_vals_Y <- rbinom(n = n_obs, size = 1, prob = prob_Y) # re-sample to obtain Y
    Y             <- binary_vals_Y                              # write binary values over previous continuous values
  }
  
  # record X, Y and Z_i
  dataset[, 'Y']   <- Y
  dataset[, 'X']   <- X
  
  # Low, low subgroup
  dataset[, 'Z1']  <- Z1
  dataset[, 'Z2']  <- Z2
  dataset[, 'Z3']  <- Z3
  dataset[, 'Z4']  <- Z4
  dataset[, 'Z5']  <- Z5
  dataset[, 'Z6']  <- Z6
  dataset[, 'Z7']  <- Z7
  dataset[, 'Z8']  <- Z8
  
  # low, high subgroup
  dataset[, 'Z9']  <- Z9
  dataset[, 'Z10'] <- Z10
  dataset[, 'Z11'] <- Z11
  dataset[, 'Z12'] <- Z12
  dataset[, 'Z13'] <- Z13
  dataset[, 'Z14'] <- Z14
  dataset[, 'Z15'] <- Z15
  dataset[, 'Z16'] <- Z16
  
  # high, low subgroup
  dataset[, 'Z17'] <- Z17
  dataset[, 'Z18'] <- Z18
  dataset[, 'Z19'] <- Z19
  dataset[, 'Z20'] <- Z20
  dataset[, 'Z21'] <- Z21
  dataset[, 'Z22'] <- Z22
  dataset[, 'Z23'] <- Z23
  dataset[, 'Z24'] <- Z24
  
  # high, high subgroup
  dataset[, 'Z25'] <- Z25
  dataset[, 'Z26'] <- Z26
  dataset[, 'Z27'] <- Z27
  dataset[, 'Z28'] <- Z28
  dataset[, 'Z29'] <- Z29
  dataset[, 'Z30'] <- Z30
  dataset[, 'Z31'] <- Z31
  dataset[, 'Z32'] <- Z32
  
  # censor covariates as appropriate
  dataset <- dataset[, !(names(dataset) %in% vars_to_make_unmeasured)]
  
  return (dataset)
}



# ----- Missingness mechanisms and missingness handling -----



# See: https://www.rdocumentation.org/packages/mice/versions/3.17.0/topics/mice
apply_naive_MI <- function(data = NULL, num_datasets = NULL, repetitions = NULL, imp_method = NULL) {
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
  
  return (imp)
}

vars_selected_string_to_binary <- function(vars_selected = NULL, var_names = NULL) {
  vars_binary        <- rep(0, times = length(var_names))
  names(vars_binary) <- var_names
  
  for (var in vars_selected) {
    vars_binary[var] <- 1.0
  }
  
  return (vars_binary)
}

vars_selected_MI_means_to_binary <- function(MI_means = NULL) {
  vars_binary        <- rep(0, times = length(MI_means))
  
  for (i in c(1:length(MI_means))) {
    if (MI_means[i] >= 0.5) {
      vars_binary[i] <- 1.0
    }
  }
  
  # include intercept here
  vars_binary        <- c(1.0, vars_binary)
  names(vars_binary) <- c("(Intercept)", names(MI_means))
  
  return (vars_binary)
}

run_naive_MI_simulation <- function(n_scenario = NULL,
                                      n_obs      = NULL,
                                      n_rep      = NULL,
                                      
                                      Z_correlation     = NULL,
                                      Z_subgroups       = NULL,
                                      target_r_sq_X     = NULL,
                                      target_r_sq_Y     = NULL,
                                      causal            = NULL,
                                      
                                      binary_X          = NULL,
                                      binary_Y          = NULL,
                                      binary_Z          = NULL,
                                      
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
                                     
                                     binary_X          = binary_X,
                                     binary_Y          = binary_Y,
                                     binary_Z          = binary_Z,
                                     
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
    
    # apply naive MI
    # imp_method = norm -> MI using Bayesian linear regression
    # see: section "Details" of https://www.rdocumentation.org/packages/mice/versions/3.17.0/topics/mice
    handled_MNAR_dataset <- apply_naive_MI(data         = MNAR_dataset,
                                             num_datasets = 3,
                                             repetitions  = 20,
                                             imp_method   = "norm")
    handled_MCAR_dataset <- apply_naive_MI(data         = MCAR_dataset,
                                             num_datasets = 3,
                                             repetitions  = 20,
                                             imp_method   = "norm")
    
    handled_MNAR_list_of_datasets <- list(complete(handled_MNAR_dataset, action = 1),
                                          complete(handled_MNAR_dataset, action = 2),
                                          complete(handled_MNAR_dataset, action = 3))

    handled_MCAR_list_of_datasets <- list(complete(handled_MCAR_dataset, action = 1),
                                          complete(handled_MCAR_dataset, action = 2),
                                          complete(handled_MCAR_dataset, action = 3))
    
    # record sample sizes before and after missingness handling is applied
    # NB: for naive MI, all imputed datasets will have the same sample size
    sample_size_table["None", "complete_cases"]             <- dim(FULL_dataset)[1]
    sample_size_table["None", "sample_size_after_handling"] <- dim(FULL_dataset)[1]
    sample_size_table["MNAR", "complete_cases"]             <- dim(MNAR_dataset[complete.cases(MNAR_dataset), ])[1]
    sample_size_table["MNAR", "sample_size_after_handling"] <- dim(handled_MNAR_dataset$data)[1]
    sample_size_table["MCAR", "complete_cases"]             <- dim(MCAR_dataset[complete.cases(MCAR_dataset), ])[1]
    sample_size_table["MCAR", "sample_size_after_handling"] <- dim(handled_MCAR_dataset$data)[1]

    # cut up versions of thje data as needed
    X_handled_MNAR_list_of_datasets <- lapply(handled_MNAR_list_of_datasets, subset, select=-c(Y))
    Z_handled_MNAR_list_of_datasets <- lapply(handled_MNAR_list_of_datasets, subset, select=-c(Y, X))
    Y_handled_MNAR_list_of_datasets <- lapply(handled_MNAR_list_of_datasets, subset, select=c(Y))
    X_handled_MNAR_list_of_datasets <- lapply(handled_MNAR_list_of_datasets, subset, select=c(X))
    
    X_handled_MCAR_list_of_datasets <- lapply(handled_MCAR_list_of_datasets, subset, select=-c(Y))
    Z_handled_MCAR_list_of_datasets <- lapply(handled_MCAR_list_of_datasets, subset, select=-c(Y, X))
    Y_handled_MCAR_list_of_datasets <- lapply(handled_MCAR_list_of_datasets, subset, select=c(Y))
    X_handled_MCAR_list_of_datasets <- lapply(handled_MCAR_list_of_datasets, subset, select=c(X))

    # MI LASSO
    cv_lasso_MNAR_list_of_models <- lapply()

    # MI LASSO_x

    # MI LASSO_UNION

    # record cov selection as normal ???
    
  } # end repetitions loop
} # function




