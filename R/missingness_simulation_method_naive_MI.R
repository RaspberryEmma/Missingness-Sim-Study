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
                                             num_datasets = 5,
                                             repetitions  = 20,
                                             imp_method   = "norm")
    handled_MCAR_dataset <- apply_naive_MI(data         = MCAR_dataset,
                                             num_datasets = 5,
                                             repetitions  = 20,
                                             imp_method   = "norm")
    
    handled_MNAR_list_of_datasets <- list(complete(handled_MNAR_dataset, action = 1),
                                          complete(handled_MNAR_dataset, action = 2),
                                          complete(handled_MNAR_dataset, action = 3),
                                          complete(handled_MNAR_dataset, action = 4),
                                          complete(handled_MNAR_dataset, action = 5))

    handled_MCAR_list_of_datasets <- list(complete(handled_MCAR_dataset, action = 1),
                                          complete(handled_MCAR_dataset, action = 2),
                                          complete(handled_MCAR_dataset, action = 3),
                                          complete(handled_MCAR_dataset, action = 4),
                                          complete(handled_MCAR_dataset, action = 5))
    
    # record sample sizes before and after missingness handling is applied
    # NB: for naive MI, all imputed datasets will have the same sample size
    sample_size_table["None", "complete_cases"]             <- dim(FULL_dataset)[1]
    sample_size_table["None", "sample_size_after_handling"] <- dim(FULL_dataset)[1]
    sample_size_table["MNAR", "complete_cases"]             <- dim(MNAR_dataset[complete.cases(MNAR_dataset), ])[1]
    sample_size_table["MNAR", "sample_size_after_handling"] <- dim(handled_MNAR_dataset$data)[1]
    sample_size_table["MCAR", "complete_cases"]             <- dim(MCAR_dataset[complete.cases(MCAR_dataset), ])[1]
    sample_size_table["MCAR", "sample_size_after_handling"] <- dim(handled_MCAR_dataset$data)[1]

    # cut up versions of the data as needed
    X_handled_MNAR_list_of_datasets <- lapply(handled_MNAR_list_of_datasets, subset, select=-c(Y))
    Z_handled_MNAR_list_of_datasets <- lapply(handled_MNAR_list_of_datasets, subset, select=-c(Y, X))
    Y_handled_MNAR_list_of_columns  <- lapply(handled_MNAR_list_of_datasets, subset, select=c(Y))
    X_handled_MNAR_list_of_columns  <- lapply(handled_MNAR_list_of_datasets, subset, select=c(X))
    
    X_handled_MCAR_list_of_datasets <- lapply(handled_MCAR_list_of_datasets, subset, select=-c(Y))
    Z_handled_MCAR_list_of_datasets <- lapply(handled_MCAR_list_of_datasets, subset, select=-c(Y, X))
    Y_handled_MCAR_list_of_columns  <- lapply(handled_MCAR_list_of_datasets, subset, select=c(Y))
    X_handled_MCAR_list_of_columns  <- lapply(handled_MCAR_list_of_datasets, subset, select=c(X))
    
    
    # MI LASSO
    cv_lasso_MNAR_model_1 <- cv.glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[1]]),
                                       y      = data.matrix(Y_handled_MNAR_list_of_columns[[1]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    cv_lasso_MNAR_model_2 <- cv.glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[2]]),
                                       y      = data.matrix(Y_handled_MNAR_list_of_columns[[2]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    cv_lasso_MNAR_model_3 <- cv.glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[3]]),
                                       y      = data.matrix(Y_handled_MNAR_list_of_columns[[3]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    cv_lasso_MNAR_model_4 <- cv.glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[4]]),
                                       y      = data.matrix(Y_handled_MNAR_list_of_columns[[4]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    cv_lasso_MNAR_model_5 <- cv.glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[5]]),
                                       y      = data.matrix(Y_handled_MNAR_list_of_columns[[5]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    
    lambda_1 <- cv_lasso_MNAR_model_1$lambda.min
    lambda_2 <- cv_lasso_MNAR_model_2$lambda.min
    lambda_3 <- cv_lasso_MNAR_model_3$lambda.min
    lambda_4 <- cv_lasso_MNAR_model_4$lambda.min
    lambda_5 <- cv_lasso_MNAR_model_5$lambda.min
    
    lasso_MNAR_model_1 <- glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[1]]),
                                  y      = data.matrix(Y_handled_MNAR_list_of_columns[[1]]),
                                  family = 'gaussian',
                                  alpha  = 1,
                                  lambda = lambda_1)
    lasso_MNAR_model_2 <- glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[2]]),
                                  y      = data.matrix(Y_handled_MNAR_list_of_columns[[2]]),
                                  family = 'gaussian',
                                  alpha  = 1,
                                  lambda = lambda_2)
    lasso_MNAR_model_3 <- glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[3]]),
                                  y      = data.matrix(Y_handled_MNAR_list_of_columns[[3]]),
                                  family = 'gaussian',
                                  alpha  = 1,
                                  lambda = lambda_3)
    lasso_MNAR_model_4 <- glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[4]]),
                                  y      = data.matrix(Y_handled_MNAR_list_of_columns[[4]]),
                                  family = 'gaussian',
                                  alpha  = 1,
                                  lambda = lambda_4)
    lasso_MNAR_model_5 <- glmnet(x      = data.matrix(X_handled_MNAR_list_of_datasets[[5]]),
                                  y      = data.matrix(Y_handled_MNAR_list_of_columns[[5]]),
                                  family = 'gaussian',
                                  alpha  = 1,
                                  lambda = lambda_5)
    
    lasso_MNAR_coefs_1 <- as.vector(lasso_MNAR_model_1$beta)
    lasso_MNAR_coefs_2 <- as.vector(lasso_MNAR_model_2$beta)
    lasso_MNAR_coefs_3 <- as.vector(lasso_MNAR_model_3$beta)
    lasso_MNAR_coefs_4 <- as.vector(lasso_MNAR_model_4$beta)
    lasso_MNAR_coefs_5 <- as.vector(lasso_MNAR_model_5$beta)
    
    names(lasso_MNAR_coefs_1)   <- rownames(lasso_MNAR_model_1$beta)
    names(lasso_MNAR_coefs_2)   <- rownames(lasso_MNAR_model_2$beta)
    names(lasso_MNAR_coefs_3)   <- rownames(lasso_MNAR_model_3$beta)
    names(lasso_MNAR_coefs_4)   <- rownames(lasso_MNAR_model_4$beta)
    names(lasso_MNAR_coefs_5)   <- rownames(lasso_MNAR_model_5$beta)
    
    lasso_MNAR_vars_selected_1  <- union(c('X'), names(lasso_MNAR_coefs_1[lasso_MNAR_coefs_1 != 0.0])) # always select X
    lasso_MNAR_vars_selected_2  <- union(c('X'), names(lasso_MNAR_coefs_2[lasso_MNAR_coefs_2 != 0.0])) # always select X
    lasso_MNAR_vars_selected_3  <- union(c('X'), names(lasso_MNAR_coefs_3[lasso_MNAR_coefs_3 != 0.0])) # always select X
    lasso_MNAR_vars_selected_4  <- union(c('X'), names(lasso_MNAR_coefs_4[lasso_MNAR_coefs_4 != 0.0])) # always select X
    lasso_MNAR_vars_selected_5  <- union(c('X'), names(lasso_MNAR_coefs_5[lasso_MNAR_coefs_5 != 0.0])) # always select X
    
    lasso_MNAR_vars_selected_1  <- lasso_MNAR_vars_selected_1[lasso_MNAR_vars_selected_1 != "(Intercept)"] # exclude intercept
    lasso_MNAR_vars_selected_2  <- lasso_MNAR_vars_selected_2[lasso_MNAR_vars_selected_2 != "(Intercept)"] # exclude intercept
    lasso_MNAR_vars_selected_3  <- lasso_MNAR_vars_selected_3[lasso_MNAR_vars_selected_3 != "(Intercept)"] # exclude intercept
    lasso_MNAR_vars_selected_4  <- lasso_MNAR_vars_selected_4[lasso_MNAR_vars_selected_4 != "(Intercept)"] # exclude intercept
    lasso_MNAR_vars_selected_5  <- lasso_MNAR_vars_selected_5[lasso_MNAR_vars_selected_5 != "(Intercept)"] # exclude intercept
    
    lasso_MNAR_vars_selected_1 <- vars_selected_string_to_binary(vars_selected = lasso_MNAR_vars_selected_1, var_names = var_names_except_Y)
    lasso_MNAR_vars_selected_2 <- vars_selected_string_to_binary(vars_selected = lasso_MNAR_vars_selected_2, var_names = var_names_except_Y)
    lasso_MNAR_vars_selected_3 <- vars_selected_string_to_binary(vars_selected = lasso_MNAR_vars_selected_3, var_names = var_names_except_Y)
    lasso_MNAR_vars_selected_4 <- vars_selected_string_to_binary(vars_selected = lasso_MNAR_vars_selected_4, var_names = var_names_except_Y)
    lasso_MNAR_vars_selected_5 <- vars_selected_string_to_binary(vars_selected = lasso_MNAR_vars_selected_5, var_names = var_names_except_Y)
    
    lasso_MNAR_vars_selected_mean <- colMeans(rbind(lasso_MNAR_vars_selected_1,
                                                    lasso_MNAR_vars_selected_2,
                                                    lasso_MNAR_vars_selected_3,
                                                    lasso_MNAR_vars_selected_4,
                                                    lasso_MNAR_vars_selected_5),
                                              na.rm = TRUE)
    
    lasso_MNAR_vars_selected_more_than_half <- names(lasso_MNAR_vars_selected_mean[lasso_MNAR_vars_selected_mean >= 0.5])
    
    two_step_LASSO_MNAR_model_1 <- glm(make_model_formula(vars_selected = lasso_MNAR_vars_selected_more_than_half),
                                       data   = handled_MNAR_list_of_datasets[[1]],
                                       family = "gaussian")
    two_step_LASSO_MNAR_model_2 <- glm(make_model_formula(vars_selected = lasso_MNAR_vars_selected_more_than_half),
                                       data   = handled_MNAR_list_of_datasets[[2]],
                                       family = "gaussian")
    two_step_LASSO_MNAR_model_3 <- glm(make_model_formula(vars_selected = lasso_MNAR_vars_selected_more_than_half),
                                       data   = handled_MNAR_list_of_datasets[[3]],
                                       family = "gaussian")
    two_step_LASSO_MNAR_model_4 <- glm(make_model_formula(vars_selected = lasso_MNAR_vars_selected_more_than_half),
                                       data   = handled_MNAR_list_of_datasets[[4]],
                                       family = "gaussian")
    two_step_LASSO_MNAR_model_5 <- glm(make_model_formula(vars_selected = lasso_MNAR_vars_selected_more_than_half),
                                       data   = handled_MNAR_list_of_datasets[[5]],
                                       family = "gaussian")
    
    two_step_LASSO_MNAR_causal_effect_estimate <- mean(two_step_LASSO_MNAR_model_1$coefficients['X'],
                                                       two_step_LASSO_MNAR_model_2$coefficients['X'],
                                                       two_step_LASSO_MNAR_model_3$coefficients['X'],
                                                       two_step_LASSO_MNAR_model_4$coefficients['X'],
                                                       two_step_LASSO_MNAR_model_5$coefficients['X'])
    
    cv_lasso_MCAR_model_1 <- cv.glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[1]]),
                                       y      = data.matrix(Y_handled_MCAR_list_of_columns[[1]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    cv_lasso_MCAR_model_2 <- cv.glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[2]]),
                                       y      = data.matrix(Y_handled_MCAR_list_of_columns[[2]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    cv_lasso_MCAR_model_3 <- cv.glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[3]]),
                                       y      = data.matrix(Y_handled_MCAR_list_of_columns[[3]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    cv_lasso_MCAR_model_4 <- cv.glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[4]]),
                                       y      = data.matrix(Y_handled_MCAR_list_of_columns[[4]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    cv_lasso_MCAR_model_5 <- cv.glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[5]]),
                                       y      = data.matrix(Y_handled_MCAR_list_of_columns[[5]]),
                                       family = 'gaussian',
                                       alpha  = 1)
    
    lambda_1 <- cv_lasso_MCAR_model_1$lambda.min
    lambda_2 <- cv_lasso_MCAR_model_2$lambda.min
    lambda_3 <- cv_lasso_MCAR_model_3$lambda.min
    lambda_4 <- cv_lasso_MCAR_model_4$lambda.min
    lambda_5 <- cv_lasso_MCAR_model_5$lambda.min
    
    lasso_MCAR_model_1 <- glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[1]]),
                                 y      = data.matrix(Y_handled_MCAR_list_of_columns[[1]]),
                                 family = 'gaussian',
                                 alpha  = 1,
                                 lambda = lambda_1)
    lasso_MCAR_model_2 <- glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[2]]),
                                 y      = data.matrix(Y_handled_MCAR_list_of_columns[[2]]),
                                 family = 'gaussian',
                                 alpha  = 1,
                                 lambda = lambda_2)
    lasso_MCAR_model_3 <- glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[3]]),
                                 y      = data.matrix(Y_handled_MCAR_list_of_columns[[3]]),
                                 family = 'gaussian',
                                 alpha  = 1,
                                 lambda = lambda_3)
    lasso_MCAR_model_4 <- glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[4]]),
                                 y      = data.matrix(Y_handled_MCAR_list_of_columns[[4]]),
                                 family = 'gaussian',
                                 alpha  = 1,
                                 lambda = lambda_4)
    lasso_MCAR_model_5 <- glmnet(x      = data.matrix(X_handled_MCAR_list_of_datasets[[5]]),
                                 y      = data.matrix(Y_handled_MCAR_list_of_columns[[5]]),
                                 family = 'gaussian',
                                 alpha  = 1,
                                 lambda = lambda_5)
    
    lasso_MCAR_coefs_1 <- as.vector(lasso_MCAR_model_1$beta)
    lasso_MCAR_coefs_2 <- as.vector(lasso_MCAR_model_2$beta)
    lasso_MCAR_coefs_3 <- as.vector(lasso_MCAR_model_3$beta)
    lasso_MCAR_coefs_4 <- as.vector(lasso_MCAR_model_4$beta)
    lasso_MCAR_coefs_5 <- as.vector(lasso_MCAR_model_5$beta)
    
    names(lasso_MCAR_coefs_1)   <- rownames(lasso_MCAR_model_1$beta)
    names(lasso_MCAR_coefs_2)   <- rownames(lasso_MCAR_model_2$beta)
    names(lasso_MCAR_coefs_3)   <- rownames(lasso_MCAR_model_3$beta)
    names(lasso_MCAR_coefs_4)   <- rownames(lasso_MCAR_model_4$beta)
    names(lasso_MCAR_coefs_5)   <- rownames(lasso_MCAR_model_5$beta)
    
    lasso_MCAR_vars_selected_1  <- union(c('X'), names(lasso_MCAR_coefs_1[lasso_MCAR_coefs_1 != 0.0])) # always select X
    lasso_MCAR_vars_selected_2  <- union(c('X'), names(lasso_MCAR_coefs_2[lasso_MCAR_coefs_2 != 0.0])) # always select X
    lasso_MCAR_vars_selected_3  <- union(c('X'), names(lasso_MCAR_coefs_3[lasso_MCAR_coefs_3 != 0.0])) # always select X
    lasso_MCAR_vars_selected_4  <- union(c('X'), names(lasso_MCAR_coefs_4[lasso_MCAR_coefs_4 != 0.0])) # always select X
    lasso_MCAR_vars_selected_5  <- union(c('X'), names(lasso_MCAR_coefs_5[lasso_MCAR_coefs_5 != 0.0])) # always select X
    
    lasso_MCAR_vars_selected_1  <- lasso_MCAR_vars_selected_1[lasso_MCAR_vars_selected_1 != "(Intercept)"] # exclude intercept
    lasso_MCAR_vars_selected_2  <- lasso_MCAR_vars_selected_2[lasso_MCAR_vars_selected_2 != "(Intercept)"] # exclude intercept
    lasso_MCAR_vars_selected_3  <- lasso_MCAR_vars_selected_3[lasso_MCAR_vars_selected_3 != "(Intercept)"] # exclude intercept
    lasso_MCAR_vars_selected_4  <- lasso_MCAR_vars_selected_4[lasso_MCAR_vars_selected_4 != "(Intercept)"] # exclude intercept
    lasso_MCAR_vars_selected_5  <- lasso_MCAR_vars_selected_5[lasso_MCAR_vars_selected_5 != "(Intercept)"] # exclude intercept
    
    lasso_MCAR_vars_selected_1 <- vars_selected_string_to_binary(vars_selected = lasso_MCAR_vars_selected_1, var_names = var_names_except_Y)
    lasso_MCAR_vars_selected_2 <- vars_selected_string_to_binary(vars_selected = lasso_MCAR_vars_selected_2, var_names = var_names_except_Y)
    lasso_MCAR_vars_selected_3 <- vars_selected_string_to_binary(vars_selected = lasso_MCAR_vars_selected_3, var_names = var_names_except_Y)
    lasso_MCAR_vars_selected_4 <- vars_selected_string_to_binary(vars_selected = lasso_MCAR_vars_selected_4, var_names = var_names_except_Y)
    lasso_MCAR_vars_selected_5 <- vars_selected_string_to_binary(vars_selected = lasso_MCAR_vars_selected_5, var_names = var_names_except_Y)
    
    lasso_MCAR_vars_selected_mean <- colMeans(rbind(lasso_MCAR_vars_selected_1,
                                                    lasso_MCAR_vars_selected_2,
                                                    lasso_MCAR_vars_selected_3,
                                                    lasso_MCAR_vars_selected_4,
                                                    lasso_MCAR_vars_selected_5),
                                              na.rm = TRUE)
    
    lasso_MCAR_vars_selected_more_than_half <- names(lasso_MCAR_vars_selected_mean[lasso_MCAR_vars_selected_mean >= 0.5])
    
    two_step_LASSO_MCAR_model_1 <- glm(make_model_formula(vars_selected = lasso_MCAR_vars_selected_more_than_half),
                                       data   = handled_MCAR_list_of_datasets[[1]],
                                       family = "gaussian")
    two_step_LASSO_MCAR_model_2 <- glm(make_model_formula(vars_selected = lasso_MCAR_vars_selected_more_than_half),
                                       data   = handled_MCAR_list_of_datasets[[2]],
                                       family = "gaussian")
    two_step_LASSO_MCAR_model_3 <- glm(make_model_formula(vars_selected = lasso_MCAR_vars_selected_more_than_half),
                                       data   = handled_MCAR_list_of_datasets[[3]],
                                       family = "gaussian")
    two_step_LASSO_MCAR_model_4 <- glm(make_model_formula(vars_selected = lasso_MCAR_vars_selected_more_than_half),
                                       data   = handled_MCAR_list_of_datasets[[4]],
                                       family = "gaussian")
    two_step_LASSO_MCAR_model_5 <- glm(make_model_formula(vars_selected = lasso_MCAR_vars_selected_more_than_half),
                                       data   = handled_MCAR_list_of_datasets[[5]],
                                       family = "gaussian")
    
    two_step_LASSO_MCAR_causal_effect_estimate <- mean(two_step_LASSO_MCAR_model_1$coefficients['X'],
                                                       two_step_LASSO_MCAR_model_2$coefficients['X'],
                                                       two_step_LASSO_MCAR_model_3$coefficients['X'],
                                                       two_step_LASSO_MCAR_model_4$coefficients['X'],
                                                       two_step_LASSO_MCAR_model_5$coefficients['X'])

    # MI LASSO_x
    
    cv_lasso_X_MNAR_model_1 <- cv.glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[1]]),
                                         y      = data.matrix(X_handled_MNAR_list_of_columns[[1]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    cv_lasso_X_MNAR_model_2 <- cv.glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[2]]),
                                         y      = data.matrix(X_handled_MNAR_list_of_columns[[2]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    cv_lasso_X_MNAR_model_3 <- cv.glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[3]]),
                                         y      = data.matrix(X_handled_MNAR_list_of_columns[[3]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    cv_lasso_X_MNAR_model_4 <- cv.glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[4]]),
                                         y      = data.matrix(X_handled_MNAR_list_of_columns[[4]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    cv_lasso_X_MNAR_model_5 <- cv.glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[5]]),
                                         y      = data.matrix(X_handled_MNAR_list_of_columns[[5]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    
    lambda_1 <- cv_lasso_X_MNAR_model_1$lambda.min
    lambda_2 <- cv_lasso_X_MNAR_model_2$lambda.min
    lambda_3 <- cv_lasso_X_MNAR_model_3$lambda.min
    lambda_4 <- cv_lasso_X_MNAR_model_4$lambda.min
    lambda_5 <- cv_lasso_X_MNAR_model_5$lambda.min
    
    lasso_X_MNAR_model_1 <- glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[1]]),
                                   y      = data.matrix(X_handled_MNAR_list_of_columns[[1]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_1)
    lasso_X_MNAR_model_2 <- glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[2]]),
                                   y      = data.matrix(X_handled_MNAR_list_of_columns[[2]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_2)
    lasso_X_MNAR_model_3 <- glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[3]]),
                                   y      = data.matrix(X_handled_MNAR_list_of_columns[[3]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_3)
    lasso_X_MNAR_model_4 <- glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[4]]),
                                   y      = data.matrix(X_handled_MNAR_list_of_columns[[4]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_4)
    lasso_X_MNAR_model_5 <- glmnet(x      = data.matrix(Z_handled_MNAR_list_of_datasets[[5]]),
                                   y      = data.matrix(X_handled_MNAR_list_of_columns[[5]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_5)
    
    lasso_X_MNAR_coefs_1 <- as.vector(lasso_X_MNAR_model_1$beta)
    lasso_X_MNAR_coefs_2 <- as.vector(lasso_X_MNAR_model_2$beta)
    lasso_X_MNAR_coefs_3 <- as.vector(lasso_X_MNAR_model_3$beta)
    lasso_X_MNAR_coefs_4 <- as.vector(lasso_X_MNAR_model_4$beta)
    lasso_X_MNAR_coefs_5 <- as.vector(lasso_X_MNAR_model_5$beta)
    
    names(lasso_X_MNAR_coefs_1)   <- rownames(lasso_X_MNAR_model_1$beta)
    names(lasso_X_MNAR_coefs_2)   <- rownames(lasso_X_MNAR_model_2$beta)
    names(lasso_X_MNAR_coefs_3)   <- rownames(lasso_X_MNAR_model_3$beta)
    names(lasso_X_MNAR_coefs_4)   <- rownames(lasso_X_MNAR_model_4$beta)
    names(lasso_X_MNAR_coefs_5)   <- rownames(lasso_X_MNAR_model_5$beta)
    
    lasso_X_MNAR_vars_selected_1  <- union(c('X'), names(lasso_X_MNAR_coefs_1[lasso_X_MNAR_coefs_1 != 0.0])) # always select X
    lasso_X_MNAR_vars_selected_2  <- union(c('X'), names(lasso_X_MNAR_coefs_2[lasso_X_MNAR_coefs_2 != 0.0])) # always select X
    lasso_X_MNAR_vars_selected_3  <- union(c('X'), names(lasso_X_MNAR_coefs_3[lasso_X_MNAR_coefs_3 != 0.0])) # always select X
    lasso_X_MNAR_vars_selected_4  <- union(c('X'), names(lasso_X_MNAR_coefs_4[lasso_X_MNAR_coefs_4 != 0.0])) # always select X
    lasso_X_MNAR_vars_selected_5  <- union(c('X'), names(lasso_X_MNAR_coefs_5[lasso_X_MNAR_coefs_5 != 0.0])) # always select X
    
    lasso_X_MNAR_vars_selected_1  <- lasso_X_MNAR_vars_selected_1[lasso_X_MNAR_vars_selected_1 != "(Intercept)"] # exclude intercept
    lasso_X_MNAR_vars_selected_2  <- lasso_X_MNAR_vars_selected_2[lasso_X_MNAR_vars_selected_2 != "(Intercept)"] # exclude intercept
    lasso_X_MNAR_vars_selected_3  <- lasso_X_MNAR_vars_selected_3[lasso_X_MNAR_vars_selected_3 != "(Intercept)"] # exclude intercept
    lasso_X_MNAR_vars_selected_4  <- lasso_X_MNAR_vars_selected_4[lasso_X_MNAR_vars_selected_4 != "(Intercept)"] # exclude intercept
    lasso_X_MNAR_vars_selected_5  <- lasso_X_MNAR_vars_selected_5[lasso_X_MNAR_vars_selected_5 != "(Intercept)"] # exclude intercept
    
    lasso_X_MNAR_vars_selected_1 <- vars_selected_string_to_binary(vars_selected = lasso_X_MNAR_vars_selected_1, var_names = var_names_except_Y)
    lasso_X_MNAR_vars_selected_2 <- vars_selected_string_to_binary(vars_selected = lasso_X_MNAR_vars_selected_2, var_names = var_names_except_Y)
    lasso_X_MNAR_vars_selected_3 <- vars_selected_string_to_binary(vars_selected = lasso_X_MNAR_vars_selected_3, var_names = var_names_except_Y)
    lasso_X_MNAR_vars_selected_4 <- vars_selected_string_to_binary(vars_selected = lasso_X_MNAR_vars_selected_4, var_names = var_names_except_Y)
    lasso_X_MNAR_vars_selected_5 <- vars_selected_string_to_binary(vars_selected = lasso_X_MNAR_vars_selected_5, var_names = var_names_except_Y)
    
    lasso_X_MNAR_vars_selected_mean <- colMeans(rbind(lasso_X_MNAR_vars_selected_1,
                                                      lasso_X_MNAR_vars_selected_2,
                                                      lasso_X_MNAR_vars_selected_3,
                                                      lasso_X_MNAR_vars_selected_4,
                                                      lasso_X_MNAR_vars_selected_5),
                                                na.rm = TRUE)
    
    lasso_X_MNAR_vars_selected_more_than_half <- names(lasso_X_MNAR_vars_selected_mean[lasso_X_MNAR_vars_selected_mean >= 0.5])
    
    two_step_LASSO_X_MNAR_model_1 <- glm(make_model_formula(vars_selected = lasso_X_MNAR_vars_selected_more_than_half),
                                         data   = handled_MNAR_list_of_datasets[[1]],
                                         family = "gaussian")
    two_step_LASSO_X_MNAR_model_2 <- glm(make_model_formula(vars_selected = lasso_X_MNAR_vars_selected_more_than_half),
                                         data   = handled_MNAR_list_of_datasets[[2]],
                                         family = "gaussian")
    two_step_LASSO_X_MNAR_model_3 <- glm(make_model_formula(vars_selected = lasso_X_MNAR_vars_selected_more_than_half),
                                         data   = handled_MNAR_list_of_datasets[[3]],
                                         family = "gaussian")
    two_step_LASSO_X_MNAR_model_4 <- glm(make_model_formula(vars_selected = lasso_X_MNAR_vars_selected_more_than_half),
                                         data   = handled_MNAR_list_of_datasets[[4]],
                                         family = "gaussian")
    two_step_LASSO_X_MNAR_model_5 <- glm(make_model_formula(vars_selected = lasso_X_MNAR_vars_selected_more_than_half),
                                         data   = handled_MNAR_list_of_datasets[[5]],
                                         family = "gaussian")
    
    two_step_LASSO_X_MNAR_causal_effect_estimate <- mean(two_step_LASSO_X_MNAR_model_1$coefficients['X'],
                                                         two_step_LASSO_X_MNAR_model_2$coefficients['X'],
                                                         two_step_LASSO_X_MNAR_model_3$coefficients['X'],
                                                         two_step_LASSO_X_MNAR_model_4$coefficients['X'],
                                                         two_step_LASSO_X_MNAR_model_5$coefficients['X'])
    
    cv_lasso_X_MCAR_model_1 <- cv.glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[1]]),
                                         y      = data.matrix(X_handled_MCAR_list_of_columns[[1]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    cv_lasso_X_MCAR_model_2 <- cv.glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[2]]),
                                         y      = data.matrix(X_handled_MCAR_list_of_columns[[2]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    cv_lasso_X_MCAR_model_3 <- cv.glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[3]]),
                                         y      = data.matrix(X_handled_MCAR_list_of_columns[[3]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    cv_lasso_X_MCAR_model_4 <- cv.glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[4]]),
                                         y      = data.matrix(X_handled_MCAR_list_of_columns[[4]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    cv_lasso_X_MCAR_model_5 <- cv.glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[5]]),
                                         y      = data.matrix(X_handled_MCAR_list_of_columns[[5]]),
                                         family = 'gaussian',
                                         alpha  = 1)
    
    lambda_1 <- cv_lasso_X_MCAR_model_1$lambda.min
    lambda_2 <- cv_lasso_X_MCAR_model_2$lambda.min
    lambda_3 <- cv_lasso_X_MCAR_model_3$lambda.min
    lambda_4 <- cv_lasso_X_MCAR_model_4$lambda.min
    lambda_5 <- cv_lasso_X_MCAR_model_5$lambda.min
    
    lasso_X_MCAR_model_1 <- glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[1]]),
                                   y      = data.matrix(X_handled_MCAR_list_of_columns[[1]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_1)
    lasso_X_MCAR_model_2 <- glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[2]]),
                                   y      = data.matrix(X_handled_MCAR_list_of_columns[[2]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_2)
    lasso_X_MCAR_model_3 <- glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[3]]),
                                   y      = data.matrix(X_handled_MCAR_list_of_columns[[3]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_3)
    lasso_X_MCAR_model_4 <- glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[4]]),
                                   y      = data.matrix(X_handled_MCAR_list_of_columns[[4]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_4)
    lasso_X_MCAR_model_5 <- glmnet(x      = data.matrix(Z_handled_MCAR_list_of_datasets[[5]]),
                                   y      = data.matrix(X_handled_MCAR_list_of_columns[[5]]),
                                   family = 'gaussian',
                                   alpha  = 1,
                                   lambda = lambda_5)
    
    lasso_X_MCAR_coefs_1 <- as.vector(lasso_X_MCAR_model_1$beta)
    lasso_X_MCAR_coefs_2 <- as.vector(lasso_X_MCAR_model_2$beta)
    lasso_X_MCAR_coefs_3 <- as.vector(lasso_X_MCAR_model_3$beta)
    lasso_X_MCAR_coefs_4 <- as.vector(lasso_X_MCAR_model_4$beta)
    lasso_X_MCAR_coefs_5 <- as.vector(lasso_X_MCAR_model_5$beta)
    
    names(lasso_X_MCAR_coefs_1)   <- rownames(lasso_X_MCAR_model_1$beta)
    names(lasso_X_MCAR_coefs_2)   <- rownames(lasso_X_MCAR_model_2$beta)
    names(lasso_X_MCAR_coefs_3)   <- rownames(lasso_X_MCAR_model_3$beta)
    names(lasso_X_MCAR_coefs_4)   <- rownames(lasso_X_MCAR_model_4$beta)
    names(lasso_X_MCAR_coefs_5)   <- rownames(lasso_X_MCAR_model_5$beta)
    
    lasso_X_MCAR_vars_selected_1  <- union(c('X'), names(lasso_X_MCAR_coefs_1[lasso_X_MCAR_coefs_1 != 0.0])) # always select X
    lasso_X_MCAR_vars_selected_2  <- union(c('X'), names(lasso_X_MCAR_coefs_2[lasso_X_MCAR_coefs_2 != 0.0])) # always select X
    lasso_X_MCAR_vars_selected_3  <- union(c('X'), names(lasso_X_MCAR_coefs_3[lasso_X_MCAR_coefs_3 != 0.0])) # always select X
    lasso_X_MCAR_vars_selected_4  <- union(c('X'), names(lasso_X_MCAR_coefs_4[lasso_X_MCAR_coefs_4 != 0.0])) # always select X
    lasso_X_MCAR_vars_selected_5  <- union(c('X'), names(lasso_X_MCAR_coefs_5[lasso_X_MCAR_coefs_5 != 0.0])) # always select X
    
    lasso_X_MCAR_vars_selected_1  <- lasso_X_MCAR_vars_selected_1[lasso_X_MCAR_vars_selected_1 != "(Intercept)"] # exclude intercept
    lasso_X_MCAR_vars_selected_2  <- lasso_X_MCAR_vars_selected_2[lasso_X_MCAR_vars_selected_2 != "(Intercept)"] # exclude intercept
    lasso_X_MCAR_vars_selected_3  <- lasso_X_MCAR_vars_selected_3[lasso_X_MCAR_vars_selected_3 != "(Intercept)"] # exclude intercept
    lasso_X_MCAR_vars_selected_4  <- lasso_X_MCAR_vars_selected_4[lasso_X_MCAR_vars_selected_4 != "(Intercept)"] # exclude intercept
    lasso_X_MCAR_vars_selected_5  <- lasso_X_MCAR_vars_selected_5[lasso_X_MCAR_vars_selected_5 != "(Intercept)"] # exclude intercept
    
    lasso_X_MCAR_vars_selected_1 <- vars_selected_string_to_binary(vars_selected = lasso_X_MCAR_vars_selected_1, var_names = var_names_except_Y)
    lasso_X_MCAR_vars_selected_2 <- vars_selected_string_to_binary(vars_selected = lasso_X_MCAR_vars_selected_2, var_names = var_names_except_Y)
    lasso_X_MCAR_vars_selected_3 <- vars_selected_string_to_binary(vars_selected = lasso_X_MCAR_vars_selected_3, var_names = var_names_except_Y)
    lasso_X_MCAR_vars_selected_4 <- vars_selected_string_to_binary(vars_selected = lasso_X_MCAR_vars_selected_4, var_names = var_names_except_Y)
    lasso_X_MCAR_vars_selected_5 <- vars_selected_string_to_binary(vars_selected = lasso_X_MCAR_vars_selected_5, var_names = var_names_except_Y)
    
    lasso_X_MCAR_vars_selected_mean <- colMeans(rbind(lasso_X_MCAR_vars_selected_1,
                                                      lasso_X_MCAR_vars_selected_2,
                                                      lasso_X_MCAR_vars_selected_3,
                                                      lasso_X_MCAR_vars_selected_4,
                                                      lasso_X_MCAR_vars_selected_5),
                                                na.rm = TRUE)
    
    lasso_X_MCAR_vars_selected_more_than_half <- names(lasso_X_MCAR_vars_selected_mean[lasso_X_MCAR_vars_selected_mean >= 0.5])
    
    two_step_LASSO_X_MCAR_model_1 <- glm(make_model_formula(vars_selected = lasso_X_MCAR_vars_selected_more_than_half),
                                         data   = handled_MCAR_list_of_datasets[[1]],
                                         family = "gaussian")
    two_step_LASSO_X_MCAR_model_2 <- glm(make_model_formula(vars_selected = lasso_X_MCAR_vars_selected_more_than_half),
                                         data   = handled_MCAR_list_of_datasets[[2]],
                                         family = "gaussian")
    two_step_LASSO_X_MCAR_model_3 <- glm(make_model_formula(vars_selected = lasso_X_MCAR_vars_selected_more_than_half),
                                         data   = handled_MCAR_list_of_datasets[[3]],
                                         family = "gaussian")
    two_step_LASSO_X_MCAR_model_4 <- glm(make_model_formula(vars_selected = lasso_X_MCAR_vars_selected_more_than_half),
                                         data   = handled_MCAR_list_of_datasets[[4]],
                                         family = "gaussian")
    two_step_LASSO_X_MCAR_model_5 <- glm(make_model_formula(vars_selected = lasso_X_MCAR_vars_selected_more_than_half),
                                         data   = handled_MCAR_list_of_datasets[[5]],
                                         family = "gaussian")
    
    two_step_LASSO_X_MCAR_causal_effect_estimate <- mean(two_step_LASSO_X_MCAR_model_1$coefficients['X'],
                                                         two_step_LASSO_X_MCAR_model_2$coefficients['X'],
                                                         two_step_LASSO_X_MCAR_model_3$coefficients['X'],
                                                         two_step_LASSO_X_MCAR_model_4$coefficients['X'],
                                                         two_step_LASSO_X_MCAR_model_5$coefficients['X'])

    # MI LASSO_UNION
    
    lasso_union_MNAR_vars_selected_1 <- lasso_MNAR_vars_selected_1 + lasso_X_MNAR_vars_selected_1
    lasso_union_MNAR_vars_selected_1[lasso_union_MNAR_vars_selected_1 == 2] <- 1
    
    lasso_union_MNAR_vars_selected_2 <- lasso_MNAR_vars_selected_2 + lasso_X_MNAR_vars_selected_2
    lasso_union_MNAR_vars_selected_2[lasso_union_MNAR_vars_selected_2 == 2] <- 1
    
    lasso_union_MNAR_vars_selected_3 <- lasso_MNAR_vars_selected_3 + lasso_X_MNAR_vars_selected_3
    lasso_union_MNAR_vars_selected_3[lasso_union_MNAR_vars_selected_3 == 2] <- 1
    
    lasso_union_MNAR_vars_selected_4 <- lasso_MNAR_vars_selected_4 + lasso_X_MNAR_vars_selected_4
    lasso_union_MNAR_vars_selected_4[lasso_union_MNAR_vars_selected_4 == 2] <- 1
    
    lasso_union_MNAR_vars_selected_5 <- lasso_MNAR_vars_selected_5 + lasso_X_MNAR_vars_selected_5
    lasso_union_MNAR_vars_selected_5[lasso_union_MNAR_vars_selected_5 == 2] <- 1
    
    lasso_union_MNAR_vars_selected_mean <- colMeans(rbind(lasso_union_MNAR_vars_selected_1,
                                                          lasso_union_MNAR_vars_selected_2,
                                                          lasso_union_MNAR_vars_selected_3,
                                                          lasso_union_MNAR_vars_selected_4,
                                                          lasso_union_MNAR_vars_selected_5),
                                                    na.rm = TRUE)
    
    lasso_union_MNAR_vars_selected_more_than_half <- names(lasso_union_MNAR_vars_selected_mean[lasso_union_MNAR_vars_selected_mean >= 0.5])
    
    two_step_LASSO_union_MNAR_model_1 <- glm(make_model_formula(vars_selected = lasso_union_MNAR_vars_selected_more_than_half),
                                             data   = handled_MNAR_list_of_datasets[[1]],
                                             family = "gaussian")
    two_step_LASSO_union_MNAR_model_2 <- glm(make_model_formula(vars_selected = lasso_union_MNAR_vars_selected_more_than_half),
                                             data   = handled_MNAR_list_of_datasets[[2]],
                                             family = "gaussian")
    two_step_LASSO_union_MNAR_model_3 <- glm(make_model_formula(vars_selected = lasso_union_MNAR_vars_selected_more_than_half),
                                             data   = handled_MNAR_list_of_datasets[[3]],
                                             family = "gaussian")
    two_step_LASSO_union_MNAR_model_4 <- glm(make_model_formula(vars_selected = lasso_union_MNAR_vars_selected_more_than_half),
                                             data   = handled_MNAR_list_of_datasets[[4]],
                                             family = "gaussian")
    two_step_LASSO_union_MNAR_model_5 <- glm(make_model_formula(vars_selected = lasso_union_MNAR_vars_selected_more_than_half),
                                             data   = handled_MNAR_list_of_datasets[[5]],
                                             family = "gaussian")
    
    two_step_LASSO_union_MNAR_causal_effect_estimate <- mean(two_step_LASSO_union_MNAR_model_1$coefficients['X'],
                                                             two_step_LASSO_union_MNAR_model_2$coefficients['X'],
                                                             two_step_LASSO_union_MNAR_model_3$coefficients['X'],
                                                             two_step_LASSO_union_MNAR_model_4$coefficients['X'],
                                                             two_step_LASSO_union_MNAR_model_5$coefficients['X'])
    
    lasso_union_MCAR_vars_selected_1 <- lasso_MCAR_vars_selected_1 + lasso_X_MCAR_vars_selected_1
    lasso_union_MCAR_vars_selected_1[lasso_union_MCAR_vars_selected_1 == 2] <- 1
    
    lasso_union_MCAR_vars_selected_2 <- lasso_MCAR_vars_selected_2 + lasso_X_MCAR_vars_selected_2
    lasso_union_MCAR_vars_selected_2[lasso_union_MCAR_vars_selected_2 == 2] <- 1
    
    lasso_union_MCAR_vars_selected_3 <- lasso_MCAR_vars_selected_3 + lasso_X_MCAR_vars_selected_3
    lasso_union_MCAR_vars_selected_3[lasso_union_MCAR_vars_selected_3 == 2] <- 1
    
    lasso_union_MCAR_vars_selected_4 <- lasso_MCAR_vars_selected_4 + lasso_X_MCAR_vars_selected_4
    lasso_union_MCAR_vars_selected_4[lasso_union_MCAR_vars_selected_4 == 2] <- 1
    
    lasso_union_MCAR_vars_selected_5 <- lasso_MCAR_vars_selected_5 + lasso_X_MCAR_vars_selected_5
    lasso_union_MCAR_vars_selected_5[lasso_union_MCAR_vars_selected_5 == 2] <- 1
    
    lasso_union_MCAR_vars_selected_mean <- colMeans(rbind(lasso_union_MCAR_vars_selected_1,
                                                          lasso_union_MCAR_vars_selected_2,
                                                          lasso_union_MCAR_vars_selected_3,
                                                          lasso_union_MCAR_vars_selected_4,
                                                          lasso_union_MCAR_vars_selected_5),
                                                    na.rm = TRUE)
    
    lasso_union_MCAR_vars_selected_more_than_half <- names(lasso_union_MCAR_vars_selected_mean[lasso_union_MCAR_vars_selected_mean >= 0.5])
    
    two_step_LASSO_union_MCAR_model_1 <- glm(make_model_formula(vars_selected = lasso_union_MCAR_vars_selected_more_than_half),
                                             data   = handled_MCAR_list_of_datasets[[1]],
                                             family = "gaussian")
    two_step_LASSO_union_MCAR_model_2 <- glm(make_model_formula(vars_selected = lasso_union_MCAR_vars_selected_more_than_half),
                                             data   = handled_MCAR_list_of_datasets[[2]],
                                             family = "gaussian")
    two_step_LASSO_union_MCAR_model_3 <- glm(make_model_formula(vars_selected = lasso_union_MCAR_vars_selected_more_than_half),
                                             data   = handled_MCAR_list_of_datasets[[3]],
                                             family = "gaussian")
    two_step_LASSO_union_MCAR_model_4 <- glm(make_model_formula(vars_selected = lasso_union_MCAR_vars_selected_more_than_half),
                                             data   = handled_MCAR_list_of_datasets[[4]],
                                             family = "gaussian")
    two_step_LASSO_union_MCAR_model_5 <- glm(make_model_formula(vars_selected = lasso_union_MCAR_vars_selected_more_than_half),
                                             data   = handled_MCAR_list_of_datasets[[5]],
                                             family = "gaussian")
    
    two_step_LASSO_union_MCAR_causal_effect_estimate <- mean(two_step_LASSO_union_MCAR_model_1$coefficients['X'],
                                                             two_step_LASSO_union_MCAR_model_2$coefficients['X'],
                                                             two_step_LASSO_union_MCAR_model_3$coefficients['X'],
                                                             two_step_LASSO_union_MCAR_model_4$coefficients['X'],
                                                             two_step_LASSO_union_MCAR_model_5$coefficients['X'])
    
    # ----- Record covariate selection -----
    
    MNAR_cov_selection["two_step_lasso", , repetition]       <- vars_selected_string_to_binary(vars_selected = lasso_MNAR_vars_selected_more_than_half,       var_names = var_names_except_Y_with_intercept)
    MNAR_cov_selection["two_step_lasso_X", , repetition]     <- vars_selected_string_to_binary(vars_selected = lasso_X_MNAR_vars_selected_more_than_half,     var_names = var_names_except_Y_with_intercept)
    MNAR_cov_selection["two_step_lasso_union", , repetition] <- vars_selected_string_to_binary(vars_selected = lasso_union_MNAR_vars_selected_more_than_half, var_names = var_names_except_Y_with_intercept)
    
    MCAR_cov_selection["two_step_lasso", , repetition]       <- vars_selected_string_to_binary(vars_selected = lasso_MCAR_vars_selected_more_than_half,       var_names = var_names_except_Y_with_intercept)
    MCAR_cov_selection["two_step_lasso_X", , repetition]     <- vars_selected_string_to_binary(vars_selected = lasso_X_MCAR_vars_selected_more_than_half,     var_names = var_names_except_Y_with_intercept)
    MCAR_cov_selection["two_step_lasso_union", , repetition] <- vars_selected_string_to_binary(vars_selected = lasso_union_MCAR_vars_selected_more_than_half, var_names = var_names_except_Y_with_intercept)
    
    # ----- Record results -----
    
    MNAR_results["two_step_lasso", "causal_true_value", repetition]      <- causal
    MNAR_results["two_step_lasso", "causal_estimate", repetition]        <- unname(two_step_LASSO_MNAR_causal_effect_estimate)
    MNAR_results["two_step_lasso", "causal_bias", repetition]            <- (unname(two_step_LASSO_MNAR_causal_effect_estimate) - causal)
    MNAR_results["two_step_lasso", "causal_bias_proportion", repetition] <- ((unname(two_step_LASSO_MNAR_causal_effect_estimate) - causal)/causal)
    MNAR_results["two_step_lasso", "causal_coverage", repetition]        <- NaN
    MNAR_results["two_step_lasso", "open_paths", repetition]             <- num_total_conf
    MNAR_results["two_step_lasso", "blocked_paths", repetition]          <- length(lasso_MNAR_vars_selected_more_than_half[lasso_MNAR_vars_selected_more_than_half != 'X'])
    MNAR_results["two_step_lasso", "proportion_paths", repetition]       <- length(lasso_MNAR_vars_selected_more_than_half[lasso_MNAR_vars_selected_more_than_half != 'X']) / num_total_conf
    MNAR_results["two_step_lasso", "empirical_SE", repetition]           <- NaN
    MNAR_results["two_step_lasso", "model_SE", repetition]               <- NaN
    
    MNAR_results["two_step_lasso_X", "causal_true_value", repetition]      <- causal
    MNAR_results["two_step_lasso_X", "causal_estimate", repetition]        <- unname(two_step_LASSO_X_MNAR_causal_effect_estimate)
    MNAR_results["two_step_lasso_X", "causal_bias", repetition]            <- (unname(two_step_LASSO_X_MNAR_causal_effect_estimate) - causal)
    MNAR_results["two_step_lasso_X", "causal_bias_proportion", repetition] <- ((unname(two_step_LASSO_X_MNAR_causal_effect_estimate) - causal)/causal)
    MNAR_results["two_step_lasso_X", "causal_coverage", repetition]        <- NaN
    MNAR_results["two_step_lasso_X", "open_paths", repetition]             <- num_total_conf
    MNAR_results["two_step_lasso_X", "blocked_paths", repetition]          <- length(lasso_X_MNAR_vars_selected_more_than_half[lasso_X_MNAR_vars_selected_more_than_half != 'X'])
    MNAR_results["two_step_lasso_X", "proportion_paths", repetition]       <- length(lasso_X_MNAR_vars_selected_more_than_half[lasso_X_MNAR_vars_selected_more_than_half != 'X']) / num_total_conf
    MNAR_results["two_step_lasso_X", "empirical_SE", repetition]           <- NaN
    MNAR_results["two_step_lasso_X", "model_SE", repetition]               <- NaN
    
    MNAR_results["two_step_lasso_union", "causal_true_value", repetition]      <- causal
    MNAR_results["two_step_lasso_union", "causal_estimate", repetition]        <- unname(two_step_LASSO_union_MNAR_causal_effect_estimate)
    MNAR_results["two_step_lasso_union", "causal_bias", repetition]            <- (unname(two_step_LASSO_union_MNAR_causal_effect_estimate) - causal)
    MNAR_results["two_step_lasso_union", "causal_bias_proportion", repetition] <- ((unname(two_step_LASSO_union_MNAR_causal_effect_estimate) - causal)/causal)
    MNAR_results["two_step_lasso_union", "causal_coverage", repetition]        <- NaN
    MNAR_results["two_step_lasso_union", "open_paths", repetition]             <- num_total_conf
    MNAR_results["two_step_lasso_union", "blocked_paths", repetition]          <- length(lasso_union_MNAR_vars_selected_more_than_half[lasso_union_MNAR_vars_selected_more_than_half != 'X'])
    MNAR_results["two_step_lasso_union", "proportion_paths", repetition]       <- length(lasso_union_MNAR_vars_selected_more_than_half[lasso_union_MNAR_vars_selected_more_than_half != 'X']) / num_total_conf
    MNAR_results["two_step_lasso_union", "empirical_SE", repetition]           <- NaN
    MNAR_results["two_step_lasso_union", "model_SE", repetition]               <- NaN
    
    MCAR_results["two_step_lasso", "causal_true_value", repetition]      <- causal
    MCAR_results["two_step_lasso", "causal_estimate", repetition]        <- unname(two_step_LASSO_MCAR_causal_effect_estimate)
    MCAR_results["two_step_lasso", "causal_bias", repetition]            <- (unname(two_step_LASSO_MCAR_causal_effect_estimate) - causal)
    MCAR_results["two_step_lasso", "causal_bias_proportion", repetition] <- ((unname(two_step_LASSO_MCAR_causal_effect_estimate) - causal)/causal)
    MCAR_results["two_step_lasso", "causal_coverage", repetition]        <- NaN
    MCAR_results["two_step_lasso", "open_paths", repetition]             <- num_total_conf
    MCAR_results["two_step_lasso", "blocked_paths", repetition]          <- length(lasso_MCAR_vars_selected_more_than_half[lasso_MCAR_vars_selected_more_than_half != 'X'])
    MCAR_results["two_step_lasso", "proportion_paths", repetition]       <- length(lasso_MCAR_vars_selected_more_than_half[lasso_MCAR_vars_selected_more_than_half != 'X']) / num_total_conf
    MCAR_results["two_step_lasso", "empirical_SE", repetition]           <- NaN
    MCAR_results["two_step_lasso", "model_SE", repetition]               <- NaN
    
    MCAR_results["two_step_lasso_X", "causal_true_value", repetition]      <- causal
    MCAR_results["two_step_lasso_X", "causal_estimate", repetition]        <- unname(two_step_LASSO_X_MCAR_causal_effect_estimate)
    MCAR_results["two_step_lasso_X", "causal_bias", repetition]            <- (unname(two_step_LASSO_X_MCAR_causal_effect_estimate) - causal)
    MCAR_results["two_step_lasso_X", "causal_bias_proportion", repetition] <- ((unname(two_step_LASSO_X_MCAR_causal_effect_estimate) - causal)/causal)
    MCAR_results["two_step_lasso_X", "causal_coverage", repetition]        <- NaN
    MCAR_results["two_step_lasso_X", "open_paths", repetition]             <- num_total_conf
    MCAR_results["two_step_lasso_X", "blocked_paths", repetition]          <- length(lasso_X_MCAR_vars_selected_more_than_half[lasso_X_MCAR_vars_selected_more_than_half != 'X'])
    MCAR_results["two_step_lasso_X", "proportion_paths", repetition]       <- length(lasso_X_MCAR_vars_selected_more_than_half[lasso_X_MCAR_vars_selected_more_than_half != 'X']) / num_total_conf
    MCAR_results["two_step_lasso_X", "empirical_SE", repetition]           <- NaN
    MCAR_results["two_step_lasso_X", "model_SE", repetition]               <- NaN
    
    MCAR_results["two_step_lasso_union", "causal_true_value", repetition]      <- causal
    MCAR_results["two_step_lasso_union", "causal_estimate", repetition]        <- unname(two_step_LASSO_union_MCAR_causal_effect_estimate)
    MCAR_results["two_step_lasso_union", "causal_bias", repetition]            <- (unname(two_step_LASSO_union_MCAR_causal_effect_estimate) - causal)
    MCAR_results["two_step_lasso_union", "causal_bias_proportion", repetition] <- ((unname(two_step_LASSO_union_MCAR_causal_effect_estimate) - causal)/causal)
    MCAR_results["two_step_lasso_union", "causal_coverage", repetition]        <- NaN
    MCAR_results["two_step_lasso_union", "open_paths", repetition]             <- num_total_conf
    MCAR_results["two_step_lasso_union", "blocked_paths", repetition]          <- length(lasso_union_MCAR_vars_selected_more_than_half[lasso_union_MCAR_vars_selected_more_than_half != 'X'])
    MCAR_results["two_step_lasso_union", "proportion_paths", repetition]       <- length(lasso_union_MCAR_vars_selected_more_than_half[lasso_union_MCAR_vars_selected_more_than_half != 'X']) / num_total_conf
    MCAR_results["two_step_lasso_union", "empirical_SE", repetition]           <- NaN
    MCAR_results["two_step_lasso_union", "model_SE", repetition]               <- NaN
    
  } # end repetitions loop
  
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
               handled_MNAR_list_of_datasets[[1]],
               handled_MCAR_list_of_datasets[[1]]
  ))
} # function


