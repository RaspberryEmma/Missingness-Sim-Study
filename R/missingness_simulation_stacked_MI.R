# ****************************************
# Missingness Simulation Study
#
# Simulation procedure for confounder-handling
# with missingness considerations
#
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 21/11/2025
# ****************************************



# ----- Helper functions -----

# Helper function for recording coefficients fitted
# Accounts for the idea that some modelling methods will exclude variables
# i.e:  Variables correct order with NaNs filling in excluded values
fill_in_blanks <- function(coefs = NULL, labels = NULL) {
  for (label in labels) {
    # if a given variable doesn't exist, create as NaN
    if (!(label %in% names(coefs))) {
      coefs[label] <- NaN
    }
  }
  
  # assert variable ordering
  coefs <- coefs[order(factor(names(coefs), levels = labels))]
  
  return (coefs)
}

# Helper function for recording covariate selection
# i.e:  Binary indicators per variable, 1 = selected, 0 = not
fill_in_cov_selection <- function(coefs = NULL, labels = NULL) {
  for (label in labels) {
    if (!(label %in% names(coefs))) {
      coefs[label] <- 0
    } else {
      coefs[label] <- 1
    }
  }
  
  # assert variable ordering
  coefs <- coefs[order(factor(names(coefs), levels = labels))]
  
  return (coefs)
}

# Create a string corresponding to a regression of Y on a given set of covariates
make_model_formula <- function(vars_selected = NULL) {
  
  if (length(vars_selected) > 0) {
    # first item does not require a '+' sign
    formula_string <- paste("Y ~ ", vars_selected[1], sep = "")
    
    for (var in vars_selected[-1]) {
      formula_string <- paste(formula_string, " + ", var, sep = "")
    }
  }
  else {
    # if no variables are selected, use constant term only
    formula_string <- "Y ~ 0"
  }
  
  return(formula_string)
}


# Create a string corresponding to a regression of X on a given set of covariates
make_X_model_formula <- function(vars_selected = NULL) {
  # remove X
  if ('X' %in% vars_selected) {
    vars_selected <- vars_selected[! vars_selected%in% c('X')]
  }
  
  if (length(vars_selected) > 0) {
    # first item does not require a '+' sign
    formula_string <- paste("X ~ ", vars_selected[1], sep = "")
    
    for (var in vars_selected[-1]) {
      formula_string <- paste(formula_string, " + ", var, sep = "")
    }
  }
  else {
    # if no variables are selected, use constant term only
    formula_string <- "X ~ 0"
  }
  
  return(formula_string)
}

# Inverse-logit transform function
# Used for simulating binary data
inverse_logit <- function(real_values = NULL) {
  probabilities <- (1)/(1 + exp(-1 * real_values))
  return (probabilities)
}


# Estimate the variance in Y before error
# Idea is that we fix this for given values of R2 and b
# We then get a compatible error term later on
determine_subgroup_var_Y <- function(num_total_conf = NULL,
                                     beta_X         = NULL,
                                     causal         = NULL,
                                     Z_correlation  = NULL,
                                     target_r_sq_Y  = NULL) {
  
  subgroup_size <- num_total_conf / 4
  beta_X_high   <- beta_X
  beta_X_low    <- beta_X / 4
  
  A <- (causal * beta_X_low)  + beta_X_high
  B <- (causal * beta_X_high) + beta_X_high
  C <- (causal * beta_X_high) + beta_X_low
  D <- (causal * beta_X_low)  + beta_X_low
  
  pairwise_combinations <- (A*B) + (A*C) + (A*D) + (B*C) + (B*D) + (C*D)
  
  LHS <- (A^2 + B^2 + C^2 + D^2) * (subgroup_size + (Z_correlation * subgroup_size * (subgroup_size - 1)))
  MID <- 2 * Z_correlation * subgroup_size^2 * pairwise_combinations
  RHS <- causal^2
  
  value <- LHS + MID + RHS
  
  return (value)
}


# Estimate the variance in the error term for Y
# Takes a given variance of Y and inverts the R2 formula
determine_subgroup_var_error_Y <- function(var_Y          = NULL,
                                           target_r_sq_Y  = NULL) {
  
  LHS   <- 1 - target_r_sq_Y
  RHS   <- var_Y / target_r_sq_Y
  value <- sqrt(LHS * RHS)
  
  # double here for bias induction
  value <- 2 * value
  
  return (value)
}


# The value for all beta coefficients used for generating X
beta_X_formula <- function(num_total_conf = NULL,
                           target_r_sq_X  = NULL,
                           Z_correlation  = NULL) {
  
  numerator   <- target_r_sq_X
  denominator <- num_total_conf * (1 - target_r_sq_X) * (1 + Z_correlation*(num_total_conf - 1))
  
  value <- sqrt(numerator / denominator)
  
  return (value)
}


# The value for the beta coefficients used for generating X
# 4 different values corresponding to the 4 subgroups
beta_X_subgroups_formula <- function(beta_X = NULL) {
  beta_X_1 <- beta_X     # HH
  beta_X_2 <- beta_X     # HL
  beta_X_3 <- beta_X / 4 # LH
  beta_X_4 <- beta_X / 4 # LL
  
  beta_Xs <- c(beta_X_1, beta_X_2, beta_X_3, beta_X_4)
  
  return (beta_Xs)
}


# The value for the beta coefficients used for generating Y
# 4 different values corresponding to the 4 subgroups
beta_Y_subgroups_formula <- function(beta_X = NULL) {
  beta_Y_1 <- beta_X     # HH
  beta_Y_2 <- beta_X / 4 # HL
  beta_Y_3 <- beta_X     # LH
  beta_Y_4 <- beta_X / 4 # LL
  
  beta_Ys <- c(beta_Y_1, beta_Y_2, beta_Y_3, beta_Y_4)
  
  # double size for bias induction
  beta_Ys <- 2 * beta_Ys
  
  return (beta_Ys)
}


estimate_within_CI <- function(model = NULL) {
  within_CI <- 0.0
  
  CI        <- confint(model, 'X', level = 0.95)
  if ((!is.na(CI[1])) && (!is.na(CI[2]))) {
    if ((causal > CI[1]) && (causal < CI[2])) {
      within_CI <- 1.0
    }
  }
  
  return (within_CI)
}



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

apply_MCAR_missingness <- function(data = NULL, vars_to_censor = NULL) {
  for (var in vars_to_censor) {
    # MCAR missingness selection vector
    MCAR_missingness <- rbinom(n_obs, size=1, prob=c(0.80))
    
    # Apply censorship to variable
    data[, var][MCAR_missingness == 1] <- NA
  }
  return (data)
}


apply_MNAR_missingness <- function(data = NULL, vars_to_censor = NULL) {
  for (var in vars_to_censor) {
    # MNAR missingness probability of censorship depends on data value
    psel <- rep(0.1, times = n_obs)
    for (i in 1:n_obs) {
      if (data[i, var] >= 3) {
        psel[i] <- 0.9
      }
      else if (data[i, var] <= 0.35) {
        psel[i] <- 0.7
      }
    }
    
    # MNAR missingness selection vector
    MNAR_missingness <- rbinom(n_obs, size=1, prob=psel)
    
    # Apply censorship to variable
    data[, var][MNAR_missingness == 1] <- NA
  }
  
  return (data)
}


# # Normal Z      -> imp_method = "pmm"     = predictive mean matching
# # Binary Z      -> imp_method = "logreg"  = logistic regression
# # Categorical Z -> imp_method = "polyreg" = polytomous regression
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
  
  imp <- mice(data,
              m      = num_datasets, # number of imputations
              maxit  = repetitions,  # number of iterations
              method = imp_method)
  
  imp_data <- complete(imp,
                       action  = "long", # stacked
                       include = FALSE)  # do not include original data
  
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
    MNAR_dataset <- apply_MNAR_missingness(FULL_dataset, vars_to_censor = vars_to_censor)
    MCAR_dataset <- apply_MNAR_missingness(FULL_dataset, vars_to_censor = vars_to_censor)
    
    # apply stacked MI
    handled_MNAR_dataset <- apply_stacked_MI(data         = MNAR_dataset,
                                             num_datasets = 5,
                                             repetitions  = 20,
                                             imp_method   = "pmm")
    handled_MCAR_dataset <- apply_stacked_MI(data         = MCAR_dataset,
                                             num_datasets = 5,
                                             repetitions  = 20,
                                             imp_method   = "pmm")
    
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
    if (binary_Y) {
      fully_adjusted_FULL_model <- glm("Y ~ .", data = FULL_dataset,         family = "binomial")
      fully_adjusted_MNAR_model <- glm("Y ~ .", data = handled_MNAR_dataset, family = "binomial")
      fully_adjusted_MCAR_model <- glm("Y ~ .", data = handled_MCAR_dataset, family = "binomial")
    } else {
      fully_adjusted_FULL_model  <- lm("Y ~ .", data = FULL_dataset)
      fully_adjusted_MNAR_model  <- lm("Y ~ .", data = handled_MNAR_dataset)
      fully_adjusted_MCAR_model  <- lm("Y ~ .", data = handled_MCAR_dataset)
    }
    
    # unadjusted model
    if (binary_Y) {
      unadjusted_FULL_model <- glm("Y ~ X", data = FULL_dataset,         family = "binomial")
      unadjusted_MNAR_model <- glm("Y ~ X", data = handled_MNAR_dataset, family = "binomial")
      unadjusted_MCAR_model <- glm("Y ~ X", data = handled_MCAR_dataset, family = "binomial")
    } else {
      unadjusted_FULL_model  <- lm("Y ~ X", data = FULL_dataset)
      unadjusted_MNAR_model  <- lm("Y ~ X", data = handled_MNAR_dataset)
      unadjusted_MCAR_model  <- lm("Y ~ X", data = handled_MCAR_dataset)
    }
    
    # LASSO (outcome Y)
    if (binary_Y) {
      
      # full dataset
      
      cv_lasso_FULL_model <- cv.glmnet(x = data.matrix(X_FULL_dataset), y = data.matrix(Y_FULL_column), family = 'binomial', alpha=1)
      lambda              <- cv_lasso_FULL_model$lambda.min
      lasso_FULL_model    <- glmnet(x = data.matrix(X_FULL_dataset), y = data.matrix(Y_FULL_column), family = 'binomial', alpha=1, lambda=lambda)
      
      lasso_FULL_coefs          <- as.vector(lasso_FULL_model$beta)
      names(lasso_FULL_coefs)   <- rownames(lasso_FULL_model$beta)
      lasso_FULL_vars_selected  <- union(c('X'), names(lasso_FULL_coefs[lasso_FULL_coefs != 0.0])) # always select X
      lasso_FULL_vars_selected  <- lasso_FULL_vars_selected[lasso_FULL_vars_selected != "(Intercept)"] # exclude intercept
      two_step_lasso_FULL_model <- glm(make_model_formula(vars_selected = lasso_FULL_vars_selected), data = FULL_dataset, family = "binomial")
      
      # MNAR data
      
      cv_lasso_MNAR_model <- cv.glmnet(x = data.matrix(X_handled_MNAR_dataset), y = data.matrix(Y_handled_MNAR_column), family = 'binomial', alpha=1)
      lambda              <- cv_lasso_MNAR_model$lambda.min
      lasso_MNAR_model    <- glmnet(x = data.matrix(X_handled_MNAR_dataset), y = data.matrix(Y_handled_MNAR_column), family = 'binomial', alpha=1, lambda=lambda)
      
      lasso_MNAR_coefs          <- as.vector(lasso_MNAR_model$beta)
      names(lasso_MNAR_coefs)   <- rownames(lasso_MNAR_model$beta)
      lasso_MNAR_vars_selected  <- union(c('X'), names(lasso_MNAR_coefs[lasso_MNAR_coefs != 0.0])) # always select X
      lasso_MNAR_vars_selected  <- lasso_MNAR_vars_selected[lasso_MNAR_vars_selected != "(Intercept)"] # exclude intercept
      two_step_lasso_MNAR_model <- glm(make_model_formula(vars_selected = lasso_MNAR_vars_selected), data = handled_MNAR_dataset, family = "binomial")
      
      # MCAR data
      
      cv_lasso_MCAR_model <- cv.glmnet(x = data.matrix(X_handled_MCAR_dataset), y = data.matrix(Y_handled_MCAR_column), family = 'binomial', alpha=1)
      lambda              <- cv_lasso_MCAR_model$lambda.min
      lasso_MCAR_model    <- glmnet(x = data.matrix(X_handled_MCAR_dataset), y = data.matrix(Y_handled_MCAR_column), family = 'binomial', alpha=1, lambda=lambda)
      
      lasso_MCAR_coefs          <- as.vector(lasso_MCAR_model$beta)
      names(lasso_MCAR_coefs)   <- rownames(lasso_MCAR_model$beta)
      lasso_MCAR_vars_selected  <- union(c('X'), names(lasso_MCAR_coefs[lasso_MCAR_coefs != 0.0])) # always select X
      lasso_MCAR_vars_selected  <- lasso_MCAR_vars_selected[lasso_MCAR_vars_selected != "(Intercept)"] # exclude intercept
      two_step_lasso_MCAR_model <- glm(make_model_formula(vars_selected = lasso_MCAR_vars_selected), data = handled_MCAR_dataset, family = "binomial")
      
    } else {
      
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
    }
    
    # LASSO_X (exposure X)
    if (binary_X) {
      
      # full dataset
      
      cv_lasso_X_FULL_model <- cv.glmnet(x = data.matrix(Z_FULL_dataset), y = data.matrix(X_FULL_column), family = 'binomial', alpha=1)
      lambda                <- cv_lasso_X_FULL_model$lambda.min
      lasso_X_FULL_model    <- glmnet(x = data.matrix(Z_FULL_dataset), y = data.matrix(X_FULL_column), family = 'binomial', alpha=1, lambda=lambda)
      
      lasso_X_FULL_coefs          <- as.vector(lasso_X_FULL_model$beta)
      names(lasso_X_FULL_coefs)   <- rownames(lasso_X_FULL_model$beta)
      lasso_X_FULL_vars_selected  <- union(c('X'), names(lasso_X_FULL_coefs[lasso_X_FULL_coefs != 0.0])) # always select X
      lasso_X_FULL_vars_selected  <- lasso_X_FULL_vars_selected[lasso_X_FULL_vars_selected != "(Intercept)"] # exclude intercept
      two_step_lasso_X_FULL_model <- glm(make_model_formula(vars_selected = lasso_X_FULL_vars_selected), data = FULL_dataset, family = "binomial")
      
      # MNAR data
      
      cv_lasso_X_MNAR_model <- cv.glmnet(x = data.matrix(Z_handled_MNAR_dataset), y = data.matrix(X_handled_MNAR_column), family = 'binomial', alpha=1)
      lambda                <- cv_lasso_X_MNAR_model$lambda.min
      lasso_X_MNAR_model    <- glmnet(x = data.matrix(Z_handled_MNAR_dataset), y = data.matrix(X_handled_MNAR_column), family = 'binomial', alpha=1, lambda=lambda)
      
      lasso_X_MNAR_coefs          <- as.vector(lasso_X_MNAR_model$beta)
      names(lasso_X_MNAR_coefs)   <- rownames(lasso_X_MNAR_model$beta)
      lasso_X_MNAR_vars_selected  <- union(c('X'), names(lasso_X_MNAR_coefs[lasso_X_MNAR_coefs != 0.0])) # always select X
      lasso_X_MNAR_vars_selected  <- lasso_X_MNAR_vars_selected[lasso_X_MNAR_vars_selected != "(Intercept)"] # exclude intercept
      two_step_lasso_X_MNAR_model <- glm(make_model_formula(vars_selected = lasso_X_MNAR_vars_selected), data = handled_MNAR_dataset, family = "binomial")
      
      # MCAR data
      
      cv_lasso_X_MCAR_model <- cv.glmnet(x = data.matrix(Z_handled_MCAR_dataset), y = data.matrix(X_handled_MCAR_column), family = 'binomial', alpha=1)
      lambda                <- cv_lasso_X_MCAR_model$lambda.min
      lasso_X_MCAR_model    <- glmnet(x = data.matrix(Z_handled_MCAR_dataset), y = data.matrix(X_handled_MCAR_column), family = 'binomial', alpha=1, lambda=lambda)
      
      lasso_X_MCAR_coefs          <- as.vector(lasso_X_MCAR_model$beta)
      names(lasso_X_MCAR_coefs)   <- rownames(lasso_X_MCAR_model$beta)
      lasso_X_MCAR_vars_selected  <- union(c('X'), names(lasso_X_MCAR_coefs[lasso_X_MCAR_coefs != 0.0])) # always select X
      lasso_X_MCAR_vars_selected  <- lasso_X_MCAR_vars_selected[lasso_X_MCAR_vars_selected != "(Intercept)"] # exclude intercept
      two_step_lasso_X_MCAR_model <- glm(make_model_formula(vars_selected = lasso_X_MCAR_vars_selected), data = handled_MCAR_dataset, family = "binomial")
      
    } else {
      
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
    }
    
    # LASSO_union (exposure X and outcome Y)
    if (binary_Y) {
      
      # full dataset
      lasso_union_FULL_vars_selected  <- union(lasso_FULL_vars_selected, lasso_X_FULL_vars_selected)
      two_step_lasso_union_FULL_model <- glm(make_model_formula(vars_selected = lasso_union_FULL_vars_selected), data = FULL_dataset, family = "binomial")
      
      # MNAR
      lasso_union_MNAR_vars_selected  <- union(lasso_MNAR_vars_selected, lasso_X_MNAR_vars_selected)
      two_step_lasso_union_MNAR_model <- glm(make_model_formula(vars_selected = lasso_union_MNAR_vars_selected), data = MNAR_dataset, family = "binomial")
      
      # MCAR
      lasso_union_MCAR_vars_selected  <- union(lasso_MCAR_vars_selected, lasso_X_MCAR_vars_selected)
      two_step_lasso_union_MCAR_model <- glm(make_model_formula(vars_selected = lasso_union_MCAR_vars_selected), data = MCAR_dataset, family = "binomial")
      
    } else {
      
      # full dataset
      lasso_union_FULL_vars_selected  <- union(lasso_FULL_vars_selected, lasso_X_FULL_vars_selected)
      two_step_lasso_union_FULL_model <- lm(make_model_formula(vars_selected = lasso_union_FULL_vars_selected), data = FULL_dataset)
      
      # MNAR
      lasso_union_MNAR_vars_selected  <- union(lasso_MNAR_vars_selected, lasso_X_MNAR_vars_selected)
      two_step_lasso_union_MNAR_model <- lm(make_model_formula(vars_selected = lasso_union_MNAR_vars_selected), data = MNAR_dataset)
      
      # MCAR
      lasso_union_MCAR_vars_selected  <- union(lasso_MCAR_vars_selected, lasso_X_MCAR_vars_selected)
      two_step_lasso_union_MCAR_model <- lm(make_model_formula(vars_selected = lasso_union_MCAR_vars_selected), data = MCAR_dataset)
      
    }
    
    
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
  
  names(TRUE_coefs) <- c("Z on X (LL)",
                         "Z on X (LH)",
                         "Z on X (HL)",
                         "Z on X (HH)",
                         
                         "Z on Y (LL)",
                         "Z on Y (LH)",
                         "Z on Y (HL)",
                         "Z on Y (HH)",
                         
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
               MCAR_results))
}


