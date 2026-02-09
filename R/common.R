# ****************************************
# Missingness Simulation Study
#
# Common functions used in multiple places
#
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 19/01/2026
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
  
  value <- 44.2
  
  return (value)
}


# The value for all beta coefficients used for generating X
beta_X_formula <- function(num_total_conf = NULL,
                           target_r_sq_X  = NULL,
                           Z_correlation  = NULL) {
  
  numerator   <- target_r_sq_X
  denominator <- num_total_conf * (1 - target_r_sq_X) * (1 + Z_correlation*(num_total_conf - 1))
  
  value <- sqrt(numerator / denominator)
  
  value <- 0.305
  
  return (value)
}


# The value for the beta coefficients used for generating X
# 4 different values corresponding to the 4 subgroups
beta_X_subgroups_formula <- function(beta_X = NULL) {
  beta_X_1 <- beta_X / 4 # low X, low Y
  beta_X_2 <- beta_X / 4 # low X, high Y
  beta_X_3 <- beta_X     # high X, low Y
  beta_X_4 <- beta_X     # high X, high Y
  
  beta_Xs <- c(beta_X_1, beta_X_2, beta_X_3, beta_X_4)
  
  return (beta_Xs)
}


# The value for the beta coefficients used for generating Y
# 4 different values corresponding to the 4 subgroups
beta_Y_subgroups_formula <- function(beta_X = NULL) {
  beta_Y_1 <- beta_X / 4 # low X, low Y
  beta_Y_2 <- beta_X     # low X, high Y
  beta_Y_3 <- beta_X / 4 # high X, low Y
  beta_Y_4 <- beta_X     # high X, high Y
  
  beta_Ys <- c(beta_Y_1, beta_Y_2, beta_Y_3, beta_Y_4)
  
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



# ----- Missingness mechanisms and missingness handling -----

apply_MCAR_missingness <- function(data = NULL, vars_to_censor = NULL) {
  # MCAR probability of censorship does not depend on data value
  psel <- 0.25
  
  for (var in vars_to_censor) {
    # MCAR selection vector
    MCAR_selection <- rbinom(n_obs, size=1, prob=psel)
    
    # Apply censorship to variable for every observation which is NOT selected
    data[, var][MCAR_selection == 0] <- NA
  }
  
  return (list(data, psel, MCAR_selection))
}


apply_MNAR_missingness <- function(data = NULL, vars_to_censor = NULL) {
  for (var in vars_to_censor) {
    # MNAR probability of censorship depends on data value
    psel <- rep(0.1, times = n_obs)
    
    for (i in 1:n_obs) {
      if (data[i, var] >= 1) {
        psel[i] <- 0.9
      }
      else if (data[i, var] <= -1.65) {
        psel[i] <- 0.7
      }
    }
    
    # MNAR selection vector
    MNAR_selection <- rbinom(n_obs, size=1, prob=psel)
    
    # Apply censorship to variable for every observation which is NOT selected
    data[, var][MNAR_selection == 0] <- NA
  }
  
  return (list(data, psel, MNAR_selection))
}


apply_MAR_missingness <- function(data = NULL, vars_to_censor = NULL) {
  # generate sub of all covariates in low, low subgroup (i.e Z_new = sum(Z1, ..., Z8))
  Z_sum_subgroup_LL <- rowSums(data[, c("Z2", "Z3", "Z4", "Z5", "Z6", "Z7", "Z8")])
  
  # find 75th percentile of this sum
  percentile_75 <- quantile(Z_sum_subgroup_LL, 0.75)
  psel <- 0.25
  
  # select into study population where value of sum is above 75th percentile
  MAR_selection <- (Z_sum_subgroup_LL >= percentile_75)
  
  # Apply censorship to variable for every observation which is NOT selected
  for (var in vars_to_censor) {
    data[, var][MAR_selection == 0] <- NA
  }
  
  return (list(data, psel, MAR_selection))
}


