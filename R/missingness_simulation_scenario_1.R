# ****************************************
# Missingness Simulation Study
#
# Simulation procedure for confounder-handling
# with missingness considerations
# Scenario 1: CCA missingness handling
#
# Emma Tarmey
#
# Started:          06/10/2025
# Most Recent Edit: 01/02/2026
# ****************************************



# ----- Tech Setup ----- 

# clear R memory
rm(list=ls())

# check all external libraries
using<-function(...) {
  libs <- unlist(list(...))
  req  <- unlist(lapply(libs, require, character.only=TRUE))
  need <- libs[req==FALSE]
  if(length(need) > 0){ 
    install.packages(need)
    lapply(need, require, character.only=TRUE)
  }
}
using("dplyr", "glmnet", "mice", "speedglm", "tidyr")

# fix wd issue
# forces wd to be the location of this file
if (Sys.getenv("RSTUDIO") == "1") {
  setwd(dirname(rstudioapi::getSourceEditorContext()$path))
}

source("common.R")



# ----- Parameters ------

n_scenario   <- 1

n_obs             <- 100000
n_rep             <- 20
Z_correlation     <- 0.1
Z_subgroups       <- 4
target_r_sq_X     <- 0.8
target_r_sq_Y     <- 0.2
causal            <- 0.5

binary_X          <- FALSE
binary_Y          <- FALSE
binary_Z          <- FALSE

num_total_conf  <- 32
num_meas_conf   <- 28
num_unmeas_conf <- 4

# missingness handling mechanism
missingness_handling <- "CCA"

# confounders to be unmeasured
#vars_to_make_unmeasured <- c()
vars_to_make_unmeasured <- c("Z1", "Z9", "Z17", "Z25")

# confounders to have missingness applied
vars_to_censor <- c("Z26")


# ----- RNG -----

# # test RNG
# set.seed(2025)

# fix RNG seed based on current scenario
seeds_df <- read.csv(file = "../data/precomputed_RNG_seeds.csv")
seed     <- seeds_df %>% filter(simulation_scenario == n_scenario)
set.seed(seed$seed)



# ------ Run simulation procedure ------

# CCA
source("missingness_simulation_method_CCA.R")
simulation_results <- run_CCA_simulation(n_scenario = n_scenario,
                                         n_obs      = n_obs,
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
                                         vars_to_censor          = vars_to_censor)

FULL_cov_selection <- simulation_results[[1]]
MNAR_cov_selection <- simulation_results[[2]]
MCAR_cov_selection <- simulation_results[[3]]

TRUE_coefs <- simulation_results[[4]]
FULL_coefs <- simulation_results[[5]]
MNAR_coefs <- simulation_results[[6]]
MCAR_coefs <- simulation_results[[7]]

FULL_results <- simulation_results[[8]]
MNAR_results <- simulation_results[[9]]
MCAR_results <- simulation_results[[10]]

sample_size_table <- simulation_results[[11]]

sample_FULL_dataset         <- simulation_results[[12]]
sample_handled_MNAR_dataset <- simulation_results[[13]]
sample_handled_MCAR_dataset <- simulation_results[[14]]


# ----- Save results -----

# take mean across all repetitions
final_FULL_cov_selection <- as.data.frame(apply(FULL_cov_selection, c(1,2), mean))
final_MNAR_cov_selection <- as.data.frame(apply(MNAR_cov_selection, c(1,2), mean))
final_MCAR_cov_selection <- as.data.frame(apply(MCAR_cov_selection, c(1,2), mean))

final_FULL_coefs <- as.data.frame(apply(FULL_coefs, c(1,2), mean))
final_MNAR_coefs <- as.data.frame(apply(MNAR_coefs, c(1,2), mean))
final_MCAR_coefs <- as.data.frame(apply(MCAR_coefs, c(1,2), mean))

final_FULL_results <- as.data.frame(apply(FULL_results, c(1,2), mean))
final_MNAR_results <- as.data.frame(apply(MNAR_results, c(1,2), mean))
final_MCAR_results <- as.data.frame(apply(MCAR_results, c(1,2), mean))

sample_size_table <- as.data.frame(sample_size_table)

# record empirical standard error separately
causal        <- 0.5
model_methods <- c("fully_adjusted", "unadjusted", "two_step_lasso", "two_step_lasso_X", "two_step_lasso_union")
for (method in model_methods) {
  # FULL
  FULL_causal_effect_estimates                <- c(FULL_results[method, "causal_estimate", ])
  final_FULL_results[ method, "empirical_SE"] <- sd(FULL_causal_effect_estimates)
  
  # MNAR
  MNAR_causal_effect_estimates                <- c(MNAR_results[method, "causal_estimate", ])
  final_MNAR_results[ method, "empirical_SE"] <- sd(MNAR_causal_effect_estimates)
  
  # MCAR
  MCAR_causal_effect_estimates                <- c(MCAR_results[method, "causal_estimate", ])
  final_MCAR_results[ method, "empirical_SE"] <- sd(MCAR_causal_effect_estimates)
}

# save simulation setup
sim_setup <- c(n_obs,
               n_rep,
               Z_correlation,
               Z_subgroups,
               target_r_sq_X,
               target_r_sq_Y,
               causal,
               binary_X,
               binary_Y,
               binary_Z,
               num_total_conf,
               num_meas_conf,
               num_unmeas_conf,
               missingness_handling)

names(sim_setup) <- c("n_obs",
                      "n_rep",
                      "Z_correlation",
                      "Z_subgroups",
                      "target_r_sq_X",
                      "target_r_sq_Y",
                      "causal",
                      "binary_X",
                      "binary_Y",
                      "binary_Z",
                      "num_total_conf",
                      "num_meas_conf",
                      "num_unmeas_conf",
                      "missingness_handling")


message("\n\n\n\n***** Simulation Setup *****")
print(sim_setup)


message("\n\n\n\n ***** Confounders *****")
print("Unmeasured confounders")
print(vars_to_make_unmeasured)
print("Confounders with missingness")
print(vars_to_censor)


message("\n\n\n ***** Sample sizes *****")
print(sample_size_table)


message("\n\n\n\n ***** True Coefficients *****")
print("Effect of confounders Z_i on exposure X")
print(TRUE_coefs[c(1:4)])
print("Effect of confounders Z_i on outcome Y")
print(TRUE_coefs[c(5:8)])
print("Causal effect of exposure X on outcome Y:")
print(TRUE_coefs[c(9)])
print("Variance in error term for outcome Y:")
print(TRUE_coefs[c(10)])


message("\n\n\n\n***** FULL *****")
message("\nCovariate selection")
print(final_FULL_cov_selection)
message("\nCoefficients")
print(final_FULL_coefs[, c(1:2)])
print(final_FULL_coefs[, c(3:10)])
print(final_FULL_coefs[, c(11:18)])
print(final_FULL_coefs[, c(19:26)])
print(final_FULL_coefs[, c(27:34)])
message("\nResults")
print(final_FULL_results[, c(1:5)])
print(final_FULL_results[, c(6:8)])
print(final_FULL_results[, c(9:10)])


message("\n\n\n\n***** MNAR *****")
message("\nCovariate selection")
print(final_MNAR_cov_selection)
message("\nCoefficients")
print(final_MNAR_coefs[, c(1:2)])
print(final_MNAR_coefs[, c(3:10)])
print(final_MNAR_coefs[, c(11:18)])
print(final_MNAR_coefs[, c(19:26)])
print(final_MNAR_coefs[, c(27:34)])
message("\nResults")
print(final_MNAR_results[, c(1:5)])
print(final_MNAR_results[, c(6:8)])
print(final_MNAR_results[, c(9:10)])


message("\n\n\n\n***** MCAR *****")
message("\nCovariate selection")
print(final_MCAR_cov_selection)
message("\nCoefficients")
print(final_MCAR_coefs[, c(1:2)])
print(final_MCAR_coefs[, c(3:10)])
print(final_MCAR_coefs[, c(11:18)])
print(final_MCAR_coefs[, c(19:26)])
print(final_MCAR_coefs[, c(27:34)])
message("\nResults")
print(final_MCAR_results[, c(1:5)])
print(final_MCAR_results[, c(6:8)])
print(final_MCAR_results[, c(9:10)])


# Save to file
id_string <- paste("missingness_sim_scenario_", n_scenario, "_missingness_method_", missingness_handling, sep='')
message(paste0("\nSaving results for ", id_string))

write.csv(final_FULL_cov_selection,   paste("../data/", id_string, "_FULL_cov_selection.csv", sep=''))
write.csv(final_MNAR_cov_selection,   paste("../data/", id_string, "_MNAR_cov_selection.csv", sep=''))
write.csv(final_MCAR_cov_selection,   paste("../data/", id_string, "_MCAR_cov_selection.csv", sep=''))

write.csv(TRUE_coefs,         paste("../data/", id_string, "_TRUE_coefs.csv", sep=''))
write.csv(final_FULL_coefs,   paste("../data/", id_string, "_FULL_coefs.csv", sep=''))
write.csv(final_MNAR_coefs,   paste("../data/", id_string, "_MNAR_coefs.csv", sep=''))
write.csv(final_MCAR_coefs,   paste("../data/", id_string, "_MCAR_coefs.csv", sep=''))

write.csv(final_FULL_results,   paste("../data/", id_string, "_FULL_results.csv", sep=''))
write.csv(final_MNAR_results,   paste("../data/", id_string, "_MNAR_results.csv", sep=''))
write.csv(final_MCAR_results,   paste("../data/", id_string, "_MCAR_results.csv", sep=''))

write.csv(sample_size_table,           paste("../data/", id_string, "_sample_size_table.csv", sep=''))
write.csv(sample_FULL_dataset,         paste("../data/", id_string, "_sample_FULL_dataset.csv", sep=''))
write.csv(sample_handled_MNAR_dataset, paste("../data/", id_string, "_sample_handled_MNAR_dataset.csv", sep=''))
write.csv(sample_handled_MCAR_dataset, paste("../data/", id_string, "_sample_handled_MCAR_dataset.csv", sep=''))


