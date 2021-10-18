setwd("/home/jimenjox/Research/Combining two cytotoxic agents with continuous dose levels in phase I-II clinical trials/R code/")

library(rjags)
library(tidyverse)
library(scales)
library(mvtnorm)
library(truncnorm)
library(data.table)
library(ggpubr)

rm(list=ls())

########################
#Scenario specification#
########################

#Seet seed based on time
time_seed = as.numeric(Sys.time())
set.seed(time_seed)

# Number of Trials to simulate
M_iter = 1000

#Choose the toxicity (DLT) scenario number and efficacy (EFF) scenario number
DLT_scenario_variables = "_DLT_sc2"
EFF_scenario_variables = "_EFF_sc1"
DLT_scenario_full_description = "Dose-toxicity scenario 2"
EFF_scenario_full_description = "Dose-efficacy scenario 1"
hypothesis = "_H0"

#Choose between "Power" (if under H1) or "Type-I error (if under H0)".
power_or_typeIerror = "Type-I error"

#######################
#User-defined function#
#######################

source("generate_bivariate_outcome.R")
source("pdlt.R")
source("peff.R")
source("twodimmtd2.R")
source("parameters_declaration.R")
source("trial_simulation.R")

path_JAGS = getwd()

#######################
#Main part of the code#
#######################

true_rho00 = eval(parse(text = paste0("true_rho00",DLT_scenario_variables)))
true_rho01 = eval(parse(text = paste0("true_rho01",DLT_scenario_variables)))
true_rho10 = eval(parse(text = paste0("true_rho10",DLT_scenario_variables)))
true_alpha3 = eval(parse(text = paste0("true_alpha3",DLT_scenario_variables)))
true_beta0 = eval(parse(text = paste0("true_beta0",DLT_scenario_variables,EFF_scenario_variables,hypothesis)))
true_beta1 = eval(parse(text = paste0("true_beta1",DLT_scenario_variables,EFF_scenario_variables,hypothesis)))
true_beta2 = eval(parse(text = paste0("true_beta2",DLT_scenario_variables,EFF_scenario_variables,hypothesis)))
true_beta3 = eval(parse(text = paste0("true_beta3",DLT_scenario_variables,EFF_scenario_variables,hypothesis)))
true_beta4 = eval(parse(text = paste0("true_beta4",DLT_scenario_variables,EFF_scenario_variables,hypothesis)))
true_beta5 = eval(parse(text = paste0("true_beta5",DLT_scenario_variables,EFF_scenario_variables,hypothesis)))

start.time <- Sys.time()

do.call("<-",list(paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                  replicate(M_iter, trial_simulation(true_rho00 = eval(parse(text = paste0("true_rho00",DLT_scenario_variables))),
                                                     true_rho01 = eval(parse(text = paste0("true_rho01",DLT_scenario_variables))),
                                                     true_rho10 = eval(parse(text = paste0("true_rho10",DLT_scenario_variables))),
                                                     true_alpha3 = eval(parse(text = paste0("true_alpha3",DLT_scenario_variables))),
                                                     true_beta0 = eval(parse(text = paste0("true_beta0",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                     true_beta1 = eval(parse(text = paste0("true_beta1",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                     true_beta2 = eval(parse(text = paste0("true_beta2",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                     true_beta3 = eval(parse(text = paste0("true_beta3",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                     true_beta4 = eval(parse(text = paste0("true_beta4",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                     true_beta5 = eval(parse(text = paste0("true_beta5",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                     path_JAGS = path_JAGS), simplify = FALSE)
                  
                  
))
end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

save.image(paste0(path_JAGS,"/environment_sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".RData"))

#--------DERIVE OPERATING CHARACTERISTICS AND SAVE THEM IN .CSV FILES#--------#

setwd("/home/jimenjox/Research/Combining two cytotoxic agents with continuous dose levels in phase I-II clinical trials/R code/")

source("get_bias_MSE.R")
source("get_power.R")
source("get_posterior_probabilities_stopping_rules.R")
source("get_pointwise_MTD_bias.R")
source("get_percent_correct_MTD_recommendation.R")
source("get_average_post_efficacy_contour.R")
source("get_proportion_correct_patient_allocation.R")
source("get_dose_combinations.R")
source("get_MTD_curve.R")
source("get_average_DLT_rate.R")
source("get_sample_size_stopping_rules.R")
source("get_optimal_dose_combination.R")
source("get_proportion_optimal_dose_combination_above_efficacy_threshold.R")

setwd("/home/jimenjox/Research/Combining two cytotoxic agents with continuous dose levels in phase I-II clinical trials/R code/intermediate_outputs/")

#optimal dose combination


write.csv(do.call("<-",list(paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_optimal_dose_combination(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                         DLT_scenario = DLT_scenario_full_description,
                                                         EFF_scenario = EFF_scenario_full_description))),
          paste0("optimal_dose_combination",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))



#Bias MTD
write.csv(do.call("<-",list(paste0("bias_MTD",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_pointwise_MTD_bias(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                   DLT_scenario = DLT_scenario_full_description,
                                                   EFF_scenario = EFF_scenario_full_description))),
          paste0("bias_MTD",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))

# % correct MTD recommendation
write.csv(do.call("<-",list(paste0("correct_recommendation_MTD",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_percent_correct_MTD_recommendation(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                                   DLT_scenario = DLT_scenario_full_description,
                                                                   EFF_scenario = EFF_scenario_full_description))),
          paste0("correct_recommendation_MTD",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))

#Power / type-I error

write.csv(do.call("<-",list(paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_power(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                      type = power_or_typeIerror,
                                      DLT_scenario = DLT_scenario_full_description,
                                      EFF_scenario = EFF_scenario_full_description))),
          paste0("power_or_typeIerror",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))


#Stopping rules
write.csv(do.call("<-",list(paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_posterior_probabilities_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                                       delta_Z1 = delta_Z1,
                                                                       delta_Z2 = delta_Z2,
                                                                       delta_E0 = delta_E0,
                                                                       DLT_scenario = DLT_scenario_full_description,
                                                                       EFF_scenario = EFF_scenario_full_description))),
          paste0("stopping_rules",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))

#Sample size stopping rules
write.csv(do.call("<-",list(paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_sample_size_stopping_rules(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis)))))),
          paste0("sample_size_stopping_rules",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))


#Proportion patient optimal dose combinations with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_proportion_optimal_dose_combination_above_efficacy_threshold(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                                                             DLT_scenario = DLT_scenario_full_description,
                                                                                             EFF_scenario = EFF_scenario_full_description))),
          paste0("optimal_dose_combination_above_threshold",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))




#Proportion patient allocation is dose combination with true P(pi_E) > theta_E
write.csv(do.call("<-",list(paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_proportion_correct_patient_allocation(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                                      DLT_scenario = DLT_scenario_full_description,
                                                                      EFF_scenario = EFF_scenario_full_description))),
          paste0("proportion_allocation",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))


#Bias and MSE model parameters
write.csv(do.call("<-",list(paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_bias_MSE(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                         DLT_scenario = DLT_scenario_full_description,
                                         EFF_scenario = EFF_scenario_full_description))),
          paste0("bias_MSE_parameters",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))


#Dose combinations to which patients were allocated
write.csv(do.call("<-",list(paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_dose_combinations(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                  DLT_scenario = DLT_scenario_full_description,
                                                  EFF_scenario = EFF_scenario_full_description))),
          paste0("dose_combinations",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))

# MTD curve
write.csv(do.call("<-",list(paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_MTD_curve(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                          DLT_scenario = DLT_scenario_full_description,
                                          EFF_scenario = EFF_scenario_full_description,
                                          true_rho00 = eval(parse(text = paste0("true_rho00",DLT_scenario_variables))),
                                          true_rho01 = eval(parse(text = paste0("true_rho01",DLT_scenario_variables))),
                                          true_rho10 = eval(parse(text = paste0("true_rho10",DLT_scenario_variables))),
                                          true_alpha3 = eval(parse(text = paste0("true_alpha3",DLT_scenario_variables)))))),
          paste0("MTD_curve",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))



#DLT rate and proportion of trials with DLT rate above theta_Z + 0.1
write.csv(do.call("<-",list(paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,hypothesis),
                            get_average_DLT_rate(df = eval(parse(text = paste0("sim_data",DLT_scenario_variables,EFF_scenario_variables,hypothesis))),
                                                 DLT_scenario = DLT_scenario_full_description,
                                                 EFF_scenario = EFF_scenario_full_description))),
          paste0("DLT_rate",DLT_scenario_variables,EFF_scenario_variables,hypothesis,".csv"))

