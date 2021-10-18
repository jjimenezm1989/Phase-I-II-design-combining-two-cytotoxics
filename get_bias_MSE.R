#' Function to calculate bias and MSE from the "df" dataset
#' @param df Dataset containing outputs (data.frame)
#' @param modeling Type of modeling approach (Independent modeling, Latent variables, FGM Copula, Clayton Copula) (character)
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)

get_bias_MSE = function(df, DLT_scenario, EFF_scenario){
  
  bias_absolute_difference_rho00 = mean(unlist(map(df, "absolute_difference_rho00")))
  mse_absolute_difference_rho00 = mean(unlist(map(df, "absolute_difference_rho00"))^2)
  
  bias_absolute_difference_rho10 = mean(unlist(map(df, "absolute_difference_rho10")))
  mse_absolute_difference_rho10 = mean(unlist(map(df, "absolute_difference_rho10"))^2)
  
  bias_absolute_difference_rho01 = mean(unlist(map(df, "absolute_difference_rho01")))
  mse_absolute_difference_rho01 = mean(unlist(map(df, "absolute_difference_rho01"))^2)
  
  bias_absolute_difference_alpha3 = mean(unlist(map(df, "absolute_difference_alpha3")))
  mse_absolute_difference_alpha3 = mean(unlist(map(df, "absolute_difference_alpha3"))^2)
  
  bias_absolute_difference_beta0 = mean(unlist(map(df, "absolute_difference_beta0")))
  mse_absolute_difference_beta0 = mean(unlist(map(df, "absolute_difference_beta0"))^2)

  bias_absolute_difference_beta1 = mean(unlist(map(df, "absolute_difference_beta1")))
  mse_absolute_difference_beta1 = mean(unlist(map(df, "absolute_difference_beta1"))^2)
  
  bias_absolute_difference_beta2 = mean(unlist(map(df, "absolute_difference_beta2")))
  mse_absolute_difference_beta2 = mean(unlist(map(df, "absolute_difference_beta2"))^2)
  
  bias_absolute_difference_beta3 = mean(unlist(map(df, "absolute_difference_beta3")))
  mse_absolute_difference_beta3 = mean(unlist(map(df, "absolute_difference_beta3"))^2)
  
  bias_absolute_difference_beta4 = mean(unlist(map(df, "absolute_difference_beta4")))
  mse_absolute_difference_beta4 = mean(unlist(map(df, "absolute_difference_beta4"))^2)
  
  bias_absolute_difference_beta5 = mean(unlist(map(df, "absolute_difference_beta5")))
  mse_absolute_difference_beta5 = mean(unlist(map(df, "absolute_difference_beta5"))^2)
  
 
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         type = c(rep("Bias",10),rep("MSE",10)),
                         parameter = rep(c("rho00", "rho10", "rho01", "alpha3", "beta0", "beta1","beta2", "beta3", "beta4", "beta5"),2),
                         value = c(bias_absolute_difference_rho00,
                                   bias_absolute_difference_rho10,
                                   bias_absolute_difference_rho01,
                                   bias_absolute_difference_alpha3,
                                   bias_absolute_difference_beta0,
                                   bias_absolute_difference_beta1,
                                   bias_absolute_difference_beta2,
                                   bias_absolute_difference_beta3,
                                   bias_absolute_difference_beta4,
                                   bias_absolute_difference_beta5,
                                   mse_absolute_difference_rho00,
                                   mse_absolute_difference_rho10,
                                   mse_absolute_difference_rho01,
                                   mse_absolute_difference_alpha3,
                                   mse_absolute_difference_beta0,
                                   mse_absolute_difference_beta1,
                                   mse_absolute_difference_beta2,
                                   mse_absolute_difference_beta3,
                                   mse_absolute_difference_beta4,
                                   mse_absolute_difference_beta5))
  
  
  output_df
  
}











