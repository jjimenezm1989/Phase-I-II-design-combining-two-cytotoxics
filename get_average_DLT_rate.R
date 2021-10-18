#' Function to calculate the average DLT rate and proportion of trials with DLT rate above theta_Z + 0.1
#' @param df Dataset containing outputs (data.frame)
#' @param modeling Type (character) of modeling approach (Independent modeling, Latent variables, FGM Copula, Clayton Copula) (character)
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)

get_average_DLT_rate = function(df, DLT_scenario, EFF_scenario){
  
  Z_by_trial = map(df,"Z")
  
  #Average DLT rate
  average_DLT_rate = mean(unlist(Z_by_trial))
  
  #Proportion of trials with DLT rate above theta_Z + 0.1
  DLT_rate_by_trial = function(x){
    ifelse(sum(as.numeric(x))/length(as.numeric(x)) > (theta_Z+0.1), 1, 0)
  }
  proportion_trial_with_DLT_rate_above_threshold = mean(unlist(lapply(Z_by_trial,DLT_rate_by_trial)))
  
  #Output
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         DLT_rate_by_trial = average_DLT_rate,
                         proportion_trial_with_DLT_rate_above_threshold = proportion_trial_with_DLT_rate_above_threshold)
  
  output_df
  
}