#' Function to calculate the pointwise MTD bias
#' @param df Dataset containing outputs (data.frame)
#' @param modeling Type (character) of modeling approach (Independent modeling, Latent variables, FGM Copula, Clayton Copula) (character)
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)

get_pointwise_MTD_bias = function(df, DLT_scenario, EFF_scenario){
  
  pointwise_MTD_bias = as_tibble(as.data.frame(do.call(rbind, map(df, "pointwise_MTD_bias"))),rownames = NA) %>% 
    group_by(x_true_MTD) %>%
    summarise(average_pointwise_MTD_bias = mean(objective_bias))
  
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         x_true_MTD = pointwise_MTD_bias$x_true_MTD,
                         average_pointwise_MTD_bias = pointwise_MTD_bias$average_pointwise_MTD_bias)
  
  output_df
}
