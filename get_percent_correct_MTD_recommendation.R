#' Function to calculate the percent of correct MTD recommentation
#' @param df Dataset containing outputs (data.frame)
#' @param modeling Type (character) of modeling approach (Independent modeling, Latent variables, FGM Copula, Clayton Copula) (character)
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)

get_percent_correct_MTD_recommendation = function(df, DLT_scenario, EFF_scenario){
  
  percent_correct_MTD_recommentation = as_tibble(as.data.frame(do.call(rbind, map(df, "percent_correct_MTD_recommentation"))),rownames = NA) %>% 
    group_by(x_true_MTD) %>%
    summarise(correct_recomm_threshold_1 = mean(correct_recomm_threshold_1),
              correct_recomm_threshold_2 = mean(correct_recomm_threshold_2))
  
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         x_true_MTD = rep(percent_correct_MTD_recommentation$x_true_MTD,2),
                         percent_correct_MTD_recommentation = c(percent_correct_MTD_recommentation$correct_recomm_threshold_1,percent_correct_MTD_recommentation$correct_recomm_threshold_2),
                         type = c(rep("Threshold 1",length(percent_correct_MTD_recommentation$x_true_MTD)),
                                  rep("Threshold 2",length(percent_correct_MTD_recommentation$x_true_MTD))))
  
  output_df
}
