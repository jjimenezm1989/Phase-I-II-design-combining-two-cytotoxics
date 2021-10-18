#' Function to obtain dose combination to which patients were allocated
#' @param df Dataset containing outputs (data.frame)
#' @param modeling Type (character) of modeling approach (Independent modeling, Latent variables, FGM Copula, Clayton Copula) (character)
#' @param DLT_scenario Toxicity scenario (numeric)
#' @param EFF_scenario Efficacy scenario (numeric)


get_dose_combinations = function(df, DLT_scenario, EFF_scenario){
  
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         x = unlist(map(df,"X")),
                         y = unlist(map(df,"Y")))
  
  return(output_df)

}