
get_MTD_curve = function(df, DLT_scenario, EFF_scenario, true_rho00, true_rho01, true_rho10, true_alpha3){
  
  post_rho00 = mean(unlist(map(df, "post_rho00")))
  post_rho01 = mean(unlist(map(df, "post_rho01")))
  post_rho10 = mean(unlist(map(df, "post_rho10")))
  post_alpha3 = mean(unlist(map(df, "post_alpha3")))
  
  
  output_df = data.frame(DLT_scenario = DLT_scenario,
                         EFF_scenario = EFF_scenario,
                         x = seq(0,by = 0.01,1),
                         y = c(y_est = twodimmtd2(post_rho00,post_rho01,post_rho10,post_alpha3,theta_Z,seq(0,by = 0.01,1)),
                               y_true = twodimmtd2(true_rho00,true_rho01,true_rho10,true_alpha3,theta_Z,seq(0,by = 0.01,1))),
                         type = c(rep("Estimated MTD",101),rep("True MTD",101)))

  return(output_df)
  
}