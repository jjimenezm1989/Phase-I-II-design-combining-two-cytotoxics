get_sample_size_stopping_rules = function(df){
  mean_sample_size_stopping_safety = mean(unlist(map(df, "sample_size_stopping_safety")))
  mean_sample_size_stopping_futility = mean(unlist(map(df, "sample_size_stopping_futility")))
  
  output = data.frame(futility = mean_sample_size_stopping_futility, 
                      safety = mean_sample_size_stopping_safety)
  
  return(output)
}