#'Function that generates correlated binary outputs
#'
#' @param rho00 Parameter rho00 from marginal toxicity model
#' @param rho10 Parameter rho10 from marginal toxicity model
#' @param rho00 Parameter rho01 from marginal toxicity model
#' @param alpha3 Parameter alpha3 from marginal toxicity model
#' @param beta0 Parameter beta0 from marginal efficacy model
#' @param beta1 Parameter beta1 from marginal efficacy model
#' @param beta2 Parameter beta2 from marginal efficacy model
#' @param beta3 Parameter beta3 from marginal efficacy model
#' @param x Dose from compound X
#' @param y Dose from compound Y

generate_bivariate_outcome = function(rho00,rho01,rho10,alpha3,beta0,beta1,beta2,beta3,beta4,beta5,corr,x,y){
  
  alpha0<-qnorm(rho00)
  alpha1<-qnorm(rho10) - qnorm(rho00)
  alpha2<-qnorm(rho01) - qnorm(rho00)
  mu_Z  = alpha0+alpha1*x+alpha2*y+alpha3*x*y
  
  mu_E = beta0 + beta1*x + beta2*y + beta3*x*y + beta4*x^2 + beta5*y^2
  
  Sigma = matrix(c(1,corr,corr,1),2,2)
  
  output1 = rmvnorm(1, mean = c(mu_Z ,mu_E), Sigma)
  
  output = as.numeric(ifelse(output1 < 0, 0, 1))
  
  return(output)
  
}