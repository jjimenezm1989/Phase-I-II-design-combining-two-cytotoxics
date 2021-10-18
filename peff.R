#'Function that calculates marginal probability of efficacy
#'
#' @param beta0 Parameter beta0 from marginal efficacy model
#' @param beta1 Parameter beta1 from marginal efficacy model
#' @param beta2 Parameter beta2 from marginal efficacy model
#' @param beta3 Parameter beta3 from marginal efficacy model


peff = function(beta0,beta1,beta2,beta3,beta4,beta5,x,y){
  
  p = pnorm(beta0 + beta1*x + beta2*y + beta3*x*y + beta4*x^2 + beta5*y^2)
  
  return(p)
  
}