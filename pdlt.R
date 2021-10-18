#'Function that calculates marginal probability of DLT
#'
#' @param rho00 Parameter rho00 from marginal toxicity model
#' @param rho10 Parameter rho10 from marginal toxicity model
#' @param rho00 Parameter rho01 from marginal toxicity model
#' @param alpha3 Parameter alpha3 from marginal toxicity model


pdlt<-function(rho00,rho01,rho10,alpha3,x,y){
  alpha0<-qnorm(rho00)
  alpha1<-qnorm(rho10) - qnorm(rho00)
  alpha2<-qnorm(rho01) - qnorm(rho00)
  p = pnorm(alpha0+alpha1*x+alpha2*y+alpha3*x*y)
  return(p)
}