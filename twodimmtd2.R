#'Function that calculates the two-dimensional MTD curve
#'
#' @param rho00 Parameter rho00 from marginal toxicity model
#' @param rho10 Parameter rho10 from marginal toxicity model
#' @param rho00 Parameter rho01 from marginal toxicity model
#' @param alpha3 Parameter alpha3 from marginal toxicity model
#' @param x Dose x

twodimmtd2 = function(rho00,rho01,rho10,alpha3,theta,x){
  alpha0<-qnorm(rho00)
  alpha1<-qnorm(rho10) - qnorm(rho00)
  alpha2<-qnorm(rho01) - qnorm(rho00)
  y = (qnorm(theta)-alpha0-alpha1*x)/(alpha2+alpha3*x)
  return(y)
}




