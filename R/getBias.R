#' Calculate a bias of a randomization sequence in convergence strategy.
#' 
#' @param R randomization sequence
#' @param eta mean shift due to selection bias
#' 
#' @export 
getBias <- function(R,eta) {
  R_cor <- R[-length(R)]  
  rho <-   1-c(0.5,cumsum(R_cor)/1:(length(R_cor)))
  good <- round(rho,6) > round(0.5,6)
  weak <- round(rho,6) < round(0.5,6)
  return((good - weak)*eta)
}
