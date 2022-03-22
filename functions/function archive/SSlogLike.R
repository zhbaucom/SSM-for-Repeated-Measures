logLik.fn <- function(y, X, Beta.Initial = Beta.Initial, par, a1 = 0, P1 = 10^8, logLikBurn = 3, P1B = NULL, trans = "square"){
  if(trans == "square") par <- par^2
  if(trans == "1/2 log") par <- exp(2*par)
  t <- ncol(y)
  logLiks <- KalmanRegEffec(y = y, X = X, Beta.Initial = Beta.Initial, sigma2.eps = (par[1]), sigma2.eta = (par[2]), a1 = a1, P1 = P1, P1B = P1B, Lik = TRUE)
  if(sum(!is.finite(logLiks)) == length(logLiks)){
    return(10^8)
  }else{
    return(sum(logLiks[is.finite(logLiks)]))
  }
}

