source("functions/KalmanRecEffic.R")
# source("functions/KalmanRegReform.R")
source("functions/SSlogLike.R")

KalmanReg <- function(y, X, Beta.Initial = NULL, a1 = 0, P1 = 10^2, sigma2.eps.initial = 1, sigma2.eta.initial = 1, logLikBurn = 3, P1B = NULL, trans = "square"){
  
  t1 <- system.time({
    opOut <- suppressWarnings(
      optim(
        par = sqrt(c(sigma2.eps.initial, sigma2.eta.initial)), 
        fn = logLik.fn, 
        y = y, X = X, Beta.Initial = Beta.Initial, a1 = a1, P1 = P1, logLikBurn = logLikBurn, P1B = P1B, trans = trans,
        method='L-BFGS-B', 
        hessian=T, 
        control=list(trace=0, REPORT=1)
      )
    )
  })
  
  
  # if(opOut$convergence != 0)warning("Liklihood did not converge for variance estimation")
  if(trans == "square") par <- opOut$par^2
  if(trans == "1/2 log") par <- exp(2*opOut$par)
  
  t2 <- system.time({
    kout <- KalmanRegEffec(y = y, X = X, sigma2.eps = par[1], sigma2.eta = par[2], a1 = a1, P1 = P1, P1B = P1B)
  })
  
  
  
  kout$convergence <- opOut$convergences
  kout$likelihoodCovergenceMessage <- opOut$message
  kout$par <- par
  kout$Optimization.Time <- t1
  kout$Filter.Time <- t2
  kout$counts <- opOut$counts
  return(kout)
}


logLik.fn <- function(y, X, Beta.Initial = Beta.Initial, par, a1 = 0, P1 = 10^8, logLikBurn = 3, P1B = NULL, trans = "square"){
  if(trans == "square") par <- par^2
  if(trans == "1/2 log") par <- exp(2*par)
  t <- ncol(y)
  logLiks <- KalmanRegEffec(
    y = y, X = X, 
    Beta.Initial = Beta.Initial, sigma2.eps = (par[1]), sigma2.eta = (par[2]), a1 = a1, P1 = P1, P1B = P1B, Lik = TRUE
  )
  if(sum(!is.finite(logLiks)) == length(logLiks)){
    return(10^8)
  }else{
    return(sum(logLiks[is.finite(logLiks)]))
  }
}



