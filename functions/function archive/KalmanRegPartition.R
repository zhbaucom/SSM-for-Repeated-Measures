source("functions/KalmanRegression.R")
# source("functions/stateSpaceSim.R")

m.det <- function(m){
  if(class(m) == "matrix")
    return(det(m))
  else
    return(m)
}

KalmanRegPart <- function(y, X, a1 = 0, Beta.Initial = NULL, P1 = 10^2, sigma2.eps.initial = 1, sigma2.eta.initial = 1, logLikBurn = 3, P1B = NULL, trans = "square", groups = 1, gs = NULL){

  n <- nrow(X)
  p <- ncol(X)
  t <- ncol(y)
  
  #Only need to supply group sizes or number of groups
  if(is.null(groups)) groups <- n/gs
  if(is.null(gs)) gs <- n/groups
  #Index for ramdomly partitioned data
  gi <- sample(rep(1:groups, gs))
  #output holders
  matrixEs <- matrix(NA, nrow = p, ncol = groups)
  arrayVar <- array(NA, dim = c(nrow = p, nrow = p, groups))
  varepsVec <- numeric(groups)
  varetaVec <- numeric(groups)
  
  for(gnum in 1:groups){
    kout <- KalmanReg(y[gi == gnum,], X[gi == gnum,], Beta.Initial = Beta.Initial, a1 = a1, P1 = P1, sigma2.eps.initial = sigma2.eps.initial, sigma2.eta.initial = sigma2.eta.initial, logLikBurn = logLikBurn, P1B = P1B, trans = trans)
    matrixEs[,gnum] <- kout$alpha.hat[(gs+1):(gs+p),t]
    arrayVar[,,gnum] <- kout$V[[t]][(gs+1):(gs+p),(gs+1):(gs+p)]
    varepsVec[gnum] <- kout$sigma2.eps
    varetaVec[gnum] <- kout$sigma2.eta
  }

  meanEstimate <- apply(matrixEs, 1, mean)
  varEstimate <- apply(arrayVar, 1:2, mean)/groups
  
  
  list(
    Estimate = meanEstimate, 
    Variance = varEstimate, 
    sigma2.eps = mean(varepsVec), 
    sigma2.eta = mean(varetaVec)
  )

}

