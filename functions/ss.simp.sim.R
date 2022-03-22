#Simulate outcome based on the simple state space model
ss.simp.sim <- function(X, B, t = 100, u0 = 0, P0 = 1, sigma2.eps = 1, sigma2.eta = 1){
  XB <- X%*%B
  n <- length(XB)
  #Initiate y and mu
  y <- matrix(NA, nrow = n, ncol = t)
  mu <- matrix(NA, nrow = n, ncol = t)
  mu[,1] <- rnorm(n, u0, sd = sqrt(P0)) + XB + rnorm(n, 0, sd = sqrt(sigma2.eta))
  y[,1] <- mu[,1] + rnorm(n,0, sd = sqrt(sigma2.eps))
  #iteritively simulate y based on the model
  for(i in 2:t){
    mu[,i] <- mu[,i-1]  + XB + rnorm(n,0, sd = sqrt(sigma2.eta))
    y[,i] <- mu[,i] + rnorm(n,0,sd = sqrt(sigma2.eps))
  }
  y
}


