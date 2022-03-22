
#Simulate outcome based on the simple state space model
ss.sim.Uneq <- function(X, B, times, u0 = 0, P0 = 1, sigma2.eps = 1, sigma2.eta = 1){
  tb <- lapply(times, diff)
  XB <- X%*%B
  n <- length(XB)
  #Initiate y and mu
  y.out <- lapply(1:n, function(x){
    T <- length(tb[[x]])+1
    y <- numeric(T)
    X.out <- matrix(NA, T, p)
    mu <- numeric(T)
    X.out[1,] <- X[x,]
    mu[1] <- rnorm(1, u0, sd = sqrt(P0)) + rnorm(1, 0, sd = sqrt(sigma2.eta))
    y[1] <- mu[1] + XB[x] + rnorm(1,0, sd = sqrt(sigma2.eps))
    
    
    #iteritively simulate y based on the model
    for(i in 2:T){
      X.out[i,] <- i*X[x,]
      mu[i] <- mu[i-1] + rnorm(1,0, sd = sqrt(tb[[x]][i-1]*sigma2.eta))
      y[i] <- mu[i]  + i * XB[x] + rnorm(1,0,sd = sqrt(sigma2.eps))
    }
    data.frame(id = x, y = y, X.out, time = times[[x]], mu = mu)
  })
  df.out <- do.call("rbind", y.out)
  colnames(df.out)[3:(2+p)] <- paste("X", 1:p, sep = "")
  df.out
}





