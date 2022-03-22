m.det <- function(m){
  if(class(m)[1] == "matrix")
    return(det(m))
  else
    return(m)
}

m.inv <- function(data, iteration) MASS::ginv(data)



KalmanRegReform2 <- function(y.long, X.long, id, time, Beta.Initial = NULL, sigma2.eps = 1, sigma2.eta = 1, a1 = 0, P1 = 1e7, P1B = NULL, Lik = FALSE){
  # logVar <- c(0.1, -20.9)
  # sigma2.eps = exp(logVar[1])
  # sigma2.eta = exp(logVar[2])
  i <- 1
  #Initialize vectors
  p <- ncol(X.long)
  n <- length(unique(id))
  n.totobs <- length(y.long)
  NS <- max(table(id))
  
  
  
  ###CONVERT FROM LONG TO SHORT DATA
  #Configure the outcome
  yl <- lapply(unique(id), function(x){
    ys <- y.long[id == x]
    c(ys, rep(NA, NS - length(ys)))
  })
  y <- do.call("rbind", yl)
  isnay <- y ^ is.na(y)
  #Configure time   
  tl <- lapply(unique(id), function(x){
    t2 <- time[id == x]
    c(t2, rep(NA, NS - length(t2)))
  })
  timeMat <- do.call("rbind", tl)
  timeDiff <- rbind(1, apply(timeMat, 1, diff))
  timeDiff[is.na(timeDiff)] <- 0
  
  NON.MISSING <- !is.na(y)
  #Configure X
  
  xl <- lapply(unique(id), function(x){
    xs <- X.long[id == x,]
    rbind(xs, matrix(0, nrow = NS-nrow(xs), ncol = p))
  })
  
  
  X <- array(NA, dim = c(n, p, NS))
  for(i in 1:length(xl))
    X[i,,] <- t(xl[[i]])
  
  
  
  
  
  # Initial KF components
  a <- matrix(NA, nrow = n+p, ncol = NS)
  Zt <- lapply(1:NS, function(x) cbind(diag(1, n, n),X[,,x]))
  v <- vector("list", NS)
  F <- F.inv <- list()
  P <- vector("list", NS)
  K <- vector("list", NS)
  L <- vector("list", NS)
  
  
  T <- diag(1, n+p, n+p)
  H <- diag(sigma2.eps, n, n)
  
  #Q is in two steps
  Qsub <- diag(0, n+p, n+p)
  Q <- lapply(1:NS, function(x){
    
    diag(append(sigma2.eta * timeDiff[x,], rep(0, p)))
    
  })
  

  SSlogLik <- rep(NA, NS)
  
  #Set starting values
  if(length(a1) == 1) a1s <- rep(a1, n) else a1s <- a1
  if(!is.null(Beta.Initial)) a[,1] <- c(a1s,Beta.Initial) else a[,1] <- c(a1s, rep(0, p))
  
  
  #Can completely specify P1 or just the diagnol
  if(class(P1) == "matrix"){
    P[1] <- P1
  }else{
    P[[1]] <- diag(P1, n+p, n+p)
    
  }
  
  Wb <- diag(rep(1,n))
  Wt <- Wb[!is.na(y[,1]),]
  y.hold <- ifelse(is.na(y[,1]), 0, y[,1])
  y.Star <- Wt %*% y.hold
  Z.Star <- Wt %*% Zt[[1]]
  H.Star <- tcrossprod(Wt %*% H, Wt) 
  
  #First step of KF
  #Set other starting values
  v[[1]] <- y.Star-Z.Star %*% a[,1]
  F[[1]] <- Z.Star %*% tcrossprod(P[[1]],Z.Star) + H.Star
  F.inv[[1]] <- m.inv(F[[1]], i)
  K[[1]] <- tcrossprod(P[[1]], Z.Star) %*% F.inv[[1]]
  L[[1]] <- T - K[[1]] %*% Z.Star
  
  SSlogLik[1] <- (log(m.det(F[[1]])) + crossprod(v[[1]], F.inv[[1]] %*% v[[1]]))
  

  i <- 2
  CalcP <- TRUE
  thresh <- .001
  ConvergedIt <- NA
  change <- FALSE
  for(i in 2:NS){
    #Accomadate missing
    Wt <- Wb[!is.na(y[,i]),]
    y.hold <- ifelse(is.na(y[,i]), 0, y[,i])
    y.Star <- Wt %*% y.hold
    Z.Star <- Wt %*% Zt[[i]]
    H.Star <- tcrossprod(Wt %*% H, Wt) 
    
    #Accomadating a single subject
    if(n == 1){
      a[,i] <- a[,i-1]  + K[[i-1]] * v[[i-1]]
    }else{
      a[,i] <- a[,i-1]  +  K[[i-1]] %*% v[[i-1]]
    }
    
    v[[i]] <- y.Star-Z.Star %*% a[,i]
    
    #P and F aren't recalculated if the variance of the betas converge
    if(CalcP){
      P[[i]] <- tcrossprod(P[[i-1]],L[[i-1]]) + Q[[i]]
      F[[i]] <- Z.Star %*% tcrossprod(P[[i]],Z.Star) + H.Star

      F.inv[[i]] <- m.inv(F[[i]], i)

      K[[i]] <- tcrossprod(P[[i]], Z.Star) %*% F.inv[[i]]
      L[[i]] <- T - K[[i]] %*% Z.Star
      
    }
    SSlogLik[i] <- (log(m.det(F[[i]])) + crossprod(v[[i]], F.inv[[i]] %*% v[[i]]))
  }
  
  if(Lik){
    lis <- SSlogLik
  }else{
    r <- matrix(NA, nrow = n+p, ncol = NS+1)
    alpha.hat <-  matrix(NA, nrow = n+p, ncol = NS)
    N <- list()
    V <- list()
    
    r[,NS+1] <- 0
    N[[NS+1]] <- matrix(0, nrow = n+p, ncol = n+p)
    
    k <- NS
    
    for(k in NS:1){
      Wt <- Wb[!is.na(y[,k]),]
      y.hold <- ifelse(is.na(y[,k]), 0, y[,k])
      y.Star <- Wt %*% y.hold
      Z.Star <- Wt %*% Zt[[k]]
      H.Star <- tcrossprod(Wt %*% H, Wt) 
      #E(alpha|y1, ..., yn)
      r[,k] <- crossprod(Z.Star, F.inv[[k]] %*% v[[k]]) + crossprod(L[[k]], r[,k+1])
      alpha.hat[,k] <- a[,k] + P[[k]] %*% r[,k]
      #Var(alpha|y1, ..., yn)
      N[[k]] <- crossprod(Z.Star, F.inv[[k]] %*% Z.Star) + crossprod(L[[k]], N[[k+1]] %*% L[[k]])
      V[[k]] <- P[[k]] - P[[k]] %*% N[[k]] %*% P[[k]]
    }

    
    lis <- list(SSlogLik = SSlogLik, alpha.hat = alpha.hat, V = V, a = a, P = P, F = F, F.inv = F.inv, v = v, K = K, L = L, ConvergedIt = ConvergedIt, sigma2.eps = sigma2.eps, sigma2.eta = sigma2.eta, P1B = P1B)
    
  }
  
  
  
  return(lis)
}




