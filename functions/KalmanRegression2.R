source("functions/KalmanRecReform2.R")


KalmanReg2 <- function(y.long, X.long, id, time, k=1, Beta.Initial = NULL, a1 = 0, P1 = 1e7, sigma2.eps.initial = 1, sigma2.eta.initial = 1, P1B = NULL, trans = "square", maxit = 10){

  
  #Partition Data if groups are more than 1
  if(k > 1){
    uid <- unique(id)
    gs <- length(uid)/k
    gi <- sample(rep(1:k, gs))
    g <- 1
    kdat <- lapply(1:k, function(g){
      gid <- uid[gi == g]
      list(y.long = y.long[id %in% gid], X.long = X.long[id %in% gid,], id = id[id %in% gid], time = time[id %in% gid])
    })
    #Join likelihood calculation for each group
    likfn <- logLik.Part
  }else{
    #Optimization if not partitioned
    kdat <- NULL
    likfn <- logLik.fn2
  }
  
  
  
  
  
  ### Optimize variance parameters using L-BFGS-B optimization
  t1 <- system.time({
    opOut <- suppressWarnings(
      optim(
        par = sqrt(c(sigma2.eps.initial, sigma2.eta.initial)), 
        fn = likfn, 
        kdat = kdat, k = k, y.long = y.long, X.long = X.long, id = id, time = time, 
        Beta.Initial = Beta.Initial, a1 = a1, P1 = P1, 
        P1B = P1B, trans = trans,
        method='L-BFGS-B', 
        hessian=T, 
        control=list(trace=0, REPORT=1, maxit = maxit)
      )
    )
  })
  
  
  
  #transform parameters to go into kalman filter
  if(trans == "square") par <- opOut$par^2
  if(trans == "1/2 log") par <- exp(2*opOut$par)
  
  
  if(k > 1){
    # KF and KS done for each partition
    t2 <- system.time({
      g <- 1
      kout1 <- lapply(1:k, function(g){
        list2env(kdat[[g]], envir = environment())
        kn <- length(unique(id))
        KNS <- max(table(id))
        p <- ncol(X.long)
        koutg <- KalmanRegReform2(
          y.long = y.long, X.long = X.long, id = id, time = time, 
          Beta.Initial = Beta.Initial, sigma2.eps = par[1], sigma2.eta = par[2], 
          a1 = a1, P1 = P1, P1B = P1B, Lik = FALSE
        )
        list(alpha.hat = koutg$alpha.hat[(kn+1):(kn+p),KNS], V = koutg$V[[KNS]][(kn+1):(kn+p), (kn+1):(kn+p)])
      })
      
      # Gather results from each partition
      
      alpha.hat.List <- lapply(1:k, function(g){
        kout1[[g]]$alpha.hat
      })
      
      V.List <- lapply(1:k, function(g){
        kout1[[g]]$V
      })
      
      # Combine estimates for parameter summary
      kout <- list(
        alpha.hat = reduce(alpha.hat.List, `+`)/k,
        V = reduce(V.List, `+`)/k^2,
        sigma2.eps = par[1], sigma2.eta = par[2]
      )
      kout$summary <- cbind(
        estimate = kout$alpha.hat,
        LCL = kout$alpha.hat + qnorm(.025) * sqrt(diag(kout$V)),
        UCL = kout$alpha.hat + qnorm(.975) * sqrt(diag(kout$V))
      )

    })
  }else{
    # If not partioned, just rereun KF KS on the full group
    t2 <- system.time({
      kout <- KalmanRegReform2(y.long, X.long, id, time, sigma2.eps = par[1], sigma2.eta = par[2], a1 = a1, P1 = P1, P1B = P1B)
    })
  }

  
  
  # Include convergence output
  kout$convergence <- opOut$convergences
  kout$likelihoodCovergenceMessage <- opOut$message
  kout$par <- par
  kout$Optimization.Time <- t1
  kout$Filter.Time <- t2
  kout$counts <- opOut$counts
  return(kout)
}






#Optimization for the full likelihood KF KS
logLik.fn2 <- function(par, kdat = NULL, k = NULL, y.long, X.long, id, time, Beta.Initial = Beta.Initial, a1 = 0, P1 = 10^8, logLikBurn = 3, P1B = NULL, trans = "square"){
  if(trans == "square") par <- par^2
  if(trans == "1/2 log") par <- exp(2*par)

  logLiks <- KalmanRegReform2(
    y.long = y.long, X.long = X.long, id = id, time = time, 
    Beta.Initial = Beta.Initial, sigma2.eps = (par[1]), sigma2.eta = (par[2]), 
    a1 = a1, P1 = P1, P1B = P1B, Lik = TRUE
  )
  sum(logLiks[is.finite(logLiks)])
}






# Optimization function for the partioned SSM
logLik.Part <- function(par, kdat, y.long, X.long, id, time, k = k, Beta.Initial = Beta.Initial, a1 = 0, P1 = 10^8, logLikBurn = 3, P1B = NULL, trans = "square"){
  if(trans == "square") par <- par^2
  if(trans == "1/2 log") par <- exp(2*par)
  
  lkout <- unlist(lapply(1:k, function(g){
    list2env(kdat[[g]], envir = environment())
    KalmanRegReform2(
      y.long = y.long, X.long = X.long, id = id, time = time, 
      Beta.Initial = Beta.Initial, sigma2.eps = par[1], sigma2.eta = par[2], 
      a1 = a1, P1 = P1, P1B = P1B, Lik = TRUE
    )
  }))
  sum(lkout[is.finite(lkout)])
}





