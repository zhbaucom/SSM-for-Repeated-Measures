source("functions/KalmanRegression2.R")
source("functions/KalmanRecReform2.R")
source("functions/ss.simp.sim.R")
source("functions/BayesKalmUneq2.R")
library(tidyverse)
library(lme4)
library(nlme)

NBSim <- function(dataType = "SSM", simXfun = NULL, ...){
  
  
  EstMat <- uclMat <- lclMat <- matrix(NA, p, 7)
  sigma2Mat <-  matrix(NA, 2, 7)
  
  
  
  
  
  #####Data Simulation
  if(is.null(simXfun)){
    X <- matrix(rnorm((p)*n), n, p)
  }else{
    X <- matrix(simXfun(p*n, ...), n, p)
  }
  
  
  
  colnames(X) <- paste("X", 1:ncol(X), sep = "")
  
  
  if(dataType == "SSM"){
    ytot <- ss.simp.sim(X, B, t = t, sigma2.eps = sigeps, sigma2.eta = sigeta, u0 = u0.true, P0 = P0.true) #Creates a simulated y
    
    y <- ytot
    
    colnames(y) <- paste("y", 1:ncol(y), sep = "")
    
    XY1 <- cbind(id = 1:n, y, X)
    
    XYs <- XY1 %>%
      data.frame() %>%
      gather("time", "y", colnames(.)[grepl("y", colnames(.))]) %>%
      mutate(time = as.numeric(substr(time, 2, nchar(time)))) 
    
    XYs[grepl("X", colnames(XYs))] <- XYs[grepl("X", colnames(XYs))] * XYs$time
    
    
    XY <- XYs %>%
      filter(!is.na(y)) %>%
      arrange(id, time) 
  }else if(dataType == "AR(1)"){
    
    XY <- lapply(1:t, function(x){
      xout <- cbind(1:n, x, (x*X) %*% B, x*X)
      colnames(xout) <- c("id","time", "XB", paste("X", 1:p, sep = ""))
      xout
    } ) %>%
      do.call("rbind", .) %>%
      as_tibble() %>%
      group_by(id) %>%
      mutate(error = arima.sim(n = n(), model = rhoS)) %>%
      mutate(y = XB + error) %>%
      select(-XB, -error) %>%
      ungroup(id) %>%
      arrange(id, time)
  }
  
  
  
  
  
  XY <- XY %>%
    group_by(id) %>%
    mutate(Keep = c(TRUE, 2:n() %in% sample(2:n(), sample(1:(length(id)-1), 1)))) %>%
    ungroup(id) %>%
    filter(Keep == TRUE) %>%
    arrange(id, time)
  
  
  
  ############FORMAT DATA FOR REGULAR KALMAN FILTER
  y.long <- XY$y
  X.long <- as.matrix(XY[,paste("X", 1:p, sep = "")])
  id <- XY$id
  time <- XY$time
  
  Beta.Initial <- coef(lm(y.long ~ X.long - 1))
  
  ##################LME#############################
  LMEform <- paste("y ~ (", paste(colnames(XY)[grepl("X", colnames(XY))], collapse = " + "), ")", sep = "")
  lmeout <- lme(as.formula(LMEform), random = ~1|id, data = XY, method = "ML")
  LMECI <- intervals(lmeout, which = "fixed")$fixed
  
  EstMat[,1] <- LMECI[-1,2]
  uclMat[,1] <- LMECI[-1,3]
  lclMat[,1] <- LMECI[-1,1]
  sigma2Mat[1, 1] <- lmeout$sigma
  
  
  #########################AR1##############################
  ARform <- paste("y ~ (", paste(colnames(XY)[grepl("X", colnames(XY))], collapse = " + "), ")", sep = "")
  ARmod <- lme(as.formula(ARform), random = ~1|id, correlation = corAR1(form = ~time), data = XY, method = "ML")
  CIout <- intervals(ARmod, which = "fixed")$fixed
  
  EstMat[,2] <- CIout[-1,2]
  uclMat[,2] <- CIout[-1,3]
  lclMat[,2] <- CIout[-1,1]
  sigma2Mat[1, 2] <- ARmod$sigma

  ##################### Run Original Kalman Filter #################
  okout <- KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7)
  okEs <- okout$alpha.hat[(n+1):(n+p),t]
  Vmat <- okout$V[[t]]
  ssV <- sqrt(diag(Vmat)[(nrow(Vmat)-p+1):nrow(Vmat)])
  UCL <-  okEs + qnorm(.975) * ssV 
  LCL <- okEs + qnorm(.025) * ssV
  
  EstMat[,3] <- okEs
  uclMat[,3] <- UCL
  lclMat[,3] <- LCL
  sigma2Mat[,3] <- c(okout$sigma2.eps, okout$sigma2.eta)
  
  
  ######################Bayesian Kalman##########################
  ######Initialize
  
  bkout <- BayesKalm.Uneq(
    y.long = y.long, X.long = X.long, id = id, time = time,
    Burn = Burn, Its = Its, 
    Beta.Initial = Beta.Initial, sigma2.beta = sigma2.beta, 
    u0 = u0, P0 = P0, 
    a0 = a0, b0 = b0, 
    c0 = c0, d0 = d0,
    silence = TRUE
  )
  
  CI <- apply(bkout$Beta[,Burn:Its], 1, quantile, c(.025, 0.975))
  
  EstMat[,4] <- rowMeans(bkout$Beta[,Burn:Its])
  uclMat[,4] <- CI[2,]
  lclMat[,4] <- CI[1,]
  sigma2Mat[,4] <- c(mean(bkout$sigma2.eps[Burn:Its]), mean(bkout$sigma2.eta[Burn:Its]))
  
  ###################### Run Original Kalman Filter ##################
  pkout2 <- KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7, k = 2)

  EstMat[,5] <- pkout2$summary[,1]
  uclMat[,5] <- pkout2$summary[,3]
  lclMat[,5] <- pkout2$summary[,2]
  sigma2Mat[,5] <- c(pkout2$sigma2.eps, pkout2$sigma2.eta)
  
  ############ Run Original Kalman Filter ###################
  pkout4 <- KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7, k = 4)

  EstMat[,6] <- pkout4$summary[,1]
  uclMat[,6] <- pkout4$summary[,3]
  lclMat[,6] <- pkout4$summary[,2]
  sigma2Mat[,6] <- c(pkout4$sigma2.eps, pkout4$sigma2.eta)
  
  
  
  ################### Run Original Kalman Filter ####################
  pkout10 <- KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7, k = 10)

  EstMat[,7] <- pkout10$summary[,1]
  uclMat[,7] <- pkout10$summary[,3]
  lclMat[,7] <- pkout10$summary[,2]
  sigma2Mat[,7] <- c(pkout10$sigma2.eps, pkout10$sigma2.eta)
  
  ############ Create Output ##########################
  colnames(EstMat) <- colnames(uclMat) <- colnames(lclMat) <- colnames(sigma2Mat) <-c("LME", "LME AR(1)", "SS Orig", "SS Bayes", "SS Part 2", "SS Part 4", "SS Part 10")
  
  list(Estimate = EstMat, UCL = uclMat, LCL = lclMat, sigma2 = sigma2Mat)
  
  
}

i <- 1
Compiler <- function(reps = 100, dataType, simXfun, ...){
  
  EstArray <- uclArray <- lclArray <- array(NA, dim = c(p, 7, reps))
  sigma2Array <- array(NA, dim = c(2, 7, reps))
  Fails <- 0
  for(i in 1:reps){
    starttime <- Sys.time()
    Again <- TRUE
    while(Again){
      # DSout <- NBdataSim()
      # list2env(DSout, envir = globalenv())
      Fcount <- 0
      SimOut <- try(NBSim(dataType = dataType, simXfun = simXfun, ...), silent = TRUE)
      if(class(SimOut) == "list"){ Again <- FALSE }else{ Fails <- Fails + 1; Fcount <- Fcount + 1; cat("__________________FAILURE__________________")}
      if(Fcount >= 5)break  
    }
    EstArray[,,i] <- SimOut$Estimate
    uclArray[,,i] <- SimOut$UCL
    lclArray[,,i] <- SimOut$LCL
    sigma2Array[,,i] <- SimOut$sigma2
   
    endtime <- Sys.time()
    
    print(paste0("Iteration ",i, sep = ""))
    print(endtime - starttime)

  }
  list(Estimate = EstArray, UCL = uclArray, LCL = lclArray, sigma2 = sigma2Array, Fails = Fails)
}













