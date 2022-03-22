source("functions/KalmanRegression2.R")
source("functions/KalmanRecReform2.R")
source("functions/ss.simp.sim.R")
source("functions/BayesKalmUneq2.R")
library(tidyverse)
library(lme4)
library(nlme)

#Sample wrapper 
mySample <- function(x){
  if(length(x) == 1) x else sample(x, 1)
}
#To simulate alpha of the
alpha.sim <- function(n = 100, t = 10, u0 = 0, P0 = 1, sigma2.eta = 1){
  mu <- matrix(NA, nrow = n, ncol = t)
  mu[,1] <- rnorm(n, u0, sd = sqrt(P0)) + rnorm(n, 0, sd = sqrt(sigma2.eta))
  #iteritively simulate y based on the model
  for(i in 2:t){
    mu[,i] <- mu[,i-1] + rnorm(n,0, sd = sqrt(sigma2.eta))
  }
  mu
}


NBSim <- function(dataType = "SSM", simXfun, ...){
  
  
  EstMat <- uclMat <- lclMat <- matrix(NA, p, 7)
  sigma2Mat <-  matrix(NA, 2, 7)
  
  #####Data Simulation
  ##### simulates comparable covariates to the NACC
  
  BaseLine <- data.frame(
    id = 1:n,
    X1 = rbinom(n, 1, prob = .25) *2 -1,
    X2 = rnorm(n),
    X3 = rbinom(n, 1, prob = .25) *2 -1,
    X4 = rnorm(n),
    X6 = rbinom(n, 1, prob = .25) *2 -1
  ) %>%
    mutate(X7 = X1*X6)
  
  X <- lapply(1:t, function(x){
    BaseLine$time <- x
    BaseLine
  }) %>%
    do.call("rbind",.) %>%
    group_by(id) %>%
    mutate(
      coff = sample(2:length(id), 1),
      X5 = (coff <= time) * 1
    ) %>%
    arrange(id) %>%
    select(-coff) %>%
    ungroup(id) 
  
  
  
  X[grepl("X", colnames(X))] <- X[grepl("X", colnames(X))] * X$time
  
  X <- X[c("id", "time", paste("X",1:p, sep = ""))]
  
  # The XB is the constant effect induced by the covariates.
  X$XB <- (as.matrix(X[grepl("X", colnames(X))]) %*% B)[,1]

  
  if(dataType == "SSM"){
    # If the data type is "SSM" we have a specific error structure
    alpha <- alpha.sim(n, t, sigma2.eta = sigeta)
    
    adf <- alpha %>%
      data.frame() %>%
      mutate(id = 1:nrow(.)) %>%
      pivot_longer(cols = starts_with("X"), names_to = "time", values_to = "alpha") %>%
      mutate(time = substr(time, 2, nchar(time)) %>% as.numeric())

    XY <- X %>%
      left_join(adf, by = c("id", "time"))
    
    # Create outcome with the alphas and random error
    XY$y <- XY$XB + XY$alpha + rnorm(n, 0, sd = sqrt(sigeps))
    
  }else if(dataType == "AR(1)"){
    # If data fallows an AR(1) data structure model
    XY <- X %>%
      group_by(id) %>%
      mutate(error = arima.sim(n = n(), model = rhoS)) %>%
      mutate(y = XB + error) %>%
      select(-error) %>%
      ungroup(id) %>%
      arrange(id, time)
  }
  
  
  
  
  # Randomly remove visits to unsure unequal visit times
  XY <- XY %>%
    group_by(id) %>%
    mutate(d1 = mySample(time[X5 > 0])) %>%
    mutate(Keep = c(TRUE, 2:n() %in% sample(2:n(), sample(1:(length(id)-1), 1)))) %>%
    mutate(Keep = ifelse(time == d1, TRUE, Keep)) %>%
    ungroup(id) %>%
    filter(Keep == TRUE) %>%
    select(-d1, - Keep, -XB)
  
  
  
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



# Gathers all simulation results and puts them into a more reportable format
Compiler <- function(reps = 100, dataType, simXfun, ...){
  
  EstArray <- uclArray <- lclArray <- array(NA, dim = c(p, 7, reps))
  sigma2Array <- array(NA, dim = c(2, 7, reps))
  Fails <- 0
  for(i in 1:reps){
    starttime <- Sys.time()
    Again <- TRUE
    while(Again){
      Fcount <- 0
      SimOut <- try(NBSim(dataType = dataType), silent = TRUE)
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













