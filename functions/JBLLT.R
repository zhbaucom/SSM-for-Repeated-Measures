library(mvtnorm)
library(tidyverse)
source("functions/BayesKalmUneq2.R")
source("functions/ffbsJoint.R")

# BayesKalmJointOutcomeCor(
#   data = ldfun, outcomes = c("TRAILA", "TRAILB"), predictors = c("time", "RACEWHITE","SEX", "APOE", "APOESEX", "EDUC", "DEC", "AgeBase", "Group"),
#   timevar = "time", id = "NACCID", numits = 10, burn = 5, numitInit = 10, burnInit = 5
# )
# data = ldfun; outcomes = c("TRAILA", "TRAILB"); predictors = c("time", "RACEWHITE","SEX", "APOE", "APOESEX", "EDUC", "DEC", "AgeBase", "Group");
# timevar = "time"; id = "NACCID"; numits = 10; burn = 5; numitInit = 10; burnInit = 5; silence = FALSE; seed = NULL
seed = NULL

BayesJointInitializer <- function(data, outcomes, predictors, timevar, id, numits = 1000, burn = 500, silence = FALSE, seed = NULL, numitInit = 1500, burnInit = 500, seed.init = NULL){
  if(!is.null(seed.init))set.seed(seed.init)
  sigma2.beta <- 20
  BayesInit <- map(outcomes, function(outcome){
    
    initmodel <- lm(as.formula(paste(outcome, "~", paste(predictors, collapse = "+"))), data = data)
    
    data$resid <- resid(initmodel)
    
    Beta.Initial <- coef(initmodel)[-1]
    
    # u0 <- coef(initmodel)[1]
    # 
    # P0 <- data %>%
    #   group_by(vars(id)) %>%
    #   filter(time == min(time)) %>%
    #   .[["resid"]] %>%
    #   var()
    # 
    
    bkout <- BayesKalm.Uneq(
      y.long = data[[outcome]], X.long = as.matrix(data[predictors]), id = data[[id]], time = data[[timevar]],
      Burn = burnInit, Its = numitInit, 
      Beta.Initial = Beta.Initial, sigma2.beta = sigma2.beta, 
      u0 = 0, P0 = 20, 
      a0 = 0.01, b0 = 0.01, c0 = 0.01, d0 = 0.01, 
      silence = silence
    )
    cat("outcome ", outcome, " initialization complete...\n")
    list(
      Beta.Initial = apply(bkout$Beta[,burnInit:numitInit], 1, mean), 
      sigma2.eta = mean(bkout$sigma2.eta[burnInit:numitInit]),
      sigma2.eps = mean(bkout$sigma2.eps[burnInit:numitInit]),
      u0 = 0,
      P0 = 20,
      BetaIts = bkout$Beta,
      sigepsIts = bkout$sigma2.eps,
      sigetaIts = bkout$sigma2.eta
    )
    
    
  })
  
  init = list(
    BayesInit = BayesInit,
    
    args = list(
      data = data, 
      outcomes = outcomes, 
      predictors = predictors, 
      timevar = timevar, 
      id = id, 
      numits = numits, 
      burn = burn, 
      silence = silence, 
      seed = seed 
    )
  )
  init
  

}





BayesJointOnly <- function(init, method = "s"){
  
  BayesInit <- init$BayesInit
  list2env(init$args, environment())
  
  # PreSpecify
  
  cs <- length(outcomes)

  
  # DataFormat ------------------------------------------------------------
  #Nest data for each id
  
  data[["time"]] <- data[[timevar]]
  data[["id"]] <- data[[id]]
  n <- length(unique(data$id))
  ndat <- data %>%
    group_by(id)%>%
    mutate(timeDiff = c(diff(time), NA)) %>%
    nest() %>%
    .[["data"]] 
  ##Define matrix of output for each id
  y <- ndat %>%
    map(~select(.x, all_of(outcomes))) %>%
    map(~t(as.matrix(.x)))
  ##Define matrix of predictors for each id
  X <- ndat %>%
    map(~select(.x, all_of(predictors))) %>%
    map(~(as.matrix(.x)))
  
  CPbig <- lapply(1:n, function(i){
    lapply(1:nrow(X[[i]]), function(x)outer(X[[i]][x,],X[[i]][x,])) %>%
      Reduce("+", .)
  }) %>%
    Reduce("+", .)
  ##create vector of time and timediff for each id
  time <- ndat %>%
    map("time")
  
  timeDiff <- ndat %>%
    map("timeDiff") %>%
    map(~.x[!is.na(.x)])
  
  
  
  
  sigma2.beta <- 20
  ## Transform for quicker inverse calculation in KF
  XtX <- data %>%
    select(all_of(predictors)) %>%
    as.matrix() %>%
    crossprod()
  
  sig2beta_XtX <- sigma2.beta * XtX
  
  ex <- eigen(sig2beta_XtX, symm = TRUE)
  
  #####################
  
  ## Prior specfication ------------------------------------------------------
  
  ##Specify number of state differences and total number of observations
  Neta <- (map_dbl(timeDiff, length) %>% sum())
  Ntotal <-  nrow(data)
  
  
  
  W = BayesInit %>%
    map_dbl("sigma2.eta") %>%
    diag()
  
  V = BayesInit %>%
    map_dbl("sigma2.eps") %>%
    diag()
  
  m0 <- BayesInit %>%
    map_dbl("u0") 
  
  C0 <- BayesInit %>%
    map_dbl("P0") %>%
    diag()
  ###Need another way to set "TT", not ncol (y)
  prior.nu.eta <- cs - 1 + 0.01
  ## prior.nu.eta <- cs + 1
  prior.Gamma.eta <- diag(0.01, cs)
  ## prior.Gamma.eta <- diag(1, cs)
  
  a0 <- 0.01
  b0 <- 0.01
  c0 <- 0.01
  d0 <- 0.01
  
  ### Beta initialize
  B.star <- BayesInit %>%
    map("Beta.Initial") %>%
    do.call("cbind", .)
  B2 <- matrix(B.star, ncol = cs)
  Beta.Initial <- numeric(length(B.star))
  
  
  
  
  ## Tracking obj
  Beta.Track <- array(NA, dim = c(ncol(X[[1]]), cs, numits))
  vcovWish <- wcovWish <- array(NA, dim = c(cs, cs, numits))
  likTrack <- numeric(numits)
  
  It <- 0
  It <- It+1
  sigZeros <- matrix(0, nrow = nrow(B2), ncol = nrow(B2))
  
  bigEmptyX <- matrix(0, nrow = cs, ncol = length(predictors) *cs)
  BigEmptySig <- matrix(0, nrow = length(predictors) * cs, ncol = length(predictors) * cs)
  
  if(!is.null(seed))set.seed(seed)
  cat("Fitting Model...\n")
  if(!silence) pb <- txtProgressBar(min = 0, max = numits, style = 3)
  for(It in 1:numits){
    y.star <- lapply(1:n, function(i){
      y[[i]]-t(X[[i]] %*% B.star)
    })
    ####### alpha
    i <- 3
    bout <- lapply(1:n, function(i){
      
      ffbs.joint(y = y.star[[i]], V = V,W = W,m0 = m0,C0 = C0, timeDiff = timeDiff[[i]])$x
    })
    
    if(method == "os"){
      ####### Sigma eta
      fp <- lapply(1:n, function(i){
        apply(bout[[i]], 1, function(z)diff(z)) / sqrt( timeDiff[[i]] )
        
      }) %>%
        do.call("rbind", .) %>%
        crossprod()
      
      
      # W <- cIRT::riwishart(prior.nu.eta + Neta, fp + prior.Gamma.eta)
      W <- cIRT::riwishart(Neta, fp)
      wcovWish[,,It] <- W
      
      ###### Sigma eps
      
      or <- lapply(1:n, function(i){
        ((y.star[[i]] - bout[[i]]))
      })  %>%
        do.call("cbind", .) %>%
        tcrossprod()
      
      
      # V <- vcovWish[,,It] <- cIRT::riwishart(prior.nu.eta + Ntotal, or + prior.Gamma.eta)
      V <- vcovWish[,,It] <- cIRT::riwishart(Ntotal, or)
      
      
      
      
      # #### Beta Post
      
      
      
      SIGepsInv <- solve(V)
      daSig <- kronecker(SIGepsInv, CPbig)
      diag(daSig) <- diag(daSig) + 1/sigma2.beta
      
      daSiginv <- solve(daSig)
      
      daB <- lapply(1:n, function(i){
        yas <- crossprod(y[[i]] - bout[[i]], SIGepsInv)
        lapply(1:ncol(yas), function(x)yas[,x]*X[[i]]) %>%
          do.call("cbind", .)
      })  %>%
        do.call("rbind", .) %>%
        colSums()
      
      
      BetaSim <- mvtnorm::rmvnorm(1, daB %*% daSiginv, sigma = daSiginv)
      
      
      B.star <- Beta.Track[,,It] <- matrix(BetaSim, ncol = cs)
      
      # likelihood Tracker
      # BigSig <- daSiginv
      # 
      # likTrack[It] <- lapply(1:n, function(i){
      #   
      #   
      #   
      #   ycent <- y.star[[i]] - bout[[i]]
      #   
      #   xiter <- X[[i]]
      #   itt <- nrow(xiter)
      #   
      #   ittrack <- numeric(cs)
      #   for(k in 1:itt){
      #     bigX <- bigEmptyX
      #     for(j in 1:cs)
      #       bigX[j,(j-1)*length(predictors) + (1:length(predictors))] <- xiter[k,]
      #     
      #     # bigX <- cbind(rbind(xiter[k,], 0), rbind(0, xiter[k,]))
      #     BetaVar <- tcrossprod(bigX %*% daSiginv, bigX)
      #     ittrack <- ittrack + dmvnorm(ycent[,k], sigma = BetaVar + W + V)
      #   }
      #   
      #   ittrack
      #   
      # }) %>% unlist() %>% sum()
      
    }
    
    if(method == "o"){
      ####### Sigma eta
      fp <- lapply(1:n, function(i){
        apply(bout[[i]], 1, function(z)diff(z)^2) / timeDiff[[i]]
        
      }) %>%
        do.call("rbind", .) %>%
        colSums()
      
      W <- diag(sapply(fp, function(x){
        1/rgamma(1, (Neta/ 2 + a0), (b0 + x/ 2))
      }))
      
      wcovWish[,,It] <- W
      
      ###### Sigma eps
      
      or <- lapply(1:n, function(i){
        ((y.star[[i]] - bout[[i]]))
      })  %>%
        do.call("cbind", .) %>%
        tcrossprod()
      
      
      # V <- vcovWish[,,It] <- cIRT::riwishart(prior.nu.eta + Ntotal, or + prior.Gamma.eta)
      V <- vcovWish[,,It] <- cIRT::riwishart(Ntotal, or)
      
      
      
      # #### Beta Post
      
      SIGepsInv <- solve(V)
      daSig <- kronecker(SIGepsInv, CPbig)
      diag(daSig) <- diag(daSig) + 1/sigma2.beta
      
      daSiginv <- solve(daSig)
      
      daB <- lapply(1:n, function(i){
        yas <- crossprod(y[[i]] - bout[[i]], SIGepsInv)
        lapply(1:ncol(yas), function(x)yas[,x]*X[[i]]) %>%
          do.call("cbind", .)
      })  %>%
        do.call("rbind", .) %>%
        colSums()
      
      
      BetaSim <- mvtnorm::rmvnorm(1, daB %*% daSiginv, sigma = daSiginv)
      
      
      B.star <- Beta.Track[,,It] <- matrix(BetaSim, ncol = cs)
      
      # likelihood Tracker
      # BigSig <- daSiginv
      # 
      # likTrack[It] <- lapply(1:n, function(i){
      #   
      #   
      #   
      #   ycent <- y.star[[i]] - bout[[i]]
      #   
      #   xiter <- X[[i]]
      #   itt <- nrow(xiter)
      #   
      #   ittrack <- numeric(cs)
      #   for(k in 1:itt){
      #     bigX <- bigEmptyX
      #     for(j in 1:cs)
      #       bigX[j,(j-1)*length(predictors) + (1:length(predictors))] <- xiter[k,]
      #     
      #     # bigX <- cbind(rbind(xiter[k,], 0), rbind(0, xiter[k,]))
      #     BetaVar <- tcrossprod(bigX %*% daSiginv, bigX)
      #     ittrack <- ittrack + dmvnorm(ycent[,k], sigma = BetaVar + W + V)
      #   }
      #   
      #   ittrack
      #   
      # }) %>% unlist() %>% sum()
    }
    
    if(method == "s"){
      ####### Sigma eta
      fp <- lapply(1:n, function(i){
        apply(bout[[i]], 1, function(z)diff(z)) / sqrt( timeDiff[[i]] )
        
      }) %>%
        do.call("rbind", .) %>%
        crossprod()
      
      
      # W <- cIRT::riwishart(prior.nu.eta + Neta, fp + prior.Gamma.eta)
      W <- cIRT::riwishart(Neta, fp)
      wcovWish[,,It] <- W
      
      ###### Sigma eps
      
      or <- lapply(1:n, function(i){
        (y.star[[i]] - bout[[i]])^2
      })  %>%
        do.call("cbind", .) %>%
        apply(1, function(x)sum((x)))
      
      
      sigma2.eps.star <- sapply(or, function(x)1/rgamma(1, ((Ntotal)/2 +c0), d0+x/2))
      V <- vcovWish[,,It] <-sigma2.eps.star * diag(cs)
      
      
      
      
      # #### Beta Post
      
      
      v.star <- lapply(1:n, function(i){
        (y[[i]]-bout[[i]])
      })
      
      B.sum <- lapply(1:n, function(i){
        (v.star[[i]] %*% X[[i]])
      }) %>%
        Reduce("+", .)
      
      # B.Big <- sigma2.beta * B.sum - tcrossprod(V, Beta.Initial)
      
      B.Big <- sigma2.beta * B.sum
      
      Sigma.Inv <- lapply(sigma2.eps.star, function(eppps)tcrossprod(ex$vectors/(ex$values + eppps)[col(ex$vectors)], ex$vectors))
      
      BetaSim <- lapply(1:cs, function(i){
        mvtnorm::rmvnorm(1, mean = crossprod(Sigma.Inv[[i]], B.Big[i,]), sigma = Sigma.Inv[[i]] * sigma2.eps.star[i] * sigma2.beta)
      }) %>%
        do.call("rbind", .) %>%
        t()
      
      
      B.star <- Beta.Track[,,It] <- BetaSim
      
      # likelihood Tracker
      
      # BigSig <- cbind(rbind(Sigma.Inv[[1]], sigZeros), rbind(sigZeros, Sigma.Inv[[2]]))
      # BigSig <- BigEmptySig
      # for(j in 1:cs)
      #   BigSig[(j-1)*cs + (1:length(predictors)),(j-1)*cs + (1:length(predictors))] <- Sigma.Inv[[j]]
      # 
      # 
      # 
      # likTrack[It] <- lapply(1:n, function(i){
      #   
      #   ycent <- y.star[[i]] - bout[[i]]
      #   
      #   xiter <- X[[i]]
      #   itt <- nrow(xiter)
      #   
      #   ittrack <- numeric(cs)
      #   for(k in 1:itt){
      #     bigX <- bigEmptyX
      #     for(j in 1:cs)
      #       bigX[j,(j-1)*length(predictors) + (1:length(predictors))] <- xiter[k,]
      #     # bigX <- cbind(rbind(xiter[i,], 0), rbind(0, xiter[i,]))
      #     BetaVar <- tcrossprod(bigX %*% BigSig, bigX)
      #     ittrack <- ittrack + dmvnorm(ycent[,k], sigma = BetaVar + W + V)
      #   }
      #   
      #   ittrack
      #   
      # }) %>% unlist() %>% sum()
      
    }
    
    
    if(!silence) setTxtProgressBar(pb, It)
    
  }
  
  list(Beta.Track = Beta.Track, Eps.Track = vcovWish, Eta.Track = wcovWish, Initial = BayesInit, Lik = likTrack)
  
}



JBLLT <- function(method = "os", data, outcomes, predictors, timevar, id, initialization = "Bayes", numits = 1000, burn = 500, silence = FALSE, seed = NULL, numitInit = 1500, burnInit = 500, seed.init = NULL){
  init = BayesJointInitializer(data = data, 
                               outcomes = outcomes, 
                               predictors = predictors, 
                               timevar = timevar, 
                               id = id, 
                               numits = numits, 
                               burn = burn, 
                               silence = silence, 
                               seed = seed, 
                               numitInit = numitInit, 
                               burnInit = burnInit, 
                               seed.init = seed.init)
  
  BayesJointOnly(init = init, method = method)
}






