nits <- 1000

library(tidyverse)
library(lme4)
library(nlme)
source("functions/KalmanRegression2.R")
source("functions/KalmanRecReform2.R")
source("functions/ss.simp.sim.R")
source("functions/BayesKalmUneq2.R")

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


#####Bayesian Kalman Filter
#Beta ~ N(Beta.Initial, sigma2.beta)
Beta.Initial <- 0
# Beta.Initial <- B
sigma2.beta <- 10
#sigma2.eta ~ IG(a0/2, b0/2)
a0 <- .01
b0 <- .01
#sigma2.ps ~ IG(c0/2, d0/2)
c0 <- .01
d0 <- .01
#mu0 ~ N(u0, P0)
u0 <- 0
P0 <- 10
########### Select Iterations and Burn ########### 
Its <- 2000
Burn <- ceiling(Its/2)

B <- readRDS("NACC/B.RDS") #Beta Vector
B <- c(B[-1], B[1])
n <- 100
t <- 10
p <- length(B)

sigeps <- 3
sigeta <- 0
dataType <- "SSM"
rhoS <-list(order = c(1, 0, 0), ar = 0, sd = 1)





bkList <- list()
for(i in 1:nits){
  BaseLine <- data.frame(
    id = 1:n,
    X1 = rbinom(n, 1, prob = .25),
    X2 = rnorm(n),
    X3 = rbinom(n, 1, prob = .25),
    X4 = rnorm(n),
    X6 = rbinom(n, 1, prob = .25)
  ) %>%
    mutate(X7 = X1*X6)
  
  X <- lapply(1:t, function(x){
    BaseLine$time <- x
    BaseLine
  }) %>%
    do.call("rbind",.) %>%
    group_by(id) %>%
    mutate(
      coff = sample(2:(length(id)-1), 1),
      X5 = (coff <= time) * 1
    ) %>%
    arrange(id) %>%
    select(-coff) %>%
    ungroup(id) 
  
  
  
  X[grepl("X", colnames(X))] <- X[grepl("X", colnames(X))] * X$time
  X <- X %>%
    group_by(id) %>%
    mutate(X5 = ifelse(X5 > 0, X5 - min(X5[X5 > 0]), X5)) %>%
    ungroup()
  
  X$X8 <- X$time
  X <- X[c("id", "time", paste("X",1:p, sep = ""))]
  
  
  X$XB <- (as.matrix(X[grepl("X", colnames(X))]) %*% B)[,1]
  
  if(dataType == "SSM"){
    alpha <- alpha.sim(n, t, sigma2.eta = sigeta)
    
    adf <- alpha %>%
      data.frame() %>%
      mutate(id = 1:nrow(.)) %>%
      pivot_longer(cols = starts_with("X"), names_to = "time", values_to = "alpha") %>%
      mutate(time = substr(time, 2, nchar(time)) %>% as.numeric())
    
    XY <- X %>%
      left_join(adf, by = c("id", "time"))
    XY$y <- XY$XB + XY$alpha + rnorm(n, 0, sd = sqrt(sigeps))
    
  }else if(dataType == "AR(1)"){
    
    XY <- X %>%
      group_by(id) %>%
      mutate(error = arima.sim(n = n(), model = rhoS)) %>%
      mutate(y = XB + error) %>%
      select(-error) %>%
      ungroup(id) %>%
      arrange(id, time)
  }
  
  
  
  
  XY <- XY %>%
    group_by(id) %>%
    mutate(
      d1 = ifelse(length(time[X5>0]) == 0, max(time), ifelse(length(time[X5>0]) == 1, time[X5>0], sample(time[X5>0], 1)))
    ) %>%
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
  
  bkList[[i]] <- BayesKalm.Uneq(
    y.long = y.long, X.long = X.long, id = id, time = time,
    Burn = Burn, Its = Its, 
    Beta.Initial = Beta.Initial, sigma2.beta = sigma2.beta, 
    u0 = u0, P0 = P0, 
    a0 = a0, b0 = b0, 
    c0 = c0, d0 = d0,
    silence = TRUE
  )
  if((i %% 10) == 0) print(i)
}


saveRDS(bkList, "Simulations/bkList.RDS")

# readRDS("Simulations/bkList.RDS") %>%
#   map("Beta") %>%
#   map(~.x[,Burn:Its]) %>%
#   map(~apply(.x, 1, mean, na.rm = TRUE))
