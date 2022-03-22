#for batch job
wo <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(wo)) wo <- 1



nobs1 <- c(20, 50, 100, 150, 200, 250, 300, 500, 750, 1000)
nobs2 <- c(50, 100, 150, 200, 250, 300, 500, 750, 1000)

# refData <- data.frame(
#   # l = rep(c("FL", "GS5", "GS10", "GS50", "Bayes"), length(nobs)),
#   # l = rep(c("GS50"), length(nobs)),
#   l = c(rep("GS100", 10), rep("GS250", 4)),
#   # nobs = sort(rep(nobs, 4))
#   nobs = c((1:10)*100, 250 * (1:4))
# )

refData <- rbind(
  data.frame(
    l = rep(c("FL", "GS5", "GS10", "Bayes"), length(nobs)),
    nobs = sort(rep(nobs1, 4))
  ),
  data.frame(
    l = "GS50",
    nobs = nobs2
  ),
  data.frame(
    l = c(rep("GS100", 10), rep("GS250", 4)),
    nobs = c((1:10)*100, 250 * (1:4))
  )
)





library(tidyverse)
source("functions/KalmanRegression2.R")
source("functions/KalmanRecReform2.R")
source("functions/ss.simp.sim.R")
source("functions/BayesKalmUneq2.R")

times <- 1000

B <- c(1, 1)

u0.true <- 0 #True mean of mu_0
P0.true <- 1 #True variance of mu_0

t <- 5
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
###True Params
sigeps <- 8
sigeta <- 3







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






DataCR <- function(){
  #####Data Simulation
  BaseLine <- data.frame(
    id = 1:n,
    X1 = rnorm(n),
    X2 = rnorm(n)
  )
  
  X <- lapply(1:t, function(x){
    BaseLine$time <- x
    BaseLine
  }) %>%
    do.call("rbind",.)
  
  
  
  X[grepl("X", colnames(X))] <- X[grepl("X", colnames(X))] * X$time
  
  X <- X[c("id", "time", paste("X",1:2, sep = ""))]
  
  X$XB <- (as.matrix(X[grepl("X", colnames(X))]) %*% B)[,1]
  
  
  alpha <- alpha.sim(n, t, sigma2.eta = sigeta)
  
  adf <- alpha %>%
    data.frame() %>%
    mutate(id = 1:nrow(.)) %>%
    pivot_longer(cols = starts_with("X"), names_to = "time", values_to = "alpha") %>%
    mutate(time = substr(time, 2, nchar(time)) %>% as.numeric())
  
  XY <- X %>%
    left_join(adf, by = c("id", "time"))
  XY$y <- XY$XB + XY$alpha + rnorm(n, 0, sd = sqrt(sigeps))
  
  
  y.long <- XY$y
  X.long <- as.matrix(XY[,paste("X", 1:2, sep = "")])
  id <- XY$id
  time <- XY$time
  
  Beta.Initial <- coef(lm(y.long ~ X.long))[-1]
  list(y.long = y.long, X.long = X.long, id = id, time = time, Beta.Initial = Beta.Initial)
  
}




l <- unique(refData$l)[wo]
nall <- refData$nobs[refData$l == l]
n <- nall[1]
for(n in nall){
  set.seed(08312021)
  SimTime <- numeric(length = times)
  for(i in 1:times){
    
    
    tryCatch(
      {

        list2env(DataCR(), envir = environment())
        SimTime[i] <- system.time({
          if(l == "FL"){
            KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7)
          }else if(l == "GS5"){
            KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7, k = n/5)
          }else if(l == "GS10"){
            KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7, k = n/10)
          }else if(l == "GS50"){
            KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7, k = n/50)
          }else if(l == "GS100"){
            KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7, k = n/100)
          }else if(l == "GS250"){
            KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial, a1 = u0, P1 = 1e7, k = n/250)
          }else if(l == "Bayes"){
            BayesKalm.Uneq(
              y.long = y.long, X.long = X.long, id = id, time = time,
              Burn = Burn, Its = Its, 
              Beta.Initial = Beta.Initial, sigma2.beta = sigma2.beta, 
              u0 = u0, P0 = P0, 
              a0 = a0, b0 = b0, 
              c0 = c0, d0 = d0,
              silence = TRUE
            )
          }
        })[1]
      },
      error = function(cond){
        SimTime[i] <<- NA
        print(paste("Error on iteration", i))
      }
      
    )
  }
  
  
  
  
  
  saveRDS(SimTime, paste("Simulations/ComputationSim4/CS", l, n, ".RDS"))
}





