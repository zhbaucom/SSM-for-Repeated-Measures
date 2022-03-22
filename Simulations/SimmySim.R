library(tidyverse)
library(lme4)
library(nlme)

#Get back to main directory
### Set working directory to correct location
cdir <- stringr::str_split(getwd(), "/")[[1]]
udir <- cdir[1:which(cdir == "State-Space-Methods")]
setwd(paste(udir, collapse = "/"))







source("functions/NBsim.R")
source("functions/BayesKalmUneq2.R") #Holds the Bayesian Kalman Regression
source("functions/KalmanRegPartition.R")
source("functions/KalmanRegression.R") #Holds the Regular Kalman Regression
# source("functions/stateSpaceSim.R") #To simulate the data
source("functions/ss.simp.sim.R") #To simulate the data



#for batch job
l <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(l)) l <- 1


#Fixed Params
SimNum <- 1000

n <- 100 #Number of subjects
B <- c(4, 2, -1) #Beta Vector
p <- length(B) #Number of Betas



u0.true <- 0 #True mean of mu_0
P0.true <- 1 #True variance of mu_0

t <- 6
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







# SIG <- data.frame(
#   sigeps = c(1, 1 , 10, 1, 3, 1),
#   sigeta = c(1, 0.1, 1, 0.5, 1, 0)
# )

SIG <- list(
  list(
    sigeps = 1,
    sigeta = 0,
    dataType = "SSM"
  ),
  list(
    sigeps = 1,
    sigeta = .1,
    dataType = "SSM"
  ),
  list(
    sigeps = 10,
    sigeta = 1,
    dataType = "SSM"
  ),
  list(
    sigeps = 1,
    sigeta = 0.5,
    dataType = "SSM"
  ),
  list(
    sigeps = 1,
    sigeta = 1,
    dataType = "SSM"
  ),
  list(
    sigeps = 3,
    sigeta = 1,
    dataType = "SSM"
  ),
  list(
    sigeps = "AR1",
    sigeta = "Small",
    dataType = "AR(1)",
    rhoS = list(order = c(1, 0, 0), ar = 0.1, sd = 1)
  ),
  list(
    sigeps = "AR1",
    sigeta = "Large",
    dataType = "AR(1)",
    rhoS = list(order = c(1, 0, 0), ar = 0.9, sd = 1)
  )
)



list2env(SIG[[l]], envir = globalenv())



Itout <- Compiler(reps = SimNum, dataType = dataType, simXfun = runif, 0, 20)



saveRDS(Itout, paste("Simulations/SimOut/ItoutPrezResults ", paste(gsub("[.]", "_",as.character(c(sigeps,sigeta))), collapse = " "), ".RDS", sep = ""))








