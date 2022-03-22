library(tidyverse)
library(lme4)
library(nlme)

#Get back to main directory
Continue <- FALSE
while(!Continue){
  cdir <- str_split(getwd(), "/")[[1]]
  if(length(cdir) == 0)errorCondition("Wrong Directory")
  if(tail(cdir, 1) != "State-Space-Methods") setwd(paste(cdir[-length(cdir)], collapse = "/")) else Continue <- TRUE
}


source("functions/NBsim.R")
source("functions/BayesKalmUneq.R") #Holds the Bayesian Kalman Regression
source("functions/KalmanRegPartition.R")
source("functions/KalmanRegression.R") #Holds the Regular Kalman Regression
# source("functions/stateSpaceSim.R") #To simulate the data
source("functions/ss.simp.sim.R") #To simulate the data



#for batch job
l <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(l)) l <- 1


#Fixed Params

SimNum <- 100


n <- 100 #Number of subjects
B <- c(4, 2, -1) #Beta Vector
p <- length(B) #Number of Betas



u0.true <- 0 #True mean of mu_0
P0.true <- 1 #True variance of mu_0

t <- 4
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







SIG <- data.frame(
  sigeps = c(1, 1 , 10, 1, 3, 1),
  sigeta = c(1, 0.1, 1, 0.5, 1, 0)
)


sigeps <- SIG$sigeps[l]
sigeta <- SIG$sigeta[l]


Itout <- Compiler(reps = SimNum, simXfun = runif, 0, 20)

saveRDS(Itout, paste("Simulations/SimOut/ItoutUT4 ", paste(gsub("[.]", "_",as.character(c(sigeps,sigeta))), collapse = " "), ".RDS", sep = ""))








