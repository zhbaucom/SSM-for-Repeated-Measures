library(tidyverse)
library(lme4)
library(nlme)

#Get back to main directory
### Set working directory to correct location
cdir <- stringr::str_split(getwd(), "/")[[1]]
udir <- cdir[1:which(cdir == "StateSpace")]
setwd(paste(udir, collapse = "/"))

source("functions/SimpleSitNBsim.R")





#for batch job
l <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(l)) l <- 1


SimNum <- 1000

n <- 50 #Number of subjects
p <- 2




u0.true <- 0 #True mean of mu_0
P0.true <- 1 #True variance of mu_0

t <- 10
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


SIG <- c("SSM", "AR(1)")

sigeps <- 5
sigeta <- 3

dataType <- SIG[l]

Itout <- Compiler(reps = SimNum, dataType = dataType)



saveRDS(Itout, paste("Simulations/SimpSimOut/ItoutNACC ", dataType, ".RDS", sep = ""))























