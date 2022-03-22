library(tidyverse)
library(lme4)
library(nlme)

#Get back to main directory
### Set working directory to correct location
cdir <- stringr::str_split(getwd(), "/")[[1]]
udir <- cdir[1:which(cdir == "StateSpace")]
setwd(paste(udir, collapse = "/"))

source("functions/NBsim3.R")





#for batch job
l <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(l)) l <- 1


#Fixed Params
SimNum <- 1000

n <- 100 #Number of subjects
B <- readRDS("NACC/B.RDS") #Beta Vector
B <- c(B[-1], B[1])
p <- length(B) #Number of Betas




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




SIG <- list(
  list(
    sigeps = 1,
    sigeta = 0,
    dataType = "SSM"
  ),
  list(
    sigeps = 3,
    sigeta = 0,
    dataType = "SSM"
  ),
  list(
    sigeps = 3,
    sigeta = 1,
    dataType = "SSM"
  ),
  list(
    sigeps = 3,
    sigeta = 2,
    dataType = "SSM"
  ),
  list(
    sigeps = 3,
    sigeta = 3,
    dataType = "SSM"
  ),
  list(
    sigeps = 30,
    sigeta = 10,
    dataType = "SSM"
  ),
  list(
    sigeps = 60,
    sigeta = 20,
    dataType = "SSM"
  ),
  list(
    sigeps = "AR1",
    sigeta = "None",
    dataType = "AR(1)",
    rhoS = list(order = c(1, 0, 0), ar = 0, sd = 1)
  ),
  list(
    sigeps = "AR1",
    sigeta = "Small",
    dataType = "AR(1)",
    rhoS = list(order = c(1, 0, 0), ar = 0.1, sd = 1)
  ),
  list(
    sigeps = "AR1",
    sigeta = "Medium",
    dataType = "AR(1)",
    rhoS = list(order = c(1, 0, 0), ar = 0.5, sd = 1)
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



saveRDS(Itout, paste("Simulations/SimOutReform/ItoutNACC4 ", paste(gsub("[.]", "_",as.character(c(sigeps,sigeta))), collapse = " "), ".RDS", sep = ""))







