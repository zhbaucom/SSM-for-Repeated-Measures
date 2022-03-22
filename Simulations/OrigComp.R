library(knitr)
library(tidyselect)
library(tidyverse)
library(microbenchmark)
library(lme4)
library(nlme)
source("functions/ss.simp.sim.R")
source("functions/KalmanRegression.R")
source("functions/KalmanRecEffic.R")

Bex <- function(kout, n, p){
  Tot <- ncol(kout$alpha.hat) 
  kout$alpha.hat[(n+1):(n+p),Tot]
}
VBex <- function(kout, n, p){
  Tot <- ncol(kout$alpha.hat) 
  diag(kout$V[[Tot]][(n+1):(n+p),(n+1):(n+p)])
}

#Back to correct Directory
Continue <- FALSE
while(!Continue){
  cdir <- str_split(getwd(), "/")[[1]]
  if(length(cdir) == 0)errorCondition("Wrong Directory")
  if(tail(cdir, 1) != "State-Space-Methods") setwd(paste(cdir[-length(cdir)], collapse = "/")) else Continue <- TRUE
}


#for batch job
l <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(l)) l <- 1



n <- 50

t <- 4
P1 <- 1
a1 <- 0
Sims <- 1000


sigma2.eps.test <- 4
sigma2.eta.test <- 1


ssorigSim <- function(mCh = 0){
  # misVec <- c(rep(1,n), sample(c(NA, 1), n*(t-1), replace = TRUE, prob = c(mr, 1-mr)))
  # misMat <- matrix(misVec, n, t)
  # X <- cbind(1,matrix(rnorm((p-1)*n), n, p-1))
  X <- matrix(rnorm((p)*n), n, p) + mCh
  # yorig <- ss.simp.sim(X, B, t, sigma2.eta = sigma2.eta.test, sigma2.eps = sigma2.eps.test, P1 = P1, a1 = a1)
  yorig <- ss.simp.sim(X, B, t, sigma2.eta = sigma2.eta.test, sigma2.eps = sigma2.eps.test, P0 = 1, u0 = 0)
  # y <- yorig * misMat
  y <- yorig
  yn <- y
  colnames(yn) <- paste("y",1:ncol(y), sep = "")
  Xn <- X 
  colnames(Xn) <- paste("X", 1:ncol(X), sep = "")
  
  XYs <- cbind(yn, Xn) %>%
    as_tibble() %>%
    gather("time", "y", colnames(.)[grepl("y", colnames(.))]) %>%
    mutate(time = substr(time, 2, nchar(time)) %>% as.numeric()) 
  
  XYs[grepl("X", colnames(XYs))] <- XYs[grepl("X", colnames(XYs))] * XYs$time
  y.long <- XYs$y
  X.long <- as.matrix(XYs[,paste("X", 1:p, sep = "")])
  
  Beta.Initial <- coef(lm(y.long ~ X.long - 1))
  
  
  kout1 <- KalmanReg(y = y, X = X, P1 = 100, sigma2.eta = sigma2.eta.test, sigma2.eps = sigma2.eps.test)
  
  kout2 <- KalmanReg(y = y, X = X, Beta.Initial = Beta.Initial, P1 = 100, sigma2.eta = sigma2.eta.test, sigma2.eps = sigma2.eps.test)
  
  kout3 <- KalmanReg(y = y, X = X, Beta.Initial = NULL, P1 = 1e7, sigma2.eta = sigma2.eta.test, sigma2.eps = sigma2.eps.test)
  
  kout4 <- KalmanReg(y = y, X = X, Beta.Initial = Beta.Initial, P1 = 1e7, sigma2.eta = sigma2.eta.test, sigma2.eps = sigma2.eps.test)
  
  
  
  EsMat <- cbind(
    Bex(kout1, n = n, p = p),
    Bex(kout2, n = n, p = p),
    Bex(kout3, n = n, p = p),
    Bex(kout4, n = n, p = p)
  )
  
  Vmat <- cbind(
    VBex(kout1, n = n, p = p),
    VBex(kout2, n = n, p = p),
    VBex(kout3, n = n, p = p),
    VBex(kout4, n = n, p = p)
  )
  
  UCL <- EsMat + qnorm(0.975) * sqrt(Vmat)
  LCL <- EsMat + qnorm(0.025) * sqrt(Vmat)
  
  CovMat <- (B < UCL) & (B > LCL)
  list(EsMat, CovMat)

}


OrigComp <- function(Sims = 100, mCh = 0){
  Fails <- 0
  for(i in 1:Sims){
    Again <- TRUE
    print(i)
    while(Again){
      # DSout <- NBdataSim()
      # list2env(DSout, envir = globalenv())
      sos <- try(ssorigSim(mCh = mCh), silent = TRUE)
      if(class(sos) == "list"){ Again <- FALSE }else{ Fails <- Fails + 1}
      if(Fails > 100) break
    }
    if(Fails > 100) break
    if(i == 1){
      EsArray <- CovArray <- array(dim = c(dim(sos[[1]]), Sims))
    }
    EsArray[,,i] <- sos[[1]]
    CovArray[,,i] <- sos[[2]]
   }
  list(Estimate = EsArray, Coverage = CovArray, B = B, Fails = Fails)
  }




  





paramList <- list(
  list(
    B = c(4,2, 10, -1),
    mCh = 0
  ),
  list(
    B = c(4,2, 10, -1),
    mCh = 20
  ),
  list(
    B = c(4,2, 10, -1),
    mCh = -20
  ),
  list(
    B = rev(c(4,2, 10, -1)),
    mCh = 0
  )
)


B <- paramList[[l]]$B
mCh <- paramList[[l]]$mCh

p <- length(B)


origs <- OrigComp(Sims = Sims, mCh = 0)

saveRDS(origs, paste("Simulations/OrigOut/Scenerio_", l, ".RDS",sep = ""))




