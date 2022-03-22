library(tidyverse)
library(lme4)
library(nlme)
source("functions/BayesKalmUneq2.R") 
source("functions/KalmanRegression2.R")
source("functions/KalmanRecReform2.R")


#Get back to main directory
### Set working directory to correct location
cdir <- stringr::str_split(getwd(), "/")[[1]]
udir <- cdir[1:which(cdir == "StateSpace")]
setwd(paste(udir, collapse = "/"))

#for batch job
l <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(l)) l <- 1


#Create New Directory
td <- format(Sys.Date(), "%m%d%Y")
if(!(dir.exists(paste("Simulations/RealDataSimQ/",td, sep = ""))))
  dir.create(paste("Simulations/RealDataSimQ/",td, sep =""))

###################### Read In Data #################

#Functions to remove missing code in tests
m1 <- function(x){
  ifelse(x %in% c(95, 96, 97, 98, -4), NA, x)
}
m2 <- function(x){
  ifelse(x %in% c(995, 996, 997, 998, -4), NA, x)
}



# Read in NACC data
longdata1 <- read.csv("NACC/Data2/investigator_nacc49.csv")

longdata <- longdata1 %>%
  #Select variables needed for this analysis
  select(
    NACCID, NACCUDSD, NACCFDYS, VISITYR, BIRTHYR, RACE, EDUC, SEX, NACCAPOE,
    ANIMALS
  ) %>%
  #Remove missing codes in tests
  mutate_at(vars(ANIMALS), m1) %>%
  #Select only those who are demented for each observation
  group_by(NACCID) %>%
  arrange(NACCID,NACCFDYS) %>%
  filter(!is.na(ANIMALS)) %>%
  #Select those who started as normal and transitioned to impairment
  filter(NACCUDSD[1] == 1 & any(NACCUDSD %in% c(3,4))) %>%
  #Derive additional variables
  mutate(
    time = NACCFDYS/365.25, 
    AGE = VISITYR - BIRTHYR,
    RACEWHITE = ifelse(RACE == 1, 1, 0),
    APOE = ifelse(NACCAPOE %in% c(2, 4), 1, ifelse(NACCAPOE == 9, NA, 0)),
    APOESEX = APOE * (SEX-1),
    Intercept = 1,
    #Do they have an impairment
    MCIoD = NACCUDSD %in% c(3,4),
    FT = min(time[MCIoD]),
    #All time past initial impairment
    DEC = time >= FT,
    #Did the subject ever go to impairment to non-impairment
    reverter = MCIoD != DEC
  ) %>%
  arrange(NACCID, AGE) %>%
  mutate(AgeBase = AGE[1]) %>%
  #Remove all reverters
  filter(!any(reverter)) %>%
  mutate(DEC = DEC * (time > 0) * (time - min(time[DEC > 0]))) %>%
  ungroup(NACCID) %>%
  mutate(A1 = AgeBase, E1 = EDUC, AgeBase = AgeBase - mean(AgeBase), EDUC = EDUC - mean(EDUC))%>%
  mutate(SEX = (2*(SEX - 1)-1) * time, RACEWHITE = (RACEWHITE * 2 -1) * time, APOE = (APOE*2-1) * time, APOESEX = (APOESEX *2 -1))

#Simulation
##Setup Data
ld2 <- longdata[,c( "NACCID", "time","Intercept", "ANIMALS", "SEX", "EDUC", "RACEWHITE", "AgeBase", "DEC","APOE", "APOESEX")] %>%
  na.omit() %>%
  arrange(NACCID, time) %>%
  group_by(NACCID) %>%
  mutate(num = length(NACCID))%>%
  filter(num > 1) %>%
  ungroup()

##Initial Params
Sims <- 100
OutArray <- array(NA, dim = c(4, 3, Sims))
sd1 <- sd(ld2$ANIMALS)/max(ld2$time)
saveRDS(sd1, "NACC/sd1.RDS")
uid <- unique(ld2$NACCID)


VarArray <- array(NA, dim = c(2, 2, Sims))
i <- 1

Its <- 2000
Burn <- floor(Its/2)

BetaMean <- rep(NA, Sims)


for(i in 1:Sims){
  
  #####SIMULATE DATA
  sid <- sample(uid, floor(length(uid)/2))
  
  
  
  # ld2$BinInt <- ld2$NACCID %in% sid
  
  ld2 <- ld2 %>%
    group_by(NACCID) %>%
    mutate(RS = rnorm(1))  %>%
    mutate(BinInt = sample(c(-1, 1), 1) * time) %>%
    ungroup()
  
  
  # ld2 <- ld2 %>%
  #   group_by(NACCID) %>%
  #   mutate(BinInt = runif(1, -1, 1)) %>%
  #   ungroup()
  
  ld2$A2 <- sd1 * ld2$BinInt + ld2$BinInt * ld2$RS+ ld2$ANIMALS
  
  
  
  Beta.Initial <- coef(lm(A2 ~ (time + SEX + EDUC + RACEWHITE + AgeBase + DEC + APOE + APOESEX + BinInt), data = ld2))[-1]
  
  tryCatch(
    error = function(cnd){
      print("Error")
    },
    {  
      ###LME
      lmeMod <- lme(A2 ~ (time +SEX + EDUC + RACEWHITE + AgeBase + DEC + APOE + APOESEX + BinInt), random =  ~1+BinInt|NACCID, data = ld2, method = "ML")
      CIout <- intervals(lmeMod, which = "fixed")$fixed[10,]
      OutArray[1,,i] <- CIout[c(2, 1, 3)]
      VarArray[1,1,i] <-lmeMod$sigma
    }
  )
  
  tryCatch(
    error = function(cnd){
      print("Error")
    },
    { 
      ###LME AR(1)
      lmeMod <- lme(A2 ~ (time +SEX + EDUC + RACEWHITE + AgeBase + DEC + APOE + APOESEX + BinInt), correlation = corAR1(form = ~time), random =  ~1+BinInt|NACCID, data = ld2, method = "ML")
      CIout <- intervals(lmeMod, which = "fixed")$fixed[10,]
      OutArray[2,,i] <- CIout[c(2, 1, 3)]
      VarArray[2,1,i] <-lmeMod$sigma
      
    }
  )
  
  tryCatch(
    error = function(cnd){
      print("Error")
    },
    {
      ###Bayesian
      bkout <- BayesKalm.Uneq(
        y.long = ld2$A2, X.long = as.matrix(ld2[,c("time", "SEX", "EDUC", "RACEWHITE", "AgeBase", "DEC", "APOE", "APOESEX", "BinInt")]), 
        id = ld2$NACCID, time = ld2$time,
        Burn = Burn, Its = Its, 
        Beta.Initial = Beta.Initial, sigma2.beta = 20, 
        u0 = 0, P0 = 100, silence = TRUE
      )
      
      
      OutArray[3,,i] <- c(mean(bkout$Beta[9,Burn:Its], na.rm = TRUE),quantile(bkout$Beta[9,Burn:Its], c(.025, .975), na.rm = TRUE))
      VarArray[1,,i] <-c(mean(bkout$sigma2.eps[Burn:Its], na.rm = T), mean(bkout$sigma2.eta[Burn:Its], na.rm = T)) 
      
      
    }
  )
  
  tryCatch(
    error = function(cnd){
      print("Error")
    },
    {
      ###SSPartition
      
      okout <- KalmanReg2(
        y.long = ld2$A2, X.long = as.matrix(ld2[,c("time","SEX", "EDUC", "RACEWHITE", "AgeBase", "DEC", "APOE", "APOESEX", "BinInt")]), 
        id = ld2$NACCID, time = ld2$time, 
        Beta.Initial = Beta.Initial, a1 = 0, P1 = 1e7, k = 15
      )
      
      
      OutArray[4,,i] <- okout$summary[9,]
      VarArray[2,,i] <-c(okout$sigma2.eps, okout$sigma2.eta)
      
      
    }
  )
  
  saveRDS(OutArray, paste("Simulations/RealDataSimQ/",td,"/NEWoutArray_",l,".RDS", sep = ""))
  saveRDS(VarArray, paste("Simulations/RealDataSimQ/",td,"/NEWVarArray_",l,".RDS", sep = ""))
  print(paste0("Iteration ",i, sep = ""))
  rep <- FALSE
}






# saveRDS(OutArray, "Simulations/RealDataSim/outArray.RDS")



