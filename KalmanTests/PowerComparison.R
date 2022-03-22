library(tidyverse)
library(lme4)
library(nlme)
source("functions/BayesKalmUneq2.R")
#Functions to remove missing code in tests
m1 <- function(x){
  ifelse(x %in% c(95, 96, 97, 98, -4), NA, x)
}
m2 <- function(x){
  ifelse(x %in% c(995, 996, 997, 998, -4), NA, x)
}


# Read in NACC data
longdata <- read.csv("NACC/Data2/investigator_nacc49.csv")

longdata <- longdata %>%
  #Select variables needed for this analysis
  select(
    NACCID, NACCUDSD, NACCFDYS, VISITYR, BIRTHYR, RACE, EDUC, SEX, NACCAPOE,
    LOGIMEM, MEMUNITS, ANIMALS, VEG, BOSTON, TRAILA, TRAILB, DIGIF, DIGIB, WAIS
  ) %>%
  #Remove missing codes in tests
  mutate_at(vars(LOGIMEM, MEMUNITS, DIGIF, DIGIB, ANIMALS, VEG, WAIS, BOSTON), m1) %>%
  mutate_at(vars(TRAILA, TRAILB), m2) %>%
  #Select only those who are demented for each observation
  group_by(NACCID) %>%
  filter(rep(all(NACCUDSD == 4), length(NACCID))) %>%
  #Derive additional variables
  mutate(
    time = NACCFDYS/365.25, 
    AGE = VISITYR - BIRTHYR,
    RACEWHITE = ifelse(RACE == 1, 1, 0),
    APOE = ifelse(NACCAPOE %in% c(2, 4), 1, ifelse(NACCAPOE == 9, NA, 0)),
    Intercept = 1
  ) %>%
  arrange(NACCID, AGE) %>%
  mutate(AgeBase = AGE[1]) %>%
  ungroup(NACCID) %>%
  mutate(AgeBase = AgeBase - mean(AgeBase), EDUC = EDUC - mean(EDUC))



Bind <- 1

ld2 <- longdata %>%
  select(ANIMALS, NACCID, time) %>%
  na.omit() %>%
  group_by(NACCID) %>%
  mutate(Keep = length(NACCID) > 1) %>%
  filter(Keep) %>%
  ungroup(NACCID)

uid <- unique(ld2$NACCID)
sid <- sample(uid, floor(length(uid)/2))



ld2$BinInt <- ld2$NACCID %in% sid * ld2$time 
ld2$A2 <- Bind * ld2$BinInt + ld2$ANIMALS  

lmemod <- lme(A2 ~ BinInt, random = ~1|NACCID, method = "ML", data = ld2)
CIout <- intervals(lmemod, which = "fixed")$fixed[2,]
Sig <- CIout[1] > 0

ld2 <- longdata %>%
  select(ANIMALS, NACCID, time) %>%
  na.omit() %>%
  group_by(NACCID) %>%
  mutate(Keep = length(NACCID) > 1) %>%
  filter(Keep) %>%
  ungroup(NACCID)



tryCatch(
  error = function(cnd){
    CIout <- NA
    CIout
  },
  {
    stop("ERROR")
    
  }
)

power <- map(seq(1, 1.05, by = .01), function(x){
  sigV <- map_lgl(1:100, function(y){


    uid <- unique(ld2$NACCID)
    sid <- sample(uid, floor(length(uid)/2))



    ld2$BinInt <- ld2$NACCID %in% sid * ld2$time
    ld2$A2 <- x * ld2$BinInt + ld2$ANIMALS


    tryCatch(
      error = function(cnd){NA},
      {
        lmemod <- lme(A2 ~ BinInt, random = ~1|NACCID, method = "ML", correlation = corAR1(form = ~time), data = ld2)
        CIout <- suppressWarnings(intervals(lmemod, which = "fixed")$fixed[2,])
        CIout[1] > 0 
      }
    )

  })
  sigV
})

map(power, mean, na.rm = TRUE)



AllB <- seq(1, 1.05, by = .01)

B <- AllB[1]
Sims <- 2
B <- 10

SimArray<- array(dim = c(Sims, 3, 3))

for(i in 1:Sims){
  
  # uid <- unique(ld2$NACCID)
  # sid <- sample(uid, floor(length(uid)/2))
  # rid <- sample(uid, floor(length(uid)/2))
  # ld2$BinInt <- (ld2$NACCID %in% sid) * ld2$time 
  # ld2$Red <- (ld2$NACCID %in% rid) * ld2$time 
  # ld2$A2 <- B * ld2$BinInt + ld2$ANIMALS  
  # 
  ld2$BinInt <- runif(nrow(ld2))
  ld2$A2 <- B * ld2$BinInt + ld2$ANIMALS  
  
  Beta.Initial <- coef(lm(A2 ~ BinInt, data = ld2))[-1]
  
  
  ###LME
  
  SimArray[i,,1] <- tryCatch(
    error = function(cnd){
      NA
    },
    {
      lmemod <- lme(A2 ~ BinInt, random = ~1|NACCID, method = "ML", data = ld2)
      suppressWarnings(intervals(lmemod, which = "fixed")$fixed[2,])
      
    }
  )
  
  
  ###LME AR(1)
  
  SimArray[i,,2] <- tryCatch(
    error = function(cnd){
      NA
    },
    {
      lmemod <- lme(A2 ~ BinInt, random = ~1|NACCID, method = "ML", correlation = corAR1(form = ~time), data = ld2)
      suppressWarnings(intervals(lmemod, which = "fixed")$fixed[2,])
      
    }
  )
  
  ###Bayesian
  
  
  bkout <- BayesKalm.Uneq(
    y.long = ld2$A2, X.long = as.matrix(ld2[,c("BinInt", "Red")]), id = ld2$NACCID, time = ld2$time,
    Its = 3000, 
    Beta.Initial = Beta.Initial, sigma2.beta = 20, 
    u0 = 0, P0 = 100,
    silence = FALSE
  )
  
  
}



plot(bkout$Beta[2,])
plot(bkout$sigma2.eps)
plot(bkout$sigma2.eta)







##### LME COVERAGE

sigV <- map(1:100, function(y){
  
  
  uid <- unique(ld2$NACCID)
  sid <- sample(uid, floor(length(uid)/2))
  
  
  
  ld2$BinInt <- ld2$NACCID %in% sid * ld2$time
  ld2$A2 <- 10 * ld2$BinInt + ld2$ANIMALS
  
  
  tryCatch(
    error = function(cnd){NA},
    {
      lmemod <- lme(A2 ~ BinInt, random = ~1|NACCID, method = "ML", data = ld2)
      CIout <- suppressWarnings(intervals(lmemod, which = "fixed")$fixed[2,])
      CIout
    }
  )
  
})


map(sigV, 2) %>%
  unlist() %>%
  mean()





###### Bayesian Coverage

sigV2 <- map(1:100, function(y){
  
  
  uid <- unique(ld2$NACCID)
  sid <- sample(uid, floor(length(uid)/2))
  
  
  
  ld2$BinInt <- ld2$NACCID %in% sid * ld2$time
  ld2$Red <- runif(nrow(ld2))
  ld2$A2 <- 10 * ld2$BinInt + ld2$ANIMALS
  
  
  tryCatch(
    error = function(cnd){NA},
    {
      bkout <- BayesKalm.Uneq(
        y.long = ld2$A2, X.long = as.matrix(ld2[,c("BinInt", "Red")]), id = ld2$NACCID, time = ld2$time,
        Its = 2000, 
        Beta.Initial = c(10, 0), sigma2.beta = 20, 
        u0 = 0, P0 = 100,
        silence = TRUE
      )
      bkout$Beta[1,]
      }
  )
  
})


sigV2 %>%
  map(~quantile(.x, c(.025, .975), na.rm = TRUE)) %>%
  map(function(x){
    x[1] < 10 & x[2] > 10
  }) %>%
  unlist() %>%
  mean()

sigV2 %>%
  map(mean, na.rm = TRUE) %>%
  unlist() %>%
  mean()









