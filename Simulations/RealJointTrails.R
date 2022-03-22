library(tidyverse)
library(cIRT)
library(mvtnorm)
library(rlang)
source("functions/BayesKalmJoint.R")



#for batch job
l <- as.numeric(Sys.getenv("SGE_TASK_ID")) 
if (is.na(l)) l <- 1


############GENERAL FUNCTIONS###################
# load("NACC/Data/RData")
# rm(list = ls()[ls() != "longdata"])
fct_case_when <- function(...) {
  args <- as.list(match.call())
  levels <- sapply(args[-1], function(f) f[[3]])  # extract RHS of formula
  levels <- levels[!is.na(levels)]
  factor(dplyr::case_when(...), levels=levels)
}
#Bins for histogram
base_bins <- function(type) {
  fun <- switch(type,
                Sturges = nclass.Sturges,
                scott = nclass.scott,
                FD = nclass.FD,
                stop("Unknown type", call. = FALSE)
  )
  
  function(x) {
    (max(x) - min(x)) / fun(x)
  }
}

#Functions to remove missing code in tests
m1 <- function(x){
  ifelse(x %in% c(95, 96, 97, 98, -4), NA, x)
}
m2 <- function(x){
  ifelse(x %in% c(995, 996, 997, 998, -4), NA, x)
}


mmap <- function(data, label, column,...){
  label <- rlang::ensym(label)
  data %>%
    mutate(!!label := map(!!!rlang::enquos(column), ...))
}

mmap2 <- function(data, label, column1, column2,...){
  label <- rlang::ensym(label)
  data %>%
    mutate(!!label := map2(!!!rlang::enquos(column1), !!!rlang::enquos(column2), ...))
}



#################DATA FORMAT######################
# Read in NACC data
longdata1 <- read.csv("NACC/Data2/investigator_nacc49.csv")


longdata <- longdata1 %>%
  #Select variables needed for this analysis
  select(
    NACCID, NACCUDSD, NACCFDYS, VISITYR, BIRTHYR, RACE, EDUC, SEX, NACCAPOE, NACCALZD,
    LOGIMEM, MEMUNITS, ANIMALS, VEG, BOSTON, TRAILA, TRAILB, DIGIF, DIGIB, WAIS
  ) %>%
  #Remove missing codes in tests
  mutate_at(vars(LOGIMEM, MEMUNITS, DIGIF, DIGIB, ANIMALS, VEG, WAIS, BOSTON), m1) %>%
  mutate_at(vars(TRAILA, TRAILB), m2) %>%
  #Select only those who are demented for each observation
  group_by(NACCID) %>%
  arrange(NACCID,NACCFDYS) %>%
  filter(
    !is.na(TRAILA), !is.na(TRAILB), !is.na(VISITYR), !is.na(BIRTHYR), !is.na(RACE), !is.na(EDUC), !is.na(SEX), !is.na(NACCAPOE)
  ) %>%
  arrange(NACCID, NACCFDYS) %>%
  #Select those who started as normal and transitioned to impairment
  filter(NACCUDSD[1] == 1 & any(NACCUDSD %in% c(3,4)), any(NACCALZD == 1)) %>%
  #Derive additional variables
  mutate(
    time = NACCFDYS/365.25, 
    AGE = VISITYR - BIRTHYR,
    AgeBase = AGE[1],
    RACEWHITE = ifelse(RACE == 1, 1, 0),
    APOE = ifelse(NACCAPOE %in% c(2, 4), 1, ifelse(NACCAPOE == 9, NA, 0)),
    SEX = SEX - 1,
    APOESEX = APOE * (SEX),
    Intercept = 1,
    #Do they have an impairment
    MCIoD = NACCUDSD %in% c(3,4),
    FT = min(time[MCIoD]),
    #All time past initial impairment
    DEC = time >= FT,
    #Did the subject ever go to impairment to non-impairment
    reverter = MCIoD != DEC,
    #Binary effect which will be added to outcomes
    Group = rbinom(1,1 ,.5)
  ) %>%
  arrange(NACCID, AGE)%>%
  #Remove all reverters
  filter(!any(reverter)) %>%
  #Select all variables that will be used for modeling
  select(NACCID, time,Intercept, SEX, EDUC, RACEWHITE, AgeBase, DEC, APOE, APOESEX, NACCUDSD, Group,
         TRAILA, TRAILB) %>%
  na.omit() %>%
  group_by(NACCID) %>%
  filter(any(DEC >0) & any(DEC == 0)) %>%
  #Create knot for those when they transitioned
  mutate(DEC = DEC * (time > 0) * (time - min(time[DEC > 0]))) %>%
  ungroup(NACCID) %>%
  mutate(A1 = AgeBase, E1 = EDUC, AgeBase = AgeBase - mean(AgeBase), EDUC = EDUC - mean(EDUC)) %>%
  mutate(
    #Make all variables interacted with time
    SEX = SEX * time, RACEWHITE = RACEWHITE * time, APOE = APOE * time, 
    APOESEX = APOESEX *time, EDUC = EDUC * time, AgeBase = AgeBase * time, 
    Group = Group * time,
    #Adjust TRAILA and TRAILB to containt "group" effect
    TRAILA = TRAILA + Group,
    TRAILB = TRAILA + Group
  ) 


bkjout <- BayesKalmJoint(
  data = longdata, outcomes = c("TRAILA", "TRAILB"), predictors = c("time", "RACEWHITE","SEX", "APOE", "APOESEX", "EDUC", "DEC", "AgeBase", "Group"),
  timevar = "time", id = "NACCID", numits = 5000, seed = 24, numitInit = 1500, burnInit = 500
)




saveRDS(bkjout, paste("Simulations/RealJointTrails/Sim_",l, ".RDS", sep = ""))
