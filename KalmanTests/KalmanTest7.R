#I used the tidyverse packages dplyr, tidyr, and ggplot2 to make the plots.
# install.packages("tidyverse") #Install the package if you haven't already.
library(tidyverse)

# source("functions/KalmanRegression.R") #Holds original Kalman Regression
source("functions/BayesKalmUneq2.R") #Holds the Bayesian Kalman Regression
source("functions/ss.simp.sim.R") #To simulate the data
# set.seed(123)
########## Data Simulation###########
n <- 100 #Number of subjects
B <- c(4, 2, -1) #Beta Vector
dataType <- "AR(1)"
rhoS <- list(order = c(1, 0, 0), ar = 0.9, sd = 1)

p <- length(B) #Number of Betas
X <- matrix(runif(n*p, 0, 20), n, p) #Randomly generate baseline variables
t <- 6

sigeps <- 1 #True sigma2.eps 
sigeta <- 1 #True sigma2.eta

u0.true <- 0 #True mean of mu_0
P0.true <- 1 #True variance of mu_0

colnames(X) <- paste("X", 1:ncol(X), sep = "")


if(dataType == "SSM"){
  ytot <- ss.simp.sim(X, B, t = t, sigma2.eps = sigeps, sigma2.eta = sigeta, u0 = u0.true, P0 = P0.true) #Creates a simulated y
  
  y <- ytot
  
  colnames(y) <- paste("y", 1:ncol(y), sep = "")
  
  XY1 <- cbind(id = 1:n, y, X)
  
  XYs <- XY1 %>%
    data.frame() %>%
    gather("time", "y", colnames(.)[grepl("y", colnames(.))]) %>%
    mutate(time = as.numeric(substr(time, 2, nchar(time)))) 
  
  XYs[grepl("X", colnames(XYs))] <- XYs[grepl("X", colnames(XYs))] * XYs$time
  XY <- XYs %>%
    filter(!is.na(y)) %>%
    arrange(id, time) 
}else if(dataType == "AR(1)"){
  XY <- lapply(1:t, function(x){
    xout <- cbind(1:n, x, (x*X) %*% B, x*X)
    colnames(xout) <- c("id","time", "XB", paste("X", 1:p, sep = ""))
    xout
  } ) %>%
    do.call("rbind", .) %>%
    as_tibble() %>%
    group_by(id) %>%
    mutate(error = arima.sim(n = n(), model = rhoS)) %>%
    mutate(y = XB + error) %>%
    select(-XB, -error) %>%
    ungroup(id) %>%
    arrange(id, time)
  y <- XY %>%
    select(id, time, y) %>%
    spread(time, y) %>%
    select(-id) %>%
    as.matrix()
  
}


XY <- XY %>%
  # group_by(id) %>%
  # mutate(Keep = c(TRUE, 2:n() %in% sample(2:n(), sample(1:(length(id)-1), 1)))) %>%
  # ungroup(id) %>%
  # filter(Keep == TRUE) %>%
  # group_by(id) %>%
  # mutate(Keep = c(TRUE, 2:n() %in% (1:sample(2:n(), 1)))) %>%
  # ungroup(id) %>%
  # filter(Keep == TRUE) %>%
  arrange(id, time)


########### PRIORS for model ##############
#Beta ~ N(Beta.Initial, sigma2.beta)
Beta.Initial <- 0
Beta.Initial <- B
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
Its <- 10000
Burn <- ceiling(Its/2)

####For long format

y.long <- XY$y
X.long <- as.matrix(XY[,paste("X", 1:p, sep = "")])
id <- XY$id
time <- XY$time


######Initialize
Beta.Initial <- coef(lm(y.long ~ X.long - 1))


###### Run Bayesian Kalman Filter ######
Bayes.Time <- system.time({
  bkout <- BayesKalm.Uneq(
    y.long = y.long, X.long = X.long, id = id, time = time,
    Burn = Burn, Its = Its, 
    Beta.Initial = Beta.Initial, sigma2.beta = sigma2.beta, 
    u0 = u0, P0 = P0, 
    a0 = a0, b0 = b0, 
    c0 = c0, d0 = d0
  )
})


# Bayesian Method Betas
cbind(
  Beta.True = B,
  Estimate = rowMeans(bkout$Beta[,Burn:Its]),
  apply(bkout$Beta[,Burn:Its], 1, quantile, c(.025, 0.975)) %>%
    t()
) 

# Bayesian Convergence Plot

bkout$Beta %>%
  t() %>%
  data.frame() %>%
  set_names(paste("B", 1:length(B), sep = "")) %>%
  mutate(ind = 1:nrow(.)) %>%
  gather("Beta", "Value", -ind) %>%
  mutate(Z = B[as.numeric(substr(Beta, 2, nchar(Beta)))]) %>%
  ggplot(aes(x = ind, y = Value, color = Beta)) + 
  geom_point(shape = 21) + 
  facet_wrap(. ~ Beta, scales = "free_y") + 
  geom_hline(aes(yintercept = Z), col = "red") +
  geom_vline(aes(xintercept = Burn)) +
  theme(legend.position = "none")




# Bayesian Method Sigmas

cbind(
  Sigma.True = c(sigeps, sigeta),
  Estimate = c(mean(bkout$sigma2.eps[Burn:Its]), mean(bkout$sigma2.eta[Burn:Its])),
  apply(data.frame(bkout$sigma2.eps, bkout$sigma2.eta)[Burn:Its,], 2, quantile, c(.025, 0.975)) %>%
    t()
) 

# Bayesian Sigma Convergence

data.frame(bkout$sigma2.eps, bkout$sigma2.eta) %>%
  set_names(c("var.eps", "var.eta")) %>%
  mutate(ind = 1:nrow(.)) %>%
  gather("Variance", "Value", var.eps, var.eta) %>%
  mutate(Z = case_when(
    Variance == "var.eps" ~ sigeps,
    Variance == "var.eta" ~ sigeta
  )) %>%
  ggplot(aes(x = ind, y = Value, color = Variance)) + 
  geom_point(shape = 21) + 
  facet_wrap(. ~ Variance, scales = "free_y") + 
  geom_hline(aes(yintercept = Z), col = "red") +
  geom_vline(aes(xintercept = Burn)) +
  theme(legend.position = "none")


# id.tag <- 0
# 
# id.tag <- id.tag + 1
# plot(time[id == id.tag], XY$mu[id == id.tag], ylim = quantile(c(bkout$mu[id.tag,,Burn:Its], XY$mu[id == id.tag]), c(0, 1)), xlim = c(1, 5), main = id.tag)
# 
# 
# 
# for(i in 1:sum(id == id.tag)){
#   points(x = rep(time[id == id.tag][i], length(bkout$mu[id.tag,i,Burn:Its])), bkout$mu[id.tag,i,Burn:Its], col = "red")
# }
# points(time[id == id.tag], XY$mu[id == id.tag])






