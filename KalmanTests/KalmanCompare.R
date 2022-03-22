
library(tidyverse)
# source("functions/NBsim.R")
# source("functions/BayesKalmUneq2.R") #Holds the Bayesian Kalman Regression
# source("functions/KalmanRegPartition.R")
source("functions/KalmanRegression.R") #Holds the Regular Kalman Regression
source("functions/KalmanRegression2.R") #Holds the Regular Kalman Regression
source("functions/KalmanRecReform.R")
source("functions/KalmanRecReform2.R")
# source("functions/stateSpaceSim.R") #To simulate the data
source("functions/ss.simp.sim.R") #To simulate the data
library(microbenchmark)

ksummary <- function(kout){
  NS <- ncol(kout$alpha.hat)
  n <- nrow(kout$F[[1]])
  Es <- kout$alpha.hat[(n+1):(n+p),NS]
  V <- diag(kout$V[[NS]])[(n+1):(n+p)]
  LCL <- Es + qnorm(0.025) * sqrt(V)
  UCL <- Es + qnorm(0.975) * sqrt(V)
  
  cbind(estimate = Es, LCL = LCL, UCL = UCL)
}



n <- 100 #Number of subjects
B <- c(4, 2, -1) #Beta Vector
p <- length(B) #Number of Betas



u0.true <- 0 #True mean of mu_0
P0.true <- 1 #True variance of mu_0

t <- 6

sigeps = 2
sigeta = 1


X <- matrix(rnorm((p)*n), n, p)
colnames(X) <- paste("X", 1:ncol(X), sep = "")
ytot <- ss.simp.sim(X, B, t = t, sigma2.eps = sigeps, sigma2.eta = sigeta, u0 = u0.true, P0 = P0.true) #Creates a simulated y
y <- ytot
colnames(y) <- paste("y", 1:ncol(y), sep = "")

####Initialize Beta
XY1 <- cbind(id = 1:n, y, X)
XYs <- XY1 %>%
  data.frame() %>%
  gather("time", "y", colnames(.)[grepl("y", colnames(.))]) %>%
  mutate(time = as.numeric(substr(time, 2, nchar(time)))) 
XYs[grepl("X", colnames(XYs))] <- XYs[grepl("X", colnames(XYs))] * XYs$time
XY <- XYs %>%
  filter(!is.na(y)) %>%
  arrange(id, time) 
XY2 <- XY %>%
  group_by(id) %>%
  mutate(Keep = c(TRUE, 2:n() %in% sample(2:n(), sample(1:(length(id)-1), 1)))) %>%
  ungroup(id) %>%
  filter(Keep == TRUE) %>%
  arrange(id, time)
y <- XY2 %>%
  select(id, time, y) %>%
  spread(time, y) %>%
  select(-id) %>%
  as.matrix()


y.long <- XY2$y
X.long <- as.matrix(XY2[,paste("X", 1:p, sep = "")])
id <- XY2$id
time <- XY2$time
Beta.Initial <- coef(lm(y.long ~ X.long - 1))

  

# 
# koutr <- KalmanRegReform(y, X, Beta.Initial = Beta.Initial)
# 
# ksummary(koutr)

# koutr2 <- KalmanRegReform2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial)
# koutr2$sigma2.eps
# koutr2$sigma2.eta
# ksummary(koutr2)


# koutreg2 <- KalmanReg2(y.long, X.long, id = id, time = time,  Beta.Initial = Beta.Initial, P1 = 1e7)
# koutreg2$sigma2.eps
# koutreg2$sigma2.eta
# ksummary(koutreg2)
# koutreg2$Optimization.Time
# koutreg2$Filter.Time
# koutreg2$counts
# 
# 
# 
koutreg3 <- KalmanReg2(y.long, X.long, id = id, time = time,  Beta.Initial = Beta.Initial, P1 = 1e7, maxit = 10)
koutreg3$sigma2.eps
koutreg3$sigma2.eta
ksummary(koutreg3)
koutreg3$Optimization.Time
koutreg3$Filter.Time
koutreg3$counts


koutreg4 <- KalmanReg2(y.long, X.long, id = id, time = time, k = 2, Beta.Initial = Beta.Initial, P1 = 1e7, maxit = 10)
koutreg4$sigma2.eps
koutreg4$sigma2.eta
# ksummary(koutreg4)
cbind(koutreg4$summary,koutreg4$summary[,1] < koutreg4$summary[,3] & koutreg4$summary[,1] > koutreg4$summary[,2] )
koutreg4$Optimization.Time
koutreg4$Filter.Time
koutreg4$counts

# 
# kouto <- KalmanReg(y, X, Beta.Initial = Beta.Initial, P1 = 1e7)
# kouto$sigma2.eps
# kouto$sigma2.eta
# ksummary(kouto)
# kouto$Optimization.Time
# kouto$Filter.Time
# kouto$counts
# 
# 
# 
# 
# 
# # microbenchmark(
# #   KalmanReg(y, X, Beta.Initial = Beta.Initial),
# #   KalmanReg2(y.long, X.long, id = id, time = time, Beta.Initial = Beta.Initial), times = 10L
# # )
# # 
# # 
# 
# 
# varSeq <- seq(0.5, 3, 0.1)
# kmat1 <- kmat2 <- matrix(nrow =length(varSeq), ncol = length(varSeq))
# 
# 
# 
# for(i in seq_along(varSeq)){
#   for(j in seq_along(varSeq)){
#     kout1 <- KalmanRegReform2(y.long, X.long, id = id, time = time,  Beta.Initial = Beta.Initial, P1 = 1e7, sigma2.eps = varSeq[i], sigma2.eta = varSeq[j])
#     kmat1[i,j] <- sum(kout1$SSlogLik[-1])
#     
#     
#     kout2 <- KalmanRegEffec(y, X, Beta.Initial = Beta.Initial, P1 = 1e7, sigma2.eps = varSeq[i], sigma2.eta = varSeq[j])
#     kmat2[i,j] <- sum(kout2$SSlogLik[-1])
#   }
# }
# 
# 
# 
# library(tidyverse)
# colnames(kmat1) <- colnames(kmat2) <- paste("i",varSeq, sep = "")
# # row.names(kmat1) <- row.names(kmat2) <- paste("j",seq_along(varSeq), sep = "")
# 
# kmat1 %>%
#   as_tibble() %>%
#   mutate("j" = varSeq) %>%
#   gather("i", "logLik", colnames(.)[grepl("i", colnames(.))]) %>%
#   mutate(i = as.numeric(substr(i, 2,nchar(i))), j = as.numeric(j)) %>%
#   mutate(rank = rank(logLik)) %>%
#   mutate(i1 = ifelse(rank == 1, i, NA), j1 = ifelse(rank == 1, j, NA)) %>%
#   ggplot(aes(x = i, y = j, fill = rank)) + 
#   geom_tile() +
#   scale_fill_gradient(low = "red", high = "blue") + 
#   geom_point(aes(x = sigeta, y = sigeps), shape = 4, size = 4) +
#   geom_point(aes(x = i1, y = j1), size = 4) +
#   geom_point(aes(x = koutreg2$sigma2.eta, y = koutreg2$sigma2.eps), shape = 24, size = 4) +
#   geom_point(aes(x = koutreg3$sigma2.eta, y = koutreg3$sigma2.eps), shape = 25, size = 4) +
#   ggtitle("Reformation")
# 
# kmat2 %>%
#   as_tibble() %>%
#   mutate("j" = varSeq) %>%
#   gather("i", "logLik", colnames(.)[grepl("i", colnames(.))]) %>%
#   mutate(i = as.numeric(substr(i, 2,nchar(i))), j = as.numeric(j)) %>%
#   mutate(rank = rank(logLik)) %>%
#   mutate(i1 = ifelse(rank == 1, i, NA), j1 = ifelse(rank == 1, j, NA)) %>%
#   ggplot(aes(x = i, y = j, fill = rank)) +
#   geom_tile() +
#   scale_fill_gradient(low = "red", high = "blue") +
#   geom_point(aes(x = sigeta, y = sigeps), shape = 4, size = 4) +
#   geom_point(aes(x = i1, y = j1), shape = 21, size = 4) +
#   geom_point(aes(x = kouto$sigma2.eta, y = kouto$sigma2.eps), shape = 24, size = 4) +
#   ggtitle("Original")
# 
# 

































