# State-Space-Methods

Functions to create and estimate state space models

## The Model

The purpose of this project is to create a model that is stringent
enough to obtain interpretable linear effects of interest, yet flexible
enough to accurately fit longitudinal neuropsychological data. We
propose the use of a specific State Space Model (SSM), the Local Linear
Trend Model (LLT), to accomplish this task. This model is meant as a
suitable replacement for commonly used mixed effect models. See 'Paper' folder for model description.

This model allows for the ability to capture the linear effects *β*
while also allowing for subject specific random variation in the *α*.
Although the application is toward cognition data, this model can be fit
to any longitudinal process where the trajectory is of importance.

## Fitting the model

To fit such a model we use 1.) a full likelihood approach, 2.) a
partitioned likelihood approach, and 3.) a Bayesian approach. For 1.)
and 2.) we can utilize the Kalman Filter and Kalman Smoother to estimate
both the latent process and linear effects. For approach 3.), we use a
Bayesian Gibb’s sampler, putting a normal prior on the linear effects
and estimating the latent process with a for filter backward sampling
algorithm. We show the Bayesian estimation process to provide the most
accurate effect estimates, maintain the highest numerical stability, and
it has the fastest computation time for samples sizes greater than 100.

## How to use the Bayesian Estimation

The Bayesian LLT can be used with the `BayesKalm.Uneq` in the functions
folder.

``` r
source('functions/BayesKalmUneq2.R')

BayesKalm.Uneq(
  y.long, X.long, id, time, 
  Burn = 500, Its = 1500, 
  Beta.Initial = 0, sigma2.beta = 10, 
  u0 = 0, P0 = 10, 
  a0 = 0.01, b0 = 0.01, 
  c0 = 0.01, d0 = 0.01, 
  silence = FALSE
)
```

The function expects the data to be in a long format. A seperate input
is supplied for the subject `id` of each row along with the `time` of
each measurement. After setting prior distributions the model can be
estimated. The posterior samples will be provided in a list which can be
used for inference.

