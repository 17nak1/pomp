rm(list=ls())
mainDir <- "C:/Users/Nazila/Downloads/stochmodel3"
setwd(mainDir)
library(pomp)
source("ModelSnippet_StochModel3.R")

Np <- 1e3
predTime <- c("2020-12-31")
modelname <- "DetModel3"

mle <- c( betaI= 0.375, theta= 0.375, iota= 62.5, beta_sd= 0, dI0= 0.525, dP0= 0.525, dT0= 0.075, dB0= 0, dI1= 0.0625, dP1= 0.4375, dT1= 0.1875, dB1= 0, qP= 0.0625, qH= 0.375, qC= 0.8125, mI= 0.0125, mC= 0.625, mV= 0.6875, sigma= 0.2, kappa= 1, gammaI= 0.2, gammaH= 0.7, gammaC= 0.6625, gammaV= 0.8875, rho= 0.125, TF= 7500, S0= 1, EQ0= 0, PQ0= 0, IQ0= 0, E0= 0, P0= 0, I0= 0, H0= 0, C0= 0, V0= 0, M0= 0 )


source("CreateModel.R")
source("CreateCovars.R")
source("CreateDataset.R")


data <- create_dataset(endTime=endTime)
covars <- create_covars(endTime=endTime)
out <- create_pomp_model(data = data, covars=covars, t0 = 0, dt=0.1)

model <- out$model

  coef(model) <- mle[c(params_mod,params_ic)]
  t <- Sys.time()
  pf <- pfilter(model, filter.mean=T, pred.mean=T,save.states=T,Np=1000);pf@loglik
  Sys.time()-t




