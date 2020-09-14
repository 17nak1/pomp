rm(list=ls())
mainDir <- "~/Git/pomp/R"
setwd(mainDir)
source("ModelSnippet_StochModel3.R")

predTime <- c("2020-12-31")
modelname <- "DetModel3"
rw_size <- 0.05
current_params <- c( betaI= 1.7088129733919104,      theta= 0.36247417077568833,      iota= 0.009072753906250272,      beta_sd= 0.722146657339815,      dI0= 0.3067998824837139,      dP0= 0.701577932458249,      dT0= 5.9680438369014155e-9,      dB0= 0.9629263473105145,      dI1= 0.1310459432396453,      dP1= 0.8823248285729485,      dT1= 1,      dB1= 0.8326569628565079,      qP= 0.05200304375822984,      qH= 0.5437813633584803,      qC= 1,      mI= 0,      mC= 1,      mV= 1,      sigma= 0.20000000000000018,      kappa= 1,      gammaI= 0.7105465943799485,      gammaH= 0.10000728795427358,      gammaC= 0.8732959116681945,      gammaV= 0.2564535075184308,      rho= 0.21797016009728593,      TF= 7.076093553614824,      S0= 1,      EQ0= 0,      PQ0= 0,      IQ0= 0,      E0= 0,      P0= 0,      I0= 0,      H0= 0,      C0= 0,      V0= 0,      M0= 0)
# current_params <-c( betaI= 0.375, theta= 0.375, iota= 62.5, beta_sd= 0, dI0= 0.525, dP0= 0.525, dT0= 0.075, dB0= 0, dI1= 0.0625, dP1= 0.4375, dT1= 0.1875, dB1= 0, qP= 0.0625, qH= 0.375, qC= 0.8125, mI= 0.0125, mC= 0.625, mV= 0.6875, sigma= 0.2, kappa= 1, gammaI= 0.2, gammaH= 0.7, gammaC= 0.6625, gammaV= 0.8875, rho= 0.125, TF= 7500, S0= 1, EQ0= 0, PQ0= 0, IQ0= 0, E0= 0, P0= 0, I0= 0, H0= 0, C0= 0, V0= 0, M0= 0 )


source("CreateModel.R")
source("CreateCovars.R")
source("CreateDataset.R")
source("DetermineRandomWalks.R")


data <- create_dataset(endTime=endTime)
covars <- create_covars(endTime=endTime)
out <- create_pomp_model(data = data, covars=covars, t0 = 0, dt=0.1)

model <- out$model
params_fixed <- c("beta_sd","dB0", "dB1","sigma","kappa")
params_fit <-params_mod[-which(params_mod %in% params_fixed)]
# coef(model) <- current_param[c(params_mod,params_ic)]
coef(model) <- c(current_params)
tt <- Sys.time()
traj.match(model,
           start=current_params,
           transform=TRUE,
           est=c(), method = 'subplex') -> sets_traj;logLik(sets_traj)
Sys.time() -tt 

# coef(model) <- c(coef(sets_traj))
# 
# mif2(model,
#      Nmif=1,
#      start=current_params,
#      transform=TRUE,
#      pars=params_fit,
#      rw.sd=rw,
#      Np=100,
#      var.factor=2,
#      cooling.type="hyperbolic",
#      cooling.fraction=0.05) -> sets_mf;sets_mf@loglik
# 
# coef(model) <- coef(sets_mf)
# tt <-Sys.time()
# pf <- pfilter(model, filter.mean=T, pred.mean=T,save.states=T,Np=1000);pf@loglik
# Sys.time() -tt

# js=read.csv("../../../Gitlab/dcp-r/pomp/samples/filterMean.csv")
# p = as.data.frame(pf@filter.mean)
# plot(pf@filter.mean[20,])
# # points(pf@filter.mean[1,], col="blue")
# points(js$S, col="red")
# 
# 
# 
