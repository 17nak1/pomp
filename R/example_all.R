rm(list=ls())
mainDir <- "~/Git/pomp/R"
setwd(mainDir)
source("ModelSnippet_StochModel3.R")

predTime <- c("2020-12-31")
modelname <- "DetModel3"
rw_size <- 0.05
current_params <- c(betaI= 0.5, theta= 0.5, iota= 49.99999999999999, beta_sd= 0, dI0= 0.30000000000000004, dP0= 0.30000000000000004, dT0= 0.30000000000000004, dB0= 0, dI1= 0.24999999999999994, dP1= 0.24999999999999994, dT1= 0.24999999999999994, dB1= 0, qP= 0.24999999999999994, qH= 0.5, qC= 0.75, mI= 0.04999999999999999, mC= 0.5, mV= 0.75, sigma= 0.2, kappa= 1, gammaI= 0.2, gammaH= 0.6000000000000001, gammaC= 0.55, gammaV= 0.55, rho= 0.5, TF= 5999.999999999995, S0= 1, EQ0= 0, PQ0= 0, IQ0= 0, E0= 0, P0= 0, I0= 0, H0= 0, C0= 0, V0= 0, M0= 0)
  # c( betaI= 0.375, theta= 0.375, iota= 62.5, beta_sd= 0, dI0= 0.525, dP0= 0.525, dT0= 0.075, dB0= 0, dI1= 0.0625, dP1= 0.4375, dT1= 0.1875, dB1= 0, qP= 0.0625, qH= 0.375, qC= 0.8125, mI= 0.0125, mC= 0.625, mV= 0.6875, sigma= 0.2, kappa= 1, gammaI= 0.2, gammaH= 0.7, gammaC= 0.6625, gammaV= 0.8875, rho= 0.125, TF= 7500, S0= 1, EQ0= 0, PQ0= 0, IQ0= 0, E0= 0, P0= 0, I0= 0, H0= 0, C0= 0, V0= 0, M0= 0 )


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
           est=, method = 'subplex') -> sets_traj;logLik(sets_traj)
Sys.time() -tt 

coef(model) <- c(coef(sets_traj))

mif2(model,
     Nmif=1,
     start=current_params,
     transform=TRUE,
     pars=params_fit,
     rw.sd=rw,
     Np=100,
     var.factor=2,
     cooling.type="hyperbolic",
     cooling.fraction=0.05) -> sets_mf;sets_mf@loglik

coef(model) <- coef(sets_mf)
tt <-Sys.time()
pf <- pfilter(model, filter.mean=T, pred.mean=T,save.states=T,Np=5000);pf@loglik
Sys.time() -tt 

js=read.csv("../../../Gitlab/dcp-r/pomp/samples/filterMean.csv")
p = as.data.frame(pf@filter.mean)
plot(pf@filter.mean[1,])
# points(pf@filter.mean[1,], col="blue")
points(js$S, col="red")



