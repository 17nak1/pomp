rm(list=ls())
mainDir <- "~/Git/pomp/R"
setwd(mainDir)
source("ModelSnippet_StochModel3.R")

predTime <- c("2020-12-31")
modelname <- "DetModel3"
rw_size <- 0.05
current_params <- c( betaI= 1.40758806343178, iota= 0.0177042779901833, beta_sd= 0, sigma= 0.2, kappa= 1, gammaI= 0.172133302616885, gammaH= 0.635123961314889, gammaC= 1.31234713108046, gammaV= 0.327003197582952, TF= 16031.6700539306, rho= 0.487576658764406, theta= 0.415700726852156, dI0= 0.485528264998708, dP0= 0.176249510892225, dT0= 0.263185388981837, dB0= 0, dI1= 0.331872295141556, dP1= 0.259732942470986, dT1= 0.630072923852895, dB1= 0, qP= 0.722056284667341, qH= 0.137405612125941, qC= 0.883271524025973, mI= 0.00550468107974399, mC= 0.197242061100997, mV= 0.685825470240576, S0= 1, EQ0= 0, PQ0= 0, IQ0= 0, E0= 0, P0= 0, I0= 0, H0= 0, C0= 0, V0= 0, M0= 0)
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



