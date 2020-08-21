rm(list=ls())
mainDir <- "~/Git/pomp/R"
setwd(mainDir)
source("ModelSnippet_StochModel3.R")

Np <- 1e3
predTime <- c("2020-12-31")
modelname <- "DetModel3"
rw_size <- 0.05
current_param <- c( betaI= 0.375, theta= 0.375, iota= 62.5, beta_sd= 0, dI0= 0.525, dP0= 0.525, dT0= 0.075, dB0= 0, dI1= 0.0625, dP1= 0.4375, dT1= 0.1875, dB1= 0, qP= 0.0625, qH= 0.375, qC= 0.8125, mI= 0.0125, mC= 0.625, mV= 0.6875, sigma= 0.2, kappa= 1, gammaI= 0.2, gammaH= 0.7, gammaC= 0.6625, gammaV= 0.8875, rho= 0.125, TF= 7500, S0= 1, EQ0= 0, PQ0= 0, IQ0= 0, E0= 0, P0= 0, I0= 0, H0= 0, C0= 0, V0= 0, M0= 0 )


source("CreateModel.R")
source("CreateCovars.R")
source("CreateDataset.R")


data <- create_dataset(endTime=endTime)
covars <- create_covars(endTime=endTime)
out <- create_pomp_model(data = data, covars=covars, t0 = 0, dt=0.1)

model <- out$model
params_fixed <- c("beta_sd","dB0", "dB1","sigma","kappa")

# coef(model) <- current_param[c(params_mod,params_ic)]
coef(model) <- c(current_params)
traj.match(model,
           start=current_params,
           transform=TRUE,
           est=params_fit, method = 'subplex') -> sets_traj;logLik(sets_traj)
tt <- Sys.time()

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
     cooling.fraction=cool.fraction) -> sets_mf;sets_mf@loglik

coef(model) <- c(sets_mf)
pf <- pfilter(model, filter.mean=T, pred.mean=T,save.states=T,Np=100);pf@loglik
Sys.time()-tt

js=read.csv("../trajMatch/oo.csv")
p = as.data.frame(pf@filter.mean)
plot(pf@filter.mean[20,])
# points(pf@filter.mean[20,], col="red")

points(js$H1, col="red")
setwd("~/Git/pomp/R")
run=1
source("DetermineRandomWalks.R")
tt <- Sys.time()
mif2(model,
     Nmif=1,
     start=mle,
     transform=TRUE,

     rw.sd=rw,
     Np=1000,
     var.factor=2,
     cooling.type="hyperbolic",
     cooling.fraction=0.05) -> sets_mf;coef(sets_mf);sets_mf@loglik
 Sys.time()-tt
 js=read.csv("oo.csv")
 p = as.data.frame(pf@filter.mean)
 plot(pf@filter.mean[2,])
 points(pf@filter.mean[2,], col="blue")

 points(js$EQ1, col="red")


