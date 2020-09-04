rm(list=ls())
mainDir <- "C:/Users/Amir/Downloads/stochmodel3"
setwd(mainDir)
source("ModelSnippet_StochModel3.R")

Np <- 1e3
predTime <- c("2020-12-31")
modelname <- "DetModel3"
rw_size <- 0.05
mle <- z


source("CreateModel.R")
source("CreateCovars.R")
source("CreateDataset.R")


data <- create_dataset(endTime=endTime)
covars <- create_covars(endTime=endTime)
out <- create_pomp_model(data = data, covars=covars, t0 = 0, dt=0.1)

model <- out$model

coef(model) <- mle[c(params_mod,params_ic)]
tt <- Sys.time()
pf <- pfilter(model, filter.mean=T, pred.mean=T,save.states=T,Np=100);pf@loglik
Sys.time()-tt
js=read.csv("oo.csv")
p = as.data.frame(pf@filter.mean)
plot(pf@filter.mean[20,])
points(pf@filter.mean[20,], col="red")

points(js$H1, col="blue")
#setwd("~/Git/pomp/R")
#run=1
#source("DetermineRandomWalks.R")
#tt <- Sys.time()
#mif2(model,
#     Nmif=1,
#     start=mle,
#     transform=TRUE,
#
#     rw.sd=rw,
#     Np=1000,
#     var.factor=2,
#     cooling.type="hyperbolic",
#     cooling.fraction=0.05) -> sets_mf;coef(sets_mf);sets_mf@loglik
# Sys.time()-tt
# js=read.csv("oo.csv")
# p = as.data.frame(pf@filter.mean)
# plot(pf@filter.mean[2,])
# points(pf@filter.mean[2,], col="blue")
#
# points(js$EQ1, col="red")
#
#
#
