rm (list=ls())

source("SetValues.R")

modeltype <- "StochasticSEIR"
# args <- commandArgs(TRUE)
# job <- as.numeric(args[1])
# run <- as.numeric(args[2])
run <- 2
job <- 1


no_cores <- c(1)#250
no_points <- c(1)#1000

no_mif <- c(3)#30
no_particles <- 100#5e3
cool_fraction <- 0.05
no_particles_pf <- 10#50e3
no_pf <- 5
rw_size <- 0.05


# Create POMP model
source("CreateDataset2.R")
source("CreateModel.R")
source("DetermineRandomWalks.R")
source("DetermineRunProperties.R")


out <- determine_run_properties(run,modeltype)
params_fixed <- out$param
params_fixed <- c("gamma","sigma",params_fixed)
params_mod_fit <- params_mod[-which(params_mod %in% params_fixed)]
params_ic_fit <- params_ic[-which(params_ic %in% params_fixed)]

ParamSetFile <- paste0("ParamSet_",modeltype,"_run",run,".csv")
full_set <- read.csv(file=ParamSetFile,header=T)
full_set <- full_set[1:no_points,]

n <- ceiling(no_points/no_cores)
ind_start <- 1 + (job-1)*n
ind_end <- job*n

dir1 <- modeltype
dir2 <- paste0(modeltype,"_run",run)
if (file.exists(file.path(mainDir,dir1))) {
  setwd(file.path(mainDir,dir1))
} else {
  dir.create(file.path(mainDir,dir1))
  setwd(file.path(mainDir,dir1))  
}
if (file.exists(file.path(mainDir,dir1,dir2))) {
  setwd(file.path(mainDir,dir1,dir2))
} else {
  dir.create(file.path(mainDir,dir1,dir2))
  setwd(file.path(mainDir,dir1,dir2))  
}

current_set <- mat.or.vec(n,dim(full_set)[2])
current_set <-as.data.frame(current_set)
names(current_set) <- names(full_set)
index <- 0
i=1
for (i in ind_start:ind_end) {
  current_params <- unlist(full_set[i,])
  
  coef(m1) <- c(current_params)
  mif2(m1,
           Nmif=no_mif,
           start=current_params,
           transform=T,
           ivps=params_ic_fit,
           pars=params_mod_fit,
           rw.sd=rw,
           Np=no_particles,
           var.factor=2,
           cooling.type="hyperbolic",
           cooling.fraction=cool_fraction
           ) -> sets_mf
  
    
  current_params <- coef(sets_mf)
  pfilter(m1,params=current_params,Np=no_particles_pf,filter.mean = T,pred.mean=T,pred.var=T,save.params = TRUE, max.fail=3000) -> ss
  loglik_pf <- replicate(n=no_pf,logLik(pfilter(m1,params=current_params,Np=1000,filter.mean = T,pred.mean=T,pred.var=T,save.states = TRUE, max.fail=3000)))
  loglik_pf_est <- logmeanexp(loglik_pf)
  loglik_pf_span <- abs(max(loglik_pf)-min(loglik_pf))
  
  index <- index + 1
  current_set[index,names(current_params)] <- current_params
  # current_set[index,"LogLik"] <- loglik_pf_est
  # current_set[index,"spanLogLik"] <- loglik_pf_span
  # 
  write.csv(current_set,file=paste0(modeltype,"_run",run,"_job",job,".csv"),row.names=F)
}


write.csv(current_set,file=paste0(modeltype,"_run",run,"_job",job,".csv"),row.names=F)

setwd(mainDir)
