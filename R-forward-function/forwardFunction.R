# package "pomp" neeeds to be installed

forward_model <- function(nsim, mainDir ,path_to_params, endTime, predTime ) {
  rm(list=ls())
  setwd(mainDir)
  library(reshape2)
  library(pomp)
  library(ggplot2)
  
  source("ModelSnippet.R")
  nstageH <- modeltype["nstageH"]
  nstageC <- modeltype["nstageC"]
  nstageV <- modeltype["nstageV"]
  
  convertDate <- function (date, startTime = "2019-12-31") {
    out <- as.numeric(as.Date(date)) - as.numeric(as.Date(startTime))
    return(out)
  }
  
  #SelectSample.R : read the parameters from RT
  set <- read.csv(file.path(path_to_params))
  mle <- unlist(set)
  
  # Upload states at time[end] for all Np points.
  last_saved_states <- as.data.frame(read.csv(file.path(paste0(mainDir,"/../results/savedStates.csv"))))
  Np <- nrow(last_saved_states)
  
  source("CreateModel.R")
  source("CreateCovars.R")
  source("CreateFutureDataset.R")
  
  full_data <- create_dataset(paste0(mainDir,"/../samples/data.csv"), endTime=endTime,predTime=predTime)
  add_data <- subset(full_data,subset = time>convertDate(endTime) & time<=convertDate(date=predTime))
  covars <- as.data.frame(read.csv(paste0(mainDir,"/../samples/covar.csv")))
  out <- create_pomp_model(data = add_data, covars=covars, t0 = convertDate(endTime), dt=0.1)
  model <- out$model
  
  unifs <- round(1+runif(nsim)*(Np-1))
  sims <- NULL
  
  for (q in 1:nsim) {
    temp <- last_saved_states[unifs[q],]
    
    S0 <- temp$S
    EQ0 <- sum(temp[paste0("EQ",seq(1,modeltype["nstageE"]))])
    PQ0 <- sum(temp[paste0("PQ",seq(1,modeltype["nstageP"]))])
    IQ0 <- sum(temp[paste0("IQ",seq(1,modeltype["nstageI"]))])
    E0 <- sum(temp[paste0("E",seq(1,modeltype["nstageE"]))])
    P0 <- sum(temp[paste0("P",seq(1,modeltype["nstageP"]))])
    I0 <- sum(temp[paste0("I",seq(1,modeltype["nstageI"]))])
    H0 <- sum(temp[paste0("H",seq(1,modeltype["nstageH"]))])
    C0 <- sum(temp[paste0("C",seq(1,modeltype["nstageC"]))])
    V0 <- sum(temp[paste0("V",seq(1,modeltype["nstageV"]))])
    M0 <- temp$M
    R0 <- temp$R
    
    ic <- c(S0,EQ0,PQ0,IQ0,E0,P0,I0,H0,C0,V0,M0)
    ic <- ic / (sum(ic)+R0)
    paramset <- mle[c(params_mod,params_ic)]     
    paramset[params_ic] <- ic               
    coef(model) <- unlist(paramset)
    x <- simulate(model,nsim=1,as.data.frame=TRUE)
    pred <- x[,c("time", "reports", "deaths", "hospital", "ICU", "ventilator")]
    pred$sim <- q
    sims <- rbind(sims,pred)
  }
  write.csv(sims, "./simulations/prediction.csv", row.names=FALSE)
  
}

nsim = 1e2
mainDir = "~/Gitlab/epidemiological-model/R-forward-function"
path_to_params = paste0(mainDir,"/../results/modelParams.csv")
endTime =c("2020-09-07")
predTime=c("2020-12-31")
forward_model(nsim, mainDir, path_to_params, endTime, predTime)
