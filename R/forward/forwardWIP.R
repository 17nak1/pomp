rm(list=ls())
mainDir <- setwd("~/Git/pomp/R/forward")
setwd(mainDir)
set <- read.csv(file.path(mainDir,"sample/sample.csv"))
RT_paramSet <- unlist(set[1,])
nsim <- 2#1e2
forward_function <- function(RT_paramSet, nsim) {
  library(pomp)
  source("ModelSnippet.R")
  predTime <- c("2020-12-31")
  nstageH <- modeltype["nstageH"]
  nstageC <- modeltype["nstageC"]
  nstageV <- modeltype["nstageV"]
  convertDate <- function (date, startTime = "2019-12-31") {
    out <- as.numeric(as.Date(date)) - as.numeric(as.Date(startTime))
    return(out)
  }
  
  # set <- read.csv(file.path(mainDir,"sample/sample.csv"))
  # RT_paramSet <- unlist(set[1,])
  # mle["beta_sd"] <- 0.01
  # mle["dB1"] <- 0.2
  
  last_saved_states <- as.data.frame(read.csv("./sample/savedStates.csv"))
  Np <- nrow(last_saved_states)
  
  
  source("CreateModel.R")
  source("CreateCovars.R")
  source("CreateDataset.R")
  
  full_data <- create_dataset(endTime=endTime,predTime=predTime)
  add_data <- subset(full_data,subset = time>convertDate(endTime) & time<=convertDate(date=predTime))
  covars <- create_covars(endTime=endTime)
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
    paramset <- RT_paramSet # here add RT params
    paramset[params_ic] <- ic
    coef(model) <- unlist(paramset)
    x <- simulate(model,nsim=1,as.data.frame=TRUE)
    pred <- x[,c("time", "reports", "deaths", "hospital", "ICU", "ventilator", "M", "R", "casesI", "casesIQ", "casesH", "deathsIIQ", "deathsCV", "S")]
    pred$sim <- q
    sims <- rbind(sims,pred)
  }
  rm(temp, S0, EQ0, PQ0, IQ0, E0, P0, I0, H0, C0, V0, M0, x, pred, unifs, q, ic)
  
  full_data$M <- cumsum(full_data$deaths)
  full_data$sim <- nsim + 1
  states <- as.data.frame(read.csv("./sample/filterMean.csv"))
  tests <- approx(covars$time,covars$tests,xout=seq(1,nrow(states)))$y
  tests[is.na(tests)] <- 1
  
  temp <- data.frame(seq(1,nrow(states),by=1),states[,"casesH"]+coef(model,"rho")*tests/(tests+coef(model,"TF"))*states[,"casesI"])
  temp$sim <- nsim+2
  colnames(temp) <- c("time","reports", "sim")
  df_reports <- subset(sims, select=c("time","reports","sim"))
  df_reports$sim <- as.numeric(df_reports$sim)
  df_reports <- rbind(df_reports,full_data[,c("time","reports","sim")],temp)
  df_reports$sim <- as.factor(df_reports$sim)
  
  temp <- data.frame(seq(1,nrow(states),by=1),states[,"deathsIIQ"]+states[,"deathsCV"])
  temp$sim <- nsim+2
  colnames(temp) <- c("time","deaths", "sim")
  df_deaths <- subset(sims, select=c("time","deaths","sim"))
  df_deaths$sim <- as.numeric(df_deaths$sim)
  df_deaths <- rbind(df_deaths,full_data[,c("time","deaths","sim")],temp)
  df_deaths$sim <- as.factor(df_deaths$sim)
  
  temp <- data.frame(seq(1,nrow(states),by=1),states[,"M"])
  temp$sim <- nsim+2
  colnames(temp) <- c("time","M", "sim")
  df_Total_deaths <- subset(sims, select=c("time","M","sim"))
  df_Total_deaths$sim <- as.numeric(df_Total_deaths$sim)
  df_Total_deaths <- rbind(df_Total_deaths,full_data[,c("time","M","sim")],temp)
  df_Total_deaths$sim <- as.factor(df_Total_deaths$sim)
  
  
  temp <- data.frame(seq(1,nrow(states),by=1),
                     rowSums(states[,c(paste0("H",seq(1,nstageH)),paste0("C",seq(1,nstageC)),paste0("V",seq(1,nstageC)))]),
                     rowSums(states[,c(paste0("C",seq(1,nstageC)),paste0("V",seq(1,nstageC)))]),
                     rowSums(states[,paste0("V",seq(1,nstageC))]))
  temp$sim <- nsim+3
  colnames(temp) <- c("time","hospital", "ICU", "ventilator", "sim")
  df_pred_hospital <- subset(sims,select=c(time,hospital,ICU,ventilator,sim))
  df_pred_hospital$sim <- as.numeric(df_pred_hospital$sim)
  df_pred_hospital <- rbind(df_pred_hospital,full_data[,c("time","hospital","ICU","ventilator", "sim")],temp)
  df_pred_hospital$sim <- as.factor(df_pred_hospital$sim)
  df_pred_hospital <- melt(df_pred_hospital,id=c("time","sim"))
  list("pred_reports" = df_reports,
       "pred_deaths"= df_deaths,
       "pred_total_deaths" = df_Total_deaths,
       "pred_hospital" = df_pred_hospital)
}

prediction = forward_function(RT_paramSet, nsim)
library(reshape2)
library(ggplot2)


theme_set(theme_bw()+theme(legend.position = "none"))
cols <- c(rep("gray", nsim),"red","blue")
p <- ggplot(prediction$pred_reports,aes(x=time,y=reports,group=sim)) +
  geom_step(aes(color=sim),alpha=0.6) +
  scale_color_manual(values= cols) +
  geom_vline(xintercept=T1, linetype=1,colour="#808080")  + coord_cartesian(ylim=c(0,1000))
pdf(file=paste0("pred_reports_DM.pdf"))
plot(p)
dev.off()
plot(p)


# Plot new reported deaths

p <- ggplot(prediction$prediction_deaths,aes(x=time,y=deaths,group=sim)) +
  geom_step(aes(color=sim),alpha=0.6) +
  scale_color_manual(values= cols) +
  geom_vline(xintercept=T1, linetype=1,colour="#808080")  + coord_cartesian(ylim=c(0, 100))
pdf(file=paste0("pred_deaths_DM.pdf"))
plot(p)
dev.off()
plot(p)

# Plot total true deaths

p <- ggplot(reports$pred_total_deaths,aes(x=time,y=M,group=sim)) +
  geom_step(aes(color=sim),alpha=0.6) +
  scale_color_manual(values= cols) +
  geom_vline(xintercept=T1, linetype=1,colour="#808080")  + coord_cartesian(ylim=c(0, 3000)) +
  ylab("Total deaths")
pdf(file=paste0("Total_deaths_DM.pdf"))
plot(p)
dev.off()
plot(p)


# Plot ICU and hospitalization
p <- ggplot(prediction$pred_hospital,aes(x=time,y=value,group=sim,col=sim)) +
  geom_step(aes(color=sim),alpha=0.6) +
  facet_grid(variable ~ .,scale="free_y")  +
  scale_color_manual(values= cols)
pdf(file=paste0("pred_hospital_DM.pdf"))
plot(p)
dev.off()
plot(p)
# 
# 
# 
# 
# # Plot susceptible
# df <- subset(sims, select=c(time,S,sim))
# df$sim <- as.numeric(df$sim)
# temp <- data.frame(seq(1,nrow(states),by=1),states["S",])
# temp$sim <- nsim+2
# colnames(temp) <- c("time","S", "sim")
# df <- rbind(df,temp)
# df$sim <- as.factor(df$sim)
# p <- ggplot(df,aes(x=time,y=S/(10e6),group=sim)) +
#   geom_step(aes(color=sim),alpha=0.6) +
#   scale_color_manual(values= cols)  + ylab("Fraction susceptible") +
#   geom_vline(xintercept=T1, linetype=1,colour="#808080") + coord_cartesian(ylim=c(0.5, 1))
# plot(p)
# rm(temp)
