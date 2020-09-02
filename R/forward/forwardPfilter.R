rm(list=ls())

mainDir <- setwd("~/Git/pomp/R/forward")
setwd(mainDir)

library(reshape2)
library(pomp)
library(ggplot2)

# modelname <- "model"
source("ModelSnippet.R")


# samplename <- "Test0"
theme_set(theme_bw()+theme(legend.position = "none"))
rw_size <- 0.05
nsim <- 1e2
predTime <- c("2020-12-31")
nstageH <- modeltype["nstageH"]
nstageC <- modeltype["nstageC"]
nstageV <- modeltype["nstageV"]
T0 = 75;
T1 = 139;
rw <- rw.sd(betaI=ifelse(time<T0,0,0),
            theta=ifelse(time<T0,0,0),
            iota=ifelse(time<T0,0,0),
            beta_sd=ifelse(time<T0,rw_size,0),
            dI0=ifelse(time<T0,0,0),
            dP0=ifelse(time<T0,0,0),
            dT0=ifelse(time<T0,0,0),
            dB0=ifelse(T0<time,ifelse(time<T1,rw_size,0),0),
            dI1=ifelse(time<T0,0,0),
            dP1=ifelse(time<T0,0,0),
            dT1=ifelse(time<T0,0,0),
            dB1=ifelse(time<T1,0,rw_size),
            qP=ifelse(time<T0,0,0),
            qH=ifelse(time<T0,0,0),
            qC=ifelse(time<T0,0,0),
            mI=ifelse(time<T0,0,0),
            mC=ifelse(time<T0,0,0),
            mV=ifelse(time<T0,0,0),
            sigma=ifelse(time<T0,0.5*0,0.5*0),
            kappa=ifelse(time<T0,0.5*0,0.5*0),
            gammaI=ifelse(time<T0,0.5*0,0.5*0),
            gammaH=ifelse(time<93,0,0),
            gammaC=ifelse(time<93,0,0),
            gammaV=ifelse(time<93,0,0),
            rho=ifelse(time<T0,0,0),
            TF=ifelse(time<35,0,0))
cols <- c(rep("gray", nsim),"red","blue")

convertDate <- function (date, startTime = "2019-12-31") {
  out <- as.numeric(as.Date(date)) - as.numeric(as.Date(startTime))
  return(out)
}

mle <- c(betaI= 2.46744881887213, theta= 0.9999999999999872, iota= 5.334242965238772e-7, beta_sd= 0, dI0= 1, dP0= 0.9999999999999925, dT0= 1, dB0= 0, dI1= 0.4585835582162743, dP1= 1, dT1= 2.7738922270259536e-13, dB1= 0, qP= 0.16871281967335705, qH= 0.9999999999999998, qC= 0.45543754894231614, mI= 0, mC= 0, mV= 1, sigma= 0.2, kappa= 1, gammaI= 0.027210002650336837, gammaH= 0.11042090059074512, gammaC= 2.0496759026297866, gammaV= 0.23293202582712966, rho= 0.9911877804418902, TF= 22.090915777310727, S0= 1, EQ0= 0, PQ0= 0, IQ0= 0, E0= 0, P0= 0, I0= 0, H0= 0, C0= 0, V0= 0, M0= 0 )
mle["beta_sd"] <- 0.01
mle["dB0"] <- 0.2
mle["dB1"] <- 0.2
source("CreateModel.R")
source("CreateCovars.R")
source("CreateDataset.R")

full_data <- create_dataset(endTime=endTime,predTime=predTime)
add_data <- subset(full_data,subset = time>convertDate(endTime) & time<=convertDate(date=predTime))
covars <- create_covars(endTime=endTime)
out <- create_pomp_model(data = add_data, covars=covars, t0 = convertDate(endTime), dt=0.1)
model <- out$model
# coef(model) <- c(mle)
# mif2(model,
#      Nmif=1,
#      start=mle,
#      transform=TRUE,
#      rw.sd=rw,
#      Np=1000,
#      var.factor=2,
#      cooling.type="hyperbolic",
#      cooling.fraction=0.05) -> sets_mf


coef(model) <- c(mle)
pf <- pfilter(model, filter.mean=T,save.states=T,Np=1000,max.fail=300)
times <- pf$times
pf_table <- pf@saved.states
Np <- ncol(pf_table[[1]])
states <- pf@filter.mean
ind <- length(times)
unifs <- round(1+runif(nsim)*(Np-1))
sims <- NULL

for (q in 1:nsim) {
  temp <- pf_table[[ind]]
  temp <- as.data.frame(temp[,unifs[q]])
  
  S0 <- temp["S",]
  EQ0 <- sum(temp[paste0("EQ",seq(1,modeltype["nstageE"])),])
  PQ0 <- sum(temp[paste0("PQ",seq(1,modeltype["nstageP"])),])
  IQ0 <- sum(temp[paste0("IQ",seq(1,modeltype["nstageI"])),])
  E0 <- sum(temp[paste0("E",seq(1,modeltype["nstageE"])),])
  P0 <- sum(temp[paste0("P",seq(1,modeltype["nstageP"])),])
  I0 <- sum(temp[paste0("I",seq(1,modeltype["nstageI"])),])
  H0 <- sum(temp[paste0("H",seq(1,modeltype["nstageH"])),])
  C0 <- sum(temp[paste0("C",seq(1,modeltype["nstageC"])),])
  V0 <- sum(temp[paste0("V",seq(1,modeltype["nstageV"])),])
  M0 <- temp["M",]
  R0 <- temp["R",]
  
  ic <- c(S0,EQ0,PQ0,IQ0,E0,P0,I0,H0,C0,V0,M0)
  ic <- ic / (sum(ic)+R0)
  paramset <- mle[c(params_mod,params_ic)]
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

# Plot new reported deaths
tests <- approx(covars$time,covars$tests,xout=seq(1,ncol(states)))$y
tests[is.na(tests)] <- 15
temp <- data.frame(seq(1,ncol(states),by=1),states["casesH",]+coef(model,"rho")*tests/(tests+coef(model,"TF"))*states["casesI",])
temp$sim <- nsim+2
colnames(temp) <- c("time","reports", "sim")

df <- subset(sims, select=c("time","reports","sim"))
df$sim <- as.numeric(df$sim)
df <- rbind(df,full_data[,c("time","reports","sim")],temp)
df$sim <- as.factor(df$sim)
p <- ggplot(df,aes(x=time,y=reports,group=sim)) +
  geom_step(aes(color=sim),alpha=0.6) +
  scale_color_manual(values= cols) +
  geom_vline(xintercept=T1, linetype=1,colour="#808080")  + coord_cartesian(ylim=c(0,1000))
plot(p)
