rm(list=ls())

source("SetValues.R")
runFilter <- T
samplename <- "Test0"
theme_set(theme_bw()+theme(legend.position = "none"))
Np <- 1e3
nsim <- 1e2
predTime <- c("2020-12-31")
nstageH <- modeltype["nstageH"]
nstageC <- modeltype["nstageC"]
nstageV <- modeltype["nstageV"]

cols <- c(rep("gray", nsim),"red","blue")

convertDate <- function (date, startTime = "2019-12-31") {
  out <- as.numeric(as.Date(date)) - as.numeric(as.Date(startTime))
  return(out)
}

modelname <- paste0("DetModel",modelnumber)
source("SelectSample.R")
mle <- c(betaI= 1.40758806343178, iota= 0.0177042779901833, beta_sd= 0, sigma= 0.2, kappa= 1, gammaI= 0.172133302616885, gammaH= 0.635123961314889, gammaC= 1.31234713108046, gammaV= 0.327003197582952, TF= 16031.6700539306, rho= 0.487576658764406, theta= 0.415700726852156, dI0= 0.485528264998708, dP0= 0.176249510892225, dT0= 0.263185388981837, dB0= 0, dI1= 0.331872295141556, dP1= 0.259732942470986, dT1= 0.630072923852895, dB1= 0, qP= 0.722056284667341, qH= 0.137405612125941, qC= 0.883271524025973, mI= 0.00550468107974399, mC= 0.197242061100997, mV= 0.685825470240576, S0= 1, EQ0= 0, PQ0= 0, IQ0= 0, E0= 0, P0= 0, I0= 0, H0= 0, C0= 0, V0= 0, M0= 0)

# params
mle["beta_sd"] <- 0.01
mle["dB0"] <- 0.2
mle["dB1"] <- 0.2
rm(params)

# if (runFilter) {
source("RunModel.R")
coef(model) <- mle[c(params_mod,params_ic)]
pf <- pfilter(model, filter.mean=T,save.states=T,Np=Np)

times <- pf$times
pf_table <- pf@saved.states
Np <- ncol(pf_table[[1]])

modelname <- paste0("DetModel",modelnumber)
full_data <- create_dataset(endTime=endTime,predTime=predTime)
covars <- create_covars(endTime)
add_data <- subset(full_data,subset = time>convertDate(endTime) & time<=convertDate(date=predTime))
out <- create_pomp_model(data = add_data, covars=covars, t0 = convertDate(endTime), dt=0.1)
model <- out$model


ind <- length(times)
unifs <- round(1+runif(nsim)*(Np-1))
sims <- NULL
q=1

for (q in 1:nsim) {
  temp <- pf_table[[ind]]
  temp <- as.data.frame(temp[,unifs[q]])
  if (modelname%in% c("DetModel3")) {
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
  }
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
rm(temp, S0, EQ0, PQ0, IQ0, E0, P0, I0, H0, C0, V0, M0, x, pred, ind, unifs, times, q, ic)

full_data$M <- cumsum(full_data$deaths)
full_data$sim <- nsim + 1


states <- pf@filter.mean

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



