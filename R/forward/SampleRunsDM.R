rm(list=ls())
setwd("~/Downloads/stochmodel3")
source("SetValues.R")

save_plot <- T
samplename <- "Test0"
predTime <- c("2020-08-31")
nsim <- 200L
nstageH <- modeltype["nstageH"]
nstageC <- modeltype["nstageC"]
nstageV <- modeltype["nstageV"]
cols <- c(rep("gray", nsim),"red")

theme_set(theme_bw()+theme(legend.position = "none"))



source("CreateModel.R")
source("CreateCovars.R")
source("CreateDataset.R")


data <- create_dataset(endTime=endTime,predTime=predTime)
covars <- create_covars(endTime=endTime)
out <- create_pomp_model(data = data, covars=covars, t0 = 0, dt=0.1)

model <- out$model
rm(out)

modelname <- paste0("DetModel",modelnumber)
source("SelectSample.R")


coef(model) <- params[c(params_mod,params_ic)]
coef(model,"beta_sd") <- 0.01
coef(model,"dB") <- 0.2



# Test transforms
# Check that the transforms work
sample1 <- coef(model)
sample2 <- coef(model, transform=TRUE)
model2 <- model
coef(model2, transform=TRUE) <- sample2
print(all.equal(coef(model), coef(model2)))
rm(sample1,sample2,model2)



x <- simulate(model,nsim=nsim,as.data.frame=TRUE)
data[,"sim"] <- nsim+1


# Check to make sure the population is not growing
df <- subset(x,sim==1)
df$tot <- 0
sts <- statenames()
for (st in sts) {
  df$tot <- df$tot+df[,st]
}
ggplot(df) + geom_line(aes(x=time,y=tot))


# Check to make sure the simulations are reasonable.
df <- subset(x, select=c(time,reports,sim))
df$sim <- as.numeric(df$sim)

df <- rbind(df,data[,c("time","reports","sim")])
df$sim <- as.factor(df$sim)
p <- ggplot(df,aes(x=time,y=reports,group=sim)) +
  geom_step(aes(color=sim),alpha=0.6) +
  scale_color_manual(values= cols) + ylab(expression("reports")) +
  geom_vline(xintercept=T0, linetype=1,colour="#808080") +
  geom_vline(xintercept=T1, linetype=1,colour="#808080") +
  coord_cartesian(ylim=c(0, 1000))
# p <- ggplot(df,aes(x=time,y=sqrt(reports),group=sim)) +
#   geom_step(aes(color=sim),alpha=0.6) +
#   scale_color_manual(values= cols) + ylab(expression(sqrt("reports"))) +
#   geom_vline(xintercept=t1, linetype=1,colour="#808080") + coord_cartesian(ylim=c(0, 100))
pdf(file=paste0("reports_",modelname,"_DM.pdf"))
plot(p)
dev.off()
plot(p)



# Plot new deaths
df <- subset(x, select=c("time","deaths","sim"))
df$sim <- as.numeric(df$sim)
df <- rbind(df,data[,c("time","deaths","sim")])
df$sim <- as.factor(df$sim)
p <- ggplot(df,aes(x=time,y=deaths,group=sim)) +
  geom_step(aes(color=sim),alpha=0.6) +
  scale_color_manual(values= cols) +
  geom_vline(xintercept=T0, linetype=1,colour="#808080") +
  geom_vline(xintercept=T1, linetype=1,colour="#808080") +
  coord_cartesian(ylim=c(0,100))
pdf(file=paste0("deaths_",modelname,"_DM.pdf"))
plot(p)
dev.off()
plot(p)




# Plot ICU and hospitalization and ventilators
df <- subset(x,select=c(time,hospital,ICU,ventilator,sim))
df$sim <- as.numeric(df$sim)
df <- rbind(df,data[,c("time","hospital","ICU","ventilator", "sim")])
df$sim <- as.factor(df$sim)
df <- melt(df,id=c("time","sim"))
p <- ggplot(df,aes(x=time,y=value,group=sim,col=sim)) +
  geom_step(aes(color=sim),alpha=0.6) +
  facet_grid(variable ~ .,scale="free_y") +
  scale_color_manual(values= cols)
pdf(file=paste0("hospital_",modelname,"_DM.pdf"))
plot(p)
dev.off()
plot(p)

library(dplyr)

# Check susceptible
df <- subset(x, select=c(time,S,sim))
group_by(df,time) %>% summarise(mean=mean(S)) -> df2
colnames(df2) <- c("time","S","n")
df2 <- data.frame(df2)

df$sim <- as.factor(df$sim)
p <- ggplot() +
  geom_step(data=df,aes(x=time,y=S/(10e6),group=sim,color=sim),alpha=0.6) +
  scale_color_manual(values= cols)  + ylab("Fraction susceptible") +
  geom_vline(xintercept=T0, linetype=1,colour="#808080") +
  geom_vline(xintercept=T1, linetype=1,colour="#808080") +
  geom_line(data=df2,aes(x=time,y=S/(10e6)))
plot(p)
detach("package:dplyr",unload=T)
pdf(file=paste0("susc_",modelname,"_DM.pdf"))
plot(p)
dev.off()
plot(p)


# Compare reporting vs cases
data <- create_dataset(endTime=endTime)
covars <- create_covars(endTime=endTime)
ggplot() + geom_step(data=covars,aes(x=time,y=tests)) +
  geom_step(data=data,aes(x=time,y=10*reports),col="red",alpha=0.5) +
  geom_line(data=covars,aes(x=time,y=4000*tests/(tests+coef(model,"TF"))),col="green",alpha=0.5)


print(params)
betaI <- coef(model)["betaI"]
gammaI <- coef(model)["gammaI"]
qP <- coef(model)["qP"]
theta <- coef(model)["theta"]
kappa <- coef(model,"kappa")
dP0 <- coef(model)["dP0"]
dI0 <- coef(model)["dI0"]
dP1 <- coef(model)["dP1"]
dI1 <- coef(model)["dI1"]



# R0 plot
df <- subset(x, subset=sim==1,select=c(time,S,tests,sim))
df$tests[is.na(df$tests)] <- 15
R0 <- NULL
Reff <- NULL
for (k in 1:nrow(df)) {
  t <- df[k,"time"]
  test <- df[k,"tests"]
  p <- test/(test+coef(model,"TF"))*coef(model,"rho")
  if (t<T0) {
    temp <- (1-p)*(theta*betaI/kappa +(1-qP)*betaI/gammaI)
  } else if (t<T0+7) {
    s <- (t-T0) / 7.0;
    ss <- 3*s*s - 2*s*s*s;
    temp <- (1-p)*(((1-ss) + ss*dP0)* theta*betaI/kappa +(1-qP)*((1-ss) + ss*dI0) *betaI/gammaI)
  } else if (t<T1) {
    temp <- (1-p)*(dP0*theta*betaI/kappa +(1-qP)*dI0*betaI/gammaI)
  } else if (t<T1+7.0) {
    s <- (t-T1) / 7.0;
    ss <- 3*s*s - 2*s*s*s;
    temp <- (1-p)*(((1-ss)*dP0 + ss*dP1)* theta*betaI/kappa +(1-qP)*((1-ss)*dI0 + ss*dI1) *betaI/gammaI)
  } else {
    temp <- (1-p)*(dP1*theta*betaI/kappa +(1-qP)*dI1*betaI/gammaI)
  }
  R0 <- c(R0,temp)
  Reff <- c(Reff,temp*df[k,"S"]/pop)
}
dat <- data.frame(time=df$time,R0,Reff)
p <- ggplot(dat) + geom_line(aes(x=time,y=R0),col="black")  +
  geom_line(aes(x=time,y=Reff),col="blue") +
  coord_cartesian(ylim=c(0, 6)) + geom_abline(slope=0,intercept=1,col="gray") +
  theme(legend.position="west")
pdf(file=paste0("R0_",modelname,"_DM.pdf"))
plot(p)
dev.off()
plot(p)

print(paste0("R0=",R0[1]))
print(paste0("R0c=",Reff[length(Reff)]))




# # Plot histograms
# a <- 40
# b <- data[which(data$time==a),"reports"]
# ggplot(subset(x,select=reports,subset=time==a)) + geom_bar(aes(x=reports)) +
#   geom_vline(xintercept=b,col="red") +
#   coord_cartesian(xlim=c(0,max(b,2*a)))
