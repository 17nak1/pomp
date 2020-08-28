set <- read.csv(file.path(mainDir,modelname,paste0(modelname,"_all.csv")))


if (samplename=="Test0") {
  params <- unlist(set[1,])
} else if (samplename=="Test1") {
  set <- subset(set,TF>600) 
  params <- unlist(set[1,])
} else if (samplename=="Test2") {
  set <- subset(set,betaI<0.5 & dI<0.5) 
  params <- unlist(set[1,])
} else if (samplename=="Test3") {
  set$RC <- NA
  for (j in 1:nrow(set)) {
    
    betaI <- set[j,"betaI"]
    gammaI <- set[j,"gammaI"]
    qP <- set[j,"qP"]
    betaP <- set[j,"betaP"]
    kappa <- set[j,"kappa"]
    dP <- set[j,"dP"]
    dI <- set[j,"dI"]
    set[j,"RC"] <-(1-dI)*betaI/gammaI*(1-qP)+(1-dP)*betaI*betaP/kappa
  }
  set <- subset(set,RC<1)
  
  # set <- subset(set,(1-dI)*betaI/gammaI*(1-qP)+(1-dP)*betaI*betaP/kappa<0.1) 
  params <- unlist(set[1,])
}
# rm(set)

