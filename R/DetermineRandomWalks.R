# Set random walk values. All fixed parameter rw.sd will eventually be removed

# if (run==1) {
  rw <- rw.sd(betaI=ifelse(time<75,rw_size,0),
              theta=ifelse(time<75,rw_size,0),
              iota=ifelse(time<75,rw_size,0),
              beta_sd=ifelse(time<75,rw_size,0),
              dI0=ifelse(time<75,0,rw_size),
              dP0=ifelse(time<75,0,rw_size),
              dT0=ifelse(time<75,0,rw_size),
              dB0=ifelse(time<75,0,rw_size),
              dI1=ifelse(time<75,0,rw_size),
              dP1=ifelse(time<75,0,rw_size),
              dT1=ifelse(time<75,0,rw_size),
              dB1=ifelse(time<75,0,rw_size),
              qP=ifelse(time<75,rw_size,rw_size),
              # qI=ifelse(time<75,rw_size,rw_size),
              qH=ifelse(time<75,rw_size,rw_size),
              qC=ifelse(time<75,rw_size,rw_size),
              # qV=ifelse(time<75,rw_size,rw_size),
              mI=ifelse(time<75,rw_size,rw_size),
              mC=ifelse(time<75,rw_size,rw_size),
              mV=ifelse(time<75,rw_size,rw_size),
              sigma=ifelse(time<75,0.5*rw_size,0.5*rw_size),
              kappa=ifelse(time<75,0.5*rw_size,0.5*rw_size),
              gammaI=ifelse(time<75,0.5*rw_size,0.5*rw_size),
              gammaH=ifelse(time<93,0,rw_size),
              gammaC=ifelse(time<93,0,rw_size),
              gammaV=ifelse(time<93,0,rw_size),
              rho=ifelse(time<75,rw_size,rw_size),
              TF=ifelse(time<35,0,rw_size))
# } else if (run==2) {
#   rw <- rw.sd(betaI=ifelse(time<t1,0,0),
#               theta=ifelse(time<t1,rw_size,0),
#               iota=ifelse(time<t1,rw_size,0),
#               beta_sd=ifelse(time<t1,rw_size,0),
#               dI=ifelse(time<t1,0,rw_size),
#               dP=ifelse(time<t1,0,rw_size),
#               dT=ifelse(time<t1,0,rw_size),
#               dB=ifelse(time<t1,0,rw_size),
#               qP=ifelse(time<t1,rw_size,rw_size),
#               qI=ifelse(time<t1,rw_size,rw_size),
#               qH=ifelse(time<t1,rw_size,rw_size),
#               qC=ifelse(time<t1,rw_size,rw_size),
#               qV=ifelse(time<t1,rw_size,rw_size),
#               sigma=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               kappa=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaI=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaH=ifelse(time<93,0,rw_size),
#               gammaC=ifelse(time<93,0,rw_size),
#               gammaV=ifelse(time<93,0,rw_size),
#               rho=ifelse(time<t1,rw_size,rw_size),
#               TF=ifelse(time<35,0,rw_size))
# } else if (run==3) {
#   rw <- rw.sd(betaI=ifelse(time<t1,rw_size,0),
#               theta=ifelse(time<t1,0,0),
#               iota=ifelse(time<t1,rw_size,0),
#               beta_sd=ifelse(time<t1,rw_size,0),
#               dI=ifelse(time<t1,0,rw_size),
#               dP=ifelse(time<t1,0,rw_size),
#               dT=ifelse(time<t1,0,rw_size),
#               dB=ifelse(time<t1,0,rw_size),
#               qP=ifelse(time<t1,rw_size,rw_size),
#               qI=ifelse(time<t1,rw_size,rw_size),
#               qH=ifelse(time<t1,rw_size,rw_size),
#               qC=ifelse(time<t1,rw_size,rw_size),
#               qV=ifelse(time<t1,rw_size,rw_size),
#               sigma=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               kappa=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaI=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaH=ifelse(time<93,0,rw_size),
#               gammaC=ifelse(time<93,0,rw_size),
#               gammaV=ifelse(time<93,0,rw_size),
#               rho=ifelse(time<t1,rw_size,rw_size),
#               TF=ifelse(time<35,0,rw_size))
# } else if (run==4) {
#   rw <- rw.sd(betaI=ifelse(time<t1,rw_size,0),
#               theta=ifelse(time<t1,rw_sd,0),
#               iota=ifelse(time<t1,0,0),
#               beta_sd=ifelse(time<t1,rw_size,0),
#               dI=ifelse(time<t1,0,rw_size),
#               dP=ifelse(time<t1,0,rw_size),
#               dT=ifelse(time<t1,0,rw_size),
#               dB=ifelse(time<t1,0,rw_size),
#               qP=ifelse(time<t1,rw_size,rw_size),
#               qI=ifelse(time<t1,rw_size,rw_size),
#               qH=ifelse(time<t1,rw_size,rw_size),
#               qC=ifelse(time<t1,rw_size,rw_size),
#               qV=ifelse(time<t1,rw_size,rw_size),
#               sigma=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               kappa=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaI=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaH=ifelse(time<93,0,rw_size),
#               gammaC=ifelse(time<93,0,rw_size),
#               gammaV=ifelse(time<93,0,rw_size),
#               rho=ifelse(time<t1,rw_size,rw_size),
#               TF=ifelse(time<35,0,rw_size))
# } else if (run==5) {
#   rw <- rw.sd(betaI=ifelse(time<t1,rw_size,0),
#               theta=ifelse(time<t1,rw_size,0),
#               iota=ifelse(time<t1,rw_size,0),
#               beta_sd=ifelse(time<t1,0,0),
#               dI=ifelse(time<t1,0,rw_size),
#               dP=ifelse(time<t1,0,rw_size),
#               dT=ifelse(time<t1,0,rw_size),
#               dB=ifelse(time<t1,0,rw_size),
#               qP=ifelse(time<t1,rw_size,rw_size),
#               qI=ifelse(time<t1,rw_size,rw_size),
#               qH=ifelse(time<t1,rw_size,rw_size),
#               qC=ifelse(time<t1,rw_size,rw_size),
#               qV=ifelse(time<t1,rw_size,rw_size),
#               sigma=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               kappa=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaI=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaH=ifelse(time<93,0,rw_size),
#               gammaC=ifelse(time<93,0,rw_size),
#               gammaV=ifelse(time<93,0,rw_size),
#               rho=ifelse(time<t1,rw_size,rw_size),
#               TF=ifelse(time<35,0,rw_size))
# } else if (run==6) {
#   rw <- rw.sd(betaI=ifelse(time<t1,rw_size,0),
#               theta=ifelse(time<t1,rw_size,0),
#               iota=ifelse(time<t1,rw_size,0),
#               beta_sd=ifelse(time<t1,rw_size,0),
#               dI=ifelse(time<t1,0,0),
#               dP=ifelse(time<t1,0,rw_size),
#               dT=ifelse(time<t1,0,rw_size),
#               dB=ifelse(time<t1,0,rw_size),
#               qP=ifelse(time<t1,rw_size,rw_size),
#               qI=ifelse(time<t1,rw_size,rw_size),
#               qH=ifelse(time<t1,rw_size,rw_size),
#               qC=ifelse(time<t1,rw_size,rw_size),
#               qV=ifelse(time<t1,rw_size,rw_size),
#               sigma=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               kappa=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaI=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaH=ifelse(time<93,0,rw_size),
#               gammaC=ifelse(time<93,0,rw_size),
#               gammaV=ifelse(time<93,0,rw_size),
#               rho=ifelse(time<t1,rw_size,rw_size),
#               TF=ifelse(time<35,0,rw_size))
# } else if (run==22) {
#   rw <- rw.sd(betaI=ifelse(time<t1,rw_size,0),
#               theta=ifelse(time<t1,rw_size,0),
#               iota=ifelse(time<t1,rw_size,0),
#               beta_sd=ifelse(time<t1,rw_size,0),
#               dI=ifelse(time<t1,0,rw_size),
#               dP=ifelse(time<t1,0,rw_size),
#               dT=ifelse(time<t1,0,rw_size),
#               dB=ifelse(time<t1,0,rw_size),
#               qP=ifelse(time<t1,rw_size,rw_size),
#               qI=ifelse(time<t1,rw_size,rw_size),
#               qH=ifelse(time<t1,rw_size,rw_size),
#               qC=ifelse(time<t1,rw_size,rw_size),
#               qV=ifelse(time<t1,rw_size,rw_size),
#               sigma=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               kappa=ifelse(time<t1,0.5*rw_size,0.5*rw_size),
#               gammaI=ifelse(time<t1,0,0),
#               gammaH=ifelse(time<93,0,rw_size),
#               gammaC=ifelse(time<93,0,rw_size),
#               gammaV=ifelse(time<93,0,rw_size),
#               rho=ifelse(time<t1,rw_size,rw_size),
#               TF=ifelse(time<35,0,rw_size))
# }


