# Set random walk values. All fixed parameter rw.sd will eventually be removed

# if (run==1) {
  rw <- rw.sd(
              betaI=ifelse(time<75,0,0),
              theta=ifelse(time<75,0,0),
              iota=ifelse(time<75,0,0),
              beta_sd=ifelse(time<75,0,0),
              dI0=ifelse(time<75,0,0),
              dP0=ifelse(time<75,0,0),
              dT0=ifelse(time<75,0,0),
              dB0=ifelse(time<75,0,rw_size),
              dI1=ifelse(time<75,0,0),
              dP1=ifelse(time<75,0,0),
              dT1=ifelse(time<75,0,0),
              dB1=ifelse(time<75,0,0),
              qP=ifelse(time<75,0,0),

              qH=ifelse(time<75,0,0),
              qC=ifelse(time<75,0,0),

              mI=ifelse(time<75,0,0),
              mC=ifelse(time<75,0,0),
              mV=ifelse(time<75,0,0),
              sigma=ifelse(time<75,0.5*0,0.5*0),
              kappa=ifelse(time<75,0.5*0,0.5*0),
              gammaI=ifelse(time<75,0.5*0,0.5*0),
              gammaH=ifelse(time<93,0,0),
              gammaC=ifelse(time<93,0,0),
              gammaV=ifelse(time<93,0,0),
              rho=ifelse(time<75,0,0),
              TF=ifelse(time<35,0,0)
              )
