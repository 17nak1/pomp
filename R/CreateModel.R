# Function creating the pomp model

create_pomp_model <- function(data, 
                              covars, 
                              t0,  
                              modeltype=c(nstageE=3L, nstageI=3L, nstageP=3L, nstageH=3L, nstageC=3L, nstageV=3L), 
                              dt=0.1) {
  
  
  nstageE <- modeltype["nstageE"]
  nstageI <- modeltype["nstageI"]
  nstageP <- modeltype["nstageP"]
  nstageH <- modeltype["nstageH"]
  nstageC <- modeltype["nstageC"]
  nstageV <- modeltype["nstageV"]
  # covars <- subset(covars, select = c("time", "pop"))


  globals <- paste0("int nstageE = ", nstageE, 
                    "; int nstageP = ", nstageP, 
                    "; int nstageI = ", nstageI, 
                    "; int nstageH = ", nstageH, 
                    "; int nstageC = ", nstageC,
                    "; int nstageV = ", nstageV,
                    "; double pop = ", 10e6, 
                    "; double T0 = ", 75, 
                    "; double T1 = ", 139, ";")

  # Create the pomp model
  COVID <- pomp(
    data = subset(data, select = c("time", "reports", "deaths", "hospital", "ICU","ventilator")),
    times = "time",
    t0 = t0,
    covar = covars,
    tcovar = "time",
    zeronames = zeronames,
    paramnames = c(params_mod, params_ic),
    dmeasure = dObs,
    rmeasure = rObs,
    rprocess = euler.sim(step.fun=rSim,delta.t = dt),
    skeleton = map(skel, delta.t = dt),
    globals = "int nstageE=3; int nstageP=3; int nstageI=3; int nstageH=3; int nstageC=3; int nstageV=3; double pop=10e6; double T0=75; double T1=139;",
    params_log = params_log,
    params_logit = params_logit,
    fromEstimationScale  = function(params, params_log, params_logit,  ...){
      params[params_log] <- exp(params[params_log])
      params[params_logit] <- plogis(params[params_logit])
      params
    },
    toEstimationScale = function(params, params_log, params_logit, ...){
      params[params_log] <- log(params[params_log])
      params[params_logit] <- qlogis(params[params_logit])
      params
    },
    initializer = rInit,
    statenames=c(statenames(nstageE, nstageP, nstageI, nstageH, nstageC, nstageV),
                 zeronames)
  )

  output <- list(model=COVID)
  print("Model created.")
  return(output)
}