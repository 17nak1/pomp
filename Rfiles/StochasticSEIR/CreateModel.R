



pomp(
  data = London_BiData,
  times="time",
  t0=1940,#with(get(paste0(name,"_BiData")),2*time[1]-time[2]),
  rprocess = euler.sim(rproc,delta.t=1/365.25),
  rmeasure=rmeas,
  covar=London_covar,
  tcovar="time",
  dmeasure=dmeas,
  zeronames=zeronames,
  initializer=initz,
  toEstimationScale=toEst,
  fromEstimationScale=fromEst,
  statenames=statenames,
  paramnames=c(params_mod,params_ic)
) -> m1


rm(rmeas,rproc,dmeas,fromEst,toEst,initz)
