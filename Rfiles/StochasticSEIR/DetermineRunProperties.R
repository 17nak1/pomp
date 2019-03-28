determine_run_properties <- function  (run=2,modeltype=1) {
  if (run %in% c(0,1)) {
    param <- NULL
    lscale <- NULL
    param_lims <- NULL
    flag_bound <- NULL
  } else if (run %in% c(2)) {
    param <- "R0"
    lscale <- 0
    param_lims <- c(0,80)
    flag_bound <- 0
  } else if (run %in% c(3)) {
    param <- "amplitude"
    lscale <- 0
    param_lims <- c(0,1)
    flag_bound <- 1
  } else if (run %in% c(4)) {
    param <- "mu"
    lscale <- 0
    param_lims <- c(0,1)
    flag_bound <- 1
  } else if (run %in% c(5)) {
    param <- "rho"
    lscale <- 0
    param_lims <- c(0.3,1)
    flag_bound <- 1
  } else if (run %in% c(6)) {
    param <- "psi"
    lscale <- 0
    param_lims <- c(0,1)
    flag_bound <- 1
  }   
 
  out <- list(param=param,lscale=lscale,param_lims=param_lims,flag_bound=flag_bound)
}