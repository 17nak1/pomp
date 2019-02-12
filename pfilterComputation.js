
// examines weights for filtering failure.
// computes conditional log likelihood and effective sample size.
// computes (if desired) prediction mean, prediction variance, filtering mean.
// it is assumed that ncol(x) == ncol(params).
// weights are used in filtering mean computation.
// if length(weights) == 1, an unweighted average is computed.
// tracks ancestry of particles if desired.
// returns all of the above in a named list.
mathLib = require('./mathLib.js')
pfilterComputations = function (states, params, Np, doPredictionMean = 0, doPredictionVariance = 0, filtmean = 0,  trackancestry = 0,  doparRS = 0, normalWeights,  toler) {
  var allFail = 0, ess, nvars, returnStates, loglik, ess
  navars = states[0].length

  // check the weights and compute sum and sum of squares
   var  w = 0, ws = 0, nlost = 0
   for (let i = 0; i < Np; i++) {
     if (normalWeights[i] > toler) {
       w += normalWeights[i];
       ws += normalWeights[i]*normalWeights[i];
     } else {      // this particle is lost
       normalWeights[i] = 0;
       nlost++;
     }
   }
   if (nlost >= Np) allFail = 1 // all particles are lost
   if (allFail) {
     loglik = Math.log(toler) // minimum log-likelihood
     ess = 0  // zero effective sample size
   } else {
     loglik = Math.log(w / Np) // mean of weights is likelihood
     ess = w * w / ws  // effective sample size
   }

// compute prediction mean
    if (doPredictionMean || doPredictionVariance) {
      for (let j = 0; j< nvars + 1;j++) {
        var sum = 0, nlost = 0
        for (let nrow =0; nrow < Np; nrow++){
          if (states[nrow][j]) {
            sum += states[nrow][j]
          } else {
            nlost++
          }
        }
        predictionMean[timeCounter][j] = sum / Np
      }
    }
// compute prediction variance
    if (doPredictionVariance) {
      var sumsq = 0, vsq
      for (let j = 0; j< nvars + 1;j++) {
        for (let nrow =0; nrow < Np; nrow++){
          if (states[nrow][j]) {
            vsq = states[nrow][j] - sum
            sumsq += vsq * vsq
          }
        }
        predictionVariance[timeCounter][j] = sumsq / (Np - 1) 
      }
    }
//  compute filter mean
    if (doFilterMean) {
      if (allFail) {   // unweighted average
        ws = 0
        for (let j = 0; j< nvars + 1;j++) {
          for (let nrow =0; nrow < Np; nrow++){
            if (states[nrow][j]) {
              ws += states[nrow][j]
            }
          }
        } 
        filterMean[j] = ws / Np
      } else {      // weighted average
        ws = 0
        for (let j = 0; j< nvars + 1;j++) {
          for (let nrow =0; nrow < Np; nrow++){
            if (states[nrow][j]) {
              ws += states[nrow][j] * weights[k]
            }
          }
        }
        filterMean[j] = ws / w
      }
    }
    // if (!allFail) {
    //   var sample = new Array(np)
    //   // resample
    //   mathLib.nosortResamp(Np, normalWeights, Np, sample, 0)
    //   for (k = 0; k < np; k++) { // copy the particles
    //     xx = ss+nvars*sample[k]
    //     for (j = 0; j < nvars; j++, st++, xx++)
    //       st[j] = xx[j]
    //     if (do_pr) {
    //       for (j = 0, xp = ps+npars*sample[k]; j < npars; j++, pt++, xp++)
    //         pt = xp
    //     }
    //     if (do_ta) xanc[k] = sample[k]+1;
    //   }
  
    // } else { // don't resample: just drop 3rd dimension in x prior to return
  
    //   PROTECT(newdim = NEW_INTEGER(2)); nprotect++;
    //   dim = INTEGER(newdim);
    //   dim[0] = nvars; dim[1] = nreps;
    //   SET_DIM(x,newdim);
    //   setrownames(x,Xnames,2);
    //   fixdimnames(x,dimnm,2);
  
    //   if (do_ta)
    //     for (k = 0; k < np; k++) xanc[k] = k+1;
    // }

   if (allFail) {
    returnStates = states
  } else {
    returnStates = statesnewstates
  } 
  var returnValues = {fail:0, loglik: loglik, ess: ess, states: returnStates, params:0, pm: predictionMean, pv: predictionVariance,fm: filterMean, ancestry:0 }  
  return returnValues
}