
// examines weights for filtering failure.
// computes conditional log likelihood and effective sample size.
// computes (if desired) prediction mean, prediction variance, filtering mean.
// it is assumed that ncol(x) == ncol(params).
// weights are used in filtering mean computation.
// if length(weights) == 1, an unweighted average is computed.
// tracks ancestry of particles if desired.
// returns all of the above in a named list.
mathLib = require('./mathLib.js')
snippet = require('./modelSnippet')

// weights from dmeasure with logScale = 0
pfilterComputations = function (states, params, Np, doPredictionMean = 0, doPredictionVariance = 0, doFiltermean = 0,  doTrackancestry = 0,  doparamResample = 0,
 NormalWeights, toler, timeLen) {

  var allFail = 0, ess, nvars, nreps, returnStates, loglik, ess
  var predictionMean = 0, predictionVariance = 0, filterMean = 0, trackancestry = 0, returnStates, returnParams
  nvars = states[0].length
  nreps = states.length
  var newstates = Array(nreps).fill(Array(nvars))
  if(doPredictionMean) {
    predictionMean = new Array(timeLen).fill(Array(nvars))
  } else {
    predictionMean = new Array(0)
  }
 
  if(doPredictionVariance) {
    predictionVariance = new Array(timeLen).fill(Array(nvars))
  } else {
    predictionVariance = new Array(0)
  }

  if(doFilterMean) {
    filterMean = new Array(timeLen).fill(Array(nvars))
  } else {
    filterMean = new Array(0)
  }
  if (doparamResample) {
    if (params.length != nreps) {
      throw "ncol(states) should be equal to ncol(params)"
    } else {
      var newParams = Array(nreps).fill(Array(params[0].length))
  }
  // check the weights and compute sum and sum of squares
  var  w = 0, ws = 0, nlost = 0
  for (let i = 0; i < nreps; i++) {
    if (weights[i] > toler) {
      w += weights[i]
      ws += weights[i] ** 2
    } else { // this particle is lost
      weights[i] = 0;
      nlost++
    }
  }

  if (nlost >= nreps) { 
    allFail = 1 // all particles are lost
  } else {
    allFail = 0
  }
  if (allFail) {
    lik = Math.log(toler) // minimum log-likelihood
    ess = 0  // zero effective sample size
  } else {
    lik = Math.log(w / Np)// mean of weights is likelihood
    ess = w * w / ws  // effective sample size
  }
  var fail = allFail
  
  // compute prediction mean
  if (doPredictionMean || doPredictionVariance) {
    for (let j = 0; j< nvars;j++) {
      var sum = 0, nlost = 0
      for (let nrow =0; nrow < nreps; nrow++){
        sum += states[nrow][j]
      }
      sum /= nreps
      predictionMean[j] = sum 
    }
  }
  
  // compute prediction variance
  if (doPredictionVariance) {
    var sumsq = 0, vsq
    for (let j = 0; j< nvars;j++) {
      for (let nrow = 0; nrow < nreps; nrow++){
          vsq = states[nrow][j] - sum
          sumsq += vsq * vsq
      }
      predictionVariance[j] = sumsq / (nreps - 1) 
    }
  }

  //  compute filter mean
  if (doFilterMean) {
    if (allFail) {   // unweighted average
      ws = 0
      for (let j = 0; j< nvars;j++) {
        for (let nrow = 0; nrow < nreps; nrow++){
          ws += states[nrow][j]
        }
      } 
      filterMean[j] = ws / nreps
    } else {      // weighted average
      ws = 0
      for (let j = 0; j< nvars;j++) {
        for (let nrow = 0; nrow < nreps; nrow++){
          ws += states[nrow][j] * weights[nrow]
        }
      }
      filterMean[j] = ws / w
    }
  }


  // resample the particles unless we have filtering failure
  if (!allFail) {
    var sample = new Array(np)
    mathLib.nosortResamp(nreps, weights, Np, sample, 0) // resample
    for (k = 0; k < np; k++) { 
      if (doparamResample) {
        // xp = ps+npars*sample[k]
        // for (j = 0; j < npars; j++, pt++, xp++)
        //   pt = xp
      }
      if (doTrackancestry) {
        ancestry[k] = sample[k] + 1
      }
    }  
  } else { // don't resample: just drop 3rd dimension in x prior to return  
    if (doTrackancestry) {
      for (k = 0; k < np; k++) {
        ancestry[k] = sample[k] + 1
      }
    }
  }

  if (allFail) {
    returnStates = states
  } else {
    returnStates = newstates
  } 

  if (allFail || !doparamResample) {
    returnParams = params
  } else {
    returnParams = newParams   
  }


  var returnValues = {fail:0, loglik: loglik, ess: ess, states: returnStates, params:0, pm: returnPredictionMean, pv: returnPredictionVariance,fm: filterMean, ancestry:0 }  
  // return returnValues
  console.log(returnValues)
}

