/*
 * Particle filter
 * A plain vanilla sequential Monte Carlo (particle filter) algorithm.
 * Resampling is performed at each observation.
 * @param Np the number of particles to use.
 * @param tol positive numeric scalar.
 * @param max.fail integer; the maximum number of filtering failures allowed.
 * @param pred.mean logical; if \code{TRUE}, the prediction means are calculated for the state variables and parameters.
 * @param pred.var logical; if \code{TRUE}, the prediction variances are calculated for the state variables and parameters.
 * @param filter.mean logical; if \code{TRUE}, the filtering means are calculated for the state variables and parameters.
 * @param filter.traj logical; if \code{TRUE}, a filtered trajectory is returned for the state variables and parameters.
 * @param save.states logical.
 *  If \code{save.states=TRUE}, the state-vector for each particle at each time is saved.
 * @references
    pfilter.R by Aaron A. King. https://github.com/kingaa/pomp/blob/abc552b335319bd36b21d3313305c89e97f48f68/R/pfilter.R
    M. S. Arulampalam, S. Maskell, N. Gordon, & T. Clapp.
    A Tutorial on Particle Filters for Online Nonlinear, Non-Gaussian Bayesian Tracking.
    IEEE Trans. Sig. Proc. 50:174--188, 2002.
 */
//conditional log liklihood(time) =log(sum(w_i, i in 0:Np) / Np)
var START = new Date()
fs = require('fs')
let fmin = require ('fmin')
let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')
var linearInterpolator = require('linear-interpolator/node_main')


//////////////////////////////////////////////////////////////////////////////////////* main function//////////////////////////////////////////////////////////////////////
function pfilterCalculation (input) {//filter.traj , save.params
  // {params:inputArr, Np:100,times:times, dt:1 / 365.25,runPredMean:1,  dataCases:dataCases, interpolPop:interpolPopulation, interpolBirth:interpolBirth}
  let START =new Date()
  let defaults = {params:-1, Np:-1, tol:1e-17, maxFail:Infinity, runPredMean:0, runPredVar:0, runFilterMean:0, runSaveStates:0, times:-1, dt:-1,
                   dataCases:0}
  for(let prop in defaults) {
    if(typeof input[prop] == 'undefined') {
      input[prop] = defaults[prop]
    }
  }
  if (input.params === -1 || input.Np === -1 || input.times === -1 || input.dt === -1) {
    throw 'Some required arguments are missed'
  }
  
  var params = input.params
  var Np = input.Np
  var toler = input.tol
  var maxFail = input.maxFail
  var dt = input.dt
  var dataCases = input.dataCases
  var interpolPop = input.interpolPop
  var interpolBirth = input.interpolBirth
  
  let [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, I_0, R_0] = params
let [t0, tdata] = [1940, 1944]
let nvars = 5
let deltaT = 14 / 365.25
let doPredictionVariance = 0, doPredictionMean = 1, doFilterMean = 0 , allFail = 0

let timeLen = dataCases.length 
let nlost = 0

let rate = [], trans = []
var particles = new Array(Np).fill(null).map(() => Array(5)),
 state =[]
let sampleNum = Array.from(Array(Np).keys())
let condLoglik = []
let stateSaved = []
let temp 
let timeCountData = 0, ws ,w , vsq, sumsq, ess, loglik = 0, lik 

let predictionMean, predictionVariance, filterMean
let states = Array(Np).fill(null).map(() => Array(nvars))
let weights, normalWeights, S, E, I, R, del_t, ST, simulateValue
let modelCases, likvalue
var st, births, pop,birthrate
if (doPredictionMean) {
  predictionMean = Array(timeLen).fill(null).map(() => Array(nvars))
}
if (doPredictionVariance) {
  predictionVariance = Array(timeLen).fill(null).map(() => Array(nvars))
}
if (doFilterMean) {
  filterMean = Array(timeLen).fill(null).map(() => Array(nvars))
}
state = snippet.initz(interpolPop(t0), S_0, E_0, I_0, R_0)
// First Np sets
var particles = new Array(Np).fill(null).map(() => [].concat(state));
temp = new Array(Np).fill(null).map(() => [].concat(state))
// Time loop
for (k = t0; k <= Number(dataCases[timeLen - 2][0]) + deltaT / 3 ; k += deltaT){//Number(dataCases[timeLen - 2][0]) + deltaT / 3
  if ( k > tdata - deltaT && k <= tdata) {
    k = tdata
  }
  weights = []; normalWeights = []
  
  if ( k > t0) {
    for (np = 0; np < Np; np++) { // copy the particles
      temp[np] = [].concat(particles[sampleNum[np]])

      temp[np][nvars - 1] = 0
    }
  } 
  
  // if (k <= tdata || k > Number(dataCases[dataCases.length - 1])) {
      steps = mathLib.numMapSteps(k, k + deltaT, dt)
  // } else {
  //     steps = mathLib.numEulerSteps(k, Number(dataCases[timeCountData + 1][0]), dt)
  // }
    del_t = (1 / steps )* deltaT

  for (let stp = 0; stp < steps; stp++) { // steps in each time interval
    st = k + stp * del_t
    pop = interpolPop(st)
    birthrate = interpolBirth(st)
    births = mathLib.rpois(birthrate * (1- snippet.rprocessVaccine(st)) * del_t )
      for (np = 0; np < Np; np++){ //calc for each particle
        trans = []
        S = temp[np][0]; E = temp[np][1]; I = temp[np][2]; R = temp[np][3]; H = temp[np][4]
        simulateValue = snippet.rprocess(params, st, del_t, [S,E,I,R,H], pop, births)
        temp[np][0] = simulateValue[0]; temp[np][1] = simulateValue[1]; temp[np][2] = simulateValue[2]; temp[np][3] = simulateValue[3]; temp[np][4] = simulateValue[4]
      }
    }
    for (np = 0; np < Np; np++){ 
      particles[np][0] = temp[np][0]
      particles[np][1] = temp[np][1]
      particles[np][2] = temp[np][2]
      particles[np][3] = temp[np][3]
      particles[np][4] = temp[np][4]
      H = temp[np][4]
  //***********weight*************
      if (k >= Number(dataCases[0][0])){
        if (stateSaved) {
          stateSaved.push(particles[np]) //[S,E,I,R,H])
        }
        modelCases = Number(dataCases[timeCountData][1])
        likvalue = snippet.dmeasure(rho, psi, H, modelCases, giveLog = 0)
        weights.push(likvalue)           
      }
    }
   
  
  //normalize
  if (k >= Number(dataCases[0][0])){
    let sumOfWeights = 0
    for (let i = 0; i < Np; i++) {
      sumOfWeights += weights[i]
    }
    for (let i = 0; i < Np; i++) {
      normalWeights[i] = weights[i] / sumOfWeights
    }
    // check the weights and compute sum and sum of squares
    w = 0, ws = 0, nlost = 0
    for (let i = 0; i < Np; i++) {
      if (weights[i] > toler) {
        w += weights[i]
        ws += weights[i] ** 2
      } else { // this particle is lost
        weights[i] = 0;
        nlost++
      }
    }
    if (nlost > maxFail) {
      throw 'execution terminated. The number of filtering failures exceeds the maximum number of filtering failures allowed. '
    }
    if (nlost >= Np) { 
      allFail = 1 // all particles are lost
    } else {
      allFail = 0
    }
    
    if (allFail) {
      lik = Math.log(toler) // minimum log-likelihood
      ess = 0  // zero effective sample size
    } else {
      ess = w * w / ws  // effective sample size
      lik = Math.log(w / Np)// mean of weights is likelihood
    }
    condLoglik[timeCountData] = [timeCountData + 1, lik]
    // the total conditional logliklihood in the time process is loglik
    loglik += lik
    mathLib.nosortResamp(Np, normalWeights, Np, sampleNum, 0)
    
    // Compute outputs
    for (let j = 0; j< nvars; j++) {
      // compute prediction mean
      if (doPredictionMean || doPredictionVariance) {
        let sum = 0, nlost = 0
        for (let nrow =0; nrow < Np; nrow++){
          if (particles[nrow][j]) {
            sum += particles[nrow][j]
          } else {
            nlost++
          }
        }
        sum /= Np
        predictionMean[timeCountData][j] = sum
      }  
      // compute prediction variance
      if (doPredictionVariance) {
        sumsq = 0
        for (let nrow = 0; nrow < Np; nrow++){
          if (particles[nrow][j]) {
            vsq = particles[nrow][j] - sum
            sumsq += Math.pow(vsq, 2)
          }
        }
        predictionVariance[timeCountData][j] = sumsq / (Np - 1) 
      }
      //  compute filter mean
      if (doFilterMean) {
        if (allFail) {   // unweighted average
          ws = 0
          for (let nrow =0; nrow < Np; nrow++){
            if (particles[nrow][j]) {
              ws += particles[nrow][j]
            }
          } 
          filterMean[timeCountData][j] = ws / Np//;console.log(ws / Np)
        } else {      // weighted average
          ws = 0
          for (let nrow =0; nrow < Np; nrow++){
            if (particles[nrow][j]) {
              ws += particles[nrow][j] * weights[nrow]
            }
          }
          filterMean[timeCountData][j] = ws / w
        }
      }
    }
    }
    timeCountData++
}//endTime

  
  console.log('loglike=',loglik)
  console.log('runing time=', new Date() - START)
  activateDownload ()
   if (input.runPredMean) {
    return predictionMean
  }

  if (input.runPredVar) {
    return predictionVariance
  }

  if (input.runFilterMean) {
    return filterMean
  }
  
  if (input.runSaveStates){
    return stateSaved
  }

}

module.exports = {
  pfilterCalculation
}


    
      


