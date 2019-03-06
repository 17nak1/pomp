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

fs = require('fs')
let fmin = require ('fmin')
let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')
const mathjs = require('mathjs') 
var linearInterpolator = require('linear-interpolator/node_main')

//////////////////////////////////////////data///////////////////////////////////////
var LondonBidata = [], LondonCovar = []
var rate = new Array(6) 
// var params = [3.165652e+01 , 3.887624e-01 , 7.305000e+01 , 1.698730e-02  ,4.566000e+01,  4.813669e-01  ,1.963092e-01 , 2.831066e-03 ,3.476483e-04 ,  2.109135e-08,9.968213e-01]
var params = [3.132490e+01 , 3.883620e-01 , 7.305000e+01 , 6.469830e-04 , 4.566000e+01 , 4.598709e-01 , 1.462546e-01 , 3.399189e-02 ,2.336327e-04 ,4.221789e-07 ,9.657741e-01 ]
var times =[1940, 1944]
var Np = 100
var nvars 
var toler = 1e-17
var nlost = 0

//* 1st data set
var London_covar = fs.readFileSync('./London_covar.csv').toString()
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonCovar.push(lines[i].split(','))
}
var dataCovar = [LondonCovar][0]
//* 2nd data set
var London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonBidata.push(lines[i].split(','))
}
var dataCases = [LondonBidata][0]

var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
var interpolPop = linearInterpolator(d1)
var interpolBirth = linearInterpolator(d2)

//////////////////////////////////////////////////////////////////////////////////////* main function//////////////////////////////////////////////////////////////////////

var [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, R_0, I_0] = params
var paramsIC = [S_0, E_0, R_0, I_0, H = 0]
var [t0, tdata] = times
var nvars = paramsIC.length
var deltaT = 14 / 365.25
var dt = 1 / 365.25
var pop , birthrate
var timeLen = dataCases.length 
var va = 0, seas

var particles = [], state =[]
var sampleNum = new Array(Np)
var doPredictionVariance = 1, doPredictionMean = 1, doFilterMean = 1 , allFail = 0
var predictionMean = Array(timeLen).fill(null).map(() => Array(nvars))
var predictionVariance = Array(timeLen).fill(null).map(() => Array(nvars))
var condLoglik = Array(timeLen).fill(null).map(() => Array(1))
var filterMean = Array(timeLen).fill(null).map(() => Array(nvars))
var timeCountData = 0, ws = 0, vsq, sumsq, ess, Loglik = 0, lik //condLoglik = []
var maxFail = 300,  stateSaved =[]
var states = Array(Np).fill(null).map(() => Array(nvars))

// if ( k === t0) {
  // var Nlog = mathLib.toLogBarycentric([state[0], state[1], state[2], state[3]],4)
  // var N = mathLib.fromLogBarycentric(Nlog, 4)
  state = snippet.initz(interpolPop(t0), S_0, E_0, I_0, R_0)
for ( i=0; i < Np; i++){
  particles[i] = [].concat(state)
}


//***********************************************TIME LOOP************************************
for (k = t0  ; k < Number(dataCases[dataCases.length - 2][0]) + deltaT / 3; k += deltaT){//Number(dataCases[dataCases.length - 2][0]) + deltaT / 3
  if ( k > tdata && k <= tdata + deltaT) {
    k = tdata + deltaT
  }
  var loglik = 0
  var lik = new Array(Np)
  var weights = [], normalWeights = []

  //****************************************PARTICLE LOOP**************************************////////////////////////////////////////////////////////////////////////////
  for (np = 0; np < Np; np++){ //calc for each particle
    var trans = new Array(6).fill(0)
    var S = particles[np][0], E = particles[np][1], I = particles[np][2], R = particles[np][3], H = particles[np][4]
      
    // transitions between classes
    if (k <= tdata || k > 1965 - deltaT ) {
      steps = mathLib.numMapSteps(k, k + deltaT, dt)
    } else {
      steps = mathLib.numEulerSteps(k, Number(dataCases[timeCountData + 1][0]), dt)
    }
    var del_t = (1 / steps )* deltaT 
    for (let stp = 0; stp < steps; stp++) { // steps in each time interval
      var st = k + stp * del_t
      st = st.toFixed(10)
      var simulateValue = snippet.rprocess(params, st, del_t, [S,E,I,R,H], interpolPop(st), interpolBirth(st))
      S = simulateValue[0]; E = simulateValue[1], I = simulateValue[2], R = simulateValue[3], H = simulateValue[4]
      
    }
    particles[np][0] = S
    particles[np][1] = E
    particles[np][2] = I
    particles[np][3] = R
    particles[np][4] = H
   
    states[np][0] = S || 0
    states[np][1] = E || 0
    states[np][2] = I || 0
    states[np][3] = R || 0
    states[np][4] = H || 0
     
    //***********RESAMPLE*************
    stateSaved.push([S,E,I,R,H])
    if (k >= Number(dataCases[0][0])){
      var modelCases = Number(dataCases[timeCountData][1])
      var likvalue = snippet.dmeasure(rho, psi, H, modelCases, giveLog = 0)
      weights.push(likvalue)
      particles[np][4] = 0
    }
  }////////////////////////////////////////////////////////////////end particle loop///////////////////////////////////////////////////////////////////////////////////////
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
    var  w = 0, ws = 0, nlost = 0
    for (let i = 0; i < Np; i++) {
      if (weights[i] > toler) {
        w += weights[i]
        ws += weights[i] ** 2
      } else { // this particle is lost
        weights[i] = 0;
        nlost++
      }
    }
    if (nlost >= Np) { 
      allFail = 1 // all particles are lost
    } else {
      allFail = 0
    }
    // console.log(k,allFail)
    if (allFail) {
      lik = Math.log(toler) // minimum log-likelihood
      ess = 0  // zero effective sample size
    } else {
      ess = w * w / ws  // effective sample size
      lik = Math.log(w / Np)// mean of weights is likelihood
    }
    condLoglik[timeCountData] = [timeCountData + 1, lik]//;console.log(condLoglik)
    Loglik += lik
    mathLib.nosortResamp(Np, normalWeights, Np, sampleNum, 0)//;console.log(k, sampleNum)
    for (np = 0; np < Np; np++) { // copy the particles
      particles[np] = [].concat(particles[sampleNum[np]])
      particles[np][nvars - 1] = 0
    }
    // Compute outputs
    for (let j = 0; j< nvars; j++) {
      // compute prediction mean
      if (doPredictionMean || doPredictionVariance) {
        var sum = 0, nlost = 0
        for (let nrow =0; nrow < Np; nrow++){
          if (states[nrow][j]) {
            sum += states[nrow][j]
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
          if (states[nrow][j]) {
            vsq = states[nrow][j] - sum
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
            if (states[nrow][j]) {
              ws += states[nrow][j]
            }
          } 
          filterMean[timeCountData][j] = ws / Np//;console.log(ws / Np)
        } else {      // weighted average
          ws = 0
          for (let nrow =0; nrow < Np; nrow++){
            if (states[nrow][j]) {
              ws += states[nrow][j] * weights[nrow]
            }
          }
          filterMean[timeCountData][j] = ws / w
        }
      }
    }
    timeCountData++ 
  }
}//endTime
console.log(Loglik)
const createCsvWriter = require('csv-writer').createArrayCsvWriter;
const csvWriter = createCsvWriter({
  header: ['S', 'E', 'I', 'R', 'H'],
  path: './predmean.csv'
})
 
csvWriter.writeRecords(predictionMean)
  .then(() => {
  console.log('...predictionMean')
})

const createCsvWriter11 = require('csv-writer').createArrayCsvWriter;
const csvWriter11 = createCsvWriter11({
  header: ['S', 'E', 'I', 'R', 'H'],
  path: './predvar.csv'
})
 
csvWriter11.writeRecords(predictionVariance)
  .then(() => {
  console.log('...predictionVar')
})

var createCsvWriter1 = require('csv-writer').createArrayCsvWriter;
var csvWriter1 = createCsvWriter1({
  header: ['S', 'E', 'I', 'R', 'H'],
  path: './filterMean.csv'
})
csvWriter1.writeRecords(filterMean)
  .then(() => {
  console.log('...filterMean')
})

var createCsvWriter2 = require('csv-writer').createArrayCsvWriter;
var csvWriter2 = createCsvWriter2({
  header: [],
  path: './condLogliklihood.csv'
}) 
csvWriter2.writeRecords(condLoglik)
  .then(() => {
  console.log('...condLoglik')
})

var createCsvWriter3 = require('csv-writer').createArrayCsvWriter;
var csvWriter3 = createCsvWriter3({
  header:['S', 'E', 'I', 'R', 'H'],
  path: './stateSaved.csv'
})
csvWriter3.writeRecords(stateSaved)
  .then(() => {
  console.log('...stateSaved')
})
  
const logMeanExp = function (x) {
  var mx = Math.max(...x)
  var s = x.map((x,i) => Math.exp(x - mx))
  var q = s.reduce((a, b) => a + b, 0)
  return mx + Math.log(q / x.length)
}

const numMapSteps = function (t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}



    
      


