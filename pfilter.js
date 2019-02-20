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

fs = require('fs')
let fmin = require ('fmin')
let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')
const mathjs = require('mathjs') 
var linearInterpolator = require('linear-interpolator/node_main')

//* Set the seed for rnorm-In R:RNGkind("L'Ecuyer-CMRG", normal.kind="Box-Muller");set.seed(1234) 
const libR = require('lib-r-math.js')
const {
  Poisson,
  rng: { MersenneTwister },
  rng: { normal: { Inversion } }
} = libR
const mt = new MersenneTwister(0)// 
const { rpois } = Poisson(new Inversion(mt))
mt.init(1234)

var {
  Normal,
  rng: {
    LecuyerCMRG,
    normal: { BoxMuller }
  }
} = libR
const ad = new LecuyerCMRG(0)
const { rnorm } = Normal(new BoxMuller(ad))

//////////////////////////////////////////data///////////////////////////////////////
var LondonBidata = [], LondonCovar = []
var rate = new Array(6) 
// var params = [2.879243e+01, 3.669313e-01, 7.305000e+01, 7.020601e-03, 4.566000e+01, 4.920423e-01, 1.709873e-01, 2.560858e-02, 1.496166e-05, 5.523655e-08, 9.743764e-01]
var params = [3.013710e+01,  3.816660e-01,  7.305000e+01,  6.151897e-03,  4.566000e+01,  4.215311e-01,  1.871425e-01,  7.947781e-02,  2.218937e-04 ,5.494138e-07,9.202997e-01  ]
var times =[1940, 1944]
var Np = 500
var nvars 
var toler = 1e-17
var nlost = 0
//begin mif2
var Nmif = 1
var rw_size = .05, rw = new Array(Nmif).fill(null).map(() => Array(Np).fill(0.05))//change it to matrix
var coolFrac = 0.5 
var rwIndex = new Array(11).fill(0)
rwIndex[0] = 1; rwIndex[1] = 1; rwIndex[2] = 1; rwIndex[3] = 1 // only states

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

//* main function

var [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, R_0, I_0] = params
var state = [S_0, E_0, R_0, I_0, H = 0]
var [t0, tdata] = times
var nvars = state.length
var deltaT = 0.03832991102 // = 2/52 
var dt = 1/ 365.25

var timeLen = dataCases.length 

var va = 0, seas
var particles = []
var sampleNum = new Array(Np)
var doPredictionVariance = 0, doPredictionMean = 1, doFilterMean = 0 , allFail = 0
var predictionMean = Array(timeLen).fill(null).map(() => Array(nvars))
var predictionVariance = Array(timeLen).fill(null).map(() => Array(nvars))
var filterMean = []
var timeCounter = 0, ws = 0, vsq, sumsq, ess, condLoglik = [], Loglik = 0, lik
var stateVec =[]
var maxFail = 300,  hvalue =[]
// var R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4] ,rho = params[5], psi = params[6]

// if ( k === t0) {
  var Nlog = mathLib.toLogBarycentric([state[0], state[1], state[2], state[3]],4)
  var N = mathLib.fromLogBarycentric(Nlog, 4)
  var m = interpolPop(t0) / (N[0] + N[1] + N[2] + N[3]);
  state[0] = Math.round(m * N[0]),
  state[1] = Math.round(m * N[1]),
  state[2] = Math.round(m * N[2]),
  state[3] = Math.round(m * N[3])
//begin mif2
var stateNoiseM = Array(nvars).fill(null).map(() => Array(Np))
var states = Array(Np).fill(null).map(() => Array(nvars))
var s = (1 - 50 * timeLen * coolFrac) / (coolFrac - 1)
var cmn = (s + 1)/ (s + t0 + (Nmif - 1) * timeLen)
  for (cnt = 0; cnt < state.length; cnt++){
    for (np = 0; np < Np; np++){
      if(rwIndex[cnt] === 11){//no noise
        stateNoiseM[cnt][np] = state[cnt] * ( 1 + cmn * rw[Nmif - 1][np] * rnorm(1))// weighted??????????????????????????
      } else {
        stateNoiseM[cnt][np] = state[cnt]
      } 
    }   
  }
particles = mathjs.transpose(stateNoiseM) 
var h = []
//***********************************************TIME LOOP************************************
for (k = t0; k < Number(dataCases[dataCases.length - 2][0]) + deltaT / 3; k += deltaT){//Number(dataCases[dataCases.length - 2][0]) + deltaT / 3
  if (k <= tdata && k > tdata - deltaT) {
    k = tdata
  }
  var pop = interpolPop(k)
  var birthrate = interpolBirth(k) 
  var tt = (k - Math.floor(k)) * 365.25
  var loglik = 0
  var lik = new Array(Np)
  var weights = [], normalWeights = []

  //****************************************PARTICLE LOOP**************************************////////////////////////////////////////////////////////////////////////////
  for (np = 0; np < Np; np++){ //calc for each particle
    var trans = new Array(6).fill(0)
    
    var S = particles[np][0], E = particles[np][1], I = particles[np][2], R = particles[np][3], H = particles[np][4]
    // if ( k === 1944) {S=90986; E=688;I= 500; R=2352428.10446462;  H=319039}
    
    if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
      seas = 1 + amplitude * 0.2411 / 0.7589
    } else {
      seas = 1 - amplitude
    }                 
    var beta = R0 * (gamma + mu) * (sigma + mu) * seas / sigma //seasonal transmission rate
    var foi = beta * I / pop
    rate[0] = foi//force of infection
    rate[1] = mu// natural S death
    rate[2] = sigma// rate of ending of latent stage
    rate[3] = mu// natural E death
    rate[4] = gamma// recovery
    rate[5] = mu// natural I death       
    // transitions between classes
    steps = mathLib.numMapSteps(k, k + deltaT, dt)

    var del_t = (1 / steps )* deltaT
    for (let stp = 0; stp < steps; stp++) { // steps in each time interval
      var births = rpois(1, birthrate * (1 - va) * del_t )// Poisson births
      mathLib.reulermultinom(2, Math.round(S), 0, del_t, 0, rate, trans)
      mathLib.reulermultinom(2, Math.round(E), 2, del_t, 2, rate, trans)
      mathLib.reulermultinom(2, Math.round(I), 4, del_t, 4, rate, trans)
      S += (births - trans[0] - trans[1])
      E += (trans[0] - trans[2] - trans[3]) 
      I += (trans[2] - trans[4] - trans[5]) 
      R = pop - S - E - I
      H += trans[4] 
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
    // if ( k === 1944) {S=90986; E=688;I= 500; R=2352428.10446462;  H=319039}
    if (k >=  tdata) {hvalue.push([S,E,I,R,H])
      var rho = params[5], psi = params[6]
      var mn = rho * H
      var v = mn * (1.0 - rho + psi * psi * mn)
      var tol = 1.0e-18
      var modelCases = Number(dataCases[timeCounter][1])
      var likvalue = snippet.dmeasure(rho, psi, H, modelCases, giveLog = 0)
      weights.push(likvalue)
    } 
    // stateVec.push([states[np]])
  }////////////////////////////////////////////////////////////////end particle loop///////////////////////////////////////////////////////////////////////////////////////
if(k == 1944 || k == 1944.03832991102 || k==1944.0766598220398){console.log(k, S, E,I,R,H)}
  if (k >= tdata) {
    //normalize
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
    if (allFail) {
      lik = Math.log(toler) // minimum log-likelihood
      ess = 0  // zero effective sample size
    } else {
      ess = w * w / ws  // effective sample size
      lik = Math.log(w / Np)// mean of weights is likelihood
    }
    condLoglik.push(lik)
    Loglik += lik

    mathLib.nosortResamp(Np, normalWeights, Np, sampleNum, 0)
    for (np = 0; np < Np; np++) { // copy the particles
      particles[np] = [].concat(particles[sampleNum[np]])
      particles[np][nvars - 1] = 0
    }

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
        predictionMean[timeCounter][j] = sum
      }
  
      // compute prediction variance
      if (doPredictionVariance) {
        sumsq = 0
        for (let nrow = 0; nrow < Np; nrow++){
          if (states[nrow][j]) {
            vsq = states[nrow][j] - sum
            sumsq += vsq ** 2
          }
        }
        predictionVariance[timeCounter][j] = sumsq / (Np - 1) 
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
          filterMean[j] = ws / Np
        } else {      // weighted average
          ws = 0
          for (let nrow =0; nrow < Np; nrow++){
            if (states[nrow][j]) {
              ws += states[nrow][j] * weights[k]
            }
          }
          filterMean[j] = ws / w
        }
      }
    }
    timeCounter++
  }  
}//endTime

const createCsvWriter = require('csv-writer').createArrayCsvWriter;
const csvWriter = createCsvWriter({
  header: ['S', 'E', 'I', 'R', 'H'],
  path: './predmean.csv'
})
 
csvWriter.writeRecords(predictionMean)
  .then(() => {
  console.log('...predictionMean')
})

var createCsvWriter2 = require('csv-writer').createArrayCsvWriter;
var csvWriter2 = createCsvWriter2({
  header: ['S', 'E', 'I', 'R', 'H'],
  path: './states.csv'
})
csvWriter2.writeRecords(states)
  .then(() => {
  console.log('...states')
})

var createCsvWriter2 = require('csv-writer').createArrayCsvWriter;
var csvWriter2 = createCsvWriter2({
  header: ['S', 'E', 'I', 'R', 'H'],
  path: './predvar.csv'
})
 
csvWriter2.writeRecords(predictionVariance)
  .then(() => {
  console.log('...predvar')
})

var createCsvWriter2 = require('csv-writer').createArrayCsvWriter;
var csvWriter2 = createCsvWriter2({
  header:['S', 'E', 'I', 'R', 'H'],
  path: './hvalue.csv'
})
 
csvWriter2.writeRecords(hvalue)
  .then(() => {
  console.log('...hvalue')
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



    
      


