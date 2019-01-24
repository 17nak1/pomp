var trans = new Array(6).fill(3)
var rate = new Array(6)
fs = require('fs')
let fmin = require ('fmin')
var mathLib = require('./mathLib')
var linearInterpolator = require('linear-interpolator/node_main')
const libR = require('lib-r-math.js')
const {
  Poisson,
  rng: { MersenneTwister },
  rng: { normal: { Inversion } }
} = libR
const mt = new MersenneTwister(0)// need reference so we "reset" PRNG
const { rpois } = Poisson(new Inversion(mt))
mt.init(1234)

//* Set the seed for rnorm-In R:RNGkind("L'Ecuyer-CMRG", normal.kind="Box-Muller");set.seed(1234) 
//or normal.kind="Box-Muller" or "BuggyKindermanRamage" or "Inversion" or "KindermanRamage"
var {
  Normal,
  rng: {
    LecuyerCMRG,
    normal: { BoxMuller }
  }
} = libR
const ad = new LecuyerCMRG(1234)
const { rnorm } = Normal(new BoxMuller(ad))
// data
var LondonBidata, LondonCovar, 
params = [49.05178717, 0.490067873611598, 73.05, 0.348209790270606, 45.66, 0.999999958052764, 1.07168303292765, 0.974365234, 4.835449219, 0.066357422, 0.927490234, 0.032927404, 2.46E-06, 0.967057826, 1.23E-05, 1940, 1944]

var rw_size = .05, rw = new Array(11).fill(rw_size), 
rwIdx = [0, 1, 2, 3, 4, 5, 6, 11, 12, 13, 14]
var N = [params[11], params[12], params[13], params[14]], Np = 3
var time =1944, Nmif = 1
// indx[1] = 1; indx[3] = 1; indx[5] = 1; indx[6] = 1

//* 1st data set
var London_covar = fs.readFileSync('./London_covar.csv').toString()
var LondonCovar = []
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonCovar.push(lines[i].split(','))
}
var dataCovar = [LondonCovar][0]
//* 2nd data set
LondonBidata = []
var London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonBidata.push(lines[i].split(','))
}
var dataCases = [LondonBidata][0]
//* main function****************************************************************

var [R0, amplitude, gamma, mu, sigma, rho, psi, alpha, iota, sigmaSE, cohort, S_0, E_0, R_0, I_0, t0, t1] = params
var estim = []
var place = []
var indx = new Array(12)
var coolFrac = 0.5
var rw_size = 0.05, delT = 0.03832991102// = 2/52 
var timeLen = dataCases.length;
for (let i = 0; i < params.length - 1; i++) {
  if (indx[i] === 1) {
    place.push(i)
    estim.push(params[i])
  }
}
var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
var interpolPop = linearInterpolator(d1)
var interpolBirth = linearInterpolator(d2)
  //* ODE function and ...
  
function poly (params, time, N) {
  
  var dt = 0.1
  var va = 0, seas, dy = new Array(6).fill(0)

  // var s = (1 - 50 * timeLen * coolFrac) / ( coolFrac - 1)
  // var cmn = (s + 1)/ (s + time + (Nmif - 1) * timeLen)
  // for (j = 0; j < rw.length; j++) {
  //   params[rwIdx[j]] += cmn * rw[j] * rnorm(1)
  // }
  // console.log(params)
  var R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4],
  alpha = params[7], iota = params[8], sigmaSE = params[9], cohort = params[10] 
  var S = N[0], E = N[1], R = N[2], I = N[3], W = 0, C = 0
  var pop = interpolPop(time)
  var birthrate = interpolBirth(time)
  // cohort effect
  if (Math.abs(time - Math.floor(time) - 251 / 365) < 0.5 * dt) {
    var br = cohort * birthrate / dt + (1 - cohort) * birthrate
  } else {
      br = (1 - cohort) * birthrate
  }
  // term-time seasonality
  var tt = (time - Math.floor(time)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - amplitude
  }                 
  
  var beta = R0 * (gamma + mu) * seas// transmission rate
  var foi = beta * Math.pow(I + iota, alpha) / pop// expected force of infection
  var dw = mathLib.rgammawn(sigmaSE, dt)// white noise (extrademographic stochasticity)
  rate[0] = foi * dw / dt// stochastic force of infection
  rate[1] = mu// natural S death
  rate[2] = sigma// rate of ending of latent stage
  rate[3] = mu// natural E death
  rate[4] = gamma// recovery
  rate[5] = mu// natural I death       
  var births = rpois(1, br * (1 - va) * dt)// Poisson births
  // transitions between classes
  var m = pop / (params[11] + params[12] + params[13] + params[14]),
  S = Math.round(m * params[11]),
  E = Math.round(m * params[12]),
  R = Math.round(m * params[13]),
  I = Math.round(m * params[14]),
  W = 0,
  C = 0,
  N = [S, E, R, I, W, C]
// console.log("ll",mathLib.reulermultinom(2,1,1,.1,1,[1,2], [0,0]))
  mathLib.reulermultinom(2, S, 0, dt, 0, rate, trans)
  mathLib.reulermultinom(2, E, 2, dt, 2, rate, trans)
  mathLib.reulermultinom(2, I, 4, dt, 4, rate, trans)
  dy[0] += births - trans[0] - trans[1]
  dy[1] += trans[0] - trans[2] - trans[3]
  dy[2] += trans[2] - trans[4] - trans[5]
  dy[3] = pop - S - E - I
  dy[4] += (dw - dt) / sigmaSE // standardized i.i.d. white noise
  dy[5] += trans[4]
  return dy
}
// console.log(trans,poly(params, 1944, N))


// ODE solver
function EulersMethod (params, covarData, delT) {
  var rho = params[5], psi = params[6], t0 = params[15], tdata = params[16]
  var steps = 1, arr2, arr = [], pop = interpolPop(t0)
  var m = pop / (params[11] + params[12] + params[13] + params[14]),
    S = Math.round(m * params[11]),
    E = Math.round(m * params[12]),
    R = Math.round(m * params[13]),
    I = Math.round(m * params[14]),
    W = 0,
    C = 0,
    N = [S, E, R, I, W, C]
  for (let k = t0; k <= 1945; k += delT) {//Number(dataCases[dataCases.length - 2][0]) + delT / 3
    N[4] = 0
    N[5] = 0
    if (k <= tdata && k > tdata - delT) {
      k = tdata
    }
    // for (let stp = 0; stp < Np; stp++) { // steps in each time interval
      arr2 = poly(params, k + delT, N)
      N = N.map((a, i) => Math.round(a + arr2[i] * 1 / steps * delT))
      // console.log(N)
    // }
    C = Math.round(N[5])
    var mn = rho * C
    var v = mn * (1.0 - rho + psi * psi * mn)
    var tol = 1.0e-18
    var cases = rnorm(1, mn, Math.sqrt(v) + tol)
    if (cases > 0) {
      cases = Math.round(cases)
    } else {
      cases = 0
    }
    if (k > tdata - 2 * delT) {
      if (k <= tdata - delT ) {
      k = tdata - delT
    }
      console.log(C)
      arr.push([k + delT , N[0], N[1], N[2], N[3], C])
    }
  }
  return arr
}
//* *********************************************************************
// function logLik (estim) {
//   var lik, loglik = 0
//   for (let i = 0; i < estim.length; i++) {
//     params[place[i]] = Math.exp(estim[i])
//   }
//   // console.log(estim,params)
//   var rho = params[5], psi = params[6]
//   var simCases = EulersMethod(params, dataCovar, delT)
//   for (let i = 0; i < simCases.length ; i++) {
//     var mn = rho * simCases[i][5]
//     var v = mn * (1.0 - rho + psi * psi * mn)
//     var tol = 1.0e-18
//     var modelCases = Number(dataCases[i][1])
//     if (modelCases > 0.0) {
//       lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
//     } else {
//       lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol
//     }
//     loglik = loglik + Math.log(lik)
//   }
//   console.log(loglik, estim)
//   return [-(loglik).toFixed(6)]
// }
EulersMethod (params, dataCovar, delT)