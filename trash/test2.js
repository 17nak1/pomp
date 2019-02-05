
var rate = new Array(6)
fs = require('fs')
let fmin = require ('fmin')
var mathLib = require('./mathLib')
const mathjs = require('mathjs') 
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
var vv =[]
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
//////////////////////////////////////////data///////////////////////////////////////
var LondonBidata, LondonCovar, 
params = [29.2655634435, 0.3944025839968, 73.05, 0.004303868154365, 45.66, 0.417917730919406, 0.530358893934836, 0.020285841729738, 3.49E-05, 8.01E-05, 0.979679198700577, 1940, 1944]
var Np = 50
var toler = 10e-8
var w=0 , ws=0
//begin mif2
var Nmif = 1
var rw_size = .05, rw = new Array(Nmif).fill(null).map(() => Array(Np).fill(0.05))//change it to matrix
var rwIn = new Array(11).fill(1)
rwIn[2] = 0; rwIn[4] = 0
var indx = new Array(11).fill(1)
indx[2] = 0; indx[4] = 0
//end
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

var [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, R_0, I_0, t0, t1] = params
var estim = []
var place = []
var coolFrac = 0.5
var delT = 0.03832991102// = 2/52 
var dt = 1/ 365.25
var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
var d3 = []
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
for (let i = 0; i < dataCases.length - 1; i++) {
  d3.push([Number(dataCases[i][0]), Number(dataCases[i][1])])
}
var interpolPop = linearInterpolator(d1)
var interpolBirth = linearInterpolator(d2)
  
// function poly (paramNoise, time, trans) {
  var cases = []
  var va = 0, seas
  var t0 = params[11], tdata = params[12] 
  params.length -= 2 
  var particles = [], sampleNum = new Array(Np).fill(0).map(Number.call, Number) , newParticles = []
  //***********************************************TIME LOOP************************************
  for (k = t0; k <= 1954 + 8*delT; k +=delT){//Number(dataCases[dataCases.length - 2][0]) + delT / 3
    if (k <= tdata && k > tdata - delT) {
      k = tdata
    }
  if ( k === t0) {
    var Nlog = mathLib.toLogBarycentric([params[7], params[8], params[9], params[10]],4)
    var N = mathLib.fromLogBarycentric(Nlog, 4)
    var m = interpolPop(k) / (N[0] + N[1] + N[2] + N[3]);
    params[7] = Math.round(m * N[0]),
    params[8] = Math.round(m * N[1]),
    params[9] = Math.round(m * N[2]),
    params[10] = Math.round(m * N[3])
  //begin mif2
  var paramNoiseM = Array(params.length).fill(null).map(() => Array(Np))

  var timeLen = dataCases.length;
  var s = (1 - 50 * timeLen * coolFrac) / (coolFrac - 1)
  var cmn = (s + 1)/ (s + k + (Nmif - 1) * timeLen)
    for (cnt = 0; cnt < params.length; cnt++){
      for (np = 0; np < Np; np++){
        if(rwIn[cnt] === 1){
          paramNoiseM[cnt][np] = params[cnt] * ( 1 + cmn * rw[Nmif - 1][np] * rnorm(1))// weighted??????????????????????????
        } else {
          paramNoiseM[cnt][np] = params[cnt]
        } 
      }   
    }
    particles = mathjs.transpose(paramNoiseM)
  }//endif ( k === t0)
  var pop = interpolPop(k)
  var birthrate = interpolBirth(k) 
  var tt = (k - Math.floor(k)) * 365.25
  var loglik = 0
  var lik = new Array(Np)
  var weights = []
  if (k > tdata){
   particles = newParticles
   newParticles = []
 }
 // if ( k >=1944 -2*delT) console.log("enter1",newParticles,k)
// console.log("enter1",particles,k)
var aa = 1,t_1 = 1944.1916495550995, t_2 = 1944.2299794661194
// if ( k >=1944 -2*delT) console.log("enter1",particles,k)//;console.log(particles[1] === particles[10])

  //****************************************PARTICLE LOOP**************************************////////////////////////////////////////////////////////////////////////////
  for (np = 0; np < particles.length; np++){//calc for each particle
    // if ( k == t_2 ) 
    // console.log(np,"enter2",particles,k)
    var trans = new Array(6).fill(0)
    var R0 = particles[np][0], amplitude = particles[np][1], gamma = particles[np][2], mu = particles[np][3], sigma = particles[np][4] 
    var S = particles[np][7], E = particles[np][8], I = particles[np][9], R = particles[np][10], H = 0
    if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
      seas = 1 + amplitude * 0.2411 / 0.7589
    } else {
      seas = 1 - amplitude
    }                 
    var beta = R0 * (gamma + mu) * (sigma + mu) * seas / sigma//seasonal transmission rate
    var foi = beta * I / pop
    rate[0] = foi//force of infection
    rate[1] = mu// natural S death
    rate[2] = sigma// rate of ending of latent stage
    rate[3] = mu// natural E death
    rate[4] = gamma// recovery
    rate[5] = mu// natural I death       
    // transitions between classes
    steps = mathLib.numEulerSteps(k, k + delT, dt)
    var del_t = 1 / steps * delT
    for (let stp = 0; stp < steps; stp++) { // steps in each time interval
      var births = rpois(1, birthrate * (1 - va) * del_t )// Poisson births
      mathLib.reulermultinom(2, Math.round(S), 0, del_t, 0, rate, trans)
      mathLib.reulermultinom(2, Math.round(E), 2, del_t, 2, rate, trans)
      mathLib.reulermultinom(2, Math.round(I), 4, del_t, 4, rate, trans)//;if ( k == 1944.3066392881592  && np == 0) console.log(trans, S,E,I,R, H)
      S += (births - trans[0] - trans[1])
      E += (trans[0] - trans[2] - trans[3]) 
      I += (trans[2] - trans[4] - trans[5]) 
      // E = (E > 0) ? (E + trans[0] - trans[2] - trans[3]) : 0
      // I = (I > 0) ? (I + trans[2] - trans[4] - trans[5]) : 0
      R = pop - S - E - I
      H += trans[4] 
    }
    // if ( k == t_2 && np == aa) console.log(trans, S,E,I,R, H)
    particles[np][7] = S
    particles[np][8] = E
    particles[np][9] = I
    particles[np][10] = R
    // if (k == t_2 && np == aa)
     // console.log("ex",particles[np])
    //***********RESAMPLE*************
    if (k >= tdata) {
      var rho = particles[np][5], psi = particles[np][6]
      var mn = rho * H
      var v = mn * (1.0 - rho + psi * psi * mn)
      var tol = 1.0e-18
      var modelCases = Number(dataCases[Math.ceil((k - tdata) / delT)][1])
      if(!isNaN(modelCases)){
        if (modelCases > 0.0) {
          lik[np] = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
        } else {
          lik[np] = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0)) + tol
        }
      } else {
        lik[np] = 1
      }
      loglik += Math.log(lik[np])
    }
  }//end particle loop/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  if (k >= tdata) {
    if (loglik != 0) {
      for(let np = 0; np < particles.length; np++) {
        weights.push(Math.log(lik[np]) / loglik)
      }
      var sampleNum = mathLib.nosortResamp(particles.length, weights, particles.length, 0)
      for (sn = 0; sn < sampleNum.length; sn++) { // copy the particles
        while (sampleNum[sn] === sampleNum[sn + 1]){
          sn++
        } 
        newParticles.push(particles[sampleNum[sn]])
      }
    } else {
      newParticles = particles
    }
  }//if (k >= tdata) 

}//endTime

console.log(newParticles)


