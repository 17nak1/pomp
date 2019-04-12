

fs = require('fs')
let fmin = require ('fmin')
let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')
let simulator = require ('./simulator.js')

//////////////////////////////////////////data///////////////////////////////////////
let dataCases = [], dataCovar = []
let params = [3.132490e+01, 3.883620e-01, 7.305000e+01, 6.469830e-04, 4.566000e+01, 4.598709e-01, 1.462546e-01, 3.399189e-02, 2.336327e-04, 4.221789e-07, 9.657741e-01 ]
let maxFail = Infinity
let Np = 10
console.log("Np",Np)
let toler = 1e-17


//* 1st data set
let London_covar = fs.readFileSync('./London_covar.csv').toString()
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length - 1; i++) {
  dataCovar.push(lines[i].split(','))
}

//* 2nd data set
let London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length - 1; i++) {
  dataCases.push(lines[i].split(','))
}


let d1 = []// read time and population from 1st data and make interpolation function
let d2 = []// read time and birthrate from 1st data and make interpolation function
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
let interpolPop = mathLib.interpolator (d1)
let interpolBirth = mathLib.interpolator (d2)
let START = new Date()


let [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, I_0, R_0] = params
let [t0, tdata] = [1940, 1944]
let nvars = 5
let deltaT = 14 / 365.25
let dt = 1 / 365.25
let doPredictionVariance = 0, doPredictionMean = 1, doFilterMean = 0 , allFail = 0

let timeLen = dataCases.length 
let nlost = 0

let rate = [], trans = []
let particles = [], state =[]
let sampleNum = Array.from(Array(Np).keys())
let condLoglik = []
let stateSaved =[]

let timeCountData = 0, ws ,w , vsq, sumsq, ess, loglik = 0, lik 

let predictionMean, predictionVariance, filterMean
let states = Array(Np).fill(null).map(() => Array(nvars))
let weights, normalWeights, S, E, I, R, del_t, ST, simulateValue
let modelCases, likvalue
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
var aa = new Array(Np).fill(null).map(() => [].concat(state))

// Time loop
for (k = t0; k <= Number(dataCases[timeLen - 2][0]) + deltaT / 3 ; k += deltaT){//Number(dataCases[timeLen - 2][0]) + deltaT / 3
  if ( k > tdata - deltaT && k <= tdata) {
    k = tdata
  }
  weights = []; normalWeights = []
  
  if ( k > t0) {
    for (np = 0; np < Np; np++) { // copy the particles
      aa[np] = [].concat(particles[sampleNum[np]])
      aa[np][nvars - 1] = 0
    }
  } 
  
  //**PARTICLE LOOP
  for (np = 0; np < Np; np++){ //calc for each particle
    trans = []

    particles[np] = simulator.simulate(aa[np], k, tdata, deltaT, dt, timeCountData, interpolPop, interpolBirth,params, Number(dataCases[dataCases.length - 1]), Number(dataCases[timeCountData + 1][0]))
    
    S = particles[np][0] 
    E = particles[np][1]
    I = particles[np][2] 
    R = particles[np][3] 
    H = particles[np][4] 

    // console.log(k)
     
    //***********weight*************
    if (k >= Number(dataCases[0][0])){
      if (stateSaved) {
        stateSaved.push(particles[np]) //[S,E,I,R,H])
      }
      modelCases = Number(dataCases[timeCountData][1])
      likvalue = snippet.dmeasure(rho, psi, H, modelCases, giveLog = 0)
      weights.push(likvalue)
      
    }
  }//  end particle loop
  
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
    timeCountData++ 
  }
}//endTime

console.log(loglik)
const createCsvWriter = require('csv-writer').createArrayCsvWriter;
const csvWriter = createCsvWriter({
  header: ['S', 'E', 'I', 'R', 'H'],
  path: './predmean.csv'
})
 
csvWriter.writeRecords(predictionMean)
  .then(() => {
  console.log('...predictionMean')
})

// const createCsvWriter11 = require('csv-writer').createArrayCsvWriter;
// const csvWriter11 = createCsvWriter11({
//   header: ['S', 'E', 'I', 'R', 'H'],
//   path: './predvar.csv'
// })
// csvWriter11.writeRecords(predictionVariance)
//   .then(() => {
//   console.log('...predictionVar')
// })

// var createCsvWriter1 = require('csv-writer').createArrayCsvWriter;
// var csvWriter1 = createCsvWriter1({
//   header: ['S', 'E', 'I', 'R', 'H'],
//   path: './filterMean.csv'
// })
// csvWriter1.writeRecords(filterMean)
//   .then(() => {
//   console.log('...filterMean')
// })

// var createCsvWriter2 = require('csv-writer').createArrayCsvWriter;
// var csvWriter2 = createCsvWriter2({
//   header: [],
//   path: './condLogliklihood.csv'
// }) 
// csvWriter2.writeRecords(condLoglik)
//   .then(() => {
//   console.log('...condLoglik')
// })

// var createCsvWriter3 = require('csv-writer').createArrayCsvWriter;
// var csvWriter3 = createCsvWriter3({
//   header:['S', 'E', 'I', 'R', 'H'],
//   path: './stateSaved.csv'
// })
// csvWriter3.writeRecords(stateSaved)
//   .then(() => {
//   console.log('...stateSaved')
// })
  
console.log('running time:',new Date() - START)



    
      


