
const fs = require('fs');
const { trajMatch } = require('./src/trajMatch.js');
// const create_dataset = require('../library/CreateDataset.js');
// const create_covars = require('../library/CreateCovars.js');
const snippet = require('../library/modelSnippet.js');
const { coef } = require("../library/helpers");
const sobolSeq =require('../library/generate-sobol/sobolSeq.js');
const { pfilter } = require('../pfilter/src/pfilter.js');
const { mif2 } = require('../mif2/src/mif2.js');

const globals = { nstageE: 3, nstageP: 3, nstageI: 3, nstageH: 3, nstageC: 3, nstageV: 3, pop: 3e6, T0: 75, T1: 139, T2: 212, TC: 21};
SobolNumberOfPoints = 10;
let lowerBounds = {betaI: 0, theta: 0, iota: Math.log((globals.TC-6)/((globals.TC-6)-1)), beta_sd: 0,
  dI0: 0, dP0: 0, dT0: 0, dB0: 0,
  dI1: 0, dP1: 0, dT1: 0, dB1: 0,
  qP: 0, qH: 0, qC: 0.5, mI: 0, mC: 0, mV: 0.5,
  sigma: 1/5, kappa: 1/1, gammaI: 1/5, gammaH: 1/5, gammaC: 1/10, gammaV: 1/10, rho: 0, TF: 4e3,
  S0: 1,EQ0: 0,PQ0: 0,IQ0: 0,E0: 0,P0: 0,I0: 0,H0: 0,C0: 0,V0: 0,M0: 0}

let upperBounds = {betaI: 1, theta: 1, iota: Math.log((globals.TC-6)/((globals.TC-6)-1)), beta_sd: 0,
  dI0: 0.6, dP0: 0.6, dT0: 0.6, dB0: 0,
  dI1: 0.5, dP1: 0.5, dT1: 0.5, dB1: 0,
  qP: 0.5, qH: 1, qC: 1, mI: 0.1, mC: 1, mV: 1,
  sigma: 1/5, kappa: 1/1, gammaI: 1/5, gammaH: 1/1, gammaC: 1/1, gammaV: 1/1,rho: 1, TF: 8e3,
  S0: 1,EQ0: 0,PQ0: 0,IQ0: 0,E0: 0,P0: 0,I0: 0,H0: 0,C0: 0,V0: 0,M0: 0}

let sobolSet = sobolSeq.sobolDesign( lowerBounds, upperBounds, SobolNumberOfPoints);

paramsFixed = ["beta_sd","dB0", "dB1", "dB2", "sigma","kappa","gammaI","iota"];
selectedParams = [...snippet.paramsMod, ...snippet.paramsIc];
let paramsFit = snippet.paramsMod;
for (let i = 0; i < paramsFixed.length; i++) {
  paramsFit = paramsFit.filter(e => e !== paramsFixed[i])
}

// read all rows and chaeck the border time and convert the selected colnames
// let file;
// file = fs.readFileSync('./DetModel_Toronto_all.csv').toString();
// let lines = file.split(/\r\n|\n/);

// let paramsetData = [];
// let paramsetHeader = lines[0].replace(/['"]+/g, '').split(',');
// for (let i = 1; i < lines.length; i++) {
//   let temp = lines[i].split(',');
//   if(temp.length > 1) {
//     let tempParamset =	{};
//     for(let j = 0; j < temp.length; j++){
//       tempParamset[paramsetHeader[j]] = Number(temp[j]);
//     }
//     paramsetData.push(tempParamset);
//   }
// }

let dataCases = [];
let dataCasesTimes = [];
let dataCovar = [];
let dataCovarTimes = [];

// 1st data set; read all rows and delete last one if it is ['']
let temp;
let data;
file = fs.readFileSync('../samples/covar.csv').toString();
lines = file.split(/\r\n|\n/);
let covarHeader = lines[0].replace(/['"]+/g, '').split(',');
covarHeader.shift();
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => (x));
    dataCovarTimes.push(Number(temp[0]));
    data = {};
    for(let j = 0; j < temp.length - 1; j++){
      data[covarHeader[j]] = temp[j + 1];
    }
    dataCovar.push(data)
  }
}

//* 2nd data set
file = fs.readFileSync('../samples/data.csv').toString()
lines = file.split(/\r\n|\n/);
let dataHeader = lines[0].replace(/['"]+/g, '').split(',');
dataHeader.shift();
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCasesTimes.push(temp[0]);
    data = {};
    for(let j = 0; j < temp.length - 1; j++){
      data[dataHeader[j]] = temp[j + 1];
    }
    dataCases.push(data)
  }
}

let current_params = [{ betaI: 0.61702834021533,  iota: 0.0689928714869514,  beta_sd: 0,  sigma: 0.2,  kappa: 1,  gammaI: 0.2,  gammaH: 0.0751112327637673,  gammaC: 1.01453232294609,  gammaV: 0.340505807671711,  TF: 19.3892866984778,  rho: 0.00654555126593593,  theta: 0.765522585934224,  dI0: 0.453813588651899,  dP0: 0.0608391515618878,  dT0: 0.164963537069865,  dB0: 0,  dI1: 0.539588845693196,  dP1: 0.0521235391641814,  dT1: 0.100341602811159,  dB1: 0,  dI2: 0.9937651045604,  dP2: 0.519916004211288,  dT2: 0.311360480290506,  dB2: 0,  qP: 0.00102436696897488,  qH: 0.997347577623437,  qC: 0.992015571875051,  mI: 0.000430326186262243,  mC: 0.62011387256283,  mV: 0.149268369336004,  S0: 1,  EQ0: 0,  PQ0: 0,  IQ0: 0,  E0: 0,  P0: 0,  I0: 0,  H0: 0,  C0: 0,  V0: 0,  M0: 0,  LogLik: -4446.1799412075 }]
const pompData = {
  data :  dataCases,
  times:  dataCasesTimes,
  t0: 0,
  skeletonDetail:  { type:"map", deltaT: 0.1 },
  rprocessDetail:  { type:"euler_sim", deltaT: 0.1 },
  covar: dataCovar,
  tcovar: dataCovarTimes,
  zeronames: snippet.zeronames,
  statenames: snippet.statenames,
  paramnames: [...snippet.paramsMod, ...snippet.paramsIc],
  globals: globals,
};
t = new Date()
let tm = trajMatch(current_params[0],{object: pompData, est: [], transform: true, method: "subplex"});
// console.log('finished.',tm, tm.loglik, new Date() - t);

// current_params[0]['beta_sd'] = 0.01;
// current_params[0]['dB0'] = 0.2;
// current_params[0]['dB1'] = 0.2;
let pf = pfilter(current_params[0],{object: pompData, Np: 1000,filterMean: true, saveStates: true, maxFail: 3000})
console.log(pf.loglik)
let h = [];  
let finalParams = [pf.params];                              /* Log params */
let keys = Object.keys(finalParams[0]);
for ( let i = 0; i < keys.length; i++) {
  h.push({id:keys[i], title:keys[i]})
}
// const createCsvWriter = require('csv-writer').createObjectCsvWriter;
// const csvWriter = createCsvWriter({
//   path: './finalParams.csv',
//   header: h
// });
// csvWriter.writeRecords(finalParams)      
// .then(() => {
//   console.log('...Done');
// });

// let headerStates = [];
// for (let i = 0; i < Object.keys(pf.filterMean[0]).length; i++) {
//   headerStates.push({id: Object.keys(pf.filterMean[0])[i], title: Object.keys(pf.filterMean[0])[i]})
// }

// const csvWriterFilter = createCsvWriter({   /* Log filterMean */
//   path: './filterMean.csv',
//   header: headerStates
// });
// csvWriterFilter.writeRecords(pf.filterMean)      
// .then(() => {
//   console.log('...Done');
// });

// const csvStates = createCsvWriter({          /* Log the last saved states (at the last time step) */
//   path: './savedStates.csv',
//   header: headerStates
// });
// csvStates.writeRecords(pf.saveStates)        /* pf.saveStates is the last saved  states */
// .then(() => {
//   console.log('...Done');
// })
  


// let mf = mif2(paramsetData[0],
//   {object: pompData,
//     Nmif: 1,
//     transform: true,
//     rw_sd: snippet.determineRW(),
//     Np: 1000,
//     varFactor: 2,
//     coolingType: "hyperbolic",
//     coolingFraction: 0.05
//   }
// )
// console.log((coef(mf),new Date() - t)/1000, mf.loglik);