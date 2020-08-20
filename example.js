/**
 * Complete example using mif2 and replicating pfilter
 *
 */
let rootDir ='.'
const fs = require('fs');
const { mif2 } = require('./mif2/src/mif2.js');
const { pfilter } = require('./pfilter/src/pfilter.js');
const snippet = require('./library/ModelSnippetCOVID3.js');
const { coef } = require("./mif2/src/mif2Helpers.js");

let dataCases = [];
let dataCasesTimes = [];
let dataCovar = [];
let dataCovarTimes = [];
let currentParams = []; 

// 1st data set; read all rows and delete last one if it is ['']
let temp, file;
let data;
file = fs.readFileSync(rootDir+'/samples/covars.csv').toString();
let lines = file.split(/\r\n|\n/);
let dataCovar_name = lines[0].replace(/['"]+/g, '').split(',');
dataCovar_name.shift();
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCovarTimes.push(temp[0]);
    data = {};
    for(let j = 0; j < temp.length - 1; j++){
      data[dataCovar_name[j]] = temp[j + 1];
    }
    dataCovar.push(data)
  }
}

//* 2nd data set
file = fs.readFileSync(rootDir+'/samples/data.csv').toString()
lines = file.split(/\r\n|\n/);
let dataCases_name = lines[0].replace(/['"]+/g, '').split(',');
dataCases_name.shift();
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCasesTimes.push(temp[0]);
    data = {};
    for(let j = 0; j < temp.length - 1; j++){
      data[dataCases_name[j]] = temp[j + 1];
    }
    dataCases.push(data)
  }
}

//* 3nd data set and names
file = fs.readFileSync(rootDir+'/samples/initial_parameters.csv').toString()
lines = file.split(/\r\n|\n/);
let currentParams_name = lines[0].replace(/['"]+/g, '').split(',');
for (let i = 1; i < lines.length ; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    data = {};
    for(let j = 0; j < temp.length; j++){
      data[currentParams_name[j]] = temp[j];
    }
    currentParams.push(data)
  }
}

currentParams = {
  betaI:1.40758806343178,
  iota:0.0177042779901833,
  beta_sd:0,
  sigma:0.2,
  kappa:1,
  gammaI:0.172133302616885,
  gammaH:0.635123961314889,
  gammaC:1.31234713108046,
  gammaV:0.327003197582952,
  TF:16031.6700539306,
  rho:0.487576658764406,
  theta:0.415700726852156,
  dI0:0.485528264998708,
  dP0:0.176249510892225,
  dT0:0.263185388981837,
  dB0:0,
  dI1:0.331872295141556,
  dP1:0.259732942470986,
  dT1:0.630072923852895,
  dB1:0,
  qP:0.722056284667341,
  qH:0.137405612125941,
  qC:0.883271524025973,
  mI:0.00550468107974399,
  mC:0.197242061100997,
  mV:0.685825470240576,
  S0:1,
  EQ0:0,
  PQ0:0,
  IQ0:0,
  E0:0,
  P0:0,
  I0:0,
  H0:0,
  C0:0,
  V0:0,
  M0:0}
  
let globals = { nstageE: 3, nstageP: 3, nstageI: 3, nstageH: 3, nstageC: 3, nstageV: 3, pop: 10e6, T0: 75, T1: 139 };

let cool_fraction = 0.05;
const pompData = {
  data :  dataCases,
  times:  dataCasesTimes,
  t0: 0,
  rprocessDetail:  { type:"euler_sim", deltaT: 0.1 },
  covar: dataCovar,
  tcovar: dataCovarTimes,
  zeronames: snippet.zeronames,
  statenames: snippet.statenames,
  paramnames: [...snippet.paramsMod, ...snippet.paramsIc],
  covarnames: dataCovar_name,
  obsnames: dataCases_name,
  globals: globals,
};
let t = new Date();
// let mf = mif2(currentParams,
//   {object: pompData,
//     Nmif: 1,
//     transform: true,
//     rw_sd: snippet.determineRW(1),
//     Np: 200,
//     varFactor: 2,
//     coolingType: "hyperbolic",
//     coolingFraction: cool_fraction
//   }
// )
// console.log((new Date() - t)/1000, mf.loglik);

t = new Date()
let pf = pfilter(currentParams,{object: pompData, params: currentParams, Np: 2,predMean: true, filterMean: true, saveStates: true, maxFail: 3000})
console.log((new Date() - t)/1000, pf.loglik);

const createCsvWriter = require('csv-writer').createObjectCsvWriter;
let header = [];
for (let i = 0; i < Object.keys(pf.predMean[0]).length; i++) {
  header.push({id: Object.keys(pf.predMean[0])[i], title: Object.keys(pf.predMean[0])[i]})
}

const csvWriter = createCsvWriter({
    path: './oo.csv',
    header:header,
});
 
 
csvWriter.writeRecords(pf.predMean)
.then(() => {
    console.log('...Done');
});

  
