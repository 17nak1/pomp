
const fs = require('fs');
const { trajMatch } = require('./src/trajMatch.js');

const create_dataset = require('../library/CreateDataset.js');
const create_covars = require('../library/CreateCovars.js');
const snippet = require('../library/modelSnippetCOVID3.js');
const { coef } = require("../library/helpers");
const sobolSeq =require('../library/generate-sobol/sobolSeq.js');
const { pfilter } = require('../pfilter/src/pfilter.js');
const { mif2 } = require('../mif2/src/mif2.js');

let endTime = "2020-07-16";
SobolNumberOfPoints = 10;
let lowerBounds = {betaI: 0, theta: 0, iota: 0, beta_sd: 0,
  dI0: 0, dP0: 0, dT0: 0, dB0: 0,
  dI1: 0, dP1: 0, dT1: 0, dB1: 0,
  qP: 0, qH: 0, qC: 0.5, mI: 0, mC: 0, mV: 0.5,
  sigma: 1/5, kappa: 1/1, gammaI: 1/5, gammaH: 1/5, gammaC: 1/10, gammaV: 1/10, rho: 0, TF: 4e3,
  S0: 1,EQ0: 0,PQ0: 0,IQ0: 0,E0: 0,P0: 0,I0: 0,H0: 0,C0: 0,V0: 0,M0: 0}

let upperBounds = {betaI: 1, theta: 1, iota: 100, beta_sd: 0,
  dI0: 0.6, dP0: 0.6, dT0: 0.6, dB0: 0,
  dI1: 0.5, dP1: 0.5, dT1: 0.5, dB1: 0,
  qP: 0.5, qH: 1, qC: 1, mI: 0.1, mC: 1, mV: 1,
  sigma: 1/5, kappa: 1/1, gammaI: 1/5, gammaH: 1/1, gammaC: 1/1, gammaV: 1/1,rho: 1, TF: 8e3,
  S0: 1,EQ0: 0,PQ0: 0,IQ0: 0,E0: 0,P0: 0,I0: 0,H0: 0,C0: 0,V0: 0,M0: 0}

let sobolSet = sobolSeq.sobolDesign( lowerBounds, upperBounds, SobolNumberOfPoints);

paramsFixed = ["beta_sd","dB0", "dB1","sigma","kappa"];
selectedParams = [...snippet.paramsMod, ...snippet.paramsIc];
let paramsFit = snippet.paramsMod;
for (let i = 0; i < paramsFixed.length; i++) {
  paramsFit = paramsFit.filter(e => e !== paramsFixed[i])
}

// read all rows and chaeck the border time and convert the selected colnames
let file;
file = fs.readFileSync('../samples/ParamSet_run1.csv').toString();
let lines = file.split(/\r\n|\n/);

let paramsetData = [];
let paramsetHeader = lines[0].replace(/['"]+/g, '').split(',');
for (let i = 1; i < lines.length; i++) {
  let temp = lines[i].split(',');
  if(temp.length > 1) {
    let tempParamset =	{};
    for(let j = 0; j < temp.length; j++){
      tempParamset[paramsetHeader[j]] = Number(temp[j]);
    }
    paramsetData.push(tempParamset);
  }
}


// Generate covars, data and pomp object
data = create_dataset('../samples/ON.csv','../samples/covidtesting.csv', endTime)
covars = create_covars('../samples/covidtesting.csv',endTime)
let t1 = 75;
let t2 = 139;
globals = { nstageE: 3, nstageP: 3, nstageI: 3, nstageH: 3, nstageC: 3, nstageV: 3, pop: 10e6, T0: 75, T1: 139 };
console.log(data.length)
let dataHeader = data.shift();
let dataCases = data.map((x,i,arr) => { 
  a = {};
  a[dataHeader[1]] = x[1];
  a[dataHeader[2]] = x[2];
  a[dataHeader[3]] = x[3];
  a[dataHeader[4]] = x[4];
  a[dataHeader[5]] = x[5];
  return a;});
let dataCasesTimes = data.map((x,i,arr) => x[0]);
dataHeader.shift();

let covarHeader = covars.shift();
let dataCovar = covars.map((x,i,arr) => { 
  a = {};
  a[covarHeader[1]] = x[1];
  return a;});
let dataCovarTimes = covars.map((x,i,arr) => x[0]);
covarHeader.shift();

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
  covarnames: covarHeader,
  obsnames: dataHeader,
  globals: globals,
};
current_params = [{betaI: 0.0326576783487093,
theta: 1.1102230246251565e-16,
iota: 104.83105214413663,
beta_sd: 0,
dI0: 1,
dP0: 0.0561651857350145,
dT0: 0.9999999999999996,
dB0: 0,
dI1: 0.004495992802272308,
dP1: 0.8295330118538808,
dT1: 0.9999935133062721,
dB1: 0,
qP: 0.018757642630112537,
qH: 1,
qC: 0.9999999999999996,
mI: 0.9999999999999913,
mC: 0,
mV: 0,
sigma: 0.2,
kappa: 1,
gammaI: 0.003894693431270531,
gammaH: 0.02123614427836244,
gammaC: 0.15669834006885247,
gammaV: 0.02037978756358223,
rho: 0.551373266677974,
TF: 886.5017868735129,
S0: 1,
EQ0: 0,
PQ0: 0,
IQ0: 0,
E0: 0,
P0: 0,
I0: 0,
H0: 0,
C0: 0,
V0: 0,
M0: 0 }]
let tm = trajMatch(current_params[0],{object: pompData, est: [], transform: true, method: "subplex"})

console.log('finished.',coef(tm), tm.value);
t = new Date()

// let pf = pfilter(paramsetData[0],{object: pompData, params: paramsetData[0], Np: 100,filterMean: true, maxFail: 3000})
// console.log((new Date() - t)/1000, pf.loglik);
// const createCsvWriter = require('csv-writer').createObjectCsvWriter;
// let header = [];
// for (let i = 0; i < Object.keys(pf.filterMean[0]).length; i++) {
//   header.push({id: Object.keys(pf.filterMean[0])[i], title: Object.keys(pf.filterMean[0])[i]})
// }

// const csvWriter = createCsvWriter({
//     path: './oo.csv',
//     header:header,
// });
 
// csvWriter.writeRecords(pf.filterMean)
// .then(() => {
//     console.log('...Done');
// });

// // let mf = mif2(paramsetData[0],
// //   {object: pompData,
// //     Nmif: 1,
// //     transform: true,
// //     rw_sd: snippet.determineRW(1),
// //     Np: 1000,
// //     varFactor: 2,
// //     coolingType: "hyperbolic",
// //     coolingFraction: 0.05
// //   }
// // )
// // console.log((coef(mf),new Date() - t)/1000, mf.loglik);