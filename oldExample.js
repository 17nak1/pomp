/**
 * Complete example using old example in Git/mif2 and Git/pfilter
 * Np=10000, Nmif=1 => loglik_mif ~ 3400-3500, loglik_pfilter ~ 5500-5800
 */
let rootDir ='.'
const fs = require('fs');
const { mif2 } = require('./mif2/src/mif2.js');
const { pfilter } = require('./pfilter/src/pfilter.js');
const snippet = require('./library/modelSnippet.js');
const { coef } = require("./mif2/src/mif2Helpers.js");

let dataCases = [];
let dataCasesTimes = [];
let dataCovar = [];
let dataCovarTimes = [];
let currentParams = []; 

// 1st data set; read all rows and delete last one if it is ['']
let temp, file;
let data;
file = fs.readFileSync(rootDir+'/samples/London_covar.csv').toString();
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
file = fs.readFileSync(rootDir+'/samples/London_BiDataMain.csv').toString()
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

currentParams = { R0: 31.3249, amplitude: 0.388362, gamma: 73.05, mu: 0.000647, sigma: 45.75, rho: 0.46, psi: 0.146, S_0: 0.034, E_0: 0.000234, I_0: 4.22e-7, R_0: 0.966 }

let cool_fraction = 0.05;
const pompData = {
  data :  dataCases,
  times:  dataCasesTimes,
  t0: 1940,
  rprocessDetail:  { type:"euler_sim", deltaT: 1/365.25 },
  covar: dataCovar,
  tcovar: dataCovarTimes,
  zeronames: snippet.zeronames,
  statenames: snippet.statenames,
  paramnames: [...snippet.paramsMod, ...snippet.paramsIc],
  covarnames: dataCovar_name,
  obsnames: dataCases_name,
};
let t = new Date();
// let mf = mif2(currentParams,
//   {object: pompData,
//     Nmif: 1,
//     transform: true,
//     rw_sd: snippet.determineRW("R0"),
//     Np: 10000,
//     varFactor: 2,
//     coolingType: "hyperbolic",
//     coolingFraction: cool_fraction
//   }
// )
// console.log((new Date() - t)/1000, mf.loglik, coef(mf));

t = new Date()
let pf = pfilter(currentParams,{object: pompData, params: currentParams, Np: 10000, predMean: true, filterMean: true, saveStates: true, maxFail: 3000})
console.log((new Date() - t)/1000, pf.loglik)

  
