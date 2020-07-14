/**
 * 
 *
 */
// let mif2 = require('./pfilter.js');
let snippet = require('../library/modelSnippet.js');
let fs = require('fs');
let { pfilter } = require('./pfilter.js');
let pomp = require('../library/pomp.js');

let dataCases = [];
let dataCasesTimes = [];
let dataCovar = [];
let dataCovarTimes = [];
let currentParams = []; 
let rootDir = '..'
// 1st data set; read all rows and delete last one if it is ['']
let temp, file;
file = fs.readFileSync(rootDir+'/samples/London_covar.csv').toString();
let lines = file.split(/\r\n|\n/);
let dataCovar_name = lines[0].replace(/['"]+/g, '').split(',');
dataCovar_name.shift();
for (let i = 1; i < lines.length; i++) {
  temp = lines[i].split(',');
  if(temp.length > 1) {
    temp = temp.map(x => Number(x));
    dataCovarTimes.push(temp[0]);
    dataCovar.push(temp.slice(1));
  }
}

//* 2nd data set
let data;
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

currentParams = currentParams[0]//Only for this example, we need loop over currentParams

const mypomp = new pomp({
  data :  dataCases,
  times:  dataCasesTimes,
  t0: 1940,
  rprocess :  { type:"euler_sim", stepFunction: snippet.rprocess, deltaT: 1/365.25 },
  rmeasure: snippet.rmeas,
  covar: dataCovar,
  tcovar: dataCovarTimes,
  dmeasure: snippet.dmeasure,
  zeronames: snippet.zeronames,
  initializer: snippet.initz,
  toEstimationScale: snippet.toEst,
  fromEstimationScale: snippet.fromEst,
  statenames: snippet.statenames,
  paramnames: snippet.paramnames,
  coef: currentParams,
  covarnames: dataCovar_name,
  obsnames: dataCases_name,
});
let t = new Date()
let pf = pfilter({object: mypomp, params: currentParams, Np: 200, filterMean: true, predMean: true, maxFail: 3000})

console.log(new Date() - t, pf.loglik)

let createCsvWriter = require('csv-writer').createArrayCsvWriter;
let csvWriter = createCsvWriter({
  header: [],
  path: rootDir+'/samples/predmean.csv',
})
csvWriter.writeRecords( pf.predMean)
  
