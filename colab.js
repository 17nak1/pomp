
const fs = require("fs");
const csv = require('csvtojson')

const process = require("process");
const snippet = require("./library/modelSnippet");


// var { sobolSet } = require('/content/gdrive/My Drive/pomp/library/generate-sobol/generateSobol.js');
// /** Main program entry point */
// async function start(workerFn) {
 
  require('dcp-client').initSync(process.argv);
  const compute = require('dcp/compute');
  const wallet = require('dcp/wallet');
  const dcpCli = require('dcp/dcp-cli');
  

//   let dataCovar = [];
//   let dataCovarTimes = [];
//   let dataCovar_name = [];
//   let data;
//   data = await csv().fromFile('./samples/London_covar.csv');
//   if(!data) throw new Error( "Cannot read the first data set");
//   dataCovar_name = Object.keys(data[0]);
//   dataCovar_name.shift();
//   for (let i = 0; i < data.length; i++) {
//     dataCovarTimes.push(Number(data[i].time));
//     delete data[i].time;
//     temp = Object.values(data[i]).map(x => Number(x));
//     dataCovar.push(temp);
//   }
  
//   let dataCases = [];
//   let dataCasesTimes = [];
//   let dataCases_name = [];
  
//   data = await csv().fromFile('./samples/London_BiDataMain.csv');
//   if(!data) throw new Error( "Cannot read the second data set");
//   dataCases_name = Object.keys(data[0]);
//   dataCases_name.shift();
//   for (let i = 0; i < data.length; i++) {
//     dataCasesTimes.push(Number(data[i].time));
//     delete data[i].time;
//     Object.keys(data[i]).forEach(function(key) {
//       data[i][key] = Number(data[i][key]);
//     });
//   }
//   dataCases = data;

//   let currentParams = [];
//   let currentParams_name = [];
  
//   data = await csv().fromFile('./samples/initial_parameters.csv');
//   if(!data) throw new Error( "Cannot read the second data set");
//   currentParams_name = Object.keys(data[0]);
//   for (let i = 0; i < data.length; i++) {
//     Object.keys(data[i]).forEach(function(key) {
//       data[i][key] = Number(data[i][key]);
//     });
//   }
//   currentParams = data;
  
/** Main program entry point */
async function start(workerFn) {
  const identityKeystore = await dcpCli.getIdentityKeystore();
  wallet.addId(identityKeystore);
  const accountKeystore = await dcpCli.getAccountKeystore();
  let dataCovar = [];
  let dataCovarTimes = [];
  let covarHeader = [];
  let data;
   
  data = await csv().fromFile('./samples/London_covar.csv');
  if(!data) throw new Error( "Cannot read the first data set");
  covarHeader = Object.keys(data[0]);
  covarHeader.shift();
  for (let i = 0; i < data.length; i++) {
    dataCovarTimes.push(Number(data[i].time));
    delete data[i].time;
  }
  dataCovar = data;
  
  let dataCases = [];
  let dataCasesTimes = [];
  let dataHeader = [];
  
  data = await csv().fromFile('./samples/London_BiDataMain.csv');
  if(!data) throw new Error( "Cannot read the second data set");
  dataHeader = Object.keys(data[0]);
  dataHeader.shift();
  for (let i = 0; i < data.length; i++) {
    dataCasesTimes.push(Number(data[i].time));
    delete data[i].time;
  }
  dataCases = data;
  
  let currentParams = [];
  
  data = await csv().fromFile('./samples/initial_parameters.csv');
  if(!data) throw new Error( "Cannot read the second data set");
  // currentParams_name = Object.keys(data[0]);
  for (let i = 0; i < data.length; i++) {
    Object.keys(data[i]).forEach(function(key) {
      data[i][key] = Number(data[i][key]);
    });
  }
  currentParams = data;

  globals = { nstageE: 3, nstageP: 3, nstageI: 3, nstageH: 3, nstageC: 3, nstageV: 3, pop: 10e6, T0: 75, T1: 139 };
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
  

  let dataPfilter = [{
    object: pompData,
    Np: 100,
    filterMean: true,
    saveStates: true,
    maxFail: 3000,
    replicate: 1
  }];
 //START DEPLOYING JOBS
  console.log("Deploying job for pfilter...");
  job = compute.for([currentParams[0]], workerFn, dataPfilter);
  
  job.public = { name: 'test-pfilter' };

  job.on('accepted', () => {
    console.log("Job accepted for pfilter", job.id);
  });

  job.on('status', (status) => {
    console.log("Status update in pfilter:", status);
  });

  job.on('result', function(res) {
      console.log( [res.result.params,res.result.loglik]);
    });
  await job.exec(compute.marketValue, accountKeystore);
  
};


async function run(){
    let workerBundle = fs.readFileSync('./pfilter/www/js/worker-bundle.js', 'utf8');
    const workerFn = `(async (...args) => {
      ${workerBundle}
      return await self.workerfn(...args);
    })`;

    await start(workerFn);
}

run();



