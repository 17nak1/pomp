
const fs = require("fs");
const csv = require('csvtojson')

const process = require("process");
const snippet = require("./library/modelSnippet");
const { startFilteringOnNode } = require('./pfilter/src/startFilteringOnNode.js');

/** Main program entry point */
async function start(workerFn) {
 
  require('dcp-client').initSync(process.argv);
  const compute = require('dcp/compute');
  // const wallet = require('dcp/wallet');
  // const dcpCli = require('dcp/dcp-cli');
  // const identityKeystore = await dcpCli.getIdentityKeystore();
  // wallet.addId(identityKeystore);
  // const accountKeystore = await dcpCli.getAccountKeystore();

  let dataCovar = [];
  let dataCovarTimes = [];
  let dataCovar_name = [];
  let data;
  data = await csv().fromFile('./samples/London_covar.csv');
  if(!data) throw new Error( "Cannot read the first data set");
  dataCovar_name = Object.keys(data[0]);
  dataCovar_name.shift();
  for (let i = 0; i < data.length; i++) {
    dataCovarTimes.push(Number(data[i].time));
    delete data[i].time;
    temp = Object.values(data[i]).map(x => Number(x));
    dataCovar.push(temp);
  }
  
  let dataCases = [];
  let dataCasesTimes = [];
  let dataCases_name = [];
  
  data = await csv().fromFile('./samples/London_BiDataMain.csv');
  if(!data) throw new Error( "Cannot read the second data set");
  dataCases_name = Object.keys(data[0]);
  dataCases_name.shift();
  for (let i = 0; i < data.length; i++) {
    dataCasesTimes.push(Number(data[i].time));
    delete data[i].time;
    Object.keys(data[i]).forEach(function(key) {
      data[i][key] = Number(data[i][key]);
    });
  }
  dataCases = data;

  let currentParams = [];
  let currentParams_name = [];
  
  data = await csv().fromFile('./samples/initial_parameters.csv');
  if(!data) throw new Error( "Cannot read the second data set");
  currentParams_name = Object.keys(data[0]);
  for (let i = 0; i < data.length; i++) {
    Object.keys(data[i]).forEach(function(key) {
      data[i][key] = Number(data[i][key]);
    });
  }
  currentParams = data;
  
  const runComputeFor = async function(paramSet, dataMif, dataPfilter) {
    console.log("Deploying job...");
    
    let timeRun = new Date();
    let job = compute.for(paramSet, workerFn, ...dataMif);
    job.public = { name: 'Mif2' } 
    job.on('console', (msg) => console.log("Got console event:", msg));
    job.on('uncaughtException', (error) => console.error(error));

    job.on('accepted', () => {
      console.log("Job accepted");
    });

    job.on('status', (status) => {
      console.log("Got a status update:", status);
    });
   
    // job.on('result', (ev) => console.log("Got a result:", ev));
    job.on('complete', async function() {
      await job.results.fetch();
      let result = job.results.values();// array of objects
      console.log("On complete", result);
     
      startFilteringOnNode(result, dataPfilter, timeRun)
    });
    
    await job.exec(compute.marketValue());
  }
  

  currentParams = [{R0: 30.34720726, amplitude: 0.480957796, gamma: 73.05, mu: 0.00743129, sigma: 45.66,psi: 0.393478325, rho: 0.46174302, S_0: 0.886372769, E_0: 0.003366557, I_0: 0.00000347, R_0: 0.110257201}];
  
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
  
  
  if (typeof mf === "undefined") mf = {};
  let dataMif = [[{
    object: pompData,
    Nmif: 2,
    transform: true,
    rw_sd: snippet.determineRW("R0"),
    Np: 500,
    varFactor: 2,
    coolingType: "hyperbolic",
    coolingFraction: cool_fraction,
    snapShotStart: mf.snapShotStart,
    paramMatrix: mf.paramMatrix,
    convRec: mf.convRec,
    _indices: mf._indices,
  }]];

  let dataPfilter = [{
    object: pompData,
    Np: 5000,
    filterMean: true,
    predMean: true,
    maxFail: 3000,
    replicate: 2
  }];
  runComputeFor(currentParams, dataMif, dataPfilter);
}

const main  = function() {
  fs.readFile('./mif2/www/js/worker-bundle.js',"utf8", (err, workerBundle) => {
    const workerFn = `(async (...args) => {
      ${workerBundle}
      return await self.workerfn(...args);
    })`;
    start(workerFn);
  })
}
main();

// node index_node.js  --scheduler="http://scheduler.nazila.office.kingsds.network"  --default-bank-account-file=~/.dcp/default.keystore.nazila --identity-file=~/.dcp/id.keystore.nazila

