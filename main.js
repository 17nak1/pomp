/**
 * Adapting with the new model the following will be affected:
 * - trajMatch/generateSobol
 * - library/modelSnippet/
 * - library/modelSnippet(determineRW)
 * - RuncomputeRorTM( noise parameters)
 * - epi-fit (paramsFixed, globals, data, covar)
 */
const fs = require("fs");
const csv = require('csvtojson');
const process = require("process");
const snippet = require("./library/modelSnippet");
const { sobolSet } = require('./trajMatch/src/generateSobol');
const { runComputeForTM } = require('./trajMatch/src/runComputeForTM.js');
const { startMifOnNode } = require("./mif2/src/startMifOnNode");
const { startFilteringOnNode } = require('./pfilter/src/startFilteringOnNode.js');


/** Main program entry point */
async function start(workerFn, options) {
  if (typeof options.pathToData === undefined) throw new Error("Path to data is undefined!")
  if (typeof options.pathToCovar === undefined) throw new Error("Path to covar is undefined!")
  if (!options.pop) throw new Error("Population is undefined!")
  if (!options.TC) throw new Error("accurate episode date of first case!")
  let pop = options.pop;
  let TC = options.TC;
  let T0 = options.T0 ? options.T0 : 75;
  let T1 = options.T1 ? options.T1 : 139;
  let T2 = options.T2 ? options.T2 : 212;
  
  await require('dcp-client').init(process.argv);
  const compute = require('dcp/compute');
  const dcpCli = require('dcp/dcp-cli');

  dcpCli.base().argv;

  let dataCovar = [];
  let dataCovarTimes = [];
  let covarHeader = [];
  let data;
  data = await csv().fromFile(pathToCovar);
  if(!data) throw new Error( "Cannot read the covar set");
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
  
  data = await csv().fromFile(pathToData);
  if(!data) throw new Error( "Cannot read the data set");
  dataHeader = Object.keys(data[0]);
  dataHeader.shift();
  for (let i = 0; i < data.length; i++) {
    dataCasesTimes.push(Number(data[i].time));
    delete data[i].time;
  }
  dataCases = data;
  
  /* Determine fixed and fitted parameters in trajMatch */
  paramsFixed = ["beta_sd","dB0", "dB1", "dB2", "sigma","kappa","gammaI","iota"];
  let paramsFit = snippet.paramsMod;
  for (let i = 0; i < paramsFixed.length; i++) {
    paramsFit = paramsFit.filter(e => e !== paramsFixed[i])
  }
  /**TC is the accurate episode date of first case */
  const globals = { nstageE: 3, nstageP: 3, nstageI: 3, nstageH: 3, nstageC: 3, nstageV: 3, pop: pop, T0: T0, T1: T1, T2: T2, TC: TC};
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
  let SobolNumberOfPoints = 1000//1e6     
  let dataTrajMatch = [{
    object: pompData, 
    est: paramsFit,
    transform: true,
    method: "subplex"}
  ];

  let dataMif = [{
    object: pompData,
    Nmif: 50,
    transform: true,
    rw_sd: snippet.determineRW(),
    Np: 8000,
    varFactor: 2,
    coolingType: "hyperbolic",
    coolingFraction: 0.05
  }];

  let dataPfilter = [{
    object: pompData,
    Np: 60e3,
    filterMean: true,
    saveStates: true,
    maxFail: 3000,
    replicate: 1
  }];

  
  //START DEPLOYING JOBS
  let bank = await require('dcp/wallet').getId('default');
  let jsonString;
  try {
    jsonString = fs.readFileSync('./COVID.json');
  } catch (error) {
    jsonString = JSON.stringify([{}]);
  }
  const jobResume = JSON.parse(jsonString);
  let isResumingPfilter = false;
  let jobResumePfilter;
  let isResumingMIF = false;
  let jobResumeMIF;
  let isResumingTM = false;
  let jobResumeTM;
  for (let i = 0; i < jobResume.length; i++) {
    if (jobResume[i].key === "COVID-pfilter") {
      isResumingPfilter = !!jobResume[i].jobId;
      jobResumePfilter = jobResume[i];
    } else if (jobResume[i].key === "COVID-MIF2") {
      isResumingMIF = !!jobResume[i].jobId;
      jobResumeMIF = jobResume[i];
    } else if (jobResume[i].key === "COVID-TM") {
      isResumingTM = !!jobResume[i].jobId;
      jobResumeTM = jobResume[i];
    }
  }
  let job;
  if (isResumingPfilter) {
    startFilteringOnNode([],  dataPfilter, jobResumePfilter);
  } else {
    if (isResumingMIF) {
      startMifOnNode([], dataMif, dataPfilter, jobResumeMIF);
    } else {
      if (isResumingTM) {
        job = await compute.resume(jobResumeTM.jobId, bank);
        console.warn(`Job ${jobResumeTM.key} with jobId = ${job.id} is resuming!`);
      } else {
        console.log("Deploying job for TM...");
        job = compute.for(sobolSet(SobolNumberOfPoints, globals), workerFn, dataTrajMatch);
        job.public = { name: 'COVID-TM' };
      }
      runComputeForTM(job, dataMif, dataPfilter, SobolNumberOfPoints);
    } 
  }  
}

const main  = function(options) {
  fs.readFile('./trajMatch/www/js/worker-bundle.js',"utf8", (err, workerBundle) => {
    const workerFn = `(async (...args) => {
      ${workerBundle}
      return await self.workerfn(...args);
    })`;
    start(workerFn, options);
  })
}

module.exports = main

// node epi-fit.js  --scheduler="http://scheduler.nazila.office.kingsds.network"  --default-bank-account-file=~/.dcp/default.nazila.keystore --identity-file=~/.dcp/id.nazila.keystore
