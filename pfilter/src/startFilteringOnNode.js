const mathLib= require('../../library/mathLib.js');
const fs = require("fs");
const { utility } = require('../../library/utility.js')

exports.startFilteringOnNode = function (paramSet, data, jobResume) {
  
  fs.readFile('./pfilter/www/js/worker-bundle.js',"utf8", (err, workerBundle) => {
    const workerFn = `(async (...args) => {
      ${workerBundle}
      return await self.workerfn(...args);
    })`;
    
    runComputeFor(workerFn);
  })
  
  const runComputeFor = async function(workerFn) {
    const compute = require('dcp/compute');
    let bank = await require('dcp/wallet').get();

    let job;
    if (jobResume) {
      job = await compute.resume(jobResume.jobId, bank);
      console.warn(`Job ${jobResume.key} with jobId = ${job.id} is resuming!`);
    } else {
      let repeatedParamSet = [];   
      for (let j = 0; j < data[0].replicate; j++) {      /** If we need to replicate pfilter */
        repeatedParamSet.push(paramSet.params);
      }
      console.log("Deploying job for pfilter...");
      job = compute.for(repeatedParamSet, workerFn, data);
      job.public = { name: 'COVID-pfilter' };
    }
    
    job.on('console', (msg) => console.log("Got console event:", msg));
    job.on('uncaughtException', (error) => console.error(error));

    job.on('accepted', () => {
      console.log("Job accepted for pfilter", job.id);
      utility("COVID-pfilter", job.id);
    });

    job.on('status', (status) => {
      console.log("Status update in pfilter:", status);
    });
    
    // job.on('result', function(res) {
    //   console.log( res.result.loglik)
    // });

    //'complete'
    const resultHandle = await job.exec(compute.marketValue(), bank);
    let pf = Array.from(resultHandle);
    pf = pf[0];                                /* In this model there is only one point for filtering */
    let finalParams = [pf.params];
  
    let h = [];                                /* Log params */
    let keys = Object.keys(finalParams[0]);
    for ( let i = 0; i < keys.length; i++) {
      h.push({id:keys[i], title:keys[i]})
    }
    const createCsvWriter = require('csv-writer').createObjectCsvWriter;
    const csvWriter = createCsvWriter({
      path: './results/modelParams.csv',
      header: h
    });
    csvWriter.writeRecords(finalParams)      
    .then(() => {
      console.log('...Done');
    });

    let headerStates = [];
    for (let i = 0; i < Object.keys(pf.filterMean[0]).length; i++) {
      headerStates.push({id: Object.keys(pf.filterMean[0])[i], title: Object.keys(pf.filterMean[0])[i]})
    }
    
    const csvWriterFilter = createCsvWriter({   /* Log filterMean */
      path: './results/filterMean.csv',
      header: headerStates
    });
    csvWriterFilter.writeRecords(pf.filterMean)      
    .then(() => {
      console.log('...Done');
    });
    
    const csvStates = createCsvWriter({          /* Log the last saved states (at the last time step) */
      path: './results/savedStates.csv',
      header: headerStates
    });
    csvStates.writeRecords(pf.saveStates)        /* pf.saveStates is the last saved  states */
    .then(() => {
      console.log('...Done');
      fs.writeFileSync('./COVID.json', JSON.stringify([{}]));
    })
  }
}