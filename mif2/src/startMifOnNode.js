
const fs = require("fs");
const mathLib = require('../../library/mathLib.js');
const { startFilteringOnNode } = require('../../pfilter/src/startFilteringOnNode.js');
const { utility } = require('../../library/utility')


exports.startMifOnNode = function (paramSet, dataMif, dataPfilter, jobResume) {
  fs.readFile('./mif2/www/js/worker-bundle.js',"utf8", (err, workerBundle) => {
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
      console.log("Deploying job for MIF2...");
      job = compute.for(paramSet, workerFn, dataMif);
      job.public = { name: 'COVID-MIF2' };
    }

    job.on('console', (msg) => console.log("Got console event:", msg));
    job.on('uncaughtException', (error) => console.error(error));

    job.on('accepted', () => {
      console.log("Job accepted for MIF2", job.id);
      utility("COVID-MIF2", job.id)
    });

    job.on('status', (status) => {
      console.log("Status update in MIF2:", status);
    });

    
    const resultHandle = await job.exec(compute.marketValue(), bank);
    let mifResults = Array.from(resultHandle);
      mifResults.sort(function(a, b) {                /* Sort and choose the best point to calculate in pfilter */
        return b.loglik - a.loglik;
      });
      console.log("MIF2 on complete", mifResults[0]);
     
      startFilteringOnNode(mifResults[0],  dataPfilter);           
  }
}