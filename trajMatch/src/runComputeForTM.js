const { startMifOnNode } = require('../../mif2/src/startMifOnNode.js');
const fs = require("fs");
const { utility } = require("../../library/utility.js");

exports.runComputeForTM = async function(job, dataMif, dataPfilter, SobolNumberOfPoints) {
  const compute = require('dcp/compute');
  let bank = await require('dcp/wallet').get();
  const FETCH_BATCH_SIZE = 1000;

  job.on('console', (msg) => {
    console.log("Got console event:", msg);
  });
  job.on('uncaughtException', (error) => console.error(error));

  job.on('accepted', () => {
    console.log("Job accepted for TM", job.id);
    utility("COVID-TM", job.id)
  });

  job.on('status', (status) => {
    console.log("Status update in TM:", status);
  });
 
  /* Sort the results based on loglik */   
  job.collateResults =false;
  const resultHandle = await job.exec(compute.marketValue(), bank);
  let trajMatchResults = [];
  for (let i = 0; i * FETCH_BATCH_SIZE < SobolNumberOfPoints; i++) {
    await resultHandle.fetch({                                     /* Only fetch limited tasks from the scheduler */
      start: i * FETCH_BATCH_SIZE,
      end: Math.min((i + 1) * FETCH_BATCH_SIZE - 1, SobolNumberOfPoints - 1)
    });
    let result = Array.from(resultHandle);
    resultHandle.reset();                                         /* Reset fetched tasks from the scheduler */
    for (let j = 0; j < result.length; j++) {
      if(result[j].loglik < -0) {
        trajMatchResults.push(result[j]);
      }
    }
    trajMatchResults.sort(function(a, b) {
      return b.loglik - a.loglik;
    });

    trajMatchResults = trajMatchResults.slice(0,101);              /* Only save the best 100 points*/
  }

  for (let i = 0; i < trajMatchResults.length; i++) {              /* Add random number to the noise parameters */
    for(key in trajMatchResults[0].params) {
      if(key === 'beta_sd' || key === 'dB0' || key === 'dB1' || key === 'dB2') {
        trajMatchResults[i].params[key] = Math.random();
      } 
    }
  }
  console.log("Best TrajMatch set", trajMatchResults[0])
  startMifOnNode(trajMatchResults, dataMif, dataPfilter);
}