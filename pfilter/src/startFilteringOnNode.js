const mathLib= require('../../library/mathLib.js');
const fs = require("fs");

exports.startFilteringOnNode = function (paramSet, data, runTime) {
  fs.readFile('./pfilter/www/js/worker-bundle.js',"utf8", (err, workerBundle) => {
    const workerFn = `(async (...args) => {
      ${workerBundle}
      return await self.workerfn(...args);
    })`;
    
    runComputeFor(workerFn);
  })
  
  const runComputeFor = async function(workerFn) {
    const compute = require('dcp/compute');
    const wallet = require('dcp/wallet');

    let repeatedParamSet = [];
    console.log("Deploying job...");
    for (let i = 0; i < paramSet.length; i++) {
      for (let j = 0; j < data[0].replicate; j++) {
        repeatedParamSet.push(paramSet[i]);
      }
    }
    
    let job = compute.for(repeatedParamSet, workerFn, data);
    job.public = { name: 'Pfilter' } 
    job.on('console', (msg) => console.log("Got console event:", msg));
    job.on('uncaughtException', (error) => console.error(error));

    job.on('accepted', () => {
      console.log("Job accepted");
    });

    job.on('status', (status) => {
      console.log("Got a status update:", status);
    });
    
    // job.on('result', function(res) {
    //   console.log("accBefore", res.result.predMean)
    // });

    job.on('complete', function(res) {
      console.log(" time", new Date() - runTime)
      let accumulatedResults = [];
      let results = res.values();
      let final;
      if(results) {
        for (let i = 0; i < results.length; i++) {
          final = {loglik: results[i].loglik};
          Object.assign(final, results[i].params);
          accumulatedResults.push(final);
        }
      }

      if(typeof accumulatedResults != 'undefined') {

        let iter = 1;
        let finalResult = [];
        let slice;
        for( let k = 0; k < accumulatedResults.length; k++) {
          slice = [];
          while(k < iter * data[0].replicate) {
            slice.push(accumulatedResults[k].loglik);
            k++;
            if (typeof accumulatedResults[k] === "undefined") break;
          }
          iter++;
          finalResult.push(Object.assign(accumulatedResults[k -1], {loglik: mathLib.logMeanExp(slice)}));
        }
        console.log("FinalArray",finalResult);
        let h = [];
        let keys = Object.keys(finalResult[0]);
        for ( let i = 0; i < keys.length; i++) {
          h.push({id:keys[i], title:keys[i]})
        }
        const createCsvWriter = require('csv-writer').createObjectCsvWriter;
        const csvWriter = createCsvWriter({
          path: './samples/finalResults.csv',
          header: h
        });

        csvWriter.writeRecords(finalResult)      
        .then(() => {
          console.log('...Done');
        });
      }
    });
    await job.exec(compute.marketValue() );
  }
}