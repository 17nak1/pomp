let mathLib= require('../../library/mathLib.js');

exports.startFilteringOnWeb = function (paramSet, data, t) {
  fetch('/js/worker-bundle-pfilter.js')
  .then(response => response.text())
  .then((workerBundle) => {
    // self.workerfn gets set by the workerBundle (see worker-src/index.js)
    const workerFn = `(async (...args) => {
      ${workerBundle}
      return await self.workerfn(...args);
    })`;
  
    runComputeFor(workerFn);
  })
  
  const runComputeFor = async function(workerFn) {
    let repeatedParamSet = [];
    console.log("Deploying job...");
    for (let i = 0; i < paramSet.length; i++) {
      for (let j = 0; j < data[0].replicate; j++) {
        repeatedParamSet.push(paramSet[i]);
      }
    }
    console.log(repeatedParamSet.length)
    // console.log( "re", repeatedParamSet.length, repeatedParamSet)
    let job = dcp.compute.for(repeatedParamSet, workerFn, data);
    job.public = { name: 'Pfilter' } 
    job.on('console', (msg) => console.log("Got console event:", msg));
    job.on('uncaughtException', (error) => console.error(error));

    job.on('accepted', () => {
      console.log("Job accepted");
    });

    // job.on('status', (status) => {
    //   console.log("Got a status update:", status);
    // });
    
    // job.on('result', function(res) {
    //   console.log("accBefore", accumulatedResults)
    //   // if(typeof res !== 'undefined') {
    //   //   let results = res.result;
    //   //   if(results) {
    //   //     // console.log("acc", res.sliceNumber)
    //   //     let final = {loglik: results.loglik, sliceNumber: res.sliceNumber};
    //   //     Object.assign(final, results.params);
    //   //     console.log("resFinal",final)
    //   //     accumulatedResults.push(final);
    //   //     console.log("resAcc", accumulatedResults)
    //   //   }
    //   // }
    // });

    job.on('complete', function(res) {
      console.log(" time", new Date() - t)
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
          finalResult.push([Object.assign(accumulatedResults[k -1], {loglik: mathLib.logMeanExp(slice)})]);
        }
        console.log("FinalArray",finalResult);
      }
    });
    await job.exec(dcp.compute.marketValue);
  }
  
}
