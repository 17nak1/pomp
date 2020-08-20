/**
 * Example of using pfilter.js
 *
 */

const snippet = require("../../library/modelSnippet");
let mathLib= require('../../library/mathLib.js');

// Loading the page
document.addEventListener('DOMContentLoaded', () => {
  fetch('js/worker-bundle.js')
  .then(response => response.text())
  .then((workerBundle) => {
    // self.workerfn gets set by the workerBundle (see worker-src/index.js)
    const workerFn = `(async (...args) => {
      ${workerBundle}
      return await self.workerfn(...args);
    })`;
    start(workerFn);
  });
})

function start (workerFn) {
  let accumulatedResults = [];
  let dataCases = [];
  let dataCasesTimes = [];
  let dataCovar = [];
  let dataCovarTimes = [];
  let currentParams = [];
  let dataCovar_name;
  let dataCases_name;
  let temp, file;

  let repeat = 2;
  document.getElementById('file1-upload').onchange = function () {
    file = this.files[0]
    dataCovar = [];
    dataCovarTimes = [];
    dataCovar_name = [];
    let reader = new FileReader();
    reader.onload = function () {
      let lines = this.result.split(/\r\n|\n/);
      dataCovar_name = lines[0].replace(/['"]+/g, '').split(',');
      dataCovar_name.shift();
      for (let i = 1; i < lines.length; i++) {
        temp = lines[i].split(',');
        if(temp.length > 1) {
          temp = temp.map(x => Number(x));
          dataCovarTimes.push(temp[0]);
          dataCovar.push(temp.slice(1));
        }
      }
    }
    reader.readAsText(file)
  }
 
  document.getElementById('file2-upload').onchange = function () {
    file = this.files[0];
    dataCases = [];
    dataCasesTimes = [];
    dataCases_name = [];
    let data;
    let reader = new FileReader();
    reader.onload = function () {
      let lines = this.result.split(/\r\n|\n/);
      dataCases_name = lines[0].replace(/['"]+/g, '').split(',');
      dataCases_name.shift();
      for (let i = 1; i < lines.length; i++) {
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
    }
    reader.readAsText(file)
  }

  document.getElementById('file3-upload').onchange = function () {
    file = this.files[0];
    currentParams = [];
    let currentParams_name = [];
    let data;
    let reader = new FileReader();
    reader.onload = function () {
      let lines = this.result.split(/\r\n|\n/);
      currentParams_name = lines[0].replace(/['"]+/g, '').split(',');
      for (let i = 1; i < lines.length; i++) {
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
    }
    
    reader.readAsText(file)
  }

  const runComputeFor = async function( paramSet, ...args) {
    console.log("Deploying job...");
    
    let job = dcp.compute.for(paramSet, workerFn, args);
    job.public = { name: 'Pfilter' } 
    job.on('console', (msg) => console.log("Got console event:", msg));
    job.on('uncaughtException', (error) => console.error(error));

    job.on('accepted', () => {
      console.log("Job accepted");
    });

    job.on('status', (status) => {
      console.log("Got a status update:", status);
    });
  
    job.on('result', function(res) {
      
      if(typeof res !== 'undefined') {
        let results = res.result;
        if(results) {
          let final = {loglik: results.loglik, sliceNumber: res.sliceNumber};
          Object.assign(final, results.params);
          accumulatedResults.push(final);
         }
      }
    });

    
    job.on('complete', function() {
      if(typeof accumulatedResults != 'undefined') {
        accumulatedResults.sort(mathLib.sortObjects("sliceNumber"));
        let iter = 1;
        let finalResult = [];
        let slice;
        for( let k = 0; k < accumulatedResults.length; k++) {
          slice = [];
          while(accumulatedResults[k].sliceNumber < iter * repeat) {
            slice.push(accumulatedResults[k].loglik);
            k++;
            if (typeof accumulatedResults[k] === "undefined") break;
          }
          iter++;
          delete accumulatedResults[k -1].sliceNumber;
          finalResult.push([Object.assign(accumulatedResults[k -1], {loglik: mathLib.logMeanExp(slice)})]);
        }
        console.log("FinalArray",finalResult);

      }
    });

    await job.exec(dcp.compute.marketValue);
  }

  let runButton = document.getElementById('go-button');
  runButton.onclick = function () {

    currentParams = [
      {R0: 30.34720726, amplitude: 0.480957796, gamma: 73.05, mu: 0.00743129, sigma: 45.66,psi: 0.393478325, rho: 0.46174302, S_0: 0.886372769, E_0: 0.003366557, I_0: 0.00000347, R_0: 0.110257201},
      {R0: 30.34720726, amplitude: 0.480957796, gamma: 73.05, mu: 0.00743129, sigma: 45.66,psi: 0.393478325, rho: 0.46174302, S_0: 0.886372769, E_0: 0.003366557, I_0: 0.00000347, R_0: 0.110257201}]
    
    const pompData = {
      data :  dataCases,
      times:  dataCasesTimes,
      t0: 1940,
      rprocessDetail:  { type:"euler_sim", deltaT: 1/365.25 },
      covar: dataCovar,
      tcovar: dataCovarTimes,
      zeronames: snippet.zeronames,
      statenames: snippet.statenames,
      paramnames: snippet.paramnames,
      covarnames: dataCovar_name,
      obsnames: dataCases_name,
    };
    
    let data = {
      object: pompData,
      Np: 200,
      filterMean: true,
      predMean: true,
      maxFail: 3000
    };
    runComputeFor(currentParams, data)
  };
}







  
