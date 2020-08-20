const snippet = require("../../library/modelSnippet");
const { startFilteringOnWeb } = require('../../pfilter/src/startFilteringOnWeb.js');

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
  let dataCases = [];
  let dataCasesTimes = [];
  let dataCovar = [];
  let dataCovarTimes = [];
  let currentParams = [];
  let dataCovar_name;
  let dataCases_name;
  let temp, file;
  /**
   * Note: We pass covar data as an array of array and in the interpolator we match them by their name
   * (they are in the same order) and there it returns an object.
   * */
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

  const runComputeFor = async function( paramSet, dataMif, dataPfilter) {
    console.log("Deploying job...");
    let t = new Date();
    let job = dcp.compute.for(paramSet, workerFn, ...dataMif);
    job.public = { name: 'Mif2' } 
    job.on('console', (msg) => console.log("Got console event:", msg));
    job.on('uncaughtException', (error) => console.error(error));

    job.on('accepted', () => {
      console.log("Job accepted");
    });

    job.on('status', (status) => {
      console.log("Got a status update:", status);
    });
   
    job.on('result', (ev) => console.log("Got a result:", ev));
    job.on('complete', async function() {
      await job.results.fetch();
      let result = job.results.values();
      console.log("On complete", result);// array of objects
     
      startFilteringOnWeb(result, dataPfilter, t)
    });
    await job.exec(dcp.compute.marketValue);
  }  

  let runButton = document.getElementById('go-button');
  runButton.onclick = function () {
    
    currentParams = [{R0: 30.34720726, amplitude: 0.480957796, gamma: 73.05, mu: 0.00743129, sigma: 45.66,psi: 0.393478325, rho: 0.46174302, S_0: 0.886372769, E_0: 0.003366557, I_0: 0.00000347, R_0: 0.110257201},
      {R0: 30.34720726, amplitude: 0.480957796, gamma: 73.05, mu: 0.00743129, sigma: 45.66,psi: 0.393478325, rho: 0.46174302, S_0: 0.886372769, E_0: 0.003366557, I_0: 0.00000347, R_0: 0.110257201}];
    let params_ic_fit = [];
    let params_mod_fit = ["R0", "amplitude", "mu", "rho", "psi"];
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
      ivps: params_ic_fit,
      pars: params_mod_fit,
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
      replicate: 5
    }];
    runComputeFor(currentParams, dataMif, dataPfilter)
  };
}