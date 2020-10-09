/**
 *  @file             pfilter.js   
 *                    Particle filtering.    
 *                         
 *  @references       https://github.com/kingaa/pomp                  
 *
 *  @author           Nazila Akhavan, nazila@kingsds.network
 *  @date             June 2020
 */
const { initState } = require("./initState.js");
const { rprocessInternal } = require("../../library/rprocessInternal.js");
const { dmeasureInternal } = require("../../library/dmeasureInternal.js");
const { pfilter_computations } = require("../../library/pfilterComputations.js");
const snippet = require('../../library/modelSnippet.js');
const pomp = require('../../library/pomp.js');

/**
 * A plain vanilla sequential Monte Carlo (particle filter) algorithm.
 * Resampling is performed at each observation.
 * @param {object} paramSet            Input parameters.
 * @param {object} args 
 * @param {number} args.Np             The number of particles to use.
 * @param {number} args.tol            Particles with log likelihood below tol are considered to be "lost". A filtering failure
 *                                     occurs when, at some time point, all particles are lost. When all particles are lost,
 *                                     the conditional log likelihood at that time point is set to 0.
 * @param {number} args.maxFail        The maximum number of filtering failures allowed. If the number of filtering failures 
 *                                     exceeds this number, execution will terminate with an error.
 * @param {boolean} args.predMean      If TRUE, the prediction means are calculated for the state variables and parameters.
 * @param {boolean} args.predVar       If TRUE, the prediction variances are calculated for the state variables and parameters.
 * @param {boolean} args.filterMean    If TRUE, the filtering means are calculated for the state variables and parameters.
 * @param {boolean} args.filterTraj    Not translated.
 * @param {boolean} args.saveStates    The state-vector for each particle at each time is saved.
 * @param {boolean} args.saveParams    The params-vector for each particle at each time is saved.
 * @param {object} args.object         An object of class POMP.
 * 
 * @returns {object}
 *  @param {number} logLik         The estimated log likelihood.
 *  @param {array} condLogLik      An array of estimated conditional log likelihoods at each time point.
 *  @param {array} effSampleSize   An array containing the effective number of particles at each time point.
 *  @param {array} predMean       The mean of the approximate prediction distribution.
 *  @param {array} predVar        The variance of the approximate prediction distribution.
 *  @param {array} filterMean     The mean of the filtering distribution.
 *  @param {boolean} filterTraj   Returns false. Not translated.
 *  @param {array} savedStates    Retrieve list of saved states.
 *  @param {array} savedParams    Retrieve list of saved params.
 *  @param {number} nfail         The number of filtering failures encountered.
 *    
 */
exports.pfilter = function (paramSet, args) {
   
  let params = paramSet.params? paramSet.params : paramSet;
  const pompData = Object.assign(args.object,{
    rprocess: snippet.rprocess,
    rmeasure: snippet.rmeasure,
    dmeasure: snippet.dmeasure,
    initializer: snippet.initializer,
    params: params
  });
  
  const pompObject = new pomp(pompData);
  
  let inputParamSet = Object.assign({}, params);
  let Np = args.Np;
  let tol = args.tol ? args.tol : 1e-17;
  let maxFail = args.maxFail ? args.maxFail :  Infinity;
  let predMean = args.predMean ? args.predMean :  false;
  let predVar = args.predVar ? args.predVar :  false;
  let filterMean = args.filterMean ? args.filterMean :  false;
  let filterTraj = args.filterTraj ? args.filterTraj :  false;
  let saveStates = args.saveStates ? args.saveStates : false;
  let saveParams = args.saveParams ? args.saveParams : false;
  
  if (Object.keys(params).length === 0) {
    throw new Error("In pfilterInternal: params must be specified");
  } 
  let onePar = false;
  
  let times = [pompObject.t0, ...pompObject.times];
  let ntimes = times.length - 1;
  
  if (typeof Np === "function" || Np === undefined || Np <= 0) {
    throw new Error(`Number of particles should be a positive number. ${Np} is not translated`);
  }
  if (!Array.isArray(params) || params.every(x => !Array.isArray(x))) { //there is only one parameter vector
    onePar = true;
    params = [params]; //as.matrix(params)
  }
  
  let x = initState(pompObject, params, Np);
  let xparticles, pparticles, pedigree
  // set up storage for saving samples from filtering distributions
  if (saveStates || filterTraj) {
    xparticles =  new Array(ntimes);
  }
  
  if (saveParams) {
    pparticles = new Array(ntimes);;
  } else {
    pparticles = [];
  }
  if (filterTraj) {
    pedigree = new Array(ntimes + 1);
  }

  let loglik = new Array(ntimes);
  let effSampleSize = new Array(ntimes).fill(0);
  let nfail = 0;
  
  // set up storage for prediction means, variances, etc.
  let predm;
  if (predMean) {
    predm = new Array(ntimes).fill(null);
  } else {
    predm = [];
  }

  let predv;
  if (predVar) {
    predv = new Array(ntimes).fill(null);
  } else {
    predv = [];
  }
  
  let filtm;
  if (filterMean) {
    filtm  = new Array(ntimes).fill(null);
  } else {
    filtm = [];
  }

  if (filterTraj) {
    throw new Error ('filterTraj is not implemented')
  }
  
  // main time loop
  let X;
  for (let nt = 0; nt < ntimes; nt++) {
    if (typeof progress === 'function') progress();
    // try {
      X = rprocessInternal(pompObject, x, [times[nt],times[nt + 1]], params, 1)
    // } catch (e) {
    //   throw new Error(`In pfilterInternal: Process simulation error: ${e}`)
    // }
  
    if (predVar) { // check for nonfinite state variables and parameters
      allFinite = X.every(a => a.every(x =>isFinite(x)))
      if (!allFinite) {  // state variables
        throw new Error("In pfilter.js: non-finite state variable(s): ");
      }
    }

    let weights = [];
    try {
      weights = dmeasureInternal(
        pompObject,
        y = pompObject.data[nt],
        X,
        times[nt + 1],
        params,
        log = false
      ); 
    } catch (error) {
      throw new Error(`In pfilter.js: error in calculation of weights: ${error}`);
    }
    
    let allFinite = weights.map(w => isFinite(w)).reduce((a, b) => a & b, 1);
    if (!allFinite) {
      throw new Error("In dmeasure: weights returns non-finite value");
    }
    
    /** compute prediction mean, prediction variance, filtering mean,
    * effective sample size, log-likelihood
    * also do resampling if filtering has not failed
    */ 
    let xx;
    try {
      xx  = pfilter_computations(
        X,
        params,
        Np = Np,
        0, //rw_sd
        predMean,
        predVar ,
        filterMean ,
        trackancestry = filterTraj,
        onePar,
        weights = weights,
        tol = tol
      );
    } catch (error) {
      throw new Error(`particle-filter error: ${error}`);
    } 

    let allFail = xx.fail;
    loglik[nt] = xx.loglik;
    effSampleSize[nt] = xx.ess;
    x = xx.states;
    params = xx.params;

    if (predMean) {
      predm[nt] = xx.fm;
    }

    if (predVar)
      predv[nt] = xx.pv;

    if (filterMean)
      filtm[nt] = xx.fm;

    if (filterTraj)
      pedigree[[nt]] = xx.ancestry;

    if (allFail) { // all particles are lost
      nfail = nfail + 1;
      if (nfail > maxFail)
        throw new Error("In mif2Pfilter: too many filtering failures")
    }

    if (saveStates || filterTraj) {
      xparticles[nt] = x;
    }

    if (saveParams) {
      pparticles[nt] = params;
    } 
  } //end of main loop 

  if (filterTraj) { // select a single trajectory
    throw new Error("filterTraj is not translated")
  }

  if (!saveStates) xparticles = [];

  if (nfail > 0) {
    console.log("warning! filtering failure occurred.");
  }
  
  return {
    params: inputParamSet,
    predMean: predm,
    predVar: predv,
    filterMean: filtm,
    filterTraj: false,//filt.t
    paramMatrix: [[]],
    effSamplesize: effSampleSize,
    condLoglik: loglik,
    savedStates: xparticles,
    savedParams: pparticles,
    Np: Np,
    tol: tol,
    nfail: nfail,
    loglik: loglik.reduce((a,b) => a + b, 0)
  }
}




