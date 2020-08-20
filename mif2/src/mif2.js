/**
 *  @file             mif2.js   
 *                    Maximum likelihood by Iterated Filtering.    
 *                         
 *  @references       https://github.com/kingaa/pomp                  
 *
 *  @author           Nazila Akhavan, nazila@kingsds.network
 *  @date             June 2020
 */
const { cooling, partrans, coef } = require("./mif2Helpers.js");
const { mif2Pfilter } = require("./mif2Pfilter.js");
const snippet = require("../../library/modelSnippetCOVID3.js");
const pomp = require('../../library/pomp.js');
/**
/**
 * An iterated filtering algorithm for estimating the parameters of a partially-observed Markov process.
 * Running mif2 causes the algorithm to perform a specified number of particle-filter iterations.
 * At each iteration, the particle filter is performed on a perturbed version of the model, in which the
 *  parameters to be estimated are subjected to random
 *  perturbations at each observation.
 * This extra variability effectively smooths the likelihood surface and combats particle depletion by 
 * introducing diversity into particle population.
 * As the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
 * The algorithm is presented and justified in Ionides et al. (2015).
 * 
 * @param {array} params                Input parameters.
 * @param {object} args 
 * @param {object} args.object          An object of class POMP.
 * @param {number} args.Nmif            Number of of filtering iterations.
 * @param {number} snapShotStart        The current number of Nmif; it is < Nmif and defined to break the loop.
 * @param {array} convRec               The estimated parameters and log likelihood in each mif iteration.
 * @param {boolean} args.transform      If parmeters need to transform.
 * @param {object} args.rw_sd_param     Parametrs with no random walk to be used in rw_sd function that Specifies
 *                                      the magnitude of the random-walk perturbations that will be applied to some or all model parameters.
 * @param {number} args.Np              The number of particles to use.
 * @param {number} args.coolingFraction Cooling Rate
 * @param {string} args.coolingType     specifications for the cooling schedule,
 *                                       i.e., the manner with which the intensity of the parameter perturbations
 *                                       is reduced with successive filtering iterations.
 * @param {number} args.tol             Tolerance
 * @param {number} args.maxFail         Maximum number of allowed fail.
 * @param {array} args._paramMatrix     An array of Np arrays of args.start.
 * @param {number} args._ndone
 * @param {number} args._indices
 * 
 * @returns {object}
 *  @param {number} logLik         The estimated log likelihood.
 *  @param {array} cond.logLik     The estimated conditional log likelihood.
 *  @param {array} convRec         The estimated parameters and log likelihood in each mif iteration.
 *  @param {array} paramMatrix     An array of Np arrays of final parameters.
 *  @param {number} effSamplesize  Effective sample size.
 *  @param {number} indices         
 *  @param {number} Np             The number of particles.
 *  @param {number} tol            Tolerance 
 *  @param {number} nfail          Number of times that pfilter failed. 
 */
exports.mif2 = function (params, args) {
  const pompData = Object.assign(args.object,{
    rprocess: snippet.rprocess,
    rmeasure: snippet.rmeasure,
    dmeasure: snippet.dmeasure,
    initializer: snippet.initializer,
    toEstimationScale: snippet.toEstimationScale,
    fromEstimationScale: snippet.fromEstimationScale,
  });

  const pompObject = new pomp(pompData);
  let Nmif = args.Nmif;
  let snapShotStart = args.snapShotStart;
  let convRec = args.convRec;
  let transform = args.transform? args.transform: FALSE;
  let rw_sd = eval(args.rw_sd);
  let Np = args.Np;
  let coolingType = args.coolingType ? args.coolingType: ["hyperbolic", "geometric"];
  let coolingFraction = args.coolingFraction;
  let tol = args.tol? args.tol: 1e-17;
  let maxFail = args.maxFail? args.maxFail: Infinity;
  let _paramMatrix = args._paramMatrix? args._paramMatrix: null;
  let _ndone = args._ndone? args._ndone: 0;
  let _indices = args._indices? args._indices: 0;
  let paramMatrix = args.paramMatrix;
  if (Array.isArray(Nmif) || !isFinite(Nmif) || Nmif < 1)
    throw new Error("Nmif must be a positive integer.");
  
  Nmif = parseInt(Nmif);
  let ntimes = pompObject.times.length;
  if (typeof Np === "function" || Np === undefined || Np <= 0) {
    throw new Error(`Number of particles should be a positive number. ${Np} is not translated`);
  }
  
  if (!rw_sd) throw new Error("rw_sd function must be specified!");
  // pkern.sd: it produces the matrix of params rw in ntime 
  let rw_sd_matrix = []; 
  for (let i = 0; i < ntimes; i++) {
    rw_sd_matrix.push(rw_sd(pompObject.times[i]))
  }
  
  if (Array.isArray(coolingFraction) || !isFinite(coolingFraction) || coolingFraction <= 0 || coolingFraction > 1)
    throw new Error(`coolingFraction must be in (0,1]. coolingFraction = ${coolingFraction}`);
  let coolingFn = cooling(coolingType, coolingFraction, ntimes);
  if (!convRec)
    convRec = new Array(Nmif + 1).fill(null).map(a => new Object());
  
  if (!snapShotStart) {
    snapShotStart = 0;
    let start = params;
    if (_paramMatrix === null) {
      if (!start) start = coef(pompObject);
    } else { 
      throw new Error("Not translated, if paramMatrix is supplied");  
    }
    
    if (Object.keys(start).length === 0 || !Object.values(start).every(element => {return typeof element === 'number'}))
      throw new Error("parameters must be specified as a named numeric vector.")
    
    if (_paramMatrix === null) {
      paramMatrix = new Array(Np).fill(null).map(a => start);
    }
    
    convRec[0] = Object.assign({loglik: null, nfail: null}, start);
    if (transform)
      paramMatrix = partrans(pompObject, paramMatrix, dir="toEstimationScale");
  }
  
  let pfp;
  let n
  // one mif iteration
  // try {
  //   pfp = mif2Pfilter(
  //     object = pompObject,
  //     params = paramMatrix,
  //     Np = Np,
  //     mifiter = _ndone + n + 1,
  //     coolingFn,
  //     rw_sd =  rw_sd_matrix,
  //     tol = tol,
  //     maxFail,
  //     transform,
  //     _indices = _indices,
  //   )
  
  // } catch (error) {
  //   throw new Error(`In mif2: Iterate the filtering stopped: ${error}`)
  // }
  // if (!pfp) return {};
  
  // paramMatrix = pfp.paramMatrix;
  // convRec[n].loglik = pfp.loglik;
  // convRec[n].nfail = pfp.nfail;

  // convRec[n + 1].loglik = null;
  // convRec[n + 1].nfail = null;
  // Object.assign(convRec[n + 1], coef(pfp));
  // _indices = pfp.indices;
  // End of one mif iteration
  
  while( snapShotStart < Nmif) {
    n = snapShotStart;
    if (typeof progress === 'function') progress();
    // try {
      pfp = mif2Pfilter(
        object = pompObject,
        params = paramMatrix,
        Np = Np,
        mifiter = _ndone + n + 1,
        coolingFn,
        rw_sd =  rw_sd_matrix,
        tol = tol,
        maxFail,
        transform,
        _indices = _indices,
      )
      
    // } catch (error) {
    //   throw new Error(`In mif2: Iterate the filtering stopped: ${error}`)
    // }
    
    if (!pfp) return {};
    paramMatrix = pfp.paramMatrix;
    convRec[n].loglik = pfp.loglik;
    convRec[n].nfail = pfp.nfail;

    convRec[n + 1].loglik = null;
    convRec[n + 1].nfail = null;
    Object.assign(convRec[n + 1], coef(pfp));
    
    _indices = pfp.indices;
    
  // } else {// object for the snapshot
  //   result = {
  //     // ...args,
  //     snapShotStart: snapShotStart + 1,
  //     paramMatrix: paramMatrix,
  //     convRec: convRec,
  //     _indices: _indices,
  //   }
    snapShotStart++;  
  }
  
  if (transform)
    pfp.paramMatrix = partrans(pompObject,paramMatrix,dir="fromEstimationScale");
  
  return{
    loglik: pfp.loglik,
    params: pfp.params
  }
}

