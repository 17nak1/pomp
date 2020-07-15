/**
 *  @file             mif2.js   
 *                    Maximum likelihood by Iterated Filtering.    
 *                         
 *  @references       https://github.com/kingaa/pomp                  
 *
 *  @author           Nazila Akhavan, nazila@kingsds.network
 *  @date             June 2020
 */
const mathLib = require("../library/mathLib.js");
const { randwalk_perturbation } = require("./mif2cRW.js");
const { rprocessInternal } = require("../library/rprocessInternal.js");
const { dmeasureInternal } = require("../library/dmeasureInternal.js");
const { pfilter_computations } = require("../library/pfilterComputations.js");
const { cooling, partrans, coef } = require("./mif2Helpers.js");
/**
 * An iterated filtering algorithm for estimating the parameters of a partially-observed Markov process.
 * Running mif2 causes the algorithm to perform a specified number of particle-filter iterations.
 * At each iteration, the particle filter is performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
 * This extra variability effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.
 * As the iterations progress, the magnitude of the perturbations is diminished according to a user-specified cooling schedule.
 * The algorithm is presented and justified in Ionides et al. (2015).
 * 
 * @param {object} args 
 * @param {object} args.object          An object of class POMP.
 * @param {number} args.Nmif            Number of of filtering iterations.
 * @param {array} args.start            Input parameters.
 * @param {boolean} args.transform      If parmeters need to transform.
 * @param {object} args.rw_sd           Specification of the magnitude of the random-walk perturbations that will be applied to some or all model parameters.
 * @param {number} args.Np              The number of particles to use.
 * @param {number} args.coolingFraction Cooling Rate
 * @param {string} args.coolingType     specifications for the cooling schedule,
 *                                       i.e., the manner with which the intensity of the parameter perturbations is reduced with successive filtering iterations.
 
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
 *  @param {number} Np 
 *  @param {number} tol 
 *  @param {number} nfail 
 */
exports.mif2 = function (args) {
  let pomp = args.object;
  let Nmif = args.Nmif;
  let start = args.start;
  let transform = args.transform? args.transform: FALSE;
  let rw_sd = args.rw_sd;
  let Np = args.Np;
  let coolingType = args.coolingType ? args.coolingType: ["hyperbolic", "geometric"];
  let coolingFraction = args.coolingFraction;
  let tol = args.tol? args.tol: 1e-17;
  let maxFail = args.maxFail? args.maxFail: Infinity;
  let _paramMatrix = args._paramMatrix? args._paramMatrix: null;
  let _ndone = args._ndone? args._ndone: 0;
  let _indices = args._indices? args._indices: 0;
  if (Array.isArray(Nmif) || !isFinite(Nmif) || Nmif < 1)
    throw new Error("Nmif must be a positive integer.");
    
  Nmif = parseInt(Nmif);
  if (_paramMatrix === null) {
    if (!start) start = coef(pomp);
  } else { 
    throw new Error("Not translated, if paramMatrix is supplied");  
  }
  
  if (Object.keys(start).length === 0 || !Object.values(start).every(element => {return typeof element === 'number'}))
    throw new Error("parameters must be specified as a named numeric vector.")

  let ntimes = pomp.times.length;
  if (typeof Np === "function" || Np === undefined || Np <= 0) {
    throw new Error(`Number of particles should be a positive number. ${Np} is not translated`);
  }

  if (!rw_sd) throw new Error("rw_sd function must be specified!");
  // pkern.sd: it produces the matrix of params rw in ntime 
  let rw_sd_matrix = [];
  for (let i = 0; i < ntimes; i++) {
    rw_sd_matrix.push(rw_sd(pomp.times[i]))
  }

  if (Array.isArray(coolingFraction) || !isFinite(coolingFraction) || coolingFraction <= 0 || coolingFraction > 1)
    throw new Error(`coolingFraction must be in (0,1]. coolingFraction = ${coolingFraction}`);
  let coolingFn = cooling(coolingType, coolingFraction, ntimes);
  let paramMatrix;
  if (_paramMatrix === null) {
    paramMatrix = new Array(Np).fill(null).map(a => start);
  } else {
    paramMatrix = _paramMatrix;
  }

  convRec = new Array(Nmif + 1).fill(null).map(a => new Object());
  convRec[0] = Object.assign({loglik: null, nfail: null}, start)
  if (transform)
    paramMatrix = partrans(pomp, paramMatrix, dir="toEstimationScale");
  
  // Iterate the filtering main loop 
  let pfp;
  for (let n = 0; n < Nmif; n++) {
    try {
      pfp = mif2Pfilter(
        object=pomp,
        params=paramMatrix,
        Np=Np,
        mifiter=_ndone + n + 1,
        coolingFn,
        rw_sd= rw_sd_matrix,
        tol=tol,
        maxFail,
        transform,
        _indices=_indices,
      )
      
    } catch (error) {
      console.error(`Iterate the filtering stoped: ${error}`)
    }
    paramMatrix = pfp.paramMatrix;
    convRec[n].loglik = pfp.loglik;
    convRec[n].nfail = pfp.nfail;

    convRec[n + 1].loglik = null;
    convRec[n + 1].nfail = null;
    Object.assign(convRec[n + 1], coef(pfp));
    
    _indices = pfp.indices;
  }
  
  if (transform)
    pfp.paramMatrix = partrans(pomp,paramMatrix,dir="fromEstimationScale");
  
  return{
    ...pfp,
    Nmif: Nmif,
    rw_sd: rw_sd,
    coolingType: coolingType,
    coolingFraction: coolingFraction,
    transform: transform,
    convRec: convRec
  }
}
/**
 * Mif version of particle filter algorithm.
 * @param {object} object              An object of class POMP
 * @param {object} args.params         Input parameters.
 * @param {number} args.Np             The number of particles to use.
 * @param {number} mifiter 
 * @param {function} coolingFn 
 * @param {array} rw_sd                Array of objects of the parameters that has random walk. 
 * @param {number} args.tol            Tolerance. 
 * @param {number} maxFail 
 * @param {boolean} transform 
 * @param {number} _indices 
 * @returns {object}
 */
const mif2Pfilter = function (object, params, Np, mifiter, coolingFn, rw_sd,
  tol = 1e-17, maxFail = Inf, transform, _indices = 0)
{
  if ((Array.isArray(tol) && tol.length !== 1) || tol === Infinity || tol < 0)
    throw new Error(`${tol} should be a small positive number in mif2Pfilter.`);

  let do_ta = !!_indices.length ;
  if (do_ta)
    throw new Error(` ${_indices} has improper length in mif2Pfilter.`);

  let times = [object.t0, ...object.times];
  let ntimes = times.length - 1;
  let loglik = new Array(ntimes);
  let effSampleSize = Number(ntimes);
  let nfail = 0;
  let alpha;
  let x = [];

  for (let nt = 0; nt < ntimes; nt++) {//ntimes
    alpha = coolingFn(nt + 1,mifiter).alpha;
    pmg = Object.assign({}, rw_sd[nt]);
    Object.keys(pmg).map(key => pmg[key] *= alpha);
    params = randwalk_perturbation(params, pmg);
    
    if (transform)
      tparams = partrans(object, params, dir="fromEstimationScale");
    
    if (nt === 0) {
      //get initial states
      let initparams = transform ? tparams : params;
      for (let i = 0; i < params.length; i++) {
        x.push(object.initializer(initparams[i], object.interpolator(object.t0)));
      }
    } 
    
    let X = [];
    try {
      X = rprocessInternal(
        object,
        xstart = x,
        [times[nt],times[nt + 1]],
        transform ? tparams : params,
        offset=1
      );
    } catch (error) {
      console.error(`In mif2.js: process simulation error: ${error}`);
    }
    
    let weights = [];
    try {
      weights = dmeasureInternal(
        object,
        y = object.data[nt],
        X,
        times[nt + 1],
        transform ? tparams : params,
        log = false
      ); 
    } catch (error) {
      console.error(`In mif2.js: error in calculation of weights: ${error}`);
    }
    
    let allFinite = weights.map(w => isFinite(w)).reduce((a, b) => a & b, 1);
    if (!allFinite) {
      throw new Error("In dmeasure: weights returns non-finite value");
    }

    // compute weighted mean at last timestep
    if (nt === ntimes - 1) {
      if (weights.map(w => w>0).reduce((a, b) => a || b, 0)) {
        // replace and fill object.params instead of coef(object). This is the same thing.
        object.params = mathLib.mean(params, w = weights);
        object.params = coef(object, transform = transform);
      } else {
        console.warn("filtering failure at last filter iteration, using unweighted mean for 'coef' ");
        object.params = mathLib.mean(params);
        object.params = coef(object, transform = transform);
      }
    }

    // compute effective sample size, log-likelihood
    // also do resampling if filtering has not failed
    let xx;
    try {
      xx  = pfilter_computations(
        X,
        params,
        Np = Np,
        0, //rw_sd
        predmean = false,
        predvar = false,
        filtmean = false,
        trackancestry = do_ta,
        onepar = false,
        weights = weights,
        tol = tol
      );
    } catch (error) {
      console.error(`particle-filter error: ${error}`) 
    }

    let allFail = xx.fail;
    loglik[nt] = xx.loglik;
    effSampleSize[nt] = xx.ess;
    if (do_ta) {
      _indices = _indices[xx.ancestry];
    }

    x = xx.states;
    params = xx.params;// should be in toScale.

    if (allFail) { // all particles are lost
      nfail = nfail + 1;
      if (nfail > maxFail)
        throw new Error("In mif2Pfilter: too many filtering failures")
    }    
  } // end of nt loop   
    

  if (nfail > 0) {
    console.log("warning! filtering failure occurred.");
  }

  return {
    ...object,
    paramMatrix: params,
    effSamplesize: effSampleSize,
    condLoglik: loglik,
    indices: _indices,
    Np: Np,
    tol: tol,
    nfail: Number(nfail),
    loglik: loglik.reduce((a,b) => a+b, 0)
  }
}


