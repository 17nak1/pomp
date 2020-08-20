/**
 *  @file             mif2Pfilter.js   
 *                    Mif version of particle filter algorithm.    
 *                         
 *  @references       https://github.com/kingaa/pomp                  
 *
 *  @author           Nazila Akhavan, nazila@kingsds.network
 *  @date             June 2020
 */
const mathLib = require("../../library/mathLib.js");
const { randwalk_perturbation } = require("./mif2cRW.js");
const { rprocessInternal } = require("../../library/rprocessInternal.js");
const { dmeasureInternal } = require("../../library/dmeasureInternal.js");
const { pfilter_computations } = require("../../library/pfilterComputations.js");
const { partrans, coef } = require("./mif2Helpers.js");

/**
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
exports.mif2Pfilter = function (object, params, Np, mifiter, coolingFn, rw_sd, tol = 1e-17, maxFail = Inf, transform, _indices = 0) {
  
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
    if (typeof progress === 'function') progress();
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
        x.push(object.initializer(initparams[i], object.interpolator(object.t0), object.globals));
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
      throw new Error(`In mif2.js: process simulation error: ${error}`);
    }
    
    let weights = [];
    // try {
      weights = dmeasureInternal(
        object,
        y = object.data[nt],
        X,
        times[nt + 1],
        transform ? tparams : params,
        log = false
      ); 
    // } catch (error) {
    //   throw new Error(`In mif2.js: error in calculation of weights: ${error}`);
    // }
    
    let allFinite = weights.map(w => isFinite(w)).reduce((a, b) => a & b, 1);
    if (!allFinite) {
      throw new Error("Mif2: In dmeasure weights returns non-finite value");
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
      throw new Error(`particle-filter error: ${error}`);
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
    console.log("warning! filtering failure occurred in mif2.");
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
