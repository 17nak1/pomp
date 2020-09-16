
const { coef } = require("../../library/helpers");
const { trajectory } = require("./trajectoryInternal.js");
const snippet = require('../../library/modelSnippet.js');
const { minimInternal } = require("./minimInternal.js");
const { trajMatchObjfun } = require("./trajMatchObjfun.js");
const pomp = require("../../library/pomp.js");

/**
 * 
 * @param {array} params                Input parameters.
 * @param {object} args 
 * @param {object} args.object          An object of class POMP.
 
 * @param {array} args.est              Array of stings of esimating parameters' name.
 * @param {boolean} args.transform      If parmeters need to transform.
 * @param {string} args.method          Namr of the optomizing method.

 * 
 * @returns {object}
 *  @param {number} logLik         The estimated log likelihood.
 *  @param {object} params         The estimated parameters.
 */
exports.trajMatch = function (params, args) {

   /* Check if input data is string and convert them to numbers*/
  let covarkeys = Object.keys(args.object.covar[0])
  for(let k =0; k < args.object.covar.length; k++) {
    for (let j = 0; j < covarkeys.length; j++)
      if((args.object.covar[k][covarkeys[j]]) === 'NaN')  {
        args.object.covar[k][covarkeys[j]] = NaN;
      } else {
        args.object.covar[k][covarkeys[j]] = Number(args.object.covar[k][covarkeys[j]]);
      }
  }

  let datakeys = Object.keys(args.object.data[0])
  for(let k =0; k < args.object.data.length; k++) {
    for (let j = 0; j < datakeys.length; j++)
      if((args.object.data[k][datakeys[j]]) === 'NaN')  {
        args.object.data[k][datakeys[j]] = NaN;
      } else {
        args.object.data[k][datakeys[j]] = Number(args.object.data[k][datakeys[j]]);
      }
  }
  
  let object = args.object;
  let start = params;
  let est = args.est? args.est: [];
  let method = args.method? args.method : ["Nelder-Mead","subplex","SANN","BFGS", "sannbox","nloptr"];
  let transform = args.transform ? args.transform: false;

  const pompData = Object.assign(object,{
    skeleton: snippet.skeleton,
    dmeasure: snippet.dmeasure,
    initializer: snippet.initializer,
    toEstimationScale: snippet.toEstimationScale,
    fromEstimationScale: snippet.fromEstimationScale,
    params: start
  });
    
  object = new pomp(pompData);
  
  if (!start) start = coef(object);
  let objfun = trajMatchObjfun(
    object=object,
    params=start,
    est=est,
    transform=transform);
  
  let m = minimInternal(
    objfun,
    start,
    est,
    object,
    method,
    transform
  )

  object.params = m.params;                           /* fill params slot appropriately */
  object.states = trajectory(object);                 /* fill states slot appropriately */

  return {                                            /* Here only return what we need for the model and comment extras */
    params:object.params,
    loglik: -m.value,
  }
}



