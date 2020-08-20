
/**
 *  @file             initState.js   
 *                    Scales initial states using object.initializer()    
 *                         
 *  @references       https://github.com/kingaa/pomp                  
 *
 *  @author           Nazila Akhavan, nazila@kingsds.network
 *  @date             June 2020
 */

/**
 * @param {object} object     An object of class POMP.
 * @param {object} params     An object of initial parameters.
 * @param {number} nsim       Number of simulations (Np).
 * @return {array}
 *  Array of Np arrays; each contains population value of states;
 */
exports.initState = function (object, params, nsim) {
  if(params === undefined) params = object.coef;
  if(isNaN(nsim)  || nsim === null) nsim = params.length;
  
  return do_init_state(object, params, nsim)
}

const do_init_state = function (object, params, ns) {
  let nrep = params.length;
  let paramsObj = params[0];
  if (ns % nrep != 0) 
    throw new Error("in 'init.state': number of desired state-vectors 'nsim' is not a multiple of ncol('params')");
  let args = object.globals;
  let initVector =  object.initializer(paramsObj, object.interpolator(object.t0), args);  
  return  new Array(ns).fill(null).map(a => initVector);
}
