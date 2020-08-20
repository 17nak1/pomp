
/**
 *  @file        dmeasureInternal.js
 *               Calculates the liklihood using snippet.dmeasure.
 *                 
 *  @author      Nazila Akhavan, nazila@kingsds.network
 *  @date        june 2020
 */

/**
 * @param {object} object       An object of class POMP
 * @param {array} y             Array of objects of data.
 * @param {array} x             Array of objects of states.
 * @param {array} times         Array of data times.
 * @param {array} params        Array of objects of parameters.
 * @param {boolean} give_log    If the liklihood should be return in log scale.
 * @returns {array}
 *  Array of 'Np' liklihood values at time t_n.
 */
exports.dmeasureInternal  = function (object, y, x, times, params, give_log) {
  ntimes = Array.isArray(times) ? times.length : 1;
  if (ntimes < 1)
    throw new Error("In 'dmeasureInternal': times is not defined");
  
  let nrepsx = x[0].length? x[0].length : x.length;; //Np
  let nrepsp = Object.keys(params).length;
  
  nreps = (nrepsp > nrepsx) ? nrepsp : nrepsx;

  if ((nreps % nrepsp != 0) || (nreps % nrepsx != 0))
    throw new Error("in 'dmeasureInternal': larger number of replicates is not a multiple of smaller"); 

  let ff = object.dmeasure;
  let F;
  if(ntimes > 1) {
    F = new Array(ntimes);
    for (k = 0; k < ntimes; k++) { // loop over times.Note:Used in trajMatch
      F[k] = ff(y[k], x[k][0], params,give_log, object.globals);
    }
  } else {
    // ntimes = 1; Note:Used in pfilter and mif2
    F = new Array(nreps);
    for (let j = 0; j < nreps; j++) { // loop over replicates
      if (nrepsp === 1) params[j] = params[0];
      F[j] = ff(y, x[j], params[j],give_log, object.globals);
    }
  }

  return F;
}