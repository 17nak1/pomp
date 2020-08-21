/**
 *  @file             Helpers.js   
 *                    depentent functions.  
 *                         
 *  @references       https://github.com/kingaa/pomp                  
 *
 *  @author           Nazila Akhavan, nazila@kingsds.network
 *  @date             June 2020
 */
/**
 * Calculates alpha and gamma based on cooling type.
 * @param {string} type 
 * @param {number} fraction 
 * @param {number} ntimes
 * @returns {object}
 *  @params {number} alpha 
 *  @params {number} gamma 
 */
const cooling = function(type, fraction, ntimes) {
  switch(type){
  case "geometric":
    let factor = Math.pow(fraction, (1/50));
    return function (nt, m) {
      let alpha = Math.pow(factor, (nt / ntimes + m - 1));
      return {alpha: alpha, gamma: alpha ** 2};
    }
    
  case "hyperbolic":
    if (fraction < 1) {
      let scal = (50 * ntimes * fraction - 1) / (1 - fraction);
      return function (nt, m) {
        let alpha = (1 + scal) / (scal + nt + ntimes * (m - 1));
        return {alpha: alpha, gamma: alpha ** 2};
      }
    } else {
      return function (nt, m) {
        return {alpha: 1, gamma: 1};
      }
    }
  }
}

/**
 * Rescale the parameters.
 * @param {object} pomp 
 * @param {array} params  An array of objects of parameters.
 * @param {string} dir    "fromEstimationScale"/ "toEstimationScale"
 * @returns {array} 
 *  An array of objects of transformed parameters
 */
const partrans = function (pomp, params, dir = ["fromEstimationScale","toEstimationScale"]) {
  if (!Array.isArray(params)) params = [params];
  let transParam = [].concat(params);
  switch(dir){
  case "fromEstimationScale":
    for (let i = 0 ; i < params.length; i++) {
      transParam[i] = pomp.fromEstimationScale(params[i]);
    }
    break;
    
  case "toEstimationScale":
    for (let i = 0 ; i < params.length; i++) {
      transParam[i] = pomp.toEstimationScale(params[i]);
    }
    break;
  } 
  
  if (transParam.length === 1) return transParam[0];
  return transParam;
}
/**
 * 
 * @param {object}  POMP         An object of class POMP.
 * @param {boolean} transform    If pomp.params need to be transformed
 * @returns {object}
 *  pomp.params
 */
const coef = function (object, transform = false) {
  if (Object.keys(object).length === 0) {
    console.error("In coef(): object should be specified");
    return false;
  }

  if (Object.keys(object.params).length > 0) {
    if (transform) {
      params = partrans(object, [object.params], dir="fromEstimationScale");
    } else {
      params = [object.params];
    }    
    return params[0];
  } else {
    return 0;
  }
}

module.exports = {cooling, partrans, coef}
