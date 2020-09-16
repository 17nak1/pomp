/**
 *  @file             mif2cRW.js   
 *                    Randwalk_perturbation adds random normal value to the parameters.    
 *                         
 *  @references       https://github.com/kingaa/pomp                  
 *
 *  @author           Nazila Akhavan, nazila@kingsds.network
 *  @date             June 2020
 */
const { qnorm } = require('lib-r-math.js/dist/src/lib/normal/qnorm');

/** 
 * @param {matrix} params  Initial parameters befor transform.
 * @param {array} rw_sd    An array of objects with random walk values. 
 */
exports.randwalk_perturbation = function (params, rw_sd) { 
  let rn;
  for(let i = 0; i < params.length; i++) {
    Object.keys(rw_sd).forEach(key => {
      if (rw_sd[key] !== null) {
        rn = normalRand()
        while (!isFinite(rn)) {
          rn = normalRand()  
        }
        params[i][key] += rw_sd[key] * rn;
      } ;
    });
  }

  return params;
}


// Normal RNG using INVERSION method
const normalRand = function () {
  let BIG = 134217728; /* 2^27 */
	u1 = Math.random();
  u1 = BIG * u1 + Math.random();
  
	return qnorm(u1 / BIG, 0.0, 1.0, 1, 0);
}
