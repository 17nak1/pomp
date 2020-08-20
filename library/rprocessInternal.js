/**
 *  @file        rprocessInternal.js
 *               Simultes values of states at time t_{n+1}.
 *                 
 *  @author       Nazila Akhavan, nazila@kingsds.network
 *  @date        june 2020
 */
const { euler_model_simulator } = require("./euler.js");

/**
 * @param {object} object       An object of class POMP
 * @param {array} xstart        Array of objects of states.
 * @param {array} times         Array of [t_n, t_{n+1}].
 * @param {array} params        Array of objects of parameters.
 * @param {number} offset 
 * @returns {array}
 *  Array of 'Np' objects of states at time t_{n+1}.
 */
exports.rprocessInternal  = function (object, xstart, times, params, offset = 0) {
  
  let ntimes = times.length;
  if (ntimes < 2) {
    throw new Error("in 'rprocess': length(times) < 2: with no transitions, there is no work to do.");
  }

  let off = Number(offset) ;
  if ((off < 0)||(off >= ntimes))
    throw new Error(`illegal 'offset' value ${off}`);
  
  let nrepsx = xstart.length;  //Np
  let nreps = params.length;

  if (nrepsx > nreps) {		// more ICs than parameters
    if (params.length === 1) params = new Array(nrepsx).fill(null).map(a => params[0]);
  } else if (nrepsx < nreps) {	// more parameters than ICs
    throw new Error("More parameters than ICs is not translated!")
  }
  // extract the process function. NOTE: only Euler translated (type = 3).
  let type = object.rprocessDetail.type === "euler_sim" ? 3: 0;
  let X, fn,deltat;
  switch (type) {
    case 1: // one-step simulator
      fn = object.rprocess;
      deltat = 1.0;
      X = euler_model_simulator(fn, xstart, times, params, deltat, type, object);
      break;

    case 2: case 3: // discrete-time and Euler
      fn = object.rprocess;
      deltat = Number(object.rprocessDetail.deltaT);
      X = euler_model_simulator(fn, xstart, times, params, deltat, type, object);
      break;  

    case 4: // Gillespie's method
      throw new Error("in 'rprocess': Gillespie's method is not translated")
      
    case 0: default:
      throw new Error("'rprocess' is undefined. Note: only 'euler_sim' (discrete-time Euler) method is translated");
  }
 
  return X; 
}  
