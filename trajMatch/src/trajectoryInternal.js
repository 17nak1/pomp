const { initState } = require("./initState.js");
const { iterateMap } = require("./trajectory.js");

exports.trajectory = function (object, params, times, t0, asDataFrame, args){

  let ep = 'in \"trajectory\": ';

  if (times === undefined || !Array.isArray(times))
    times = object.times;

  if (!Array.isArray(times) || times.length == 0)
    throw new Error(ep + '\"times\" is empty, there is no work to do');

  let isAscending = times.slice(1).map((e,i) => e > times[i]).every(x => x);
  if (!isAscending)
    throw new Error(ep + '\"times\" must be an increasing sequence of times');

  if (t0 === undefined)
    t0 = object.t0;

  if (t0 > times[0])  
    throw new Error(ep + '\"times\" the zero-time \"t0\"  must occur no later than the first observation');

  let ntimes = times.length;

  if (params === undefined) params = object.params;
  
  if (Object.keys(params).length === 0) 
    throw new Error(ep + '\"params\" must be supplied');
  

  params = [params]//as.matrix(params)
  let nrep = params.length;
  let paramnames = Object.keys(params[0]);
  if (paramnames.length <= 0)
    throw new Error(ep + '\"params\" must have rownames');

  let x0 = initState(object, params);
  let nvar = x0[0].length;
  statenames = Object.keys(x0[0]);

  let type = object.skeletonDetail.type;          // map or vectorfield?

  let x;
  if (type === "map") {
    try {
      x = iterateMap(object, times, t0, x0, params);
    } catch (error) {
      throw new Error(` ${ep} in map iterator: ${error}`)
    }  

  } else if (type === "vectorfield") {
    throw new Error("vectorfield is not translated.");
  
  } else {
    throw new Erro(`${ep} deterministic skeleton has not been properly specified`);
  }

  return x
}