const { partrans } = require("../../library/helpers");
const subplex = require('../../library/subplex/subplex.js');

exports.minimInternal = function(objfun, start, est, object, method, transform, lower = null, upper = null, lb = lower, ub = upper, args) {
  
  let ep = "In minimInternal:";
  if (Object.keys(start).length < 1)
    throw new Error(ep +"start must be supplied")
  let guess = {};
  if (transform) {
    start = partrans(object,start,dir="toEstimationScale")
    if (Object.keys(start).some(a => {a === null}))
      throw new Error("'est' must refer to parameters named in partrans(object,start,dir=toEstimationScale")
  } else {
    if (Object.keys(start).some(a => {a === null}))
      throw new Error("'est' must refer to parameters named in start.")
  }
  if (est.length > 0) {
    Object.keys(start).forEach(key => {
      if (est.includes(key)) guess[key] = start[key];
    })
  }

  let val;
  if (est.length === 0) {
    
    val = objfun(guess);
    conv = NaN;
    evals = [1,0];
    msg = "no optimization performed";
  } else {

    if (method == 'subplex') {
      
      subplex.f = objfun;
      subplex.x0 = Object.values(guess);
      subplex.tol = 2.220446e-16;
      opt = subplex.run();

    } else if (method=="sannbox") {
      throw new Error ("Method 'sannbox' is not translated");
    } else if (method=="nloptr") {
      throw new Error ("Method 'nloptr' is not translated");
    } else {
      throw new Error ("Only 'subplex' is translated. You need define method ='subplex'");
      // opt <- optim(par=guess,fn=objfun,method=method,control=opts)
    }

    msg = opt[3];
    if (method == "nloptr") {
      throw new Error ("Method 'nloptr' is not translated");
    } else {

      val = opt[1];
      for(let i = 0 ; i < est.length; i++)
        start[est[i]] = opt[0][i];
      conv = 0;//opt.convergence;
      evals = opt[2];
    }
  }

  if (transform) start = partrans(object, start, dir = "fromEstimationScale");

  return {
    params : start,
    est : est,
    transform : transform,
    value : val,
    convergence : parseInt(conv),
    evals : evals,
    msg : msg
  }
}
