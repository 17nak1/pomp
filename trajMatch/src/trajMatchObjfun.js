const { coef, partrans } = require("./helpers");
const { dmeasureInternal } = require("../../library/dmeasureInternal.js");
const { trajectory } = require("./trajectoryInternal.js");

exports.trajMatchObjfun  = function (object, params, est, transform = false, args) {
  return tmofInternal(
    object=object,
    params=params,
    est=est,
    transform=transform,
    args
  )
}

const tmofInternal = function (object, params, est, transform, args) {

  ep = "in traj.match.objfun : "
  let tempParams = [];
  if (!est) est = [];

  if (Object.keys(params).length === 0) params = coef(object);
  
  if (transform) {
    tempParams = partrans(object, params, dir = "toEstimationScale");
  }
    
  let parEstIdx = est;

  if (parEstIdx.some(a => {return a === NaN}))
    throw new Error(ep + "est does not match with parameters")

  return  (par,parX) => {
    let d;
    if (parEstIdx.length > 0) {
      for (let i = 0; i < parEstIdx.length; i++) {
        tempParams[parEstIdx[i]]= parX[i+1];
      }
    } 
    
    if (transform) tparams = partrans(object, tempParams, dir="fromEstimationScale");
    
    let x = trajectory(
      object,
      params = transform? tparams : tempParams,
      args
    )
    d = dmeasureInternal(
      object,
      y=object.data,
      x,
      times = object.times,
      params = transform? tparams : tempParams,
      log = true
    )
    return -d.reduce((a, b) => a + b, 0);
  }
}
