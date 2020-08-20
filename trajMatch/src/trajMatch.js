
const { coef } = require("./helpers");
const { trajectory } = require("./trajectoryInternal.js");
const snippet = require('../../library/modelSnippetCOVID3.js');
const { minimInternal } = require("./minimInternal.js");
const { trajMatchObjfun } = require("./trajMatchObjfun.js");
const pomp = require("../../library/pomp.js");

exports.trajMatch = function (params, args) {
  let object = args.object;
  let start = params;
  let est = args.est? args.est: [];
  let method = args.method? args.method : c("Nelder-Mead","subplex","SANN","BFGS", "sannbox","nloptr");
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

  // fill params slot appropriately
  object.params = m.params

  // fill states slot appropriately
  x = trajectory(object)
  object.states = x;

  return Object.assign(object,
    {transform: transform,
    est: est,
    value: -m.value,
    evals: m.evals,
    convergence: m.convergence,
    msg: m.msg
  })
}



