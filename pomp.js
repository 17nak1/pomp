function (data, times, t0, ..., rprocess, dprocess, rmeasure, 
  dmeasure, measurementModel, skeleton, initializer, rprior, 
  dprior, params, covar, tcovar, obsnames, statenames, paramnames, 
  covarnames, zeronames, PACKAGE, fromEstimationScale, toEstimationScale, 
  globals, cdir, cfile, shlibArgs) 
{
var data , c1,c2
var covar, tcovar, obsnames, statenames, paramnames , covarnames, zeronames
var fromEstimationScale, toEstimationScale,rmeasure,dmeasure

  let ep = 'in ‘pomp’: '
  if (typeof data === 'undefined') {
    throw ep + ' "data" is a required argument.'
  }
  //nargs():When used inside a function body, nargs returns the number of arguments supplied to that function, including positional arguments left blank.
  // if (nargs() == 1) 
  //   return(data)
  // //is() Functions to test inheritance relationships between an object and a class or between two classes (extends).
  // if (!is(data, "data.frame") && !is(data, "pomp")) 
  //   stop(ep, sQuote("data"), " must be a data frame or an object of class ", 
  //     sQuote("pomp"), call. = FALSE)
  c1 = typeof fromEstimationScale === 'undefined'
  c2 = typeof vtoEstimationScale === 'undefined'
  if (!c1 !== !c2) {
    throw ep + 'if one of ‘toEstimationScale’ or ‘fromEstimationScale’ is supplied, then so must the other.'
  }
  c1 = typeof covar === 'undefined'
  c2 = typeof tcovar === 'undefined'
  if (!c1 !== !c2) { 
    throw ep + 'if one of ‘covar’ or ‘tcovar’ is supplied, then so must the other.'
  }
  if (typeof measurementModel !== 'undefined') {
    if(!(typeof dmeasure === 'undefined' || dmeasure === null) || !(typeof rmeasure === 'undefined' || rmeasure === null)) {
      console.log('Warning :' + ep + 'specifying ‘measurement.model’ overrides specification of ‘rmeasure’ and ‘dmeasure’.')
    }
    // mm = measform2pomp(measurementModel)
    // rmeasure =  mm.rmeasure
    // dmeasure =  mm.dmeasure
  }
  skelType = 'undefined'
  skelmapDelta_t = 1
  if (typeof skeleton === 'undefined') {
    skeleton = null
  } else if (skeleton === null) {
    skelType = 'remove'
  // } else if (is(skeleton, "safecall")) {
  //   skeleton <- skeleton@call
    map = function(f, deltaT = 1) {
      skelType <<- "map"
      skelmapDelta_t <<- Nember(deltaT)
      if (skelmapDelta_t <= 0) {
        throw 'in ‘map’, ‘delta.t’ must be positive.'
      }
      return f
    }
     vectorfield = function(f) {
      skelType <<- "vectorfield"
      return f
    }
    flist = list(map, vectorfield)
  //   skeleton <- eval(skeleton, envir = flist, enclos = parent.frame())
  } else {
    throw '‘skeleton’ must be specified as either a  ‘vectorfield’ or a ‘map’.'
  }
  construct_pomp(data = data, times = times, t0 = t0, ..., 
    rprocess = rprocess, dprocess = dprocess, rmeasure = rmeasure, 
    dmeasure = dmeasure, initializer = initializer, skelType = skelType, 
    skelmapDelta_t = skelmapDelta_t, skeleton = skeleton, 
    rprior = rprior, dprior = dprior, params = params, covar = covar, 
    tcovar = tcovar, obsnames = obsnames, statenames = statenames, 
    paramnames = paramnames, covarnames = covarnames, zeronames = zeronames, 
    PACKAGE = PACKAGE, fromEstimationScale = fromEstimationScale, 
    toEstimationScale = toEstimationScale, globals = globals, 
    cdir = cdir, cfile = cfile, shlibArgs = shlibArgs)
    
}
  