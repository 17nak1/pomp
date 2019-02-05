
var mathLib = require('./mathLib')
pfilterComputations = function (x, params, Np, predmean = 0, predvar = 0, filtermean = 0, trackancestry = 0, doparRS = 0, weights, tol)
{
//   var x, params, Np, predmean = 0, predvar = 0, filtmean = 0, trackancestry = 0, doparRS = 0, weights, tol = 1e-17
//   x = [[97497.5866824985, 1373.89351938552,  877.113877358443 , 3073346.86162676],
// [98299.5866824985,  1356.89351938552,  896.113877358443, 3072542.86162676],
// [96063.5866824985,  1486.89351938552,  925.113877358443,  3074619.86162676],
// [97373.5866824985,  1326.89351938552,  876.113877358443,  3073518.86162676],
// [97184.5866824985,  1401.89351938552,  885.113877358443,  3073623.86162676]]

//   params = [[29.2655634435, 0.3944025839968, 73.05, 0.004303868154365, 45.66, 0.417917730919406, 0.530358893934836]]
//   Np = 5
//   weights = [0.19963968806878257,0.20165353289993293,0.19678460595724773,0.20090831859789957,0.20101385447613734 ]
  var nprotect = 0
  var pm = null, pv = null, fm = null, anc = null
  var ess, fail, loglik
  var newstates = null, newparams = null
  var retval = new Array(9)
  // const char *dimnm[2] = {"variable","rep"}
  var xpm = 0, xpv = 0, xfm = 0, xw, xx = 0, xp = 0//pointer
  var xanc = 0//pointer
  var dim = [], dimP, newdim, Xnames, Pnames
  var dimAdd//pointer
  var np
  var nvars, npars = 0, nreps, nlost
  var do_pm, do_pv, do_fm, do_ta, do_pr, all_fail = 0
  var sum = 0, sumsq = 0, vsq, ws, w, toler
  var j, k

  dim[0] = x.length 
  dim[1] = x[0].length
  nvars = dim[1]//number of states
  nreps = dim[0]//states repeats
  xx = [].concat(x)//REAL(x)
  dim[0] = params.length//params
  dim[1] = params[0].length 
  npars = dim[0]
  // if (nreps % dim[0] != 0)
  //   errorcall(null,"ncol('states') should be a multiple of ncol('params')"); // # nocov : (((paramsCol * int = statesCol)))
  // np = *(INTEGER(AS_INTEGER(Np))); // number of particles to resample
  np = Number(Np)

  do_pm = predmean // calculate prediction means?
  do_pv = predvar// calculate prediction variances?
  do_fm = filtmean// calculate filtering means?
  do_ta = trackancestry// track ancestry?
  do_pr = doparRS// Do we need to do parameter resampling?

  if (do_pr) {
    if (dim[1] != nreps)// #states repeats should be equal to # params repeats
      throw "ncol('states') should be equal to ncol('params')" // # nocov
  }

  var Pess = 1 //= NEW_NUMERIC(1)// effective sample size
  var loglik = 1 //= NEW_NUMERIC(1)// log likelihood
  var fail = 1// = NEW_LOGICAL(1)// particle failure?

  xw = [].concat(weights)
  toler = tol   // failure tolerance

  // check the weights and compute sum and sum of squares
  for (k = 0, w = 0, ws = 0, nlost = 0; k < nreps; k++) {
    if (xw[k] > toler) {
      w += xw[k]
      ws += Math.pow(xw[k],2)
    } else {      // this particle is lost
      xw[k] = 0
      nlost++
    }
  }
  if (nlost >= nreps) {
    all_fail = 1 // all particles are lost
  }
  if (all_fail) {
    loglik = Math.log(toler) // minimum log-likelihood
    ess = 0 // zero effective sample size
  } else {
    loglik = Math.log(w / nreps) // mean of weights is likelihood
    ess = Math.pow(w, 2) / ws  // effective sample size
  }
  fail = all_fail

  if (do_pm || do_pv) {
    pm = Number(nvars)//NEW_NUMERIC
    xpm = new Array(pm)//REAL(pm);
  }

  if (do_pv) {
    pv = Number(nvars)
    xpv = new Array(pv)//REAL(pv);
  }

  if (do_fm) {
    fm = Number(nvars)
    xfm = new Array(fm)//REAL(fm);
  }

  if (do_ta) {
    anc = parseInt(np)//NEW_INTEGER
    xanc = new Array(anc)//INTEGER(anc);
  }

  for (j = 0; j < nvars; j++) { // state variables

    // compute prediction mean
    if (do_pm || do_pv) {
      for (k = 0, sum = 0; k < nreps; k++) {
        sum += xx[k][j]
      }
      sum /= nreps
      xpm[j] = sum
    }

    // compute prediction variance
    if (do_pv) {
      for (k = 0, sumsq = 0; k < nreps; k++) {
        vsq = xx[k][j] - sum//xx[j + k * nvars] - sum
        sumsq += vsq * vsq
      }
      xpv[j] = sumsq / (nreps - 1)
    }
    //  compute filter mean
    if (do_fm) {
      if (all_fail) { // unweighted average
        for (k = 0, ws = 0; k < nreps; k++) {
          ws += xx[k][j]
        }
        xfm[j] = ws / nreps
      } else { // weighted average
        for (k = 0, ws = 0; k < nreps; k++) {
          ws += xx[k][j] * xw[k]
        }
        xfm[j] = ws / w
      }
    }
  }
  if (!all_fail) { // resample the particles unless we have filtering failure
    var xdim = new Array(2)
    var sample = new Array(np)
    var ss, st, ps, pt 

    // create storage for new states
    xdim[0] = nvars; xdim[1] = np;
    newstates = new Array(xdim[1]).fill(0).map(() => Array(xdim[0]))
    // setrownames(newstates,Xnames,2);
    // fixdimnames(newstates,dimnm,2);
    // ss = [].concat(x)
    st = newstates

    // create storage for new parameters
    if (do_pr) {
      xdim[0] = npars; xdim[1] = np;
      newparams = new Array(xdim[1]).fill(0).map(() => Array(xdim[0]))
      // setrownames(newparams,Pnames,2);
      // fixdimnames(newparams,dimnm,2);
      ps = params
      pt = newparams
    }

    // resample
    mathLib.nosortResamp(nreps, weights, np, sample, 0)
    for (k = 0; k < np; k++) { // copy the particles
      for (j = 0; j < nvars; j++) {
        st[k][j] = xx[sample[k]][j]
      } 
      if (do_pr) {
        for (j = 0; j < npars; j++) {
          pt[k][j] = xp[sample[k][j]]
        }
      }
      if (do_ta) {
        xanc[k] = sample[k] + 1
      }
    }
  } else { // don't resample: just drop 3rd dimension in x prior to return

    newdim = Number(2)
    dim = Number(newdim)
    dim[0] = nvars; dim[1] = nreps;
    // SET_DIM(x,newdim);
    // setrownames(x,Xnames,2);
    // fixdimnames(x,dimnm,2);

    if (do_ta)
      for (k = 0; k < np; k++) {
        xanc[k] = k + 1
      }
  }
  // return st
  
  retval[0] = fail
  retval[1] = loglik
  retval[2] = ess
  if (all_fail) {
    retval[3] = x
  } else {
    retval[3] = newstates
  }

  if (all_fail || !do_pr) {
    retval[4] = params
  } else {
    retval[4] = newparams
  }

  if (do_pm) {
    retval[5] = pm
  }
  if (do_pv) {
    retval[6] = pv
  }
  if (do_fm) {
    retval[7] = fm
  }
  if (do_ta) {
    retval[8] = anc
  }

  return retval
}


