
exports.iterateMap = function (object, times, t0, x0, params) {

  let deltat = object.skeletonDetail.deltaT;
  let t = t0;

  let nvars = Object.keys(x0[0]).length;
  let nreps = x0.length;

  let npars = Object.keys(params[0]).length;

  if (nreps !== params.length)
    throw new Error("in 'trajectory': dimension mismatch between 'x0' and 'params'"); // # nocov

  
  let ntimes = times.length;

  let Snames = Object.keys(x0[0]);
  let Pnames = Object.keys(params[0]);
  let Cnames = object.covarnames;

  let zeronames = object.zeronames;
  let nzeros = zeronames.length;
  
  let X = [];
  for(let i = 0; i < ntimes; i++){
    a = [{...(x0[0])}];
    X.push(a);
  }
  // set up the computations
  let ff = object.skeleton;
  let XX = iterate_map_native(times, params, deltat, t, x0, nreps, ff, object);

  return XX;
}

const iterate_map_native = function (times, p, deltat, t, x, nreps, ff, object) {
  let ntimes = times.length;
  let X = new Array(ntimes).fill(0).map( a => []);
  for (let i = 0; i < ntimes; i++) {
    for (let j = 0; j < nreps; j++) {
     X[i].push({});
    }
  }
  let args = object.globals;
  for (let k = 0; k < ntimes; k++) { //loop over object.times
  
    for (let i = 0; i < object.zeronames.length; i++) {
      for (j = 0; j < nreps; j++) x[j][object.zeronames[i]] = 0; // zeroStates=0 at the beginning of each time
    }

    let nsteps = numMapSteps(t, times[k], deltat);
    for (let h = 0; h < nsteps; h++) {
      let interpolatorObj = object.interpolator(t);
      let Xs = new Array(nreps);
      for (let j = 0; j < nreps; j++) {
        Xs[j] = ff(x[j], p[j], t, deltat, interpolatorObj,args);
      }
      x = Xs;
      t += deltat;
    }
    for (let j = 0; j < nreps; j++) {
      Object.assign(X[k][j],x[j])
    }
    if (nsteps === 0) X[k] = [].concat(x); 
  }
  return X;
}

const numMapSteps = function (t1, t2, dt) {
  if (typeof progress === 'function') progress();
  let DOUBLE_EPS = 10e-8
  let tol = Math.sqrt(DOUBLE_EPS)
  let nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}