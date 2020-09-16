
/**
 *  @file        MathLib.js
 *               Library of depentent functions.  
 *
 *  @author       Nazila Akhavan, nazila@kingsds.network
 *  @date        March 2019
 */

let mathLib = {}

let erf = require('math-erf');
let { rbinom } = require('./rbinom');
let libUnif = require('lib-r-math.js');
const { rng: { MersenneTwister, timeseed }} = libUnif
let U = new MersenneTwister(0);

const {
  Gamma,
  rng: {
    LecuyerCMRG,
    normal: { BoxMuller }
  }
} = libUnif
const lc = new LecuyerCMRG(1234);
const { rgamma } = Gamma(new BoxMuller(lc));

mathLib.rnorm = function (mu = 0, sd = 1) {
  let val = Math.sqrt(-2 * Math.log(rng())) * Math.cos(2 * pi * rng());
  return val * sd + mu
}

// Normal  distribution function
mathLib.pnorm = function (x, mu = 0, sd = 1, lower_tail = true, give_log = false) {
  if (sd < 0) {
    return NaN;
  }
  let ans = 1 / 2 * (1 + erf((x - mu) / sd / Math.sqrt(2)));
  if (!lower_tail) {
    ans = 1 - ans;
  }
  if (give_log) {
    ans = Math.log(ans);
  }
  return ans;
}

/** 
 * Density function for the Poisson distribution with parameter lambda.
 *  x : vector of (non-negative integer) quantiles.
 *  lambda : vector of (non-negative) means.
 *  give_log : output is in log scale.
 */
mathLib.dpois = function (x, lambda, give_log = 1) {
  let ans, logAns, total = 0;
  if (isNaN(x) || isNaN(lambda) || lambda < 0) {
    return NaN;
  }
  if (!Number.isInteger(x)) {
    return 0;
  }
  if (x < 0 || !isFinite(x)) {
    return 0;
  }
  x = Math.round(x);
  ans = -lambda + x * Math.log(lambda);
  for (let i = 1; i <= x; i++) {
    total += Math.log(i);
  }
  logAns = ans - total;
  logAns = (give_log) ? logAns : Math.exp(logAns);
  return logAns;
}

// Resampling function
mathLib.nosortResamp = function (nw, w, np, p, offset) {
  for (j = 1; j < nw; j++) {
   w[j] += w[j-1];
 }
  if (w[nw - 1] <= 0) {
    throw "in 'systematic_resampling': non-positive sum of weight"
  }
  let du = w[nw - 1] / np
  let u = -du * U.unif_rand()//Math.random()
  let i = 0;
  for (let j = 0; j < np; j++) {
    u += du;
    while ((u > w[i]) && (i < nw - 1)) i++;//looking for the low weight
    p[j] = i;
  }
  
  if (offset){// add offset if needed
    for (j = 0; j < np; j++) p[j] += offset;
  }
}

// The Euler-multinomial distributions
mathLib.reulermultinom = function (m = 1, size, rateAdd, dt, transAdd, rate, trans) {
  let p = 0
  let j, k
  if ((size < 0) || (dt < 0) || (Math.floor(size + 0.5) !== size)) {
    for (k = 0; k < m; k++) trans[k + transAdd] = NaN
    return 0
  }
  for (k = 0; k < m; k++) {
    if (rate[k + rateAdd] < 0.0) {
      for (j = 0; j < m; j++) trans[j + transAdd] = NaN
      return 0
    }
    p += rate[k + rateAdd]// total event rate
  }
  if (p > 0) {
    size = rbinom(size, 1 - Math.exp(-p * dt))// total number of events
    if (!(isFinite(size)))
      throw 'result of binomial draw is not finite.'
    m -= 1
    for (k = 0; k < m; k++) {
      if (rate[k + rateAdd] > p) p = rate[k + rateAdd]
      trans[k + transAdd] = ((size > 0) && (p > 0)) ? rbinom(size, rate[k + rateAdd] / p) : 0
      if (!(isFinite(size) && isFinite(p) && isFinite(rate[k + rateAdd]) && isFinite(trans[k + transAdd]))) {
        throw 'result of binomial draw is not finite.'
      }
      size -= trans[k + transAdd]
      p -= rate[k + rateAdd]
    }
    trans[m + transAdd] = size
  } else {
    for (k = 0; k < m; k++) trans[k + transAdd] = 0
  }
}

mathLib.sign = function (x, signal) {
  if (isNaN(x))
      return x
  return signal ? Math.abs(x) : -Math.abs(x);
}


mathLib.expRand = function (uniformRand) {
    let q = [
        0.6931471805599453,
        0.9333736875190459,
        0.9888777961838675,
        0.9984959252914960,
        0.9998292811061389,
        0.9999833164100727,
        0.9999985691438767,
        0.9999998906925558,
        0.9999999924734159,
        0.9999999995283275,
        0.9999999999728814,
        0.9999999999985598,
        0.9999999999999289,
        0.9999999999999968,
        0.9999999999999999,
        1.0000000000000000
    ];
    let a = 0.;
    let u = uniformRand();
    while (u <= 0. || u >= 1.)
        u = uniformRand();
    while (true) {
        u += u;
        if (u > 1.)
            break;
        a += q[0];
    }
    u -= 1.;
    if (u <= q[0])
      return a + u
    let i = 0;
    let ustar = uniformRand();
    let umin = ustar;
    do {
        ustar = uniformRand();
        if (umin > ustar)
            umin = ustar;
        i++;
    } while (u > q[i]);
    return a + umin * q[0];
}

mathLib.fromLogBarycentric = function (xN) {
  var sum = 0;
  for (let i = 0; i < xN.length; i++) {
    xN[i] = Math.exp(xN[i]);
    sum += xN[i];
  }
  for (let i = 0; i < xN.length; i++) {
    xN[i] = xN[i] / sum;
  }
  return xN;
}

mathLib.toLogBarycentric = function (xN) {
  let sum = 0;
  for (let i = 0; i < xN.length; i++) {
    sum += xN[i];
  }
  for (let i = 0; i < xN.length; i++) {
    xN[i] = Math.log(xN[i] / sum);
  }
  return xN;
}

mathLib.logit = function (p) {
  return Math.log(p / (1 - p));
}

mathLib.plogis = function (x) {
  return (1+ Math.tanh(x/2))/2;
}

mathLib.expit = function (x) {
  return 1 / (1 + Math.exp(-x));
}

mathLib.qlogis = mathLib.logit;

/**
 * Gamma white noise process with intensity sigma. to find more about rgamma(n, shape, _scale, norm) check the following link.
 * https://github.com/R-js/libRmath.js/blob/306df9aa9a17c323a0dfd7418a4c859270606f7c/src/lib/gamma/index.ts#L142
 * 
 */
const _rgamma = function (n, shape, rate, scale) {
  let _scale = gammaNormalizeParams(rate, scale);
  if (_scale !== undefined) {
      return rgamma(n, shape, _scale);
  }
  console.warn('Cannot normalize to [scale]');
}

const gammaNormalizeParams = function (rate, scale) {
  if (scale === undefined && rate === undefined) {
    return 1;
  }
  if (scale !== undefined && rate !== undefined) {
    if (Math.abs(scale * rate - 1) >= 1e-16) {
      console.warn('Both scale:%d and rate:%d are defined but scale <> 1/rate');
      return undefined;
    }
    return scale;
  }
  if (scale !== undefined && rate === undefined) {
    return scale;
  }
  if (scale === undefined && rate !== undefined) {
    return 1 / rate;
  }
  throw new Error('unreachable code, you cant be here!');
}

mathLib.rgammawn = function (sigma, dt) {
  let sigmasq = Math.pow(sigma, 2);
  return (sigmasq > 0) ?  _rgamma(1, dt / sigmasq, sigmasq) : dt;
}

/**
 * The Log-Mean-Exp Trick avoiding over- and under-flow in doing so. It can optionally return an estimate of the standard error in this quantity.
 */
mathLib.logMeanExp = function (x) {
  var mx = Math.max(...x)
  var s = x.map((x,i) => Math.exp(x - mx))
  var q = s.reduce((a, b) => a + b, 0)
  return mx + Math.log(q / x.length)
}

mathLib.mean = function (x , w = 0) {
  let nrow = x.length;
  let ncol = Object.keys(x[0]).length;
  let mean = {};
  let temp = 0;
  if (w === 0)
    w = new Array(nrow).fill(1);
  let sumw = w.reduce((a,b) => a+b, 0);
  for (let i = 0; i < ncol; i ++) {
    temp = 0;
    for (let j = 0; j < nrow; j++) {
      if (w[j] !== 0) temp += Object.values(x[j])[i] * w[j];
    }
    mean[Object.keys(x[0])[i]] = temp / sumw; 
  }  
  return mean;
}

module.exports = mathLib;
