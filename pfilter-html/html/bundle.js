require=(function(){function r(e,n,t){function o(i,f){if(!n[i]){if(!e[i]){var c="function"==typeof require&&require;if(!f&&c)return c(i,!0);if(u)return u(i,!0);var a=new Error("Cannot find module '"+i+"'");throw a.code="MODULE_NOT_FOUND",a}var p=n[i]={exports:{}};e[i][0].call(p.exports,function(r){var n=e[i][1][r];return o(n||r)},p,p.exports,r,e,n,t)}return n[i].exports}for(var u="function"==typeof require&&require,i=0;i<t.length;i++)o(t[i]);return o}return r})()({1:[function(require,module,exports){

var mathLib = {}
let exp = 2.718281828
let pi = 3.141592654
var erf = require('math-erf')
var rbinom = require('./rbinom')

mathLib.pnorm = function (x, mu = 0, sd = 1, lower_tail = true, give_log = false) {
  if (sd < 0) {
    return NaN
  }
  let ans = 1 / 2 * (1 + erf((x - mu) / sd / Math.sqrt(2)))
  if (!lower_tail) {
    ans = 1 - ans
  }
  if (give_log) {
    ans = Math.log(ans)
  }
  return ans
}


mathLib.numEulerSteps = function(t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  if (t1 >= t2) {
    dt = 0
    nstep = 0
  } else if (t1 + dt >= t2) {
    dt = t2 - t1
    nstep = 1
  } else {
    nstep = Math.ceil((t2 - t1) / dt /(1 + tol))
    dt = (t2 - t1) / nstep
  }
  return nstep
}

mathLib.numMapSteps = function (t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}

mathLib.nosortResamp = function (nw, w, np, p, offset) {
  for (j = 1; j < nw; j++) {
   w[j] += w[j-1]
 }
  if (w[nw - 1] <= 0) {
    throw "in 'systematic_resampling': non-positive sum of weight"
  }
  var du = w[nw - 1] / np
  var u = -du * Math.random()// U.unif_rand()

  for (j = 0, i = 0; j < np; j++) {
    u += du
    while ((u > w[i]) && (i < nw - 1)) i++;//looking for the low weight
    p[j] = i
  }
  if (offset){// add offset if needed
    for (j = 0; j < np; j++) p[j] += offset
  }
}


mathLib.reulermultinom = function (m = 1, size, rateAdd, dt, transAdd, rate, trans) {
  var p = 0
  var j, k
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
    size = rbinom.rbinomOne(size, 1 - Math.exp(-p * dt))// total number of events
    if (!(isFinite(size)))
      throw 'result of binomial draw is not finite.'
    m -= 1
    for (k = 0; k < m; k++) {
      if (rate[k + rateAdd] > p) p = rate[k + rateAdd]
      trans[k + transAdd] = ((size > 0) && (p > 0)) ? rbinom.rbinomOne(size, rate[k + rateAdd] / p) : 0
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

mathLib.rpois = function (lambda = 1) {
  var k = 0; p = 1; l= Math.exp(-lambda)
  while (p > l) { 
    k += 1
    p = p * Math.random()
  }
  return k-1
}

mathLib.interpolator = function (points) {
  var first, n = points.length - 1,
    interpolated,
    leftExtrapolated,
    rightExtrapolated;

  if (points.length === 0) {
    return function () {
      return 0
    }
  }

  if (points.length === 1) {
    return function () {
      return points[0][1]
    }
  }

  points = points.sort(function (a, b) {
    return a[0] - b[0]
  })
  first = points[0]

  leftExtrapolated = function (x) {
    var a = points[0], b = points[1];
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  interpolated = function (x, a, b) {
    return a[1] + (x - a[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  rightExtrapolated = function (x) {
    var a = points[n - 1], b = points[n];
    return b[1] + (x - b[0]) * (b[1] - a[1]) / (b[0] - a[0])
  }

  return function (x) {
    var i
    if (x <= first[0]) {
      return leftExtrapolated(x)
    }
    for (i = 0; i < n; i += 1) {
      if (x > points[i][0] && x <= points[i + 1][0]) {
        return interpolated(x, points[i], points[i + 1])
      }
    }
    return rightExtrapolated(x);
  }
}

module.exports = mathLib;

},{"./rbinom":3,"math-erf":8}],2:[function(require,module,exports){

snippet = {}
let mathLib = require('./mathLib')
//* Set the seed for rnorm-In R:RNGkind("L'Ecuyer-CMRG", normal.kind="Box-Muller");set.seed(1234) 
// const libR = require('lib-r-math.js')
// const {
//   Poisson,
//   rng: { MersenneTwister },
//   rng: { normal: { Inversion } }
// } = libR
// const mt = new MersenneTwister(0)// 
// const { rpois } = Poisson(new Inversion(mt))
// mt.init(1234)
snippet.rprocess = function (params, t, del_t, [S,E,I,R,H], pop, births) {
  var seas, beta, beta0, foi, R0, tt, va
  var trans = new Array(6).fill(0)
  var rate = new Array(6) 
  var deltaT = 14 / 365.25
  var dt = 1 / 365.25 
  
  R0 = params[0], amplitude = params[1], gamma = params[2], mu = params[3], sigma = params[4] 
  beta0 = R0 * (gamma + mu) * (sigma + mu) / sigma
  
  va = snippet.rprocessVaccine(t)
  tt = (t - Math.floor(t)) * 365.25
  if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
    seas = 1 + amplitude * 0.2411 / 0.7589
  } else {
    seas = 1 - amplitude
  }                 
  beta = R0 * (gamma + mu) * (sigma + mu) * seas / sigma //seasonal transmission rate
  foi = beta * I / pop
  rate[0] = foi//force of infection
  rate[1] = mu// natural S death
  rate[2] = sigma// rate of ending of latent stage
  rate[3] = mu// natural E death
  rate[4] = gamma// recovery
  rate[5] = mu// natural I death 
   
  // births = mathLib.rpois(birthrate * (1 - va) * del_t )// Poisson births
  mathLib.reulermultinom(2, Math.round(S), 0, del_t, 0, rate, trans)
  mathLib.reulermultinom(2, Math.round(E), 2, del_t, 2, rate, trans)
  mathLib.reulermultinom(2, Math.round(I), 4, del_t, 4, rate, trans)
  S += (births - trans[0] - trans[1])
  E += (trans[0] - trans[2] - trans[3]) 
  I += (trans[2] - trans[4] - trans[5]) 
  R = pop - S - E - I
  H += trans[4] 
  return [S, E, I, R, H]
}
snippet.rprocessVaccine = function(t) {
  var vaccineRate
  if (t < 1968)
    vaccineRate = 0
  else if (t >= 1968 && t <= 1969)
    vaccineRate = 0.33
  else if (t >= 1969 && t <= 1970)
    vaccineRate = 0.46
  else if (t >= 1970 && t <= 1971)
    vaccineRate = 0.51
  else if (t >= 1971 && t <= 1972)
    vaccineRate = 0.53
  else if (t >= 1972 && t <= 1973)
    vaccineRate = 0.52
  else if (t >= 1973 && t <= 1974)
    vaccineRate = 0.46
  else if (t >= 1974 && t <= 1975)
    vaccineRate = 0.46
  else if (t >= 1975 && t <= 1976)
    vaccineRate = 0.48
  else if (t >= 1976 && t <= 1977)
    vaccineRate = 0.48
  else if (t >= 1977 && t <= 1978)
    vaccineRate = 0.51
  else if (t >= 1978 && t <= 1979)
    vaccineRate = 0.53;
  else if (t >= 1979 && t <= 1980)
    vaccineRate = 0.55;
  else if (t >= 1980 && t <= 1981)
    vaccineRate = 0.58;
  else if (t >= 1981 && t <= 1982)
    vaccineRate = 0.60
  else if (t >= 1982 && t <= 1983)
    vaccineRate = 0.63
  else if (t >= 1983 && t <= 1984)
    vaccineRate = 0.68
  else if (t >= 1984 && t <= 1985)
    vaccineRate = 0.71
  else if (t >= 1985 && t <= 1988)
    vaccineRate = 0.76
  else if (t >= 1988 && t <= 1989)
    vaccineRate = 0.814
  else if (t >= 1989 && t <= 1990)
    vaccineRate = 0.9488
  else if (t >= 1990 && t <= 1991)
    vaccineRate = 0.9818
  else if (t >= 1991 && t <= 1992)
    vaccineRate = 0.90
  else if (t >= 1992 && t <= 1993)
    vaccineRate = 0.92
  else if (t >= 1993 && t <= 1994)
    vaccineRate = 0.91
  else if (t >= 1994 && t <= 1995)
    vaccineRate = 0.91
  else if (t >= 1995 && t <= 1996)
    vaccineRate = 0.92
  else if (t >= 1996 && t <= 1997)
    vaccineRate = 0.92
  else if (t >= 1997 && t <= 1998)
    vaccineRate = 0.91
  else if (t >= 1998 && t <= 1999)
    vaccineRate = 0.88
  else if (t >= 1999 && t <= 2000)
    vaccineRate = 0.88
  else if (t >= 2000 && t <= 2001)
    vaccineRate = 0.87
  else if (t >= 2001 && t <= 2002)
    vaccineRate = 0.84
  else if (t >= 2002 && t <= 2003)
    vaccineRate = 0.82
  else if (t >= 2003 && t <= 2004)
    vaccineRate = 0.80
  else if (t >= 2004 && t <= 2005)
    vaccineRate = 0.81
  else if (t >= 2005 && t <= 2006)
    vaccineRate = 0.84
  else if (t >= 2006 && t <= 2007)
    vaccineRate = 0.85
  else if (t >= 2007 && t <= 2008)
    vaccineRate = 0.85
  else if (t >= 2008 && t <= 2009)
    vaccineRate = 0.85
  else if (t >= 2009 && t <= 2010)
    vaccineRate = 0.88
  else
    vaccineRate = 0.89
  return vaccineRate
}

snippet.initz = function(pop, S_0, E_0, I_0, R_0, H) {
  var m = pop / (S_0 + E_0 + R_0 + I_0),
    S = Math.round(m * S_0),
    E = Math.round(m * E_0),
    I = Math.round(m * I_0),
    R = Math.round(m * R_0),
    H = 0
  return [S, E, I, R, H]
}

snippet.dmeasure = function (rho, psi, H, dCases, giveLog = 1) {
  var lik
  var mn = rho * H
  var v = mn * (1.0 - rho + psi * psi * mn)
  var tol = 1.0e-18
  var modelCases = Number(dCases)
  if(!isNaN(modelCases)){
    if (modelCases > 0.0) {
      lik = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
    } else {
      lik = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol)) + tol
    }
    if (giveLog) lik = Math.log(lik)
  } else {
    lik = (giveLog) ? 0 : 1;
  }
  return lik
}

snippet.rmeasure = function (H, rho, psi) {
  var mn = rho * H
  var v = mn * (1.0 - rho + psi * psi * mn)
  var tol = 1.0e-18
  var cases = mathLib.rnorm(mn, Math.sqrt(v) + tol)
  if (cases > 0) {
    cases = Math.round(cases)
  } else {
    cases = 0
  }
  return cases
}

module.exports = snippet

},{"./mathLib":1}],3:[function(require,module,exports){
/*
 * References              Jacob K.F. Bogers  info@mail.jacob-bogers.com
 *                         https://github.com/R-js/libRmath.js/blob/master/src/lib/binomial/rbinom.ts
 */
rbinom = {}
rbinom. rbinomOne = function (size, pp) {
    var c = 0;
    var fm = 0;
    var npq = 0;
    var p1 = 0;
    var p2 = 0;
    var p3 = 0;
    var p4 = 0;
    var qn = 0;
    var xl = 0;
    var xll = 0;
    var xlr = 0;
    var xm = 0;
    var xr = 0;
    var psave = -1.0;
    var nsave = -1;
    var m = 0;
    var f;
    var f1;
    var f2;
    var u;
    var v;
    var w;
    var w2;
    var x;
    var x1;
    var x2;
    var z;
    var z2;
    var p;
    var q;
    var np;
    var g;
    var r;
    var al;
    var alv;
    var amaxp;
    var ffm;
    var ynorm;
    var i;
    var ix = 0;
    var k;
    var n;
    if (!isFinite(size)){
        throw "Input values should be finite"
    }
    r = Math.round(size)
    if (r !== size)
        return NaN
    if (!isFinite(pp) ||
        r < 0 ||
        pp < 0 ||
        pp > 1) {
        return NaN
    }
    if (r === 0 || pp === 0)
        return 0;
    if (pp === 1)
        return r;
    if (r >= Number.MAX_SAFE_INTEGER) {
        throw 'Evade overflow' + r + '> MAX_SAFE_INTEGER'
    }
    n = Math.trunc(r);
    p = Math.min(pp, 1 - pp);
    q = 1 - p;
    np = n * p;
    r = p / q;
    g = r * (n + 1);
    var gotoL_np_small = false;
    if (pp !== psave || n !== nsave) {
        psave = pp;
        nsave = n;
        if (np < 30.0) {
            qn = rbinom.R_pow_di(q, n);
            gotoL_np_small = true;
        } else {
            ffm = np + p;
            m = Math.trunc(ffm);
            fm = m;
            npq = np * q;
            p1 = Math.trunc(2.195 * Math.sqrt(npq) - 4.6 * q) + 0.5;
            xm = fm + 0.5;
            xl = xm - p1;
            xr = xm + p1;
            c = 0.134 + 20.5 / (15.3 + fm);
            al = (ffm - xl) / (ffm - xl * p);
            xll = al * (1.0 + 0.5 * al);
            al = (xr - ffm) / (xr * q);
            xlr = al * (1.0 + 0.5 * al);
            p2 = p1 * (1.0 + c + c);
            p3 = p2 + c / xll;
            p4 = p3 + c / xlr;
        }
    } else if (n === nsave) {
        if (np < 30.0)
            gotoL_np_small = true;
    }
    var gotoFinis = false;
    while (true && !gotoL_np_small) {
        u = Math.random() * p4;
        v = Math.random();
        if (u <= p1) {
            ix = Math.trunc(xm - p1 * v + u);
            gotoFinis = true;
            break;
        }
        if (u <= p2) {
            x = xl + (u - p1) / c;
            v = v * c + 1.0 - Math.abs(xm - x) / p1;
            if (v > 1.0 || v <= 0)
                continue;
            ix = Math.trunc(x);
        } else {
            if (u > p3) {
                ix = Math.trunc(xr - Math.log(v) / xlr);
                if (ix > n)
                    continue;
                v = v * (u - p3) * xlr;
            } else {
                ix = Math.trunc(xl + Math.log(v) / xll);
                if (ix < 0)
                    continue;
                v = v * (u - p2) * xll;
            }
        }
        k = Math.abs(ix - m);
        if (k <= 20 || k >= npq / 2 - 1) {
            f = 1.0;
            if (m < ix) {
                for (i = m + 1; i <= ix; i++)
                    f *= g / i - r;
            } else if (m !== ix) {
                for (i = ix + 1; i <= m; i++)
                    f /= g / i - r;
            }
            if (v <= f) {
                gotoFinis = true;
                break
            }
        } else {
            amaxp = k / npq * ((k * (k / 3 + 0.625) + 0.1666666666666) / npq + 0.5);
            ynorm = -k * k / (2.0 * npq);
            alv = Math.log(v);
            if (alv < ynorm - amaxp) {
                gotoFinis = true;
                break;
            }
            if (alv <= ynorm + amaxp) {
                x1 = ix + 1;
                f1 = fm + 1.0;
                z = n + 1 - fm;
                w = n - ix + 1.0;
                z2 = z * z;
                x2 = x1 * x1;
                f2 = f1 * f1;
                w2 = w * w;
                if (alv <=
                    xm * Math.log(f1 / x1) +
                        (n - m + 0.5) * Math.log(z / w) +
                        (ix - m) * Math.log(w * p / (x1 * q)) +
                        (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / f2) / f2) / f2) / f2) /
                            f1 /
                            166320.0 +
                        (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / z2) / z2) / z2) / z2) /
                            z /
                            166320.0 +
                        (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / x2) / x2) / x2) / x2) /
                            x1 /
                            166320.0 +
                        (13860.0 - (462.0 - (132.0 - (99.0 - 140.0 / w2) / w2) / w2) / w2) /
                            w /
                            166320) {
                    gotoFinis = true;
                    break;
                }
            }
        }
    }
    if (!gotoFinis) {
        while (true) {
            ix = 0;
            f = qn;
            u = Math.random()//Math.random();
            while (true) {
                if (u < f) {
                    gotoFinis = true;
                    break;
                }
                if (ix > 110)
                    break;
                u -= f;
                ix++;
                f *= g / ix - r;
            }
            if (gotoFinis) {
                break;
            }
        }
    }
    if (psave > 0.5) {
        ix = n - ix;
    }
    return ix;
}

rbinom.R_pow_di = function (x, n) {
    var pow = 1.0;
    if (Number.isNaN(x))
        return x;
    if (n !== 0) {
        if (!Number.isFinite(x))
            return R_pow(x, n);
        if (n < 0) {
            n = -n;
            x = 1 / x;
        }
        while (true) {
            if (n & 1)
                pow *= x;
            if ((n >>= 1))
                x *= x;
            else
                break;
        }
    }
    return pow;
}
module.exports = rbinom

},{}],4:[function(require,module,exports){

let mathLib = require('./mathLib.js')
let snippet = require('./modelSnippet.js')

simulator = {}
simulator.simulate = function (particles, k, tdata, deltaT, dt, timeCountData, interpolPop, interpolBirth,params, t1,t2 ) {
  var st, S, E, I, R, H 
  S = particles[0]; E = particles[1]; I = particles[2]; R = particles[3]; H = particles[4]
  // transitions between classes
  if (k <= tdata || k > t1) {
    steps = mathLib.numMapSteps(k, k + deltaT, dt)
  } else {
    steps = mathLib.numEulerSteps(k, t2, dt)
  }
  del_t = (1 / steps )* deltaT 
  for (let stp = 0; stp < steps; stp++) { // steps in each time interval
    st = k + stp * del_t
    simulateValue = snippet.rprocess(params, st, del_t, [S,E,I,R,H], interpolPop(st), interpolBirth(st))
    S = simulateValue[0]; E = simulateValue[1]; I = simulateValue[2]; R = simulateValue[3]; H = simulateValue[4]
  }
  return [S, E, I, R, H]
}

module.exports = simulator
},{"./mathLib.js":1,"./modelSnippet.js":2}],5:[function(require,module,exports){
'use strict';

// EXPORTS //

module.exports = Number.NEGATIVE_INFINITY;

},{}],6:[function(require,module,exports){
'use strict';

// EXPORTS //

module.exports = Number.POSITIVE_INFINITY;

},{}],7:[function(require,module,exports){
(function (global, factory) {
    typeof exports === 'object' && typeof module !== 'undefined' ? factory(exports) :
    typeof define === 'function' && define.amd ? define(['exports'], factory) :
    (factory((global.fmin = global.fmin || {})));
}(this, function (exports) { 'use strict';

    /** finds the zeros of a function, given two starting points (which must
     * have opposite signs */
    function bisect(f, a, b, parameters) {
        parameters = parameters || {};
        var maxIterations = parameters.maxIterations || 100,
            tolerance = parameters.tolerance || 1e-10,
            fA = f(a),
            fB = f(b),
            delta = b - a;

        if (fA * fB > 0) {
            throw "Initial bisect points must have opposite signs";
        }

        if (fA === 0) return a;
        if (fB === 0) return b;

        for (var i = 0; i < maxIterations; ++i) {
            delta /= 2;
            var mid = a + delta,
                fMid = f(mid);

            if (fMid * fA >= 0) {
                a = mid;
            }

            if ((Math.abs(delta) < tolerance) || (fMid === 0)) {
                return mid;
            }
        }
        return a + delta;
    }

    // need some basic operations on vectors, rather than adding a dependency,
    // just define here
    function zeros(x) { var r = new Array(x); for (var i = 0; i < x; ++i) { r[i] = 0; } return r; }
    function zerosM(x,y) { return zeros(x).map(function() { return zeros(y); }); }

    function dot(a, b) {
        var ret = 0;
        for (var i = 0; i < a.length; ++i) {
            ret += a[i] * b[i];
        }
        return ret;
    }

    function norm2(a)  {
        return Math.sqrt(dot(a, a));
    }

    function scale(ret, value, c) {
        for (var i = 0; i < value.length; ++i) {
            ret[i] = value[i] * c;
        }
    }

    function weightedSum(ret, w1, v1, w2, v2) {
        for (var j = 0; j < ret.length; ++j) {
            ret[j] = w1 * v1[j] + w2 * v2[j];
        }
    }

    /** minimizes a function using the downhill simplex method */
    function nelderMead(f, x0, parameters) {
        parameters = parameters || {};

        var maxIterations = parameters.maxIterations || x0.length * 200,
            nonZeroDelta = parameters.nonZeroDelta || 1.05,
            zeroDelta = parameters.zeroDelta || 0.001,
            minErrorDelta = parameters.minErrorDelta || 1e-6,
            minTolerance = parameters.minErrorDelta || 1e-5,
            rho = (parameters.rho !== undefined) ? parameters.rho : 1,
            chi = (parameters.chi !== undefined) ? parameters.chi : 2,
            psi = (parameters.psi !== undefined) ? parameters.psi : -0.5,
            sigma = (parameters.sigma !== undefined) ? parameters.sigma : 0.5,
            maxDiff;

        // initialize simplex.
        var N = x0.length,
            simplex = new Array(N + 1);
        simplex[0] = x0;
        simplex[0].fx = f(x0);
        simplex[0].id = 0;
        for (var i = 0; i < N; ++i) {
            var point = x0.slice();
            point[i] = point[i] ? point[i] * nonZeroDelta : zeroDelta;
            simplex[i+1] = point;
            simplex[i+1].fx = f(point);
            simplex[i+1].id = i+1;
        }

        function updateSimplex(value) {
            for (var i = 0; i < value.length; i++) {
                simplex[N][i] = value[i];
            }
            simplex[N].fx = value.fx;
        }

        var sortOrder = function(a, b) { return a.fx - b.fx; };

        var centroid = x0.slice(),
            reflected = x0.slice(),
            contracted = x0.slice(),
            expanded = x0.slice();

        for (var iteration = 0; iteration < maxIterations; ++iteration) {
            simplex.sort(sortOrder);

            if (parameters.history) {
                // copy the simplex (since later iterations will mutate) and
                // sort it to have a consistent order between iterations
                var sortedSimplex = simplex.map(function (x) {
                    var state = x.slice();
                    state.fx = x.fx;
                    state.id = x.id;
                    return state;
                });
                sortedSimplex.sort(function(a,b) { return a.id - b.id; });

                parameters.history.push({x: simplex[0].slice(),
                                         fx: simplex[0].fx,
                                         simplex: sortedSimplex});
            }

            maxDiff = 0;
            for (i = 0; i < N; ++i) {
                maxDiff = Math.max(maxDiff, Math.abs(simplex[0][i] - simplex[1][i]));
            }

            if ((Math.abs(simplex[0].fx - simplex[N].fx) < minErrorDelta) &&
                (maxDiff < minTolerance)) {
                break;
            }

            // compute the centroid of all but the worst point in the simplex
            for (i = 0; i < N; ++i) {
                centroid[i] = 0;
                for (var j = 0; j < N; ++j) {
                    centroid[i] += simplex[j][i];
                }
                centroid[i] /= N;
            }

            // reflect the worst point past the centroid  and compute loss at reflected
            // point
            var worst = simplex[N];
            weightedSum(reflected, 1+rho, centroid, -rho, worst);
            reflected.fx = f(reflected);

            // if the reflected point is the best seen, then possibly expand
            if (reflected.fx < simplex[0].fx) {
                weightedSum(expanded, 1+chi, centroid, -chi, worst);
                expanded.fx = f(expanded);
                if (expanded.fx < reflected.fx) {
                    updateSimplex(expanded);
                }  else {
                    updateSimplex(reflected);
                }
            }

            // if the reflected point is worse than the second worst, we need to
            // contract
            else if (reflected.fx >= simplex[N-1].fx) {
                var shouldReduce = false;

                if (reflected.fx > worst.fx) {
                    // do an inside contraction
                    weightedSum(contracted, 1+psi, centroid, -psi, worst);
                    contracted.fx = f(contracted);
                    if (contracted.fx < worst.fx) {
                        updateSimplex(contracted);
                    } else {
                        shouldReduce = true;
                    }
                } else {
                    // do an outside contraction
                    weightedSum(contracted, 1-psi * rho, centroid, psi*rho, worst);
                    contracted.fx = f(contracted);
                    if (contracted.fx < reflected.fx) {
                        updateSimplex(contracted);
                    } else {
                        shouldReduce = true;
                    }
                }

                if (shouldReduce) {
                    // if we don't contract here, we're done
                    if (sigma >= 1) break;

                    // do a reduction
                    for (i = 1; i < simplex.length; ++i) {
                        weightedSum(simplex[i], 1 - sigma, simplex[0], sigma, simplex[i]);
                        simplex[i].fx = f(simplex[i]);
                    }
                }
            } else {
                updateSimplex(reflected);
            }
        }

        simplex.sort(sortOrder);
        return {fx : simplex[0].fx,
                x : simplex[0]};
    }

    /// searches along line 'pk' for a point that satifies the wolfe conditions
    /// See 'Numerical Optimization' by Nocedal and Wright p59-60
    /// f : objective function
    /// pk : search direction
    /// current: object containing current gradient/loss
    /// next: output: contains next gradient/loss
    /// returns a: step size taken
    function wolfeLineSearch(f, pk, current, next, a, c1, c2) {
        var phi0 = current.fx, phiPrime0 = dot(current.fxprime, pk),
            phi = phi0, phi_old = phi0,
            phiPrime = phiPrime0,
            a0 = 0;

        a = a || 1;
        c1 = c1 || 1e-6;
        c2 = c2 || 0.1;

        function zoom(a_lo, a_high, phi_lo) {
            for (var iteration = 0; iteration < 16; ++iteration) {
                a = (a_lo + a_high)/2;
                weightedSum(next.x, 1.0, current.x, a, pk);
                phi = next.fx = f(next.x, next.fxprime);
                phiPrime = dot(next.fxprime, pk);

                if ((phi > (phi0 + c1 * a * phiPrime0)) ||
                    (phi >= phi_lo)) {
                    a_high = a;

                } else  {
                    if (Math.abs(phiPrime) <= -c2 * phiPrime0) {
                        return a;
                    }

                    if (phiPrime * (a_high - a_lo) >=0) {
                        a_high = a_lo;
                    }

                    a_lo = a;
                    phi_lo = phi;
                }
            }

            return 0;
        }

        for (var iteration = 0; iteration < 10; ++iteration) {
            weightedSum(next.x, 1.0, current.x, a, pk);
            phi = next.fx = f(next.x, next.fxprime);
            phiPrime = dot(next.fxprime, pk);
            if ((phi > (phi0 + c1 * a * phiPrime0)) ||
                (iteration && (phi >= phi_old))) {
                return zoom(a0, a, phi_old);
            }

            if (Math.abs(phiPrime) <= -c2 * phiPrime0) {
                return a;
            }

            if (phiPrime >= 0 ) {
                return zoom(a, a0, phi);
            }

            phi_old = phi;
            a0 = a;
            a *= 2;
        }

        return a;
    }

    function conjugateGradient(f, initial, params) {
        // allocate all memory up front here, keep out of the loop for perfomance
        // reasons
        var current = {x: initial.slice(), fx: 0, fxprime: initial.slice()},
            next = {x: initial.slice(), fx: 0, fxprime: initial.slice()},
            yk = initial.slice(),
            pk, temp,
            a = 1,
            maxIterations;

        params = params || {};
        maxIterations = params.maxIterations || initial.length * 20;

        current.fx = f(current.x, current.fxprime);
        pk = current.fxprime.slice();
        scale(pk, current.fxprime,-1);

        for (var i = 0; i < maxIterations; ++i) {
            a = wolfeLineSearch(f, pk, current, next, a);

            // todo: history in wrong spot?
            if (params.history) {
                params.history.push({x: current.x.slice(),
                                     fx: current.fx,
                                     fxprime: current.fxprime.slice(),
                                     alpha: a});
            }

            if (!a) {
                // faiiled to find point that satifies wolfe conditions.
                // reset direction for next iteration
                scale(pk, current.fxprime, -1);

            } else {
                // update direction using Polak–Ribiere CG method
                weightedSum(yk, 1, next.fxprime, -1, current.fxprime);

                var delta_k = dot(current.fxprime, current.fxprime),
                    beta_k = Math.max(0, dot(yk, next.fxprime) / delta_k);

                weightedSum(pk, beta_k, pk, -1, next.fxprime);

                temp = current;
                current = next;
                next = temp;
            }

            if (norm2(current.fxprime) <= 1e-5) {
                break;
            }
        }

        if (params.history) {
            params.history.push({x: current.x.slice(),
                                 fx: current.fx,
                                 fxprime: current.fxprime.slice(),
                                 alpha: a});
        }

        return current;
    }

    function gradientDescent(f, initial, params) {
        params = params || {};
        var maxIterations = params.maxIterations || initial.length * 100,
            learnRate = params.learnRate || 0.001,
            current = {x: initial.slice(), fx: 0, fxprime: initial.slice()};

        for (var i = 0; i < maxIterations; ++i) {
            current.fx = f(current.x, current.fxprime);
            if (params.history) {
                params.history.push({x: current.x.slice(),
                                     fx: current.fx,
                                     fxprime: current.fxprime.slice()});
            }

            weightedSum(current.x, 1, current.x, -learnRate, current.fxprime);
            if (norm2(current.fxprime) <= 1e-5) {
                break;
            }
        }

        return current;
    }

    function gradientDescentLineSearch(f, initial, params) {
        params = params || {};
        var current = {x: initial.slice(), fx: 0, fxprime: initial.slice()},
            next = {x: initial.slice(), fx: 0, fxprime: initial.slice()},
            maxIterations = params.maxIterations || initial.length * 100,
            learnRate = params.learnRate || 1,
            pk = initial.slice(),
            c1 = params.c1 || 1e-3,
            c2 = params.c2 || 0.1,
            temp,
            functionCalls = [];

        if (params.history) {
            // wrap the function call to track linesearch samples
            var inner = f;
            f = function(x, fxprime) {
                functionCalls.push(x.slice());
                return inner(x, fxprime);
            };
        }

        current.fx = f(current.x, current.fxprime);
        for (var i = 0; i < maxIterations; ++i) {
            scale(pk, current.fxprime, -1);
            learnRate = wolfeLineSearch(f, pk, current, next, learnRate, c1, c2);

            if (params.history) {
                params.history.push({x: current.x.slice(),
                                     fx: current.fx,
                                     fxprime: current.fxprime.slice(),
                                     functionCalls: functionCalls,
                                     learnRate: learnRate,
                                     alpha: learnRate});
                functionCalls = [];
            }


            temp = current;
            current = next;
            next = temp;

            if ((learnRate === 0) || (norm2(current.fxprime) < 1e-5)) break;
        }

        return current;
    }

    exports.bisect = bisect;
    exports.nelderMead = nelderMead;
    exports.conjugateGradient = conjugateGradient;
    exports.gradientDescent = gradientDescent;
    exports.gradientDescentLineSearch = gradientDescentLineSearch;
    exports.zeros = zeros;
    exports.zerosM = zerosM;
    exports.norm2 = norm2;
    exports.weightedSum = weightedSum;
    exports.scale = scale;

}));
},{}],8:[function(require,module,exports){
'use strict';

/**
* NOTE: the following copyright and license, as well as the long comment were part of the original implementation available as part of [FreeBSD]{@link https://svnweb.freebsd.org/base/release/9.3.0/lib/msun/src/s_erf.c?revision=268523&view=co}.
*
* The implementation follows the original, but has been modified for JavaScript.
*/

/**
* ====================================================
* Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
*
* Developed at SunPro, a Sun Microsystems, Inc. business.
* Permission to use, copy, modify, and distribute this
* software is freely granted, provided that this notice
* is preserved.
* ====================================================
*/

/**
* double erf(double x)
*                                 x
*                        2       |\
*        erf(x)  = -----------   | exp(-t*t)dt
*                     sqrt(pi)  \|
*                                0
*
*        erfc(x) = 1-erf(x)
*
*   Note that
*
*        erf(-x)  = -erf(x)
*        erfc(-x) = 2 - erfc(x)
*
* Method:
*   1. For |x| in [0, 0.84375),
*
*        erf(x)  = x + x*R(x^2)
*        erfc(x) = 1 - erf(x)           if x in [-.84375,0.25]
*                = 0.5 + ((0.5-x)-x*R)  if x in [0.25,0.84375]
*
*      where R = P/Q where P is an odd polynomial of degree 8 and Q is an odd polynomial of degree 10.
*
*        | R - (erf(x)-x)/x | <= 2**-57.90
*
*      Remark: the formula is derived by noting
*
*        erf(x) = (2/sqrt(pi))*(x - x^3/3 + x^5/10 - x^7/42 + ....)
*
*      and that
*
*        2/sqrt(pi) = 1.128379167095512573896158903121545171688
*
*      is close to one. The interval is chosen because the fix point of erf(x) is near 0.6174 (i.e., erf(x)=x when x is near 0.6174), and, by some experiment, 0.84375 is chosen to guarantee the error is less than one ulp for erf.
*
*   2. For |x| in [0.84375,1.25), let s = |x| - 1, and c = 0.84506291151 rounded to single (24 bits)
*
*        erf(x)  = sign(x) * (c + P1(s)/Q1(s))
*        erfc(x) = (1-c) - P1(s)/Q1(s) if x > 0
*                  1+(c+P1(s)/Q1(s))   if x < 0
*        |P1/Q1 - (erf(|x|)-c)| <= 2**-59.06
*
*      Remark: here we use the Taylor series expansion at x=1.
*
*        erf(1+s) = erf(1) + s*Poly(s)
*                 = 0.845.. + P1(s)/Q1(s)
*
*      That is, we use a rational approximation to approximate
*
*        erf(1+s) - (c = (single)0.84506291151)
*
*      Note that |P1/Q1|< 0.078 for x in [0.84375,1.25] where
*
*        P1(s) = degree 6 poly in s
*        Q1(s) = degree 6 poly in s
*
*   3. For x in [1.25,1/0.35(~2.857143)),
*
*        erfc(x) = (1/x)*exp(-x*x-0.5625+R1/S1)
*        erf(x)  = 1 - erfc(x)
*
*      where
*
*        R1(z) = degree 7 poly in z, (z=1/x^2)
*        S1(z) = degree 8 poly in z
*
*   4. For x in [1/0.35,28],
*
*        erfc(x) = (1/x)*exp(-x*x-0.5625+R2/S2)       if x > 0
*                = 2.0 - (1/x)*exp(-x*x-0.5625+R2/S2) if -6 < x < 0
*                = 2.0 - tiny                         if x <= -6
*        erf(x)  = sign(x)*(1.0 - erfc(x))            if x < 6, else
*        erf(x)  = sign(x)*(1.0 - tiny)
*
*      where
*
*        R2(z) = degree 6 poly in z, (z=1/x^2)
*        S2(z) = degree 7 poly in z
*
*   Note1:
*       To compute exp(-x*x-0.5625+R/S), let s be a single precision number and s := x; then
*
*        -x*x = -s*s + (s-x)*(s+x)
*        exp(-x*x-0.5626+R/S) = exp(-s*s-0.5625)*exp((s-x)*(s+x)+R/S);
*
*   Note2:
*       Here 4 and 5 make use of the asymptotic series
*
*                     exp(-x*x)
*         erfc(x) ~  ----------- * ( 1 + Poly(1/x^2) )
*                     x*sqrt(pi)
*
*       We use a rational approximation to approximate
*
*           g(s) = f(1/x^2) = log(erfc(x)*x) - x*x + 0.5625
*
*       Here is the error bound for R1/S1 and R2/S2
*
*           |R1/S1 - f(x)| < 2**(-62.57)
*           |R2/S2 - f(x)| < 2**(-61.52)
*
*   5. For inf > x >= 28,
*
*        erf(x)  = sign(x) * (1 - tiny)   (raise inexact)
*        erfc(x) = tiny*tiny              (raise underflow) if x > 0
*                = 2 - tiny               if x<0
*
*   6. Special cases:
*
*        erf(0) = 0
*        erf(inf) = 1
*        erf(-inf) = -1
*        erfc(0) = 1
*        erfc(inf) = 0
*        erfc(-inf) = 2,
*        erf(NaN) is NaN
*        erfc(NaN) is NaN
*/

// MODULES //

var evalpoly = require( 'math-evalpoly' ).factory;
var exp = require( 'math-exp' );
var setLowWord = require( 'math-float64-set-low-word' );


// CONSTANTS //

var PINF = require( 'const-pinf-float64' );
var NINF = require( 'const-ninf-float64' );

var TINY = 1e-300;
var VERY_TINY = 2.848094538889218e-306; // 0x00800000, 0x00000000

// 2**-28 = 1/(1<<28) = 1/268435456
var SMALL = 3.725290298461914e-9;

var ERX = 8.45062911510467529297e-1; // 0x3FEB0AC1, 0x60000000

// Coefficients for approximation to erf on [0, 0.84375)
var EFX = 1.28379167095512586316e-1;  // 0x3FC06EBA, 0x8214DB69
var EFX8 = 1.02703333676410069053;    // 0x3FF06EBA, 0x8214DB69
var PPC = 1.28379167095512558561e-1;  // 0x3FC06EBA, 0x8214DB68
var PP = [
	-3.25042107247001499370e-1, // 0xBFD4CD7D, 0x691CB913
	-2.84817495755985104766e-2, // 0xBF9D2A51, 0xDBD7194F
	-5.77027029648944159157e-3, // 0xBF77A291, 0x236668E4
	-2.37630166566501626084e-5  // 0xBEF8EAD6, 0x120016AC
];
var QQC = 1.0;
var QQ = [
	3.97917223959155352819e-1, // 0x3FD97779, 0xCDDADC09
	6.50222499887672944485e-2, // 0x3FB0A54C, 0x5536CEBA
	5.08130628187576562776e-3, // 0x3F74D022, 0xC4D36B0F
	1.32494738004321644526e-4, // 0x3F215DC9, 0x221C1A10
	-3.96022827877536812320e-6 // 0xBED09C43, 0x42A26120
];

// Coefficients for approximation to erf on [0.84375, 1.25)
var PAC = -2.36211856075265944077e-3; // 0xBF6359B8, 0xBEF77538
var PA = [
	4.14856118683748331666e-1,  // 0x3FDA8D00, 0xAD92B34D
	-3.72207876035701323847e-1, // 0xBFD7D240, 0xFBB8C3F1
	3.18346619901161753674e-1,  // 0x3FD45FCA, 0x805120E4
	-1.10894694282396677476e-1, // 0xBFBC6398, 0x3D3E28EC
	3.54783043256182359371e-2,  // 0x3FA22A36, 0x599795EB
	-2.16637559486879084300e-3  // 0xBF61BF38, 0x0A96073F
];
var QAC = 1.0;
var QA = [
	1.06420880400844228286e-1, // 0x3FBB3E66, 0x18EEE323
	5.40397917702171048937e-1, // 0x3FE14AF0, 0x92EB6F33
	7.18286544141962662868e-2, // 0x3FB2635C, 0xD99FE9A7
	1.26171219808761642112e-1, // 0x3FC02660, 0xE763351F
	1.36370839120290507362e-2, // 0x3F8BEDC2, 0x6B51DD1C
	1.19844998467991074170e-2  // 0x3F888B54, 0x5735151D
];

// Coefficients for approximation to erfc on [1.25, 1/0.35)
var RAC = -9.86494403484714822705e-3; // 0xBF843412, 0x600D6435
var RA = [
	-6.93858572707181764372e-1, // 0xBFE63416, 0xE4BA7360
	-1.05586262253232909814e1,  // 0xC0251E04, 0x41B0E726 
	-6.23753324503260060396e1,  // 0xC04F300A, 0xE4CBA38D
	-1.62396669462573470355e2,  // 0xC0644CB1, 0x84282266
	-1.84605092906711035994e2,  // 0xC067135C, 0xEBCCABB2
	-8.12874355063065934246e1,  // 0xC0545265, 0x57E4D2F2
	-9.81432934416914548592     // 0xC023A0EF, 0xC69AC25C
];
var SAC = 1.0;
var SA = [
	1.96512716674392571292e1,  // 0x4033A6B9, 0xBD707687
	1.37657754143519042600e2,  // 0x4061350C, 0x526AE721
	4.34565877475229228821e2,  // 0x407B290D, 0xD58A1A71
	6.45387271733267880336e2,  // 0x40842B19, 0x21EC2868
	4.29008140027567833386e2,  // 0x407AD021, 0x57700314
	1.08635005541779435134e2,  // 0x405B28A3, 0xEE48AE2C
	6.57024977031928170135,    // 0x401A47EF, 0x8E484A93
	-6.04244152148580987438e-2 // 0xBFAEEFF2, 0xEE749A62
];

// Coefficients for approximation to erfc on [1/0.35, 28]
var RBC = -9.86494292470009928597e-3; // 0xBF843412, 0x39E86F4A
var RB = [
	-7.99283237680523006574e-1, // 0xBFE993BA, 0x70C285DE
	-1.77579549177547519889e1,  // 0xC031C209, 0x555F995A
	-1.60636384855821916062e2,  // 0xC064145D, 0x43C5ED98
	-6.37566443368389627722e2,  // 0xC083EC88, 0x1375F228
	-1.02509513161107724954e3,  // 0xC0900461, 0x6A2E5992
	-4.83519191608651397019e2,  // 0xC07E384E, 0x9BDC383F
];
var SBC = 1.0;
var SB = [
	3.03380607434824582924e1, // 0x403E568B, 0x261D5190
	3.25792512996573918826e2, // 0x40745CAE, 0x221B9F0A
	1.53672958608443695994e3, // 0x409802EB, 0x189D5118
	3.19985821950859553908e3, // 0x40A8FFB7, 0x688C246A
	2.55305040643316442583e3, // 0x40A3F219, 0xCEDF3BE6
	4.74528541206955367215e2, // 0x407DA874, 0xE79FE763
	-2.24409524465858183362e1 // 0xC03670E2, 0x42712D62
];


// FUNCTIONS //

// Compile functions to evaluate polynomials based on the above coefficients...
var polyvalPP = evalpoly( PP );
var polyvalQQ = evalpoly( QQ );
var polyvalPA = evalpoly( PA );
var polyvalQA = evalpoly( QA );
var polyvalRA = evalpoly( RA );
var polyvalSA = evalpoly( SA );
var polyvalRB = evalpoly( RB );
var polyvalSB = evalpoly( SB );


// ERF //

/**
* FUNCTION: erf( x )
*	Evaluates the error function.
*
* @param {Number} x - input value
* @returns {Number} evaluated error function
*/
function erf( x ) {
	var sign;
	var ax;
	var z;
	var r;
	var s;
	var y;
	var p;
	var q;

	// Special case: NaN
	if ( x !== x ) {
		return NaN;
	}
	// Special case: +infinity
	if ( x === PINF ) {
		return 1;
	}
	// Special case: -infinity
	if ( x === NINF ) {
		return -1;
	}
	// Special case: +-0
	if ( x === 0 ) {
		return x;
	}
	if ( x < 0 ) {
		sign = true;
		ax = -x;
	} else {
		sign = false;
		ax = x;
	}
	// |x| < 0.84375
	if ( ax < 0.84375 ) {
		if ( ax < SMALL ) {
			if ( ax < VERY_TINY ) {
				// Avoid underflow:
				return 0.125 * (8.0*x + EFX8*x);
			}
			return x + EFX*x;
		}
		z = x * x;
		r = PPC + z*polyvalPP( z );
		s = QQC + z*polyvalQQ( z );
		y = r / s;
		return x + x*y;
	}
	// 0.84375 <= |x| < 1.25
	if ( ax < 1.25 ) {
		s = ax - 1;
		p = PAC + s*polyvalPA( s );
		q = QAC + s*polyvalQA( s );
		if ( sign ) {
			return -ERX - p/q;
		}
		return ERX + p/q;
	}
	// +inf > |x| >= 6
	if ( ax >= 6 ) {
		if ( sign ) {
			return TINY - 1.0; // raise inexact
		}
		return 1.0 - TINY; // raise inexact
	}
	s = 1.0 / (ax*ax);

	// |x| < 1/0.35 ~ 2.857143
	if ( ax < 2.857142857142857 ) {
		r = RAC + s*polyvalRA( s );
		s = SAC + s*polyvalSA( s );
	}
	// |x| >= 1/0.35 ~ 2.857143
	else {
		r = RBC + s*polyvalRB( s );
		s = SBC + s*polyvalSB( s );
	}
	z = setLowWord( ax, 0 ); // pseudo-single (20-bit) precision x
	r = exp( -z*z - 0.5625 ) * exp( (z-ax)*(z+ax) + r/s );
	if ( sign ) {
		return r/ax - 1;
	}
	return 1 - r/ax;
} // end FUNCTION erf()


// EXPORTS //

module.exports = erf;
},{"const-ninf-float64":5,"const-pinf-float64":6,"math-evalpoly":11,"math-exp":12,"math-float64-set-low-word":13}],9:[function(require,module,exports){
'use strict';

// EVALPOLY //

/**
* FUNCTION: evalpoly( c, x )
*	Evaluates a polynomial.
*
* @param {Number[]|Int8Array|Uint8Array|Uint8ClampedArray|Int16Array|Uint16Array|Int32Array|Uint32Array|Float32Array|Float64Array} c - polynomial coefficients sorted in ascending degree
* @param {Number} x - value at which to evaluate the polynomial
* @returns {Number} evaluated polynomial
*/
function evalpoly( c, x ) {
	var p;
	var i;
	
	i = c.length;
	if ( i < 2 || x === 0 ) {
		if ( i === 0 ) {
			return 0;
		}
		return c[ 0 ];
	}
	i -= 1;

	// Use Horner's rule (http://en.wikipedia.org/wiki/Horner's_method) to achieve efficient computation...
	p = c[ i ]*x + c[ i-1 ];
	i -= 2;
	while ( i >= 0 ) {
		p = p*x + c[ i ];
		i -= 1;
	}
	return p;
} // end FUNCTION evalpoly()


// EXPORTS //

module.exports = evalpoly;

},{}],10:[function(require,module,exports){
/* jshint evil:true */
'use strict';

// EVALPOLY FACTORY //

/**
* FUNCTION: factory( c )
*	Returns a function for evaluating a polynomial.
*
* @param {Number[]|Int8Array|Uint8Array|Uint8ClampedArray|Int16Array|Uint16Array|Int32Array|Uint32Array|Float32Array|Float64Array} c - polynomial coefficients sorted in ascending degree
* @returns {Function} function for evaluating a polynomial
*/
function factory( c ) {
	var f;
	var n;
	var m;
	var i;

	// Code generation. Start with the function definition...
	f = 'return function evalpoly(x){';

	// Create the function body...
	n = c.length;

	// If no coefficients, the function always returns 0...
	if ( n === 0 ) {
		f += 'return 0;';
	}
	// If only one coefficient, the function always returns that coefficient...
	else if ( n === 1 ) {
		f += 'return ' + c[ 0 ] + ';';
	}
	// If more than one coefficient, apply Horner's method...
	else {
		// If `x == 0`, return the first coefficient...
		f += 'if(x===0){return ' + c[ 0 ] + ';}';

		// Otherwise, evaluate the polynomial...
		f += 'return ' + c[ 0 ];
		m = n - 1;
		for ( i = 1; i < n; i++ ) {
			f += '+x*';
			if ( i < m ) {
				f += '(';
			}
			f += c[ i ];
		}
		// Close all the parentheses...
		for ( i = 0; i < m-1; i++ ) {
			f += ')';
		}
		f += ';';
	}
	// Close the function:
	f += '}';

	// Create the function in the global scope:
	return ( new Function( f ) )();

	/**
	* returns
	*	function evalpoly( x ) {
	*		if ( x === 0 ) {
	*			return c[ 0 ];
	*		}
	*		return c[0]+x*(c[1]+x*(c[2]+x*(c[3]+...+x*(c[n-2]+x*c[n-1]))));
	*	}
	*/
} // end FUNCTION factory()


// EXPORTS //

module.exports = factory;
},{}],11:[function(require,module,exports){
'use strict';

// EXPORTS //

module.exports = require( './evalpoly.js' );
module.exports.factory = require( './factory.js' );
},{"./evalpoly.js":9,"./factory.js":10}],12:[function(require,module,exports){
'use strict';

// EXPORTS //

module.exports = Math.exp;

},{}],13:[function(require,module,exports){
'use strict';

// MODULES //

var LOW = require( './low.js' );


// NOTES //

/**
* float64 (64 bits)
* f := fraction (significand/mantissa) (52 bits)
* e := exponent (11 bits)
* s := sign bit (1 bit)
*
* |-------- -------- -------- -------- -------- -------- -------- --------|
* |                                Float64                                |
* |-------- -------- -------- -------- -------- -------- -------- --------|
* |              Uint32               |               Uint32              |
* |-------- -------- -------- -------- -------- -------- -------- --------|
*
* If little endian (more significant bits last):
*                         <-- lower      higher -->
* |   f7       f6       f5       f4       f3       f2    e2 | f1 |s|  e1  |
*
* If big endian (more significant bits first):
*                         <-- higher      lower -->
* |s| e1    e2 | f1     f2       f3       f4       f5        f6      f7   |
*
*
* Note: in which Uint32 can we find the lower order bits? If LE, the first; if BE, the second.
* Refs: http://pubs.opengroup.org/onlinepubs/9629399/chap14.htm
*/


// VARIABLES //

var FLOAT64_VIEW = new Float64Array( 1 );
var UINT32_VIEW = new Uint32Array( FLOAT64_VIEW.buffer );


// SET LOW WORD //

/**
* FUNCTION: setLowWord( x, low )
*	Sets the less significant 32 bits of a double-precision floating-point number.
*
* @param {Number} x - double
* @param {Number} low - unsigned 32-bit integer to replace the lower order word of `x`
* @returns {Number} new double having the same higher order word as `x`
*/
function setLowWord( x, low ) {
	FLOAT64_VIEW[ 0 ] = x;
	UINT32_VIEW[ LOW ] = ( low >>> 0 ); // identity bit shift to ensure integer
	return FLOAT64_VIEW[ 0 ];
} // end FUNCTION setLowWord()


// EXPORTS //

module.exports = setLowWord;

},{"./low.js":14}],14:[function(require,module,exports){
'use strict';

// MODULES //

var isLittleEndian = require( 'utils-is-little-endian' );


// INDEX //

var LOW;
if ( isLittleEndian === true ) {
	LOW = 0; // first index
} else {
	LOW = 1; // second index
}


// EXPORTS //

module.exports = LOW;

},{"utils-is-little-endian":16}],15:[function(require,module,exports){
'use strict';

var ctors = {
	'uint16': Uint16Array,
	'uint8': Uint8Array
};


// EXPORTS //

module.exports = ctors;

},{}],16:[function(require,module,exports){
'use strict';

// MODULES //

var ctors = require( './ctors.js' );


// IS LITTLE ENDIAN //

/**
* FUNCTION: isLittleEndian()
*	Returns a boolean indicating if an environment is little endian.
*
* @returns {Boolean} boolean indicating if an environment is little endian
*/
function isLittleEndian() {
	var uint16_view;
	var uint8_view;

	uint16_view = new ctors[ 'uint16' ]( 1 );

	// Set the uint16 view to a value having distinguishable lower and higher order words.
	// 4660 => 0x1234 => 0x12 0x34 => '00010010 00110100' => (0x12,0x34) == (18,52)
	uint16_view[ 0 ] = 0x1234;

	// Create a uint8 view on top of the uint16 buffer:
	uint8_view = new ctors[ 'uint8' ]( uint16_view.buffer );

	// If little endian, the least significant byte will be first...
	return ( uint8_view[ 0 ] === 0x34 );
} // end FUNCTION isLittleEndian()


// EXPORTS //

module.exports = isLittleEndian();

},{"./ctors.js":15}],"pfilter":[function(require,module,exports){
var START = new Date()

let fmin = require ('fmin')
let mathLib = require('./mathLib')
let snippet = require('./modelSnippet.js')
let simulator = require ('./simulator.js')

function pfilterCalculation (input) {//filter.traj , save.params
  // {params:inputArr, Np:100,times:times, dt:1 / 365.25,runPredMean:1,  dataCases:dataCases, interpolPop:interpolPopulation, interpolBirth:interpolBirth}
  let START =new Date()
  let defaults = {params:-1, Np:-1, tol:1e-17, maxFail:Infinity, runPredMean:0, runPredVar:0, runFilterMean:0, runSaveStates:0, times:-1, dt:-1,
                   dataCases:0}
  for(let prop in defaults) {
    if(typeof input[prop] == 'undefined') {
      input[prop] = defaults[prop]
    }
  }
  if (input.params === -1 || input.Np === -1 || input.times === -1 || input.dt === -1) {
    throw 'Some required arguments are missed'
  }
  
  var params = input.params
  var Np = input.Np
  var toler = input.tol
  var maxFail = input.maxFail
  var dt = input.dt
  var dataCases = input.dataCases
  var interpolPop = input.interpolPop
  var interpolBirth = input.interpolBirth
  
  let [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, I_0, R_0] = params
let [t0, tdata] = [1940, 1944]
let nvars = 5
let deltaT = 14 / 365.25
let doPredictionVariance = 0, doPredictionMean = 1, doFilterMean = 0 , allFail = 0

let timeLen = dataCases.length 
let nlost = 0

let rate = [], trans = []
var particles = new Array(Np).fill(null).map(() => Array(5)),
 state =[]
let sampleNum = Array.from(Array(Np).keys())
let condLoglik = []
let stateSaved = []
let temp 
let timeCountData = 0, ws ,w , vsq, sumsq, ess, loglik = 0, lik 

let predictionMean, predictionVariance, filterMean
let states = Array(Np).fill(null).map(() => Array(nvars))
let weights, normalWeights, S, E, I, R, del_t, ST, simulateValue
let modelCases, likvalue
var st, births, pop,birthrate
if (doPredictionMean) {
  predictionMean = Array(timeLen).fill(null).map(() => Array(nvars))
}
if (doPredictionVariance) {
  predictionVariance = Array(timeLen).fill(null).map(() => Array(nvars))
}
if (doFilterMean) {
  filterMean = Array(timeLen).fill(null).map(() => Array(nvars))
}
state = snippet.initz(interpolPop(t0), S_0, E_0, I_0, R_0)
// First Np sets
var particles = new Array(Np).fill(null).map(() => [].concat(state));
temp = new Array(Np).fill(null).map(() => [].concat(state))
// Time loop
for (k = t0; k <= Number(dataCases[timeLen - 2][0]) + deltaT / 3 ; k += deltaT){//Number(dataCases[timeLen - 2][0]) + deltaT / 3
  if ( k > tdata - deltaT && k <= tdata) {
    k = tdata
  }
  weights = []; normalWeights = []
  
  if ( k > t0) {
    for (np = 0; np < Np; np++) { // copy the particles
      temp[np] = [].concat(particles[sampleNum[np]])

      temp[np][nvars - 1] = 0
    }
  } 
  
  // if (k <= tdata || k > Number(dataCases[dataCases.length - 1])) {
      steps = mathLib.numMapSteps(k, k + deltaT, dt)
  // } else {
  //     steps = mathLib.numEulerSteps(k, Number(dataCases[timeCountData + 1][0]), dt)
  // }
    del_t = (1 / steps )* deltaT

  for (let stp = 0; stp < steps; stp++) { // steps in each time interval
    st = k + stp * del_t
    pop = interpolPop(st)
    birthrate = interpolBirth(st)
    births = mathLib.rpois(birthrate * (1- snippet.rprocessVaccine(st)) * del_t )
      for (np = 0; np < Np; np++){ //calc for each particle
        trans = []
        S = temp[np][0]; E = temp[np][1]; I = temp[np][2]; R = temp[np][3]; H = temp[np][4]
        simulateValue = snippet.rprocess(params, st, del_t, [S,E,I,R,H], pop, births)
        temp[np][0] = simulateValue[0]; temp[np][1] = simulateValue[1]; temp[np][2] = simulateValue[2]; temp[np][3] = simulateValue[3]; temp[np][4] = simulateValue[4]
      }
    }
    for (np = 0; np < Np; np++){ 
      particles[np][0] = temp[np][0]
      particles[np][1] = temp[np][1]
      particles[np][2] = temp[np][2]
      particles[np][3] = temp[np][3]
      particles[np][4] = temp[np][4]
      H = temp[np][4]
  //***********weight*************
      if (k >= Number(dataCases[0][0])){
        if (stateSaved) {
          stateSaved.push(particles[np]) //[S,E,I,R,H])
        }
        modelCases = Number(dataCases[timeCountData][1])
        likvalue = snippet.dmeasure(rho, psi, H, modelCases, giveLog = 0)
        weights.push(likvalue)           
      }
    }
   
  
  //normalize
  if (k >= Number(dataCases[0][0])){
    let sumOfWeights = 0
    for (let i = 0; i < Np; i++) {
      sumOfWeights += weights[i]
    }
    for (let i = 0; i < Np; i++) {
      normalWeights[i] = weights[i] / sumOfWeights
    }
    // check the weights and compute sum and sum of squares
    w = 0, ws = 0, nlost = 0
    for (let i = 0; i < Np; i++) {
      if (weights[i] > toler) {
        w += weights[i]
        ws += weights[i] ** 2
      } else { // this particle is lost
        weights[i] = 0;
        nlost++
      }
    }
    if (nlost > maxFail) {
      throw 'execution terminated. The number of filtering failures exceeds the maximum number of filtering failures allowed. '
    }
    if (nlost >= Np) { 
      allFail = 1 // all particles are lost
    } else {
      allFail = 0
    }
    
    if (allFail) {
      lik = Math.log(toler) // minimum log-likelihood
      ess = 0  // zero effective sample size
    } else {
      ess = w * w / ws  // effective sample size
      lik = Math.log(w / Np)// mean of weights is likelihood
    }
    condLoglik[timeCountData] = [timeCountData + 1, lik]
    // the total conditional logliklihood in the time process is loglik
    loglik += lik
    mathLib.nosortResamp(Np, normalWeights, Np, sampleNum, 0)
    
    // Compute outputs
    for (let j = 0; j< nvars; j++) {
      // compute prediction mean
      if (doPredictionMean || doPredictionVariance) {
        let sum = 0, nlost = 0
        for (let nrow =0; nrow < Np; nrow++){
          if (particles[nrow][j]) {
            sum += particles[nrow][j]
          } else {
            nlost++
          }
        }
        sum /= Np
        predictionMean[timeCountData][j] = sum
      }  
      // compute prediction variance
      if (doPredictionVariance) {
        sumsq = 0
        for (let nrow = 0; nrow < Np; nrow++){
          if (particles[nrow][j]) {
            vsq = particles[nrow][j] - sum
            sumsq += Math.pow(vsq, 2)
          }
        }
        predictionVariance[timeCountData][j] = sumsq / (Np - 1) 
      }
      //  compute filter mean
      if (doFilterMean) {
        if (allFail) {   // unweighted average
          ws = 0
          for (let nrow =0; nrow < Np; nrow++){
            if (particles[nrow][j]) {
              ws += particles[nrow][j]
            }
          } 
          filterMean[timeCountData][j] = ws / Np//;console.log(ws / Np)
        } else {      // weighted average
          ws = 0
          for (let nrow =0; nrow < Np; nrow++){
            if (particles[nrow][j]) {
              ws += particles[nrow][j] * weights[nrow]
            }
          }
          filterMean[timeCountData][j] = ws / w
        }
      }
    }
    }
    timeCountData++
}//endTime

  console.log('loglike=',loglik)
  console.log('runing time=', (new Date() - START) / 1000)
  activateDownload ()
   if (input.runPredMean) {
    return predictionMean
  }

  if (input.runPredVar) {
    return predictionVariance
  }

  if (input.runFilterMean) {
    return filterMean
  }
  
  if (input.runSaveStates){
    return stateSaved
  }

}

module.exports = {
  pfilterCalculation
}
},{"./mathLib":1,"./modelSnippet.js":2,"./simulator.js":4,"fmin":7}]},{},[]);
