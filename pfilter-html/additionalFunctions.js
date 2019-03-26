// var {
//   Normal,
//   Gamma,
//   rng: {
//     LecuyerCMRG,
//     normal: { BoxMuller }
//   }
// } = libR

// const ad = new LecuyerCMRG(1234)
// const { rnorm, pnorm } = Normal(new BoxMuller(ad))
// console.log(rnorm(), pnorm(0))//-0.7129350418967081

// const lc = new LecuyerCMRG(1234)
// const { rgamma } = Gamma(new BoxMuller(lc))
// console.log(rgamma(1, 1, 0.5))//1.2294789082038762


const numMapSteps = function (t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}

mathLib.fromLogBarycentric = function (xN, n) {
  var sum = 0, xt = []
  for (let i = 0; i < n; i++) {
    xt.push(Math.exp(xN[i]))
    sum += xt[i]
  }
  for (let i = 0; i < n; i++) {
    xN[i] = xt[i] / sum
  }
  return xN
}

mathLib.toLogBarycentric = function (xN, n) {
  var sum = 0
  for (let i = 0; i < n; i++) {
    sum += xN[i]
  }
  for (let i = 0; i < n; i++) {
    xN[i] = Math.log(xN[i] / sum)
  }
  return xN
}

mathLib.logit = function (p) {
  return Math.log(p / (1 - p))
}

mathLib.expit = function (x) {
  return 1 / (1 + Math.exp(-x))
}

mathLib.rgammawn = function (sigma, dt) {
  var sigmasq
  sigmasq = Math.pow(sigma, 2)
  return (sigmasq > 0) ? rgamma(1, dt / sigmasq, sigmasq) : dt
}

const logMeanExp = function (x) {
  var mx = Math.max(...x)
  var s = x.map((x,i) => Math.exp(x - mx))
  var q = s.reduce((a, b) => a + b, 0)
  return mx + Math.log(q / x.length)
}

// mathLib.rnorm = function (mu = 0, sd = 1) {
//   var val = Math.sqrt(-2.0 * Math.log(U.unif_rand())) * Math.cos(2.0 * pi * U.unif_rand())
//   return val * sd + mu
// }
mathLib.numMapSteps = function (t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}