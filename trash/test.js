
var rate = new Array(6)
fs = require('fs')
let fmin = require ('fmin')
var mathLib = require('./mathLib')
const mathjs = require('mathjs') 
var linearInterpolator = require('linear-interpolator/node_main')
const libR = require('lib-r-math.js')
const {
  Poisson,
  rng: { MersenneTwister },
  rng: { normal: { Inversion } }
} = libR
const mt = new MersenneTwister(0)// need reference so we "reset" PRNG
const { rpois } = Poisson(new Inversion(mt))
mt.init(1234)
var vv =[]
//* Set the seed for rnorm-In R:RNGkind("L'Ecuyer-CMRG", normal.kind="Box-Muller");set.seed(1234) 
//or normal.kind="Box-Muller" or "BuggyKindermanRamage" or "Inversion" or "KindermanRamage"
var {
  Normal,
  rng: {
    LecuyerCMRG,
    normal: { BoxMuller }
  }
} = libR
const ad = new LecuyerCMRG(1234)
const { rnorm } = Normal(new BoxMuller(ad))
//////////////////////////////////////////data///////////////////////////////////////
var LondonBidata, LondonCovar, 
params = [29.2655634435, 0.3944025839968, 73.05, 0.004303868154365, 45.66, 0.417917730919406, 0.530358893934836, 0.020285841729738, 3.49E-05, 8.01E-05, 0.979679198700577, 1940, 1944]
var Np = 10, nreps = Np
var nvars = params.length
var toler = 10e-8
var w=0 , ws=0
//begin mif2
var Nmif = 1
var rw_size = .05, rw = new Array(Nmif).fill(null).map(() => Array(Np).fill(0.05))//change it to matrix
var rwIn = new Array(11).fill(1)
rwIn[2] = 0; rwIn[4] = 0
var indx = new Array(11).fill(1)
indx[2] = 0; indx[4] = 0
//end
//* 1st data set
var London_covar = fs.readFileSync('./London_covar.csv').toString()
var LondonCovar = []
var lines = London_covar.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonCovar.push(lines[i].split(','))
}
var dataCovar = [LondonCovar][0]
//* 2nd data set
LondonBidata = []
var London_BiData = fs.readFileSync('./London_BiData.csv').toString()
var lines = London_BiData.split('\n')
for (let i = 1; i < lines.length; i++) {
  LondonBidata.push(lines[i].split(','))
}
var dataCases = [LondonBidata][0]
//* main function****************************************************************

var [R0, amplitude, gamma, mu, sigma, rho, psi, S_0, E_0, R_0, I_0] = params
var timeSet = [t0, t1]
var estim = []
var place = []
var coolFrac = 0.5
var delT = 0.03832991102// = 2/52 
var dt = 1/ 365.25
var d1 = []// read time and population from 1st data and make interpolation function
var d2 = []// read time and birthrate from 1st data and make interpolation function
var d3 = []
for (let i = 0; i < dataCovar.length - 1; i++) {
  d1.push([Number(dataCovar[i][0]), Number(dataCovar[i][1])])
  d2.push([Number(dataCovar[i][0]), Number(dataCovar[i][2])])
}
for (let i = 0; i < dataCases.length - 1; i++) {
  d3.push([Number(dataCases[i][0]), Number(dataCases[i][1])])
}
var interpolPop = linearInterpolator(d1)
var interpolBirth = linearInterpolator(d2)
  
// function poly (paramNoise, time, trans) {
  var cases = []
  var va = 0, seas
  var t0 = timeSet[0], tdata = timeSet[1] 
  var particles = []
  //***********************************************TIME LOOP************************************
  for (k = t0; k <= 1945+3*delT; k +=delT){//Number(dataCases[dataCases.length - 2][0]) + delT / 3
    if (k <= tdata && k > tdata - delT) {
      k = tdata
    }
  if ( k === t0) {
    var Nlog = mathLib.toLogBarycentric([params[7], params[8], params[9], params[10]],4)
    var N = mathLib.fromLogBarycentric(Nlog, 4)
    var m = interpolPop(k) / (N[0] + N[1] + N[2] + N[3]);
    params[7] = Math.round(m * N[0]),
    params[8] = Math.round(m * N[1]),
    params[9] = Math.round(m * N[2]),
    params[10] = Math.round(m * N[3])
  //begin mif2
  var paramNoiseM = Array(params.length).fill(null).map(() => Array(Np))
  var timeLen = dataCases.length;
  var s = (1 - 50 * timeLen * coolFrac) / (coolFrac - 1)
  var cmn = (s + 1)/ (s + k + (Nmif - 1) * timeLen)
    for (cnt = 0; cnt < params.length; cnt++){
      for (np = 0; np < Np; np++){
        if(rwIn[cnt] === 1){
          paramNoiseM[cnt][np] = params[cnt] * ( 1 + cmn * rw[Nmif - 1][np] * rnorm(1))// weighted??????????????????????????
        } else {
          paramNoiseM[cnt][np] = params[cnt]
        } 
      }   
    }
    particles = mathjs.transpose(paramNoiseM)
    var nvars = particles[0].length; nreps = particles.length
    var npars = params.length
  }//endif ( k === t0)
  var pop = interpolPop(k)
  var birthrate = interpolBirth(k) 
  var tt = (k - Math.floor(k)) * 365.25
  var loglik = 0
  var lik = new Array(Np)
  var weights = []
  var sample = new Array(Np)
var aa = 1,t_1 = 1944.1916495550995, t_2 = 1944.2299794661194

// if ( k == t_2) 
  // console.log("enter1",particles[aa][9],k)
  //****************************************PARTICLE LOOP**************************************//
  for (np = 0; np < Np; np++){//calc for each particle
    
     // console.log(np, "enter2",particles[aa][9],k)

    var trans = new Array(6).fill(0)
    var R0 = particles[np][0], amplitude = particles[np][1], gamma = particles[np][2], mu = particles[np][3], sigma = particles[np][4] 
    var S = particles[np][7], E = particles[np][8], I = particles[np][9], R = particles[np][10], H = 0
    if ((tt >= 7 && tt <= 100) || (tt >= 115 && tt <= 199) || (tt >= 252 && tt <= 300) || (tt >= 308 && tt <= 356)) {
      seas = 1 + amplitude * 0.2411 / 0.7589
    } else {
      seas = 1 - amplitude
    }                 
    var beta = R0 * (gamma + mu) * (sigma + mu) * seas / sigma//seasonal transmission rate
    var foi = beta * I / pop
    rate[0] = foi//force of infection
    rate[1] = mu// natural S death
    rate[2] = sigma// rate of ending of latent stage
    rate[3] = mu// natural E death
    rate[4] = gamma// recovery
    rate[5] = mu// natural I death       
    // transitions between classes
    steps = mathLib.numEulerSteps(k, k + delT, dt)
    var del_t = 1 / steps * delT
    for (let stp = 0; stp < steps; stp++) { // steps in each time interval
      var births = rpois(1, birthrate * (1 - va) * del_t )// Poisson births
      mathLib.reulermultinom(2, Math.round(S), 0, del_t, 0, rate, trans)
      mathLib.reulermultinom(2, Math.round(E), 2, del_t, 2, rate, trans)
      mathLib.reulermultinom(2, Math.round(I), 4, del_t, 4, rate, trans)//;if ( k == 1944.3066392881592  && np == 0) console.log(trans, S,E,I,R, H)
      S += (births - trans[0] - trans[1])
      E += (trans[0] - trans[2] - trans[3]) 
      I += (trans[2] - trans[4] - trans[5]) 
      // E = (E > 0) ? (E + trans[0] - trans[2] - trans[3]) : 0
      // I = (I > 0) ? (I + trans[2] - trans[4] - trans[5]) : 0
      R = pop - S - E - I
      H += trans[4] 
    }
    // if ( k == t_2 && np == aa) console.log(trans, S,E,I,R, H)
    particles[np][7] = S
    particles[np][8] = E
    particles[np][9] = I
    particles[np][10] = R
    // if ( k == t_2 && np == aa) 
      // console.log("ex",particles[aa][9])
    
    if (k >= tdata) {
      var rho = particles[np][5], psi = particles[np][6]
      var mn = rho * H
      var v = mn * (1.0 - rho + psi * psi * mn)
      var tol = 1.0e-18
      var modelCases = Number(dataCases[Math.ceil((k - tdata) / delT)][1])
      if(!isNaN(modelCases)){
        if (modelCases > 0.0) {
          lik[np] = mathLib.pnorm(modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0) - mathLib.pnorm(modelCases - 0.5, mn, Math.sqrt(v) + tol, 1, 0) + tol
        } else {
          lik[np] = mathLib.pnorm((modelCases + 0.5, mn, Math.sqrt(v) + tol, 1, 0)) + tol
        }
      } else {
        lik[np] = 1
      }
      loglik += Math.log(lik[np])
    }
  }//end particle loop


  if (loglik != 0) {
    for(w = 0, ws = 0, np = 0; np < nreps; np++) {//make weights vector
      weights.push(Math.log(lik[np]) / loglik)

      if (weights[np] > toler) {
        w += weights[np]
        ws += Math.pow(weights[np], 2)
      } else {
        weights[np] = 0
        nlost++
      }
    }
    if (nlost >= nreps) {
      var allFail = 1
    }
    
   //*****************************************RESAMPLE*******************************************************

    if (!all_fail) { // resample the particles unless we have filtering failure
    var xdim = new Array(2)
    var sample = new Array(Np)
    var ss = 0, st = 0, ps = 0, pt = 0 

    // create storage for new states
    // xdim[0] = nvars; xdim[1] = Np
    var newstates = Array(Np).fill(Array(nvars)); nprotect++;
    ss = particles
    st = newstates

    // create storage for new parameters
    if (do_pr) {
      // xdim[0] = npars; xdim[1] = Np;
      var newparams = Array(npars).fill(Array(Np)); nprotect++;
      ps = params
      pt = newparams
    }

    // resample
    nosort_resamp(nreps,weights,Np,sample,0);
    for (k = 0; k < np; k++) { // copy the particles
      for (j = 0, xx = ss + nvars * sample[k]; j < nvars; j++, st++, xx++)
        newstates[st] = particles[xx];
      if (do_pr) {
        for (j = 0, xp = ps+npars*sample[k]; j < npars; j++, pt++, xp++)
          *pt = *xp;
      }
      if (do_ta) xanc[k] = sample[k]+1;
    }

  } else { // don't resample: just drop 3rd dimension in x prior to return

    PROTECT(newdim = NEW_INTEGER(2)); nprotect++;
    dim = INTEGER(newdim);
    dim[0] = nvars; dim[1] = nreps;
    SET_DIM(x,newdim);
    setrownames(x,Xnames,2);
    fixdimnames(x,dimnm,2);

    if (do_ta)
      for (k = 0; k < np; k++) xanc[k] = k+1;
  }
    // if (allFail) {

    // } else {
    //   loglikC = Math.log(w / nreps); // mean of weights is likelihood
    //   // ess = w * w / ws; // effective sample size
    // }
    // console.log(weights)
  // mathLib.nosortResamp(nreps, weights, Np, sample, 0)
  //   for (np = 0; np < Np; np++) { // copy the particles
  //     particles = particles + Np * sample[k]
  //     for (j = 0; j < nvars; j++, st++, xx++){
  //       *st = *xx;
  //     }
  //   }
  //   } else { // don't resample: just drop 3rd dimension in x prior to return






      
  //   }
    // if ( sampleNum[np] !== sample[np]) {
      //   particles[np] = particles[sample[np]]
      // }
}//endTime

















  //  w = 0, ws = 0
  // for(let np = 0, nlost = 0; np < Np; np++) {
  //   if (weights[np] > toler) {
  //     w += weights[np]
  //     ws += Math.pow(weights[np], 2)
  //   } else {
  //     weights[np] = 0
  //     nlost++
  //   }
  //   // console.log(w,ws)
  // }
  // // if ( w != 0) console.log(Math.log(w/Np))
  
// console.log(particles) 
const logMeanExp = function (x) {
  var mx = Math.max(...x)
  var s = x.map((x,i) => Math.exp(x - mx))
  var q = s.reduce((a, b) => a + b, 0)
  return mx + Math.log(q / x.length)
}

const numMapSteps = function (t1, t2, dt) {
  var DOUBLE_EPS = 10e-8
  var tol = Math.sqrt(DOUBLE_EPS)
  var nstep
  // nstep will be the number of discrete-time steps to take in going from t1 to t2.
  nstep = Math.floor((t2 - t1) / dt /(1 - tol))
  return (nstep > 0) ? nstep : 0
}

// const eulerModelSimulator = function (func, xstart, times, params,
//    deltat, rprocmode method, SEXP accumvars, SEXP covar, SEXP args, SEXP gnsi)
// {
//   if (deltat <= 0) throw 'error'
//   X = ret_array(nvars,nreps,ntimes,Snames)
//     // dim = INTEGER(GET_DIM(xstart)); nvars = dim[0]; nreps = dim[1];
//   var npars = params.length
//   var ntimes = times.length
//   for (step = 1, xs = REAL(X), xt = xs + nvars * nreps, time = times, t = time[0]; step < ntimes; step++, xs = xt, xt += nvars * nreps) {
//     if (t > time[step]) throw 'error'
//     // set accumulator variables to zero
//     for (j = 0; j < nreps; j++)
//       for (i = 0; i < nzeros; i++)
//         xt[zidx[i] + nvars * j] = 0.0;

//     // determine size and number of time-steps
//     switch (method) {
//     case onestep: default:  // one step
//       dt = time[step] - t;
//       nstep = (dt > 0) ? 1 : 0
//       break;
//     case discrete:      // fixed step
//       dt = deltat;
//       nstep = num_map_steps(t, time[step], dt)
//       break;
//     case euler:     // Euler method
//       dt = deltat;
//       nstep = numEulerSteps(t, time[step], dtAdd)
//       break;
//     }


//   }
// }
 
  

    
      

