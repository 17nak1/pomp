const libR = require('lib-r-math.js')
const {
Poisson,
rng: { MersenneTwister },
rng: { normal: { Inversion } }
} = libR
const mt = new MersenneTwister(0)//
const { rpois } = Poisson(new Inversion(mt))
mt.init(1234)

const libR2 = require('lib-r-math.js')
const {
 Binomial,
 rng: { KnuthTAOCP2002 }
} = libR2
const kn = new KnuthTAOCP2002(1234)
const { rnorm, rbinom } = Binomial(kn)


 function reulermultinom(m = 1, size, rateAdd, dt, transAdd, rate, trans) {
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
       size = rbinom(1, size, 1 - Math.exp(-p * dt)) // total number of events
      if (!(isFinite(size)))
        throw 'result of binomial draw is not finite.'
      m -= 1
      for (k = 0; k < m; k++) {
        if (rate[k + rateAdd] > p) p = rate[k + rateAdd]
        trans[k + transAdd] = ((size > 0) && (p > 0)) ? rbinom(1, size, rate[k + rateAdd] / p) : 0
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

function rprocess(t,dt,states,params,covar){
    let seas, beta, foi;
    let births, va, tt;
    let rate = [6], trans =[6];



    va = 0


    // term-time seasonality
    tt = (t-Math.floor(t))*365.25;
    if ((tt>=7&&tt<=100) || (tt>=115&&tt<=199) || (tt>=252&&tt<=300) || (tt>=308&&tt<=356))
    seas = 1.0+params.amplitude*0.2411/0.7589;
    else
    seas = 1.0-params.amplitude;

    // transmission rate
    beta = params.R0*(params.gamma+params.mu)*(params.sigma+params.mu)*seas/params.sigma;  //seasonal transmission rate
    // expected force of infection
    foi = beta*states.I/covar.pop;

    rate[0] = foi;  //         force of infection
    rate[1] = params.mu;             // natural S death
    rate[2] = params.sigma;        // rate of ending of latent stage
    rate[3] = params.mu;             // natural E death
    rate[4] = params.gamma;        // recovery
    rate[5] = params.mu;             // natural I death
    //if( t<=1944.03832991103 && t>=1944.03832991102)
    //printf(\"%f and %f \\n\", t, foi);
    // Poisson births
    births = Math.random()//rpois(1,covar.birthrate*(1-va)*dt);
    // transitions between classes
    reulermultinom(2,states.S,0,dt,0,rate,trans);
    reulermultinom(2,states.E,2,dt,2,rate,trans);
    reulermultinom(2,states.I,4,dt,4,rate,trans);
    states.S += births - trans[0] - trans[1];
    states.E += trans[0] - trans[2] - trans[3];
    states.I += trans[2] - trans[4] - trans[5];
    states.R = covar.pop - states.S - states.E - states.I;
    states.H += trans[4];           // true incidence
}

let timediff = 24;//years
let np = 10;
let dt = 1/365;
let nsteps = timediff / dt * np;
let t = 1394;
let states = {S: 45099, E: 310, I: 1, R: 1281337, H: 0};
let params = {R0: 31.3249, amplitude: 0.388362, gamma: 73.05, mu: 0.000646983, sigma: 45.66,}
let covar = {pop: 1326745.8826649827, birthrate: 86829.64443292982};
let now = Date.now();

for(let i = 1; i < nsteps; i++) {
    rprocess(t,dt,states,params,covar)
}

console.log('elapsed time=', Date.now() - now, 'ms')
console.log();