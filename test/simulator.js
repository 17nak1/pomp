
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
    S = simulateValue[0]; E = simulateValue[1], I = simulateValue[2], R = simulateValue[3], H = simulateValue[4]
  }
  return [S, E, I, R, H]
}

module.exports = simulator