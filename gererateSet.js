var dataset = []
var file = fs.readFileSync('./DeterministicSEIR_all.csv').toString()
var lines = file.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataset.push(lines[i].split(','))
}
determineRunProperties = function  (run=2,modeltype=1) {
  if (run == 0 || run ==1) {
    param = null
    lscale = null
    paramLimits = null
    flagBound = null
  } else if (run == 2) {
    param = "R0"
    lscale = 0
    paramLimits = [0,80]
    flagBound = 0
  } else if (run = 3) {
    param = "amplitude"
    lscale = 0
    paramLimits = [0,1]
    flagBound = 1
  } else if (run= 4) {
    param = "mu"
    lscale = 0
    paramLimits = [0,1]
    flagBound = 1
  } else if (run= 5) {
    param = "rho"
    lscale = 0
    paramLimits = [0,1]
    flagBound = 1
  } else if (run= 6) {
    param = "psi"
    lscale = 0
    paramLimits = [0,2]
    flagBound = 0
  }   
}
 modeltypes = "DeterministicSEIR"
runs = [2]

for (i = 0; i < length(runs); i++) {
	var run = runs[i]
	out = determineRunProperties(run)
    param = out[0]
    lscale = out[1]
    param_lims = out[2]
    flag_bound = out[3]
    k1 = Math.ceil(max(dataset[,"LogLik"]))
    dataset = subset(dataset,subset=LogLik>k1-tol&LogLik<0)

}