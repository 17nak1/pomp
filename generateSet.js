/**
 *  @file       generateSet.js        This function can generate sets as initial parameters for calculation in trajMatch. 
 *                                    User needs to provide values for 
 *                                    tolerance : Determine how far from the best liklihood is acceptable.
 *                                    number of profile : Number of points to be consider in the interval of paramLimits.
 *                                    number of points : Number of  points to be generated.
 *                                    s : The value that is used to add noice in the best set and generate more points.
 *                                    indexMult : Defines how to generate new values ('divide & multiply' or 'subtract & add')
 *                                    indexInc : Defines how to generate new values (divide or multiply or both)(subtract or add or both)
 * 
 *                                    
 *                                    Also user needs to define the following for each estimating parameter through the function 'determineRunProperties';
 *                                    paramLimits : Lower and upper bound.
 *                                    logScale : If we consider calculation in the log scale, this value equals one.
 *                                    flagBound : If the generated values should be in the interval (0,1), this value equals one.
 *                                      
 *
 *  @author     Nazila Akhavan
 *  @date       March 2019
 */
fs = require('fs')
var dataset = []
var file = fs.readFileSync('./DeterministicSEIR_all.csv').toString()
var lines = file.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataset.push(lines[i].split(','))
}
dataset.pop()

// Order inputs based on ['R0', 'amplitude', 'gamma', 'mu', 'sigma', 'rho', 'psi', 'S_0', 'E_0', 'R_0', 'I_0']
if(JSON.parse(lines[0].split(',')[9]) !== "R_0"){
  for (i = 0; i < dataset.length; i++) {
    var tem = dataset[i][9]
    dataset[i][9] = dataset[i][10]
    dataset[i][10] = tem
  }
}
// Indicies for params. eg. params[MU] instead of params[3]
const R0Index = 0
const AMPLITUDE = 1
const GAMMA = 2
const MU = 3
const SIGMA = 4
const RHO = 5
const PSI = 6
const S_0 = 7
const E_0 = 8
const R_0 = 9
const I_0 = 10
const LogLikIndex = 11


var tolerance = 50
var numberOfProfile = 50
var numberOfPoints = 2000

var indexInc = 0
var indexMult = 1
var s = 0.01

determineRunProperties = function  (run) {
  if (run == R0Index) {
    logScale = 0
    paramLimits = [0,80]
    flagBound = 0
  } else if (run == AMPLITUDE) {
    logScale = 0
    paramLimits = [0,1]
    flagBound = 1
  } else if (run == MU) {
    logScale = 0
    paramLimits = [0,1]
    flagBound = 1
  } else if (run == RHO) {
    logScale = 0
    paramLimits = [0,1]
    flagBound = 1
  } else if (run == PSI) {
    logScale = 0
    paramLimits = [0,1]
    flagBound = 0
  } 
  return [paramLimits, logScale, flagBound ]  
}

generateSets (dataset, LogLikIndex, [R0Index, AMPLITUDE, MU, RHO, PSI], tolerance, numberOfProfile, numberOfPoints, indexInc, indexMult, s)

function generateSets (dataset, LogLikIndex, paramIndexArray, tolerance, numberOfProfile, numberOfPoints, indexInc, indexMult, s) {
  for ( index = 0; index < paramIndexArray.length; index++) {
    var paramIndex = paramIndexArray[index]
    var paramLimits, logScale, flagBound
    [paramLimits, logScale, flagBound] = determineRunProperties (paramIndex)
    generate (dataset, LogLikIndex, paramIndex, tolerance, numberOfProfile, numberOfPoints, paramLimits, logScale, flagBound, indexInc, indexMult, s)
  }
}

function generate (dataset, LogLikIndex, paramIndex, tolerance, numberOfProfile, numberOfPoints, paramLimits, logScale, flagBound, indexInc, indexMult, s) {
  var Maxloglik, step = 0, ltemp = 0, temp = [], temp2 = []
  var newDataset =[], paramArray = []
  var set1 = [], paramProfile = []
  
  // Reorder dataset descending based on LogLik column and find the maximum LogLik
  dataset.sort(sortFunction)
  Maxloglik = dataset[0][LogLikIndex]
  
  // Calculate the step size for the parameter limits interval 
  if (logScale === 1) {
    if (paramLimits[0] <= 0 || paramLimits[1] <= 0) {
      throw "The lower(upper) bound for the parameter is not positive."
    }
    step = (Math.log(paramLimits[1]) - Math.log(paramLimits[0])) / (numberOfProfile - 1)
    ltemp = Math.log(paramLimits[0])
  } else {
    step = (paramLimits[1] - paramLimits[0]) / (numberOfProfile - 1)
    ltemp = paramLimits[0]
  }
  // newDataset include rows that has LogLik in [LogLik - tolerance, LogLik] from which Paramprogile will be make. 
  for (i = 0; i < dataset.length; i++ ) { 
    if (dataset[i][LogLikIndex] > Maxloglik - tolerance && dataset[i][LogLikIndex] < 0) {
      newDataset.push((dataset[i]).map(Number))
    } else {
      i = dataset.length
    }
  }
  
  for (i = 0; i < numberOfProfile; i++) {
    if (logScale === 1) {
      paramArray.push(Math.exp(ltemp))
    } else {
      paramArray.push(Number(ltemp.toFixed(8)))
    }
    ltemp += step
  }
  
  for (q = 1; q < paramArray.length; q++) {
    set1 = []
    for (j =0; j < newDataset.length; j++) {
      if (newDataset.length > 0){
        if (newDataset[j][paramIndex] >= paramArray[q - 1] && newDataset[j][paramIndex] <= paramArray[q]) {
          set1.push(newDataset[j])
        }
      }
    }
    if(set1.length > 0) {
      set1.sort(sortFunction) 
      paramProfile.push(set1[0])
    } 
  }
  
  temp = paramProfile.map(row => [].concat(row))
  temp2 = paramProfile.map(row => [].concat(row))
  
  for (q = 1; q <= Math.ceil(numberOfPoints / temp.length); q++) {
    if (indexMult === 1) {
      if (indexInc === -1) {
        nextDivide(temp2, paramIndex, s, paramProfile)
      } else if (indexInc === 1) {
        nextMultiply(temp, paramIndex, s, paramProfile)
      } else {
        if (q % 2 === 1) {
          nextDivide(temp2, paramIndex, s, paramProfile)
        } else {
          nextMultiply(temp, paramIndex, s, paramProfile)
        }
      }
    } else {
      if (indexInc === -1) {
        nextSubtract(temp2, paramIndex, s, paramProfile)
      } else if (indexInc === 1) {
        nextAdd(temp, paramIndex, s, paramProfile)
      } else {
        if (q % 2 === 1) {
          nextSubtract(temp2, paramIndex, s, paramProfile)
        } else {
          nextAdd(temp, paramIndex, s, paramProfile)
        }
      }
    }
  }
  paramProfile.splice(numberOfPoints)
  for (i = 0; i < paramProfile.length; i++) {
    paramProfile[i].pop()
  }
  if (flagBound === 1) {
    for (i = 0; i < paramProfile.length; i++) {
      if(paramProfile[i][paramIndex] > 1 - 1e-6) {
        paramProfile[i][paramIndex] = 1 - 1e-6
      } else if (paramProfile[i][paramIndex] < 1e-6) {
      paramProfile[i][paramIndex] = 1e-6
      }
    }
  } else if (flagBound === 2) {
    for (i = 0; i < paramProfile.length; i++) {
      if (paramprofile[i][paramIndex] < 1e-6)
        paramProfile[i][paramIndex] =  1e-6
    }
  }
  
  const createCsvWriter = require('csv-writer').createArrayCsvWriter;
  const csvWriter = createCsvWriter({
    header: ['R0', 'amplitude', 'gamma', 'mu', 'sigma', 'rho', 'psi', 'S_0', 'E_0', 'R_0', 'I_0'],
    path: verifyPath (paramIndex)
  })
   
  csvWriter.writeRecords(paramProfile)
    .then(() => {
    console.log(verifyPath (paramIndex) + '...done')
  })
}    

function sortFunction(a, b) {
  if (Number(a[LogLikIndex]) === Number(b[LogLikIndex])) {
    return 0
  }
  else {
    return (Number(a[LogLikIndex]) < Number(b[LogLikIndex])) ? 1 : -1;
  }
}

function nextDivide(temp, paramIndex, s, paramProfile) {
  for (i = 0; i < temp.length; i++) { 
    temp[i][paramIndex] /= (1 + s)
  }
  paramProfile.push(...[].concat(temp.map(row => [].concat(row))))
}
function nextMultiply(temp, paramIndex, s, paramProfile) {
  for (i = 0; i < temp.length; i++) {
    temp[i][paramIndex] *= (1 + s)
  }
  paramProfile.push(...[].concat(temp.map(row => [].concat(row))))
}
function nextAdd (temp, paramIndex, s, paramProfile) {
  for (i = 0; i < temp.length; i++) {
    temp[i][paramIndex] += s
  }
  paramProfile.push(...[].concat(temp.map(row => [].concat(row))))
}

function nextSubtract (temp, paramIndex, s, paramProfile) {
  for (i = 0; i < temp.length; i++) {
    temp[i][paramIndex] -= s
  }
  paramProfile.push(...[].concat(temp.map(row => [].concat(row))))
}

function verifyPath (paramIndex) {
  if ( paramIndex == R0Index) {
    return "./paramSet_R0.csv"
  } else if ( paramIndex == AMPLITUDE) {
    return './paramSet_amplitude.csv'
  } else if ( paramIndex == MU) {
    return './paramSet_mu.csv'
  } else if ( paramIndex == RHO) {
    return './paramSet_rho.csv'
  } else if ( paramIndex == PSI) {
    return './paramSet_psi.csv'
  } 
}

