fs = require('fs')
var dataset1 = []
var file = fs.readFileSync('./ParamSet_DeterministicSEIR_run22222.csv').toString()
var lines = file.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataset1.push(lines[i].split(','))
}

var dataset2 = []
var file = fs.readFileSync('./paramSet_R0.csv').toString()
var lines = file.split('\n')
for (let i = 1; i < lines.length; i++) {
  dataset2.push(lines[i].split(','))
}
for( i= 0; i < dataset1.length; i++) {
  for (j=0; j < dataset1[0].length; j++) {
    if (Math.abs(dataset2[i][j] - dataset1[i][j]) > 1e-12) {
      console.log(i,j, dataset2[i][j] - dataset1[i][j])
    }

  }
}