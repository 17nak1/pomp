const epiFit = require("./main.js");

/** epiFit provides the final set of parameters saved as `modelParams.csv` in the results folder.
 * 
 * To run the model for each city you need to provide covar.csv and data.csv as explained in the document. 
 * pop is the population of the city.
 * TC is the accurate episode date of first case.
 */

epiFit({
  pathToCovar:('./samples/covar.csv'),
  pathToData:('./samples/data.csv'),
  pop:3e6,
  TC:21                                   
})