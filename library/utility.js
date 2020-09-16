const fs = require("fs");

exports.utility = function (key,jobId ) {
  let jsonString
  try {
    jsonString = fs.readFileSync('./COVID.json'); 
  } catch (error) {
    jsonString = JSON.stringify([{}]);
  }
  
  let data = JSON.parse(jsonString);
  if (!data[0].key) data = [];
  jsonString = JSON.stringify([...data, { key: key, jobId: jobId }]);
  fs.writeFileSync('./COVID.json', jsonString);
}
