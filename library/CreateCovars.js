const fs = require('fs');
let linear = require('./everpolate').linear;

let create_covars = function(dir_covidtesting,endTime="2020-04-27",predTime=null){
	let t0 = new Date("2019-12-31");
	let tf = new Date(endTime);

	// read all rows and chaeck the border time and convert the selected colnames
	let temp, file;
	file = fs.readFileSync(dir_covidtesting).toString();
	let lines = file.split(/\r\n|\n/);
	let selected_colnames = ["Reported Date",
		"Total patients approved for testing as of Reporting Date",
		"Under Investigation"];
	let dataCovar_Index = [];	
	let dataCovar = [];
	let dataCovar_temp = lines[0].replace(/['"]+/g, '').split(',');
	for(let i = 0; i < selected_colnames.length; i++){
		dataCovar_Index.push(dataCovar_temp.indexOf(selected_colnames[i]));
	}
  let startdate = new Date("2020-02-03");
	for (let i = 1; i < lines.length; i++) {
		let temp = lines[i].split(',');
		if(temp.length > 1) {
			tempDate = new Date(temp[dataCovar_Index[0]]);
			if( tempDate > startdate && tempDate <= tf){
				let tempDayCount = Math.ceil(Math.abs(tempDate - t0) / (1000 * 3600 * 24));
			  let tempcovar =	[tempDayCount];
				for(let j = 1; j < selected_colnames.length; j++){
					tempcovar.push(Number(temp[dataCovar_Index[j]]));
				}
				dataCovar.push(tempcovar);
			}
		}
	}

	let coaverTimes = dataCovar.map(x => x[0]);
	let coaverData = dataCovar.map(x => x[1]);
	let timeLength = dataCovar[dataCovar.length-1][dataCovar_Index[0]];
	let time = Array.from(Array(timeLength), (_, i) => i + 1);

	let total_tests = linear(time, coaverTimes, coaverData);
	total_tests = total_tests.map(x => Math.ceil(x));
  
  tests = total_tests.map( (_,i,arr) => Math.max((arr[i]) - (arr[i-1]),0));
	tests = tests.map( (_,i,arr) => Math.ceil(arr[i] * 0.5 + arr[i+1] * 0.5));
	tests[timeLength - 1] = tests[timeLength - 2];

	let data = time.map( (x,i) => [x,tests[i]]);
	data.unshift(["time","tests"]);

	if (predTime != null) {
    predTime = new Date(predTime);
		let add_time = Math.ceil(Math.abs(tf - predTime) / (1000 * 3600 * 24));
		for(let i =0; i < add_time; i++)
			data.push([timeLength + i,10e3]);
	}
	
	return data;
}

module.exports = create_covars;