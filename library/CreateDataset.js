const fs = require('fs');

let create_dataset = function(endTime="2020-05-04",predTime=null){
  t0 = new Date("2019-12-31");
  tf = new Date(endTime);

  // read all rows and chaeck the border time and convert the selected colnames
	let temp, file;
	file = fs.readFileSync('../samples/ON.csv').toString();
	let lines = file.split(/\r\n|\n/);
	let dataReport = [];

  delay = 0;
  let startdate = new Date("2020-01-01");
	for (let i = 1; i < lines.length; i++) {
		let temp = lines[i].split(',');
		if(temp.length > 1) {
			tempDate = new Date(temp[0]);
			if( tempDate > startdate && tempDate <= tf - delay){
				let tempDayCount = Math.ceil(Math.abs(tempDate - t0) / (1000 * 3600 * 24));
				dataReport.push([tempDayCount, Number(temp[1])]);
			}
		}
  }
  
	let timeLength = dataReport[dataReport.length-1][0];
	let time = Array.from(Array(timeLength), (_, i) => i + 1);
  
  reports = Array.from(Array(timeLength)).fill(0);
  for(let i = 0; i < dataReport.length; i++)
    reports[dataReport[i][0] - 1] = dataReport[i][1];
  for(let i = timeLength - delay - 1; i < timeLength; i++)
    reports[i] = NaN;

  // read all rows and chaeck the border time and convert the selected colnames
	temp, file;
	file = fs.readFileSync('../samples/covidtesting.csv').toString();
	lines = file.split(/\r\n|\n/);
  let selected_colnames = ["Reported Date",
    "Deaths",
    "Number of patients hospitalized with COVID-19",
    "Number of patients in ICU with COVID-19",
    "Number of patients in ICU on a ventilator with COVID-19"];
	let dataCovar_Index = [];	
	let dataCovar = [];
	let dataCovar_temp = lines[0].replace(/['"]+/g, '').split(',');
	for(let i = 0; i < selected_colnames.length; i++){
		dataCovar_Index.push(dataCovar_temp.indexOf(selected_colnames[i]));
	}
  startdate = new Date("2020-03-17");
	for (let i = 1; i < lines.length; i++) {
		let temp = lines[i].split(',');
		if(temp.length > 1) {
			tempDate = new Date(temp[dataCovar_Index[0]]);
			if( tempDate > startdate && tempDate <= tf){
				let tempDayCount = Math.ceil(Math.abs(tempDate - t0) / (1000 * 3600 * 24));
			  let tempcovar =	[tempDayCount];
				for(let j = 1; j < selected_colnames.length; j++){
					tempcovar.push(temp[dataCovar_Index[j]] == '' ? NaN :Number(temp[dataCovar_Index[j]]));
				}
				dataCovar.push(tempcovar);
			}
		}
  }
  
  total_deaths = Array.from(Array(timeLength)).fill(0);
  hospital = Array.from(Array(timeLength)).fill(NaN);
  ICU = Array.from(Array(timeLength)).fill(NaN);
  ventilator = Array.from(Array(timeLength)).fill(NaN);
  for(let i = 0; i < dataCovar.length; i++){
    total_deaths[dataCovar[i][0] - 1] = dataCovar[i][1];
    hospital[dataCovar[i][0] - 1] = dataCovar[i][2];
    ICU[dataCovar[i][0] - 1] = dataCovar[i][3];
    ventilator[dataCovar[i][0] - 1] = dataCovar[i][4];
  }
  deaths = total_deaths.map( (_,i,arr) => (arr[i]) - (arr[i-1]));
	deaths[0] = total_deaths[0];
	
  let data = time.map( (x,i) => [x, reports[i], deaths[i], hospital[i], ICU[i], ventilator[i]]);

  for(let i = 0; i < 23 - 1; i++){
    data[i][3] = 0;
    data[i][4] = 0;
    data[i][5] = 0;
  }
  
  data.unshift(["time","reports","deaths","hospital","ICU","ventilator"]);

	if (predTime != null) {
    predTime = new Date(predTime);
		let add_time = Math.ceil(Math.abs(tf - predTime) / (1000 * 3600 * 24));
		for(let i =0; i < add_time; i++)
			data.push([timeLength + i, NaN, NaN, NaN, NaN, NaN]);
	}
	return data;
}

module.exports = create_dataset;