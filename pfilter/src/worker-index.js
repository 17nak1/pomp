const { pfilter } = require('./pfilter.js');

// Add the workerfn to the self namespace so it can be called by whatever wraps the output bundle
const workerfn = self.workerfn = async (...args) => {
  const work = (resolve, reject) => () => {
    let result;
    let date = new Date();
    try {
      if (typeof progress === 'function') progress();
      // console.log("STARTING PFILTER");
      result = pfilter(...args);
    } catch (e) {
      console.error(e);
      result =  NaN;
    }
    console.debug("Calculation time: " , new Date() - date)
    resolve(result);
  }

  return await new Promise((resolve, reject) => {
    work(resolve, reject)();
  });
}

export default workerfn;