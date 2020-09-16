const { trajMatch } = require('./trajMatch.js');

// Add the workerfn to the self namespace so it can be called by whatever wraps the output bundle
const workerfn = self.workerfn = async (...args) => {
  const work = (resolve, reject) => () => {
    let result;
    let date = Date.now();
    try {
      if (typeof progress === 'function') progress();
      result = trajMatch(...args);
    } catch (e) {
      console.error(e);
      result =  NaN;
    }
    console.debug("Calculation time in trajMatch: " , Date.now() - date, result.loglik)
    resolve(result);
  }

  return await new Promise((resolve, reject) => {
    work(resolve, reject)();
  });
}

export default workerfn;