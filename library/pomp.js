
let pomp = class Pomp {
  constructor(args) {
      if(args == undefined){
          args = {};
      }
      
      this.data = args.data || [];
      this.times = args.times || [];
      this.t0 = args.t0 || 0;
      this.rprocess = args.rprocess || {};
      this.rprocessDetail = args.rprocessDetail || {};
      this.dprocess = args.dprocess || {};
      this.rmeasure = args.rmeasure || {};
      this.dmeasure = args.dmeasure || {};
      this.measurement_model = args.measurement_model || {};
      this.skeleton = args.skeleton || {};
      this.skeletonDetail = args.skeletonDetail || {};
      this.initializer = args.initializer || {};
      this.rprior = args.rprior || {};
      this.dprior = args.dprior || {};
      this.params = args.params || {};
      this.covar = args.covar || {};
      this.tcovar = args.tcovar || {};
      this.statenames = args.statenames || {};
      this.paramnames = args.paramnames || {};
      this.zeronames = args.zeronames || {};
      this.PACKAGE = args.PACKAGE || {};
      this.fromEstimationScale = args.fromEstimationScale || {};
      this.toEstimationScale = args.toEstimationScale || {};
      this.globals = args.globals || {};
      this.cdir = args.cdir || {};
      this.cfile = args.cfile || {};
      this.shlib_args = args.shlib_args || {};

      let ep = 'in \'pomp\': ';
      if (this.data == undefined || this.data.length == 0) 
          throw new Error(ep + '\'data\' is a required argument.');

      //## return as quickly as possible if no work is to be done
      if (Object.keys(args).length == 1) return(data);

      if (!Array.isArray(this.data) && this.data.constructor.name != 'Pomp')
          throw new Error(ep + '\'data\'must be a data frame or an object of class \'pomp\'.');

      //## if one transformation is supplied, then both must be
      let c1 = this.fromEstimationScale === undefined || this.fromEstimationScale == {};
      let c2 = this.toEstimationScale === undefined || this.toEstimationScale == {};
      if (c1 ^ c2)
          throw new Error(ep + 'if one of \'toEstimationScale\', \'fromEstimationScale\' is supplied, then so must the other');

      //## if 'covar' is supplied, then so must 'tcovar'
      c1 = this.covar === undefined || this.covar == {};
      c2 = this.tcovar === undefined || this.tcovar == {};
      if (c1 ^ c2)
          throw new Error(ep + 'if one of \'covar\', \'tcovar\' is supplied, then so must the other');
      

      if (this.times === undefined || this.times === [])  throw new Error(ep + '\'times\' is a required argument');
          if(!Array.isArray(this.times)){        
          if ((Number.isInteger(this.times) && (this.times<1 || this.times > this.data[0].length)) ||
              (typeof this.times === 'string' && this.data[0].indexOf(this.times) == -1)) {
                  throw new Error(ep + 'when \'data\' is a data frame, \'times\' must identify a single column of \'data\' either by name or by index.');            
          }
          let tpos = -1;
          if (Number.isInteger(this.times)) {
              tpos = this.times;
          } else if (typeof this.times === 'string') {
              tpos = this.data[0].indexOf(this.times);
          }
          this.times = this.data.map(function(value,index) { return value[tpos]; });
          this.data.map( function(value) { return value.splice(tpos, 1); });
          this.obsnames = this.data.shift();
          this.times.shift();
      }
  }
}


pomp.prototype.interpolator = function(t){
  let n = this.covar.length;
  let covarnames = Object.keys(this.covar[0]);
  let nCovarnames = covarnames.length; 
      
  let point = {}; 
  if (n === 0) {
    for(let k = 0; k < nCovarnames; k++)
        point[covarnames[k]] = 0;
  }

  if (n === 1) {
    for(let k = 0; k < nCovarnames; k++)
        point[covarnames[k]] = +this.covar[0][covarnames[k]];
  }

  if (t <= this.tcovar[0]) {
    for(let k = 0; k < nCovarnames; k++) {
        point[covarnames[k]] = +this.covar[0][covarnames[k]] + (t - this.tcovar[0]) * (this.covar[1][covarnames[k]] - this.covar[0][covarnames[k]]) / (this.tcovar[1] -this.tcovar[0])
    }
  } else if(t >= this.tcovar[n - 1]){
    for(let k = 0; k < nCovarnames; k++) {
        point[covarnames[k]] = +this.covar[n-1][covarnames[k]] + (t - this.tcovar[n-1]) * (this.covar[n-2][covarnames[k]] - this.covar[n-1][covarnames[k]]) / (this.tcovar[n-2] -this.tcovar[n-1])
    }
  } else {
      for (let i = 0; i < n; i++) {
          if (t > this.tcovar[i] && t <= this.tcovar[i+1]) {
            for(let k = 0; k < nCovarnames; k++) {
                point[covarnames[k]] = +this.covar[i][covarnames[k]] + (t - this.tcovar[i]) * (this.covar[i+1][covarnames[k]] - this.covar[i][covarnames[k]]) / (this.tcovar[i+1] -this.tcovar[i])
            }
            return point;
          }
      }
    }
  return point;
};

module.exports = pomp;