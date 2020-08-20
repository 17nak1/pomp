let run = function(){
  let msg = '';
    if(this.scale.length === 1){
        this.scale[0] = -Math.abs(this.scale[0]);
    }
    else{
        for(let i = 0; i< this.scale.length; i++){
            this.scale[i] = Math.abs(this.scale[i]);
        }
    }
    this.n = this.x0.length ;
    this.work = new Array((this.n)*(this.n+6)+1)
    this.iwork = new Array(2*(this.n))

    this.subplx (this.f,this.n,this.tol,this.maxnfe,this.scale.dArray(),this.x0.dArray(),this.fx,this.nfe,this.work.dArray(),this.iwork.dArray(),this.iflag)
    // console.log("FX is ===========" , this)
    switch (this.iflag) {
        case -1:
          msg = 'number of function evaluations exceeds \'maxit\'';
          console.log(msg);
          break;
        case 0:
          msg = 'success! tolerance satisfied';
          console.log(msg);
          break;
        case 1:
          msg = 'limit of machine precision reached';
          console.log(msg);
          break;
        case -2:
          msg = '\'parscale\' is too small relative to \'par\'';
          console.log(msg);
          break;
        case 2: default:
          msg = 'impossible error in subplex';
          console.log(msg);
          break;
        }
    return [this.x0, this.fx, this.nfe + 1, msg]
}

module.exports = run;