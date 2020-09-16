/** generateSobol provides an array of objects of parameters with the length of 'SobolNumberOfPoints' considering lowerBounds/upperBounds.
 *  NOTE :The order of parameters in lower and upper bounds needs to be the same.
 */
const sobolSeq =require('../../library/generate-sobol/sobolSeq.js');


exports.sobolSet = function (SobolNumberOfPoints, args) {
  let iotaFixed = Math.log((args.TC-6)/((args.TC-6)-1));
  let lowerBounds = {betaI:0, theta:0, iota:iotaFixed, beta_sd:0,
      dI0:0, dP0:0, dT0:0, dB0:0,
      dI1:0, dP1:0, dT1:0, dB1:0,
      dI2:0, dP2:0, dT2:0, dB2:0,
      qP:0, qH:0, qC:0.5, mI:0, mC:0, mV:0.5,
      sigma:1/5, kappa:1/1, gammaI:1/5, gammaH:1/5, gammaC:1/10, gammaV:1/10,
      rho:0, TF:1e2,
      S0:1,EQ0:0,PQ0:0,IQ0:0,E0:0,P0:0,I0:0,H0:0,C0:0,V0:0,M0:0}
  
  let upperBounds = {betaI:2, theta:1, iota:iotaFixed, beta_sd:0,
    dI0:0.6, dP0:0.6, dT0:0.6, dB0:0,
    dI1:0.5, dP1:0.5, dT1:0.5, dB1:0,
    dI2:0.5, dP2:0.5, dT2:0.5, dB2:0,
    qP:0.5, qH:1, qC:1, mI:0.1, mC:1, mV:1,
    sigma:1/5, kappa:1/1, gammaI:1/5, gammaH:1/1, gammaC:1/1, gammaV:1/1,
    rho:1, TF:8e3,
    S0:1,EQ0:0,PQ0:0,IQ0:0,E0:0,P0:0,I0:0,H0:0,C0:0,V0:0,M0:0}
  
  return sobolSeq.sobolDesign(lowerBounds, upperBounds, SobolNumberOfPoints);
}