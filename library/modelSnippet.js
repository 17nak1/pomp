/**
 *  @file        ModelSnippet_StochModel3.js
 *               Makes the model and its dependencies.
 *               Note: To define rprocess, dprocess, dmeasure and initz
 *               you should use the same input structure that is defined here. 
 *                 
 *  @author       Nazila Akhavan, nazila@kingsds.network
 *  @date        june 2020
 */

let snippet = {}
let mathLib = require('./mathLib');
let { rpois } = require('./rpois');

snippet.skeleton = function (states, params, t, dt, covar, args) {
  let d ={};
  let SUSC  = [states.S];
  let EXPD  = [states.E1, states.E2, states.E3];
  let PRE   = [states.P1, states.P2, states.P3];
  let INFD  = [states.I1, states.I2, states.I3];
  let HOSP  = [states.H1, states.H2, states.H3];
  let CARE  = [states.C1, states.C2, states.C3];
  let VENT  = [states.V1, states.V2, states.V3];
  let DEAD  = [states.M];
  let RCVD  = [states.R];
  let EXPDQ = [states.EQ1, states.EQ2, states.EQ3];
  let PREQ  = [states.PQ1, states.PQ2, states.PQ3];
  let INFDQ = [states.IQ1, states.IQ2, states.IQ3];
  
  
  // Different transmission rates
  let dQdt;
  let TOT_PRE = 0;
  for(let i = 0; i < args.nstageP; i++) {
    TOT_PRE += PRE[i];
  }

  let TOT_INFD = 0;
  for(i = 0; i < args.nstageI; i++) {
  TOT_INFD += INFD[i];
  }
  
  let PD, lambdaI, lambdaP, lambdaPQ, lambda, lambdaQ;
  
  PD = params.rho * (covar.tests / (covar.tests + params.TF));
  
  lambdaI = params.betaI * TOT_INFD;
  lambdaP = params.betaI * params.theta * TOT_PRE * (1-PD);
  lambdaPQ = params.betaI * params.theta * TOT_PRE * PD;
  
  if (t < args.T0) {
    dQdt  = mathLib.rgammawn(params.beta_sd, dt)/dt;
    lambda = ( (lambdaI + lambdaP + params.iota) / args.pop ) * dQdt;
    lambdaQ = ( (lambdaPQ) / args.pop ) * dQdt;
  } else if (t < args.T0 + 7.0) {
    let x = (t-args.T0) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = mathLib.rgammawn(((1-ss) + ss*params.dB0)*params.beta_sd, dt)/dt;
    lambda = ( ( ((1-ss) + ss*params.dI0) * lambdaI + ((1-ss) + ss*params.dP0) * lambdaP + ((1-ss) + ss*params.dT0)*params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( ( ((1-ss) + ss*params.dP0) * lambdaPQ  ) / args.pop ) * dQdt;
  } else if (t < args.T1) {
    dQdt  = mathLib.rgammawn(params.dB0*params.beta_sd, dt)/dt;
    lambda = ( ( params.dI0 * lambdaI + params.dP0 * lambdaP + params.dT0 * params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( (  params.dP0 * lambdaPQ  ) / args.pop ) * dQdt;
  } else if (t < args.T1 +7.0) {
    let x = (t-args.T1) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = mathLib.rgammawn(((1-ss)*params.dB0 + ss*params.dB1)*params.beta_sd, dt)/dt;
    lambda = ( ( ((1-ss)*params.dI0 + ss*params.dI1) * lambdaI + ((1-ss)*params.dP0 + ss*params.dP1) * lambdaP + ((1-ss)*params.dT0 + ss*params.dT1)*params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( ( ((1-ss)*params.dP0 + ss*params.dP1) * lambdaPQ  ) / args.pop ) * dQdt;
  } else if (t < args.T2) {	
    dQdt  = mathLib.rgammawn(params.dB1*params.beta_sd, dt)/dt;	
    lambda = ( ( params.dI1 * lambdaI + params.dP1 * lambdaP + params.dT1 * params.iota ) / args.pop ) * dQdt;	
    lambdaQ = ( (  params.dP1 * lambdaPQ  ) / args.pop ) * dQdt;	
  } else if (t < args.T2 + 7.0) {	
    let x = (t-args.T2) / 7.0;	
    let ss = 3*x*x - 2*x*x*x;	
    dQdt  = mathLib.rgammawn(((1-ss)*params.dB1 + ss* params.dB2)*params.beta_sd, dt)/dt;	
    lambda = ( ( ((1-ss)*params.dI1 + ss * params.dI2) * lambdaI + ((1-ss)*params.dP1 + ss*params.dP2) * lambdaP + ((1-ss)*params.dT1 + ss*args.T2)*params.iota ) / args.pop ) * dQdt;	
    lambdaQ = ( ( ((1-ss)*params.dP1 + ss*params.dP2) * lambdaPQ  ) / args.pop ) * dQdt;
  } else {
    dQdt  = mathLib.rgammawn(params.dB2*params.beta_sd, dt)/dt;
    lambda = ( ( params.dI2 * lambdaI + params.dP2 * lambdaP + params.dT2 * params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( (  params.dP2 * lambdaPQ  ) / args.pop ) * dQdt;
  }
  
  // From class S
  let transS = new Array(2);
  let rateS = new Array(2);
  rateS[0] = lambda ;
  rateS[1] = lambdaQ;
  let rateS_tot = rateS[0] + rateS[1];
  let transS_tot = (1 - Math.exp(- rateS_tot * dt)) * SUSC[0];
  transS[0] = rateS[0] / rateS_tot * transS_tot;
  transS[1] = rateS[1] / rateS_tot * transS_tot;
  
  // From class EQ
  let transEQ = new Array(args.nstageE);
  let rateEQ = args.nstageE * params.sigma;
  for (let i = 0; i < args.nstageE; i++) {
    transEQ[i] =  (1 - Math.exp(- rateEQ * dt))* EXPDQ[i];
  }
  
  // From class PQ
  let transPQ = new Array(args.nstageP + 1);
  let ratePQ = args.nstageP * params.kappa;
  for (let i = 0; i < args.nstageP-1; i++) {
    transPQ[i] =  (1 - Math.exp(- ratePQ * dt))* PREQ[i];
  }
  let transPQIQH;
  transPQIQH = (1 - Math.exp(- ratePQ * dt))* PREQ[args.nstageP-1];
  transPQ[args.nstageP-1] = (1-params.qP) * transPQIQH;
  transPQ[args.nstageP] = params.qP * transPQIQH;
  
  // From class IQ
  let transIQ = new Array(args.nstageI + 1);
  let rateIQ = args.nstageI * params.gammaI;
  for (let i = 0; i < args.nstageI-1; i++) {
    transIQ[i] =  (1 - Math.exp(- rateIQ * dt))* INFDQ[i];
  }
  let transIQRD;
  transIQRD = (1 - Math.exp(- rateIQ * dt))* INFDQ[args.nstageI-1];
  transIQ[args.nstageI-1] = (1 - params.mI) * transIQRD;
  transIQ[args.nstageI] = params.mI * transIQRD;
    
  // From class E
  let transE = new Array(args.nstageE);
  let rateE = args.nstageE * params.sigma;
  for (let i = 0; i < args.nstageE; i++) {
    transE[i] =  (1 - Math.exp(- rateE * dt))* EXPD[i];
  }
  
  // From class P
  let transP = new Array(args.nstageP + 2);
  let rateP = args.nstageP * params.kappa;
  for (let i = 0; i < args.nstageP-1; i++) {
    transP[i] =  (1 - Math.exp(- rateP * dt))* PRE[i];
  }
  let transPIHIQ;
  transPIHIQ = (1.0 - Math.exp(- rateP * dt))* PRE[args.nstageP-1];
  transP[args.nstageP-1] = (1 - PD) * (1 - params.qP) * transPIHIQ;
  transP[args.nstageP] = params.qP * transPIHIQ;
  transP[args.nstageP+1] = PD * (1 - params.qP) * transPIHIQ;
  
  // From class I
  let transI = new Array(args.nstageI + 1);
  let rateI = args.nstageI * params.gammaI;
  for (let i = 0; i < args.nstageI-1; i++) {
    transI[i] =  (1 - Math.exp(- rateI * dt))* INFD[i];
  }
  
  let transIRD;
  transIRD = (1 - Math.exp(- rateI * dt))* INFD[args.nstageI-1];
  transI[args.nstageI-1] = (1 - params.mI) * transIRD;
  transI[args.nstageI] = params.mI * transIRD;
  
  // From class H
  let transH = new Array(args.nstageH + 1);
  let rateH = args.nstageH * params.gammaH;
  for (let i = 0; i < args.nstageH-1; i++) {
    transH[i] =  (1.0 - Math.exp(- rateH * dt))* HOSP[i];
  }
  
  let transHRC;
  transHRC = (1.0 - Math.exp(- rateH * dt))* HOSP[args.nstageH-1];
  transH[args.nstageH-1] = (1 - params.qH) * transHRC;
  transH[args.nstageH] = params.qH * transHRC;
  
  // From class C
  let transC= new Array(args.nstageC + 2);
  let rateC = args.nstageC * params.gammaC;
  for (let i = 0; i < args.nstageC-1; i++) {
    transC[i] =  (1.0 - Math.exp(- rateC * dt))* CARE[i];
  }
  let transCRVM;
  transCRVM = (1.0 - Math.exp(- rateC * dt))* CARE[args.nstageC-1];
  transC[args.nstageC-1] =  (1 - params.mC) * (1 - params.qC) * transCRVM;
  transC[args.nstageC] = params.qC * transCRVM;
  transC[args.nstageC+1] = params.mC * (1 - params.qC) * transCRVM;
  
  // From class V
  let transV = new Array(args.nstageV + 1);
  let rateV = args.nstageV * params.gammaV;
  for (let i = 0; i < args.nstageV-1; i++) {
    transV[i] =  (1 - Math.exp(- rateV * dt))* VENT[i];
  }
  let transVRD;
  transVRD = (1.0 - Math.exp(- rateV * dt))* VENT[args.nstageV-1];
  transV[args.nstageV-1] = (1 - params.mV) * transVRD;
  transV[args.nstageV] = params.mV * transVRD;
  
  // Balance the equations
  d.casesI = states.casesI;
  d.casesIQ = states.casesIQ;
  d.casesH = states.casesH;
  d.deathsIIQ = states.deathsIIQ;
  d.deathsCV = states.deathsCV;
  
  // Balance the equations
  SUSC[0] -= transS[0];
  EXPD[0] += transS[0];
  for (i = 0; i < args.nstageE; i++) {
    EXPD[i] -= transE[i];
    if (i > 0) EXPD[i] += transE[i-1];
  }
  
  PRE[0] += transE[args.nstageE-1];
  for (i = 0; i < args.nstageP; i++) {
    PRE[i] -= transP[i];
    if (i > 0) PRE[i] += transP[i-1];
  }
  
  INFD[0] += transP[args.nstageP-1];
  HOSP[0] += transP[args.nstageP];
  PRE[args.nstageP-1] -= transP[args.nstageP];
  for (i = 0; i < args.nstageI; i++) {
    INFD[i] -= transI[i];
    if(i > 0) INFD[i] += transI[i-1];
  }
  
  for (i = 0; i < args.nstageH; i++) {
    HOSP[i] -= transH[i];
    if (i > 0) HOSP[i] += transH[i-1];
  }
  
  CARE[0] += transH[args.nstageH];
  HOSP[args.nstageH-1] -= transH[args.nstageH];
  INFD[args.nstageI-1] -= transI[args.nstageI];
  for (i = 0; i < args.nstageC; i++) {
    CARE[i] -= transC[i];
    if (i > 0) CARE[i] += transC[i-1];
  }
  
  VENT[0] += transC[args.nstageC];
  CARE[args.nstageC-1] -= transC[args.nstageC] + transC[args.nstageC+1];
  for (i = 0; i < args.nstageV; i++) {
    VENT[i] -= transV[i];
    if (i > 0) VENT[i] += transV[i-1];
  }
  
  VENT[args.nstageV-1] -= transV[args.nstageV];               
  RCVD[0] += transI[args.nstageI-1] + transH[args.nstageH-1] + transC[args.nstageC-1] + transV[args.nstageV-1];
  DEAD[0] += transI[args.nstageI] + transC[args.nstageC+1] + transV[args.nstageV];
  
  SUSC[0] -= transS[1];
  EXPDQ[0] += transS[1];
  for (i = 0; i < args.nstageE; i++) {
    EXPDQ[i] -= transEQ[i];
    if (i > 0) EXPDQ[i] += transEQ[i-1];
  }
  
  PREQ[0] += transEQ[args.nstageE-1];
  for (i = 0; i < args.nstageP; i++) {
    PREQ[i] -= transPQ[i];
    if (i > 0) PREQ[i] += transPQ[i-1];
  }
  
  PREQ[args.nstageP-1] -= transPQ[args.nstageP];  
  PRE[args.nstageP-1] -= transP[args.nstageP+1];
  INFDQ[0] += transPQ[args.nstageP-1] + transP[args.nstageP+1];
  HOSP[0] += transPQ[args.nstageP];
  for (i = 0; i < args.nstageI; i++) {
    INFDQ[i] -= transIQ[i];
    if (i > 0) INFDQ[i] += transIQ[i-1];
  }
  
  INFDQ[args.nstageI-1] -= transIQ[args.nstageI];
  RCVD[0] += transIQ[args.nstageI-1];
  DEAD[0] += transIQ[args.nstageI];
  
  d.casesI += transP[args.nstageP-1];
  d.casesIQ += transPQ[args.nstageP-1] + transP[args.nstageP+1];
  d.casesH += transPQ[args.nstageP] + transP[args.nstageP];
  d.deathsIIQ += transI[args.nstageI] + transIQ[args.nstageI];
  d.deathsCV += transC[args.nstageC+1] + transV[args.nstageV];

  //&DS;
  d['S'] = SUSC[0];
  //&DE1;
  d['E1'] = EXPD[0];
  d['E2'] = EXPD[1];
  d['E3'] = EXPD[2];
  //&DP1;
  d['P1'] = PRE[0];
  d['P2'] = PRE[1];
  d['P3'] = PRE[2];
  //&DI1;
  d['I1'] = INFD[0];
  d['I2'] = INFD[1];
  d['I3'] = INFD[2];
  //&DH1;
  d['H1'] = HOSP[0];
  d['H2'] = HOSP[1];
  d['H3'] = HOSP[2];
  //&DC1;
  d['C1'] = CARE[0];
  d['C2'] = CARE[1];
  d['C3'] = CARE[2];
  //&DV1;
  d['V1'] = VENT[0];
  d['V2'] = VENT[1];
  d['V3'] = VENT[2];
  //&DM;
  d['M'] = DEAD[0];
  //&DR;
  d['R'] = RCVD[0];
  //&DEQ1;
  d['EQ1'] = EXPDQ[0];
  d['EQ2'] = EXPDQ[1];
  d['EQ3'] = EXPDQ[2];
  //&DPQ1;
  d['PQ1'] = PREQ[0];
  d['PQ2'] = PREQ[1];
  d['PQ3'] = PREQ[2];
  // &DIQ1;
  d['IQ1'] = INFDQ[0];
  d['IQ2'] = INFDQ[1];
  d['IQ3'] = INFDQ[2];

  return d;
}

snippet.rprocess = function (states, params, t, dt, covar, args) {
  
  let SUSC  = [states.S];
  let EXPD  = [states.E1, states.E2, states.E3];
  let PRE   = [states.P1, states.P2, states.P3];
  let INFD  = [states.I1, states.I2, states.I3];
  let HOSP  = [states.H1, states.H2, states.H3];
  let CARE  = [states.C1, states.C2, states.C3];
  let VENT  = [states.V1, states.V2, states.V3];
  let DEAD  = [states.M];
  let RCVD  = [states.R];
  let EXPDQ = [states.EQ1, states.EQ2, states.EQ3];
  let PREQ  = [states.PQ1, states.PQ2, states.PQ3];
  let INFDQ = [states.IQ1, states.IQ2, states.IQ3];
   
  
  // Different transmission rates

  let TOT_PRE = 0;
  for(let i = 0 ; i < args.nstageP; i++) {
    TOT_PRE += PRE[i];
  }

  let TOT_INFD = 0;
  for(let i = 0; i < args.nstageI; i++) {
    TOT_INFD += INFD[i];
  }
  
  let PD;
  
  PD = params.rho * (covar.tests / (covar.tests + params.TF));
  
  let lambdaI = params.betaI * TOT_INFD;
  let lambdaP = params.betaI * params.theta * TOT_PRE * (1-PD);
  let lambdaPQ = params.betaI * params.theta * TOT_PRE * PD;
  let dQdt, lambda, lambdaQ;
  if (t < args.T0) {
    dQdt  = mathLib.rgammawn(params.beta_sd, dt)/dt;
    lambda = ( (lambdaI + lambdaP + params.iota) / args.pop ) * dQdt;
    lambdaQ = ( (lambdaPQ) / args.pop ) * dQdt;
  } else if (t < args.T0 + 7.0) {
    let x = (t-args.T0) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = mathLib.rgammawn(((1-ss) + ss*params.dB0)*params.beta_sd, dt)/dt;
    lambda = ( ( ((1-ss) + ss*params.dI0) * lambdaI + ((1-ss) + ss*params.dP0) * lambdaP + ((1-ss) + ss*params.dT0)*params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( ( ((1-ss) + ss*params.dP0) * lambdaPQ  ) / args.pop ) * dQdt;
  } else if (t < args.T1) {
    dQdt  = mathLib.rgammawn(params.dB0*params.beta_sd, dt)/dt;
    lambda = ( ( params.dI0 * lambdaI + params.dP0 * lambdaP + params.dT0 * params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( (  params.dP0 * lambdaPQ  ) / args.pop ) * dQdt;
  } else if (t < args.T1 +7.0) {
    let x = (t-args.T1) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = mathLib.rgammawn(((1-ss)*params.dB0 + ss*params.dB1)*params.beta_sd, dt)/dt;
    lambda = ( ( ((1-ss)*params.dI0 + ss*params.dI1) * lambdaI + ((1-ss)*params.dP0 + ss*params.dP1) * lambdaP + ((1-ss)*params.dT0 + ss*params.dT1)*params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( ( ((1-ss)*params.dP0 + ss*params.dP1) * lambdaPQ  ) / args.pop ) * dQdt;
  } else if (t < args.T2) {	
    dQdt  = mathLib.rgammawn(params.dB1*params.beta_sd, dt)/dt;	
    lambda = ( ( params.dI1 * lambdaI + params.dP1 * lambdaP + params.dT1 * params.iota ) / args.pop ) * dQdt;	
    lambdaQ = ( (  params.dP1 * lambdaPQ  ) / args.pop ) * dQdt;	
  } else if (t < args.T2 + 7.0) {	
    let x = (t-args.T2) / 7.0;	
    let ss = 3*x*x - 2*x*x*x;	
    dQdt  = mathLib.rgammawn(((1-ss)*params.dB1 + ss* params.dB2)*params.beta_sd, dt)/dt;	
    lambda = ( ( ((1-ss)*params.dI1 + ss * params.dI2) * lambdaI + ((1-ss)*params.dP1 + ss*params.dP2) * lambdaP + ((1-ss)*params.dT1 + ss*args.T2)*params.iota ) / args.pop ) * dQdt;	
    lambdaQ = ( ( ((1-ss)*params.dP1 + ss*params.dP2) * lambdaPQ  ) / args.pop ) * dQdt;
  } else {
    dQdt  = mathLib.rgammawn(params.dB2*params.beta_sd, dt)/dt;
    lambda = ( ( params.dI2 * lambdaI + params.dP2 * lambdaP + params.dT2 * params.iota ) / args.pop ) * dQdt;
    lambdaQ = ( (  params.dP2 * lambdaPQ  ) / args.pop ) * dQdt;
  }
  
  // From class S
  let transS = new Array(2);
  let rateS = new Array(2);
  rateS[0] = lambda ;
  rateS[1] = lambdaQ;
  mathLib.reulermultinom(2, Math.round(SUSC[0]), 0, dt, 0, rateS, transS);
  
  // From class EQ
  let  transEQ = new Array(args.nstageE);
  let  rateEQ = [args.nstageE * params.sigma];
  for (let i = 0; i < args.nstageE; i++) {
    mathLib.reulermultinom(1, Math.round(EXPDQ[i]), 0, dt, i, rateEQ, transEQ);
  }
  
  // From class PQ
  let transPQ = new Array(args.nstageP + 1);
  let ratePQ = [args.nstageP * params.kappa];
  for (i = 0; i < args.nstageP - 1; i++) {
    mathLib.reulermultinom(1, Math.round(PREQ[i]), 0, dt, i, ratePQ, transPQ);
  }
  let ratePQIQH = new Array(2);
  ratePQIQH[0] = (1 - params.qP) * args.nstageP * params.kappa;
  ratePQIQH[1] = params.qP * args.nstageP * params.kappa;
  mathLib.reulermultinom(2, Math.round(PREQ[args.nstageP - 1]), 0, dt, args.nstageP - 1, ratePQIQH, transPQ);
  
  // From class IQ
  let transIQ = new Array(args.nstageI + 1);
  let rateIQ = [args.nstageI * params.gammaI];
  for (let i = 0; i < args.nstageI - 1; i++) {
    mathLib.reulermultinom(1, Math.round(INFDQ[i]), 0, dt, i, rateIQ, transIQ);
  }
  
  let rateIQRD= new Array(2);
  rateIQRD[0] = (1 - params.mI) * args.nstageI * params.gammaI;
  rateIQRD[1] = params.mI * args.nstageI * params.gammaI;
  mathLib.reulermultinom(2, Math.round(INFDQ[args.nstageI-1]), 0, dt, args.nstageI - 1, rateIQRD, transIQ);

  
  // From class E
  let transE = new Array(args.nstageE);
  let rateE = [args.nstageE * params.sigma];
  for (let i = 0; i < args.nstageE; i++) {
    mathLib.reulermultinom(1, Math.round(EXPD[i]), 0, dt, i, rateE, transE);
  }
  
  // From class P
  let transP = new Array(args.nstageP + 2);
  let rateP = [args.nstageP * params.kappa];
  for (i = 0; i < args.nstageP - 1; i++) {
    mathLib.reulermultinom(1, Math.round(PRE[i]), 0, dt, i, rateP, transP);  
  }
  
  let ratePIHIQ = new Array(3);
  ratePIHIQ[0] = (1 - PD) * (1 - params.qP) * args.nstageP * params.kappa;
  ratePIHIQ[1] = params.qP * args.nstageP * params.kappa;
  ratePIHIQ[2] = PD * (1 - params.qP) * args.nstageP * params.kappa;
  mathLib.reulermultinom(3, Math.round(PRE[args.nstageP - 1]), 0, dt, args.nstageP - 1, ratePIHIQ, transP);

  
  // From class I
  let transI = new Array(args.nstageI + 1);
  let rateI = [args.nstageI * params.gammaI];
  for (let i = 0; i < args.nstageI - 1; i++) {
    mathLib.reulermultinom(1, Math.round(INFD[i]), 0, dt, i, rateI, transI);
  }
  
  let rateIRD = new Array(2);
  rateIRD[0] = (1 - params.mI) * args.nstageI * params.gammaI;
  rateIRD[1] = params.mI * args.nstageI * params.gammaI;
  mathLib.reulermultinom(2, Math.round(INFD[args.nstageI - 1]), 0, dt, args.nstageI - 1, rateIRD, transI);
  
  // From class H
  let transH = new Array(args.nstageH + 1);
  let rateH = [args.nstageH * params.gammaH];
  for (let i = 0; i < args.nstageH - 1; i++) {
    mathLib.reulermultinom(1, Math.round(HOSP[i]), 0, dt, i, rateH, transH);
  }

  let rateHRC = new Array(2);
  rateHRC[0] = (1 - params.qH) * args.nstageH * params.gammaH;
  rateHRC[1] = params.qH * args.nstageH * params.gammaH;
  mathLib.reulermultinom(2, Math.round(HOSP[args.nstageH - 1]), 0, dt, args.nstageH - 1, rateHRC, transH);  
  
  // From class C
  let transC = new Array(args.nstageC + 2);
  let rateC = [args.nstageC * params.gammaC];
  for (let i = 0; i < args.nstageC - 1; i++) {
    mathLib.reulermultinom(1, Math.round(CARE[i]), 0, dt, i, rateC, transC);
  }
  let rateCRVM = new Array(3);
  rateCRVM[0] = (1 - params.mC) * args.nstageC * params.gammaC;
  rateCRVM[1] = params.qC * args.nstageC * params.gammaC;
  rateCRVM[2] = params.mC * (1-params.qC) * args.nstageC * params.gammaC;
  mathLib.reulermultinom(3, Math.round(CARE[args.nstageC - 1]), 0, dt, args.nstageC - 1, rateCRVM, transC);
  
  // From class V
  let transV = new Array(args.nstageV + 1);
  let rateV = [args.nstageV * params.gammaV];
  for (i = 0; i < args.nstageV - 1; i++) {
    mathLib.reulermultinom(1, Math.round(VENT[i]), 0, dt, i, rateV, transV);
  }
  let rateVRD = new Array(2);
  rateVRD[0] = (1 - params.mV) * args.nstageV * params.gammaV;
  rateVRD[1] = params.mV * args.nstageV * params.gammaV;
  mathLib.reulermultinom(2, Math.round(VENT[args.nstageV - 1]), 0, dt, args.nstageV - 1, rateVRD, transV);

  
  // Balance the equations
  SUSC[0] -= transS[0];
  EXPD[0] += transS[0];
  for (i = 0; i < args.nstageE; i++) {
    EXPD[i] -= transE[i];
    if (i > 0) EXPD[i] += transE[i-1];
  }

  PRE[0] += transE[args.nstageE-1];
  for (i = 0; i < args.nstageP; i++) {
    PRE[i] -= transP[i];
    if (i > 0) PRE[i] += transP[i-1];
  }

  INFD[0] += transP[args.nstageP-1];
  HOSP[0] += transP[args.nstageP];
  PRE[args.nstageP-1] -= transP[args.nstageP];
  for (i = 0; i < args.nstageI; i++) {
    INFD[i] -= transI[i];
    if(i > 0) INFD[i] += transI[i-1];
  }
  
  for (i = 0; i < args.nstageH; i++) {
    HOSP[i] -= transH[i];
    if (i > 0) HOSP[i] += transH[i-1];
  }

  CARE[0] += transH[args.nstageH];
  HOSP[args.nstageH-1] -= transH[args.nstageH];
  INFD[args.nstageI-1] -= transI[args.nstageI];
  for (i = 0; i < args.nstageC; i++) {
    CARE[i] -= transC[i];
    if (i > 0) CARE[i] += transC[i-1];
  }

  VENT[0] += transC[args.nstageC];
  CARE[args.nstageC-1] -= transC[args.nstageC] + transC[args.nstageC+1];
  for (i = 0; i < args.nstageV; i++) {
    VENT[i] -= transV[i];
    if (i > 0) VENT[i] += transV[i-1];
  }

  VENT[args.nstageV-1] -= transV[args.nstageV];               
  RCVD[0] += transI[args.nstageI-1] + transH[args.nstageH-1] + transC[args.nstageC-1] + transV[args.nstageV-1];
  DEAD[0] += transI[args.nstageI] + transC[args.nstageC+1] + transV[args.nstageV];
  
  SUSC[0] -= transS[1];
  EXPDQ[0] += transS[1];
  for (i = 0; i < args.nstageE; i++) {
    EXPDQ[i] -= transEQ[i];
    if (i > 0) EXPDQ[i] += transEQ[i-1];
  }
  
  PREQ[0] += transEQ[args.nstageE-1];
  for (i = 0; i < args.nstageP; i++) {
    PREQ[i] -= transPQ[i];
    if (i > 0) PREQ[i] += transPQ[i-1];
  }
  
  PREQ[args.nstageP-1] -= transPQ[args.nstageP];  
  PRE[args.nstageP-1] -= transP[args.nstageP+1];
  INFDQ[0] += transPQ[args.nstageP-1] + transP[args.nstageP+1];
  HOSP[0] += transPQ[args.nstageP];
  for (i = 0; i < args.nstageI; i++) {
    INFDQ[i] -= transIQ[i];
    if (i > 0) INFDQ[i] += transIQ[i-1];
  }
  
  INFDQ[args.nstageI-1] -= transIQ[args.nstageI];
  RCVD[0] += transIQ[args.nstageI-1];
  DEAD[0] += transIQ[args.nstageI];
  
  states.casesI += transP[args.nstageP-1];
  states.casesIQ += transPQ[args.nstageP-1] + transP[args.nstageP+1];
  states.casesH += transP[args.nstageP] + transPQ[args.nstageP];
  states.deathsIIQ += transI[args.nstageI] + transIQ[args.nstageI];;
  states.deathsCV += transC[args.nstageC+1] + transV[args.nstageV];
  
  [states.S] = SUSC;
  [states.E1, states.E2, states.E3] = EXPD;
  [states.P1, states.P2, states.P3] = PRE;
  [states.I1, states.I2, states.I3] = INFD;
  [states.H1, states.H2, states.H3] = HOSP;
  [states.C1, states.C2, states.C3] = CARE;
  [states.V1, states.V2, states.V3] = VENT;
  [states.M] = DEAD;
  [states.R] = RCVD;
  
  [states.EQ1, states.EQ2, states.EQ3] = EXPDQ;
  [states.PQ1, states.PQ2, states.PQ3] = PREQ;
  [states.IQ1, states.IQ2, states.IQ3] = INFDQ;
  return states;
}

snippet.initializer = function(params, covar, args) {
  let initObj = {};
  
  let fS = params.S0;
  let fEQ = params.EQ0 / args.nstageE;
  let fPQ = params.PQ0 / args.nstageP;
  let fIQ = params.IQ0 / args.nstageI;
  let fE = params.E0 / args.nstageE;
  let fP = params.P0 / args.nstageP;
  let fI = params.I0 / args.nstageI;
  let fH = params.H0 / args.nstageH;
  let fC = params.C0 / args.nstageC;
  let fV = params.V0 / args.nstageV;
  let fM = params.M0;
  let fR = 1 - fS - args.nstageE*fE - args.nstageP*fP - args.nstageI*fI - args.nstageH*fH - args.nstageC*fC - args.nstageV*fV - args.nstageE*fEQ - args.nstageP*fPQ - args.nstageI*fIQ - fM;
  
  initObj["S"] = Math.round(args.pop*fS);
  for (let i = 0; i < args.nstageE; i++) initObj[`EQ${i + 1}`] = Math.round(args.pop*fEQ);
  for (let i = 0; i < args.nstageP; i++) initObj[`PQ${i + 1}`] = Math.round(args.pop*fPQ);
  for (let i = 0; i < args.nstageI; i++) initObj[`IQ${i + 1}`] = Math.round(args.pop*fIQ);
  for (let i = 0; i < args.nstageE; i++) initObj[`E${i + 1}`] = Math.round(args.pop*fE);
  for (let i = 0; i < args.nstageP; i++) initObj[`P${i + 1}`] = Math.round(args.pop*fP);
  for (let i = 0; i < args.nstageI; i++) initObj[`I${i + 1}`] = Math.round(args.pop*fI);
  for (let i = 0; i < args.nstageH; i++) initObj[`H${i + 1}`] = Math.round(args.pop*fH);
  for (let i = 0; i < args.nstageC; i++) initObj[`C${i + 1}`] = Math.round(args.pop*fC);
  for (let i = 0; i < args.nstageV; i++) initObj[`V${i + 1}`] = Math.round(args.pop*fV);
  initObj["M"] = Math.round(args.pop*fM);
  initObj["R"] = Math.round(args.pop*fR);
  
  initObj["casesI"] = 0;
  initObj["casesIQ"] = 0;
  initObj["casesH"] = 0;
  initObj["deathsIIQ"] = 0;
  initObj["deathsCV"] = 0;
  
  return initObj;
}

snippet.dmeasure = function (data ,states, params, giveLog = 1, args) {
  let HOSP = [], CARE = [], VENT = [];
  for(let i = 0; i < args.nstageH; i++) HOSP.push(states[`H${i + 1}`]);
  for(let i = 0; i < args.nstageC; i++) CARE.push(states[`C${i + 1}`]);
  for(let i = 0; i < args.nstageV; i++) VENT.push(states[`V${i + 1}`]);
  
  let lik_reports, lik_deaths, lik_hospital, lik_ICU, lik_ventilator;
  
  if (isFinite(data.reports)) {
    let reported_cases = states.casesH + states.casesIQ;
    lik_reports = mathLib.dpois(data.reports, reported_cases + 1e-6, 1);
  } else {
    lik_reports = 0;
  }
  
  if (isFinite(data.deaths)) {
    let reported_deaths = states.deathsCV + states.deathsIIQ;
    lik_deaths = mathLib.dpois(data.deaths, reported_deaths + 1e-6, 1); 
  } else {
    lik_deaths = 0;
  }
  
  let TOT_H = 0;
  for(let i = 0; i < args.nstageH; i++) {
    TOT_H += HOSP[i];
  }

  let TOT_C = 0;
  for(let i = 0; i < args.nstageC; i++) {
    TOT_C += CARE[i];
  }

  let TOT_V = 0;
  for(let i = 0; i < args.nstageV; i++) {
    TOT_V += VENT[i];
  }

  if (isFinite(data.hospital) & isFinite(data.ICU)) {
    lik_hospital = mathLib.dpois(data.hospital - data.ICU, TOT_H + 1e-6, 1);
  } else {
    lik_hospital = 0;
  }
  if (isFinite(data.ICU) & isFinite(data.ventilator)) {
    lik_ICU = mathLib.dpois(data.ICU - data.ventilator, TOT_C + 1e-6, 1);
  } else {
    lik_ICU = 0;
  }
  if (isFinite(data.ventilator)) {
    lik_ventilator = mathLib.dpois(data.ventilator, TOT_V + 1e-6, 1);
  } else {
    lik_ventilator = 0;
  }
  
  let lik = lik_reports + lik_deaths + lik_hospital + lik_ICU + lik_ventilator;
  
  if (giveLog == 0 ) {
    lik = Math.exp(lik);
  }
  
  return lik;
}

snippet.rmeasure = function (states, params) {
  let HOSP = [], CARE = [], VENT = [];
  let results = {};
  for(let i = 0; i < args.nstageH; i++) HOSP.push(states[`H${i + 1}`]);
  for(let i = 0; i < args.nstageC; i++) HOSP.push(states[`C${i + 1}`]);
  for(let i = 0; i < args.nstageV; i++) HOSP.push(states[`V${i + 1}`]);
  
  let reported_cases = states.casesH + states.casesIQ;
  results.reports = rpois(reported_cases  + 1e-6);
  
  let reported_deaths = states.deathsCV + states.deathsIIQ;
  results.deaths = rpois(reported_deaths  + 1e-6);
  
  let TOT_H = 0;
  for(let i = 0; i < args.nstageH; i++) {
  TOT_H += HOSP[i];
  }

  let TOT_C = 0;
  for(i = 0; i < args.nstageC; i++) {
  TOT_C += CARE[i];
  }

  let TOT_V = 0;         
  for(i = 0; i < args.nstageV; i++) {
  TOT_V += VENT[i];
  }

  results.ventilator = rpois(TOT_V + 1e-6);
  results.ICU = rpois(TOT_C + 1e-6) + results.ventilator;
  results.hospital = rpois(TOT_H + 1e-6) + results.ICU;

}

let params_log = ["betaI", "iota","beta_sd", "sigma", "kappa", "gammaI", "gammaH", "gammaC", "gammaV","TF"];
let params_logit = ["rho", "theta",
                   "dI0", "dP0", "dT0", "dB0",
                   "dI1", "dP1", "dT1", "dB1",
                   "dI2", "dP2", "dT2", "dB2",
                   "qP", "qH", "qC", "mI", "mC", "mV"];
let statenamesFn = function() {
  let nstageE = 3;
  let nstageP = 3;
  let nstageI = 3;
  let nstageH = 3;
  let nstageC = 3;
  let nstageV = 3;
  let sn = ["S"]
  for (let i = 0; i < nstageE; i++) sn.push("EQ"+i);
  for (let i = 0; i < nstageP; i++) sn.push("PQ"+i);
  for (let i = 0; i < nstageI; i++) sn.push("IQ"+i);
  for (let i = 0; i < nstageE; i++) sn.push("E"+i);
  for (let i = 0; i < nstageP; i++) sn.push("P"+i);
  for (let i = 0; i < nstageI; i++) sn.push("I"+i);
  for (let i = 0; i < nstageH; i++) sn.push("H"+i);
  for (let i = 0; i < nstageC; i++) sn.push("C"+i);
  for (let i = 0; i < nstageV; i++) sn.push("V"+i);
  sn.push("M", "R");
  return sn;
}    

snippet.paramsMod = [...params_log, ...params_logit];
snippet.paramsIc = ["S0","EQ0", "PQ0", "IQ0", "E0", "P0", "I0", "H0", "C0", "V0", "M0"];
snippet.zeronames = ["casesI", "casesIQ", "casesH", "deathsIIQ", "deathsCV"];
snippet.statenames = [...statenamesFn(), ...snippet.zeronames];

snippet.toEstimationScale = function(params) {
  let estiParams = Object.assign({}, params);
  for (let i = 0; i < params_log.length; i++) {
    estiParams[params_log[i]] = Math.log(params[params_log[i]]);
  }

  for (let i = 0; i < params_logit.length; i++) {
    estiParams[params_logit[i]] = mathLib.qlogis(params[params_logit[i]]);
  }
  
  return estiParams;
}

snippet.fromEstimationScale = function(params) {
  let estiParams = Object.assign({}, params);
  for (let i = 0; i < params_log.length; i++) {
    estiParams[params_log[i]] = Math.exp(params[params_log[i]]);
  }

  for (let i = 0; i < params_logit.length; i++) {
    estiParams[params_logit[i]] = mathLib.plogis(params[params_logit[i]]);
  }
  
  return estiParams;
}

snippet.determineRW = function() {
    return ((time) => {
      let rw_size = 0.05;
      let T0 = 75;
      let T1 = 139;
      let T2 = 212;
      let d = {};
      for (let i = 0; i < Object.keys(snippet.paramsMod).length; i++) {
        d[snippet.paramsMod[i]] = 0;
      }
      d.beta_sd = time < T0 ? rw_size : 0;
      d.dB0 = T0 < time < T1 ? rw_size : 0;
      d.dB1 = T1 < time < T2 ? 0 : rw_size;
      d.dB2 = time < T2 ? 0 : rw_size;
      return d;
    }).toString();
}


module.exports = snippet
