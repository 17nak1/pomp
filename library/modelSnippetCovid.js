/**
 *  @file        modeSnippet.js
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

let nstageE = 3;
let nstageP = 3;
let nstageI = 3;
let nstageH = 3;
let nstageC = 3;
let nstageV = 3;
let pop = 10e6;
let t1 = 75;

snippet.rprocess = function (states, params, t, dt, covar) {
  
  let SUSC = [states.S];
  let DEAD = [states.M];
  let RCVD = [states.R];

  let EXPD = []; PRE = []; INFD = []; HOSP = []; CARE = []; VENT = [];
  for(let i = 0; i < nstageE; i++) EXPD.push(states[`E${i + 1}`]);
  for(let i = 0; i < nstageP; i++) PRE.push(states[`P${i + 1}`]);
  for(let i = 0; i < nstageI; i++) INFD.push(states[`I${i + 1}`]);
  for(let i = 0; i < nstageH; i++) HOSP.push(states[`H${i + 1}`]);
  for(let i = 0; i < nstageC; i++) CARE.push(states[`C${i + 1}`]);
  for(let i = 0; i < nstageV; i++) VENT.push(states[`V${i + 1}`]);
  
  let EXPDQ = [], PREQ = [], INFDQ = [];
  for(let i = 0; i < nstageE; i++) EXPDQ.push(states[`EQ${i + 1}`]);
  for(let i = 0; i < nstageP; i++) PREQ.push(states[`PQ${i + 1}`]);
  for(let i = 0; i < nstageI; i++) INFDQ.push(states[`IQ${i + 1}`]);
   
  
  // Different transmission rates

  let TOT_PRE = 0;
  for(let i = 0 ; i < nstageP; i++) {
    TOT_PRE += PRE[i];
  }

  let TOT_INFD = 0;
  for(let i = 0; i < nstageI; i++) {
    TOT_INFD += INFD[i];
  }
  
  let PD;
  if ( !isNaN(Number(covar.tests) )) {
    PD = params.rho * (covar.tests / (covar.tests + params.TF));
  } else {
    PD = params.rho * (15.0 / (15.0 + params.TF)) ;
  }
  let lambdaI = params.betaI * TOT_INFD;
  let lambdaP1 = params.betaI * params.theta * TOT_PRE * (1-PD);
  let lambdaP2 = params.betaI * params.theta * TOT_PRE * PD;
  let dQdt, lambda, lambdaQ;
  if (t < t1) {
    dQdt  = mathLib.rgammawn(params.beta_sd, dt)/dt;
    lambda = ( (lambdaI + lambdaP1 + params.iota) / pop ) * dQdt;
    lambdaQ = ( (lambdaP2) / pop ) * dQdt;
  } else if (t < t1 + 7.0) {
    let x = (t-t1) / 7.0;
    let ss = 3*x*x - 2*x*x*x;
    dQdt  = mathLib.rgammawn(((1-ss) + ss*params.dB)*params.beta_sd, dt)/dt;
    lambda = ( ( ((1-ss) + ss*params.dI) * lambdaI + ((1-ss) + ss*params.dP) * lambdaP1 + ((1-ss) + ss*params.dT)*params.iota ) / pop ) * dQdt;
    lambdaQ = ( ( ((1-ss) + ss*params.dP) * lambdaP2  ) / pop ) * dQdt;
  } else {
    dQdt  = mathLib.rgammawn(params.dB*params.beta_sd, dt)/dt;
    lambda = ( ( params.dI * lambdaI + params.dP * lambdaP1 + params.dT * params.iota ) / pop ) * dQdt;
    lambdaQ = ( ( params.dI * lambdaI + params.dP * lambdaP1 + params.dT * params.iota ) / pop ) * dQdt;
  }
  
  // From class S
  let transS = new Array(2);
  let rateS = new Array(2);
  rateS[0] = lambda ;
  rateS[1] = lambdaQ;
  // reulermultinom(2,SUSC[0], &rateS[0], dt, &transS[0]);
  mathLib.reulermultinom(2, Math.round(SUSC[0]), 0, dt, 0, rateS, transS);
  
  // From class EQ
  let  transEQ = new Array(nstageE);
  let  rateEQ = [nstageE * params.sigma];
  for (let i = 0; i < nstageE; i++) {
  // reulermultinom(1, EXPDQ[i], &rateEQ, dt, &transEQ[i]);
    mathLib.reulermultinom(1, Math.round(EXPDQ[i]), 0, dt, i, rateEQ, transEQ);
  }
  
  // From class PQ
  let transPQ = new Array(nstageP + 1);
  let ratePQ = [nstageP * params.kappa];
  for (i = 0; i < nstageP - 1; i++) {
  // reulermultinom(1, PREQ[i], &ratePQ, dt, &transPQ[i]);
    mathLib.reulermultinom(1, Math.round(PREQ[i]), 0, dt, i, ratePQ, transPQ);
  }
  let ratePQIQH = new Array(2);
  ratePQIQH[0] = (1 - params.qP) * nstageP * params.kappa;
  ratePQIQH[1] = params.qP * nstageP * params.kappa;
  // reulermultinom(2, PREQ[nstageP-1], &ratePQIQH[0], dt, &transPQ[nstageP-1]);
  mathLib.reulermultinom(2, Math.round(PREQ[nstageP - 1]), 0, dt, nstageP - 1, ratePQIQH, transPQ);
  
  // From class IQ
  let transIQ = new Array(nstageI + 1);
  let rateIQ = [nstageI * params.gammaI];
  for (let i = 0; i < nstageI - 1; i++) {
  // reulermultinom(1, INFDQ[i], &rateIQ, dt, &transIQ[i]);
    mathLib.reulermultinom(1, Math.round(INFDQ[i]), 0, dt, i, rateIQ, transIQ);
  }
  
  let rateIQRD= new Array(2);
  rateIQRD[0] = (1 - params.qI) * nstageI * params.gammaI;
  rateIQRD[1] = params.qI * nstageI * params.gammaI;
  // reulermultinom(2, INFDQ[nstageI-1], &rateIQRD[0], dt, &transIQ[nstageI-1]);
  mathLib.reulermultinom(2, Math.round(INFDQ[nstageI-1]), 0, dt, nstageI - 1, rateIQRD, transIQ);

  
  // From class E
  let transE = new Array(nstageE);
  let rateE = [nstageE * params.sigma];
  for (let i = 0; i < nstageE; i++) {
  // reulermultinom(1, EXPD[i], &rateE, dt, &transE[i]);
    mathLib.reulermultinom(1, Math.round(EXPD[i]), 0, dt, i, rateE, transE);
  }
  
  // From class P
  let transP = new Array(nstageP + 2);
  let rateP = [nstageP * params.kappa];
  for (i = 0; i < nstageP - 1; i++) {
  // reulermultinom(1, PRE[i], &rateP, dt, &transP[i]);
    mathLib.reulermultinom(1, Math.round(PRE[i]), 0, dt, i, rateP, transP);  
  }
  
  let ratePIHIQ = new Array(3);
  ratePIHIQ[0] = (1 - PD) * (1 - params.qP) * nstageP * params.kappa;
  ratePIHIQ[1] = params.qP * nstageP * params.kappa;
  ratePIHIQ[2] = PD * (1 - params.qP) * nstageP * params.kappa;
  // reulermultinom(3, PREQ[nstageP-1], &ratePIHIQ[0], dt, &transP[nstageP-1]);
  mathLib.reulermultinom(3, Math.round(PREQ[nstageP - 1]), 0, dt, nstageP - 1, ratePIHIQ, transP);

  
  // From class I
  let transI = new Array(nstageI + 1);
  let rateI = [nstageI * params.gammaI];
  for (let i = 0; i < nstageI - 1; i++) {
  // reulermultinom(1, INFD[i], &rateI, dt, &transI[i]);
    mathLib.reulermultinom(1, Math.round(INFD[i]), 0, dt, i, rateI, transI);
  }
  
  let rateIRD = new Array(2);
  rateIRD[0] = (1 - params.qI) * nstageI * params.gammaI;
  rateIRD[1] = params.qI * nstageI * params.gammaI;
  // reulermultinom(2, INFD[nstageI-1], &rateIRD[0], dt, &transI[nstageI-1]);
  mathLib.reulermultinom(2, Math.round(INFD[nstageI - 1]), 0, dt, nstageI - 1, rateIRD, transI);
  
  // From class H
  let transH = new Array(nstageH + 1);
  let rateH = [nstageH * params.gammaH];
  for (let i = 0; i < nstageH - 1; i++) {
  // reulermultinom(1, HOSP[i], &rateH, dt, &transH[i]);
    mathLib.reulermultinom(1, Math.round(HOSP[i]), 0, dt, i, rateH, transH);
  }

  let rateHRC = new Array(2);
  rateHRC[0] = (1 - params.qH) * nstageH * params.gammaH;
  rateHRC[1] = params.qH * nstageH * params.gammaH;
  // reulermultinom(2, HOSP[nstageH-1], &rateHRC[0], dt, &transH[nstageH-1]);
  mathLib.reulermultinom(2, Math.round(HOSP[nstageH - 1]), 0, dt, nstageH - 1, rateHRC, transH);  
  
  // From class C
  let transC = new Array(nstageC + 1);
  let rateC = [nstageC * params.gammaC];
  for (let i = 0; i < nstageC - 1; i++) {
  // reulermultinom(1, CARE[i], &rateC, dt, &transC[i]);
    mathLib.reulermultinom(1, Math.round(CARE[i]), 0, dt, i, rateC, transC);
  }
  let rateCRV = new Array(2);
  rateCRV[0] = (1 - params.qC) * nstageC * params.gammaC;
  rateCRV[1] = params.qC * nstageC * params.gammaC;
  // reulermultinom(2, CARE[nstageC-1], &rateCRV[0], dt, &transC[nstageC-1]);
  mathLib.reulermultinom(2, Math.round(CARE[nstageC - 1]), 0, dt, nstageC - 1, rateCRV, transC);
  
  // From class V
  let transV = new Array(nstageV + 1);
  let rateV = [nstageV * params.gammaV];
  for (i = 0; i < nstageV - 1; i++) {
  // reulermultinom(1, VENT[i], &rateV, dt, &transV[i]);
    mathLib.reulermultinom(1, Math.round(VENT[i]), 0, dt, i, rateV, transV);
  }
  let rateVRD = new Array(2);
  rateVRD[0] = (1 - params.qV) * nstageV * params.gammaV;
  rateVRD[1] = params.qV * nstageV * params.gammaV;
  // reulermultinom(2, VENT[nstageV-1], &rateVRD[0], dt, &transV[nstageV-1]);
  mathLib.reulermultinom(2, Math.round(VENT[nstageV - 1]), 0, dt, nstageV - 1, rateVRD, transV);

  
  // Balance the equations
  SUSC[0] -= transS[0];
  EXPD[0] += transS[0];
  for (let i = 0; i < nstageE; i++) EXPD[i] -= transE[i];
  for (let i = 1; i < nstageE; i++) EXPD[i] += transE[i-1];
  PRE[0] += transE[nstageE-1];
  for (let i = 0; i < nstageP; i++) PRE[i] -= transP[i];
  for (let i = 1; i < nstageP; i++) PRE[i] += transP[i-1];
  INFD[0] += transP[nstageP-1];
  HOSP[0] += transP[nstageP];
  PRE[nstageP-1] -= transP[nstageP];
  for (let i = 0; i < nstageI; i++) INFD[i] -= transI[i];
  for (let i = 1; i < nstageI; i++) INFD[i] += transI[i-1];
  for (let i = 0; i < nstageH; i++) HOSP[i] -= transH[i];
  for (let i = 1; i < nstageH; i++) HOSP[i] += transH[i-1];
  CARE[0] += transH[nstageH];
  HOSP[nstageH-1] -= transH[nstageH];
  INFD[nstageI-1] -= transI[nstageI];
  for (let i = 0; i < nstageC; i++) CARE[i] -= transC[i];
  for (let i = 1; i < nstageC; i++) CARE[i] += transC[i-1];
  VENT[0] += transC[nstageC];
  CARE[nstageC-1] -= transC[nstageC];
  for (let i = 0; i < nstageV; i++) VENT[i] -= transV[i];
  for (let i = 1; i < nstageV; i++) VENT[i] += transV[i-1];
  VENT[nstageV-1] -= transV[nstageV];               
  RCVD[0] += transI[nstageI-1] + transH[nstageH-1] + transC[nstageC-1] + transV[nstageV-1];
  DEAD[0] += transI[nstageI] + transV[nstageV];
  
  SUSC[0] -= transS[1];
  EXPDQ[0] += transS[1];
  for (let i = 0; i < nstageE; i++) EXPDQ[i] -= transEQ[i];
  for (let i = 1; i < nstageE; i++) EXPDQ[i] += transEQ[i-1];
  PREQ[0] += transEQ[nstageE-1];
  for (let i = 0; i < nstageP; i++) PREQ[i] -= transPQ[i];
  for (let i = 1; i < nstageP; i++) PREQ[i] += transPQ[i-1];
  PREQ[nstageP-1] -= transPQ[nstageP];  
  PRE[nstageP-1] -= transP[nstageP+1];
  INFDQ[0] += transPQ[nstageP-1] + transP[nstageP+1];
  HOSP[0] += transPQ[nstageP];
  for (let i = 0; i < nstageI; i++) INFDQ[i] -= transIQ[i];
  for (let i = 1; i < nstageI; i++) INFDQ[i] += transIQ[i-1];
  INFDQ[nstageI-1] -= transIQ[nstageI];
  RCVD[0] += transIQ[nstageI-1];
  DEAD[0] += transIQ[nstageI];
  
  
  states.casesI += transP[nstageP-1];
  states.casesIQ += transPQ[nstageP-1] + transP[nstageP+1];
  states.casesH += transP[nstageP] + transPQ[nstageP];
  states.deathsI += transI[nstageI];
  states.deathsV += transV[nstageV];
  
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

snippet.initializer = function(args, covar) {
  let initObj = {};
  
  let fS = args.S0;
  let fEQ = args.EQ0 / nstageE;
  let fPQ = args.PQ0 / nstageP;
  let fIQ = args.IQ0 / nstageI;
  let fE = args.E0 / nstageE;
  let fP = args.P0 / nstageP;
  let fI = args.I0 / nstageI;
  let fH = args.H0 / nstageH;
  let fC = args.C0 / nstageC;
  let fV = args.V0 / nstageV;
  let fM = args.M0;
  let fR = 1 - fS - nstageE*fE - nstageP*fP - nstageI*fI - nstageH*fH - nstageC*fC - nstageV*fV - nstageE*fEQ - nstageP*fPQ - nstageI*fIQ - fM;
  
  initObj["S"] = Math.round(pop*fS);
  for (let i = 0; i < nstageE; i++) initObj[`EQ${i + 1}`] = Math.round(pop*fEQ);
  for (let i = 0; i < nstageP; i++) initObj[`PQ${i + 1}`] = Math.round(pop*fPQ);
  for (let i = 0; i < nstageI; i++) initObj[`IQ${i + 1}`] = Math.round(pop*fIQ);
  for (let i = 0; i < nstageE; i++) initObj[`E${i + 1}`] = Math.round(pop*fE);
  for (let i = 0; i < nstageP; i++) initObj[`P${i + 1}`] = Math.round(pop*fP);
  for (let i = 0; i < nstageI; i++) initObj[`I${i + 1}`] = Math.round(pop*fI);
  for (let i = 0; i < nstageH; i++) initObj[`H${i + 1}`] = Math.round(pop*fH);
  for (let i = 0; i < nstageC; i++) initObj[`C${i + 1}`] = Math.round(pop*fC);
  for (let i = 0; i < nstageV; i++) initObj[`V${i + 1}`] = Math.round(pop*fV);
  initObj["M"] = Math.round(pop*fM);
  initObj["R"] = Math.round(pop*fR);
  
  initObj["casesI"] = 0;
  initObj["casesIQ"] = 0;
  initObj["casesH"] = 0;
  initObj["deathsI"] = 0;
  initObj["deathsV"] = 0;
  
  return initObj;
}

snippet.dmeasure = function (data ,states, params, giveLog = 1) {
  let HOSP = [], CARE = [], VENT = [];
  for(let i = 0; i < nstageH; i++) HOSP.push(states[`H${i + 1}`]);
  for(let i = 0; i < nstageC; i++) CARE.push(states[`C${i + 1}`]);
  for(let i = 0; i < nstageV; i++) VENT.push(states[`V${i + 1}`]);
  
  let lik_reports, lik_deaths, lik_hospital, lik_ICU, lik_ventilator;
  
  if (isFinite(data.reports)) {
    let reported_cases = states.casesH + states.casesIQ;
    lik_reports = mathLib.dpois(data.reports, reported_cases + 1e-6, 1);
  } else {
    lik_reports = 0;
  }
  
  if (isFinite(data.deaths)) {
    let reported_deaths = states.deathsV + states.deathsI;
    lik_deaths = mathLib.dpois(data.deaths, reported_deaths + 1e-6, 1);
  } else {
    lik_deaths = 0;
  }
  
  let TOT_H = 0;
  for(let i = 0; i < nstageH; i++) {
  TOT_H += HOSP[i];
  }

  let TOT_C = 0;
  for(let i = 0; i < nstageC; i++) {
  TOT_C += CARE[i];
  }

  let TOT_V = 0;
  for(let i = 0; i < nstageV; i++) {
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
  for(let i = 0; i < nstageH; i++) HOSP.push(states[`H${i + 1}`]);
  for(let i = 0; i < nstageC; i++) HOSP.push(states[`C${i + 1}`]);
  for(let i = 0; i < nstageV; i++) HOSP.push(states[`V${i + 1}`]);
  
  let TOT_H, TOT_C, TOT_V;
  let reported_cases = states.casesH + states.casesIQ;
  results.reports = rpois(reported_cases  + 1e-6);
  
  let reported_deaths = states.deathsV + states.deathsI;
  results.deaths = rpois(reported_deaths  + 1e-6);
  
  let TOT_H = 0;
  for(let i = 0; i < nstageH; i++) {
  TOT_H += HOSP[i];
  }

  let TOT_C = 0;
  for(i = 0; i < nstageC; i++) {
  TOT_C += CARE[i];
  }

  let TOT_V = 0;         
  for(i = 0; i < nstageV; i++) {
  TOT_V += VENT[i];
  }

  results.ventilator = rpois(TOT_V + 1e-6);
  results.ICU = rpois(TOT_C + 1e-6) + results.ventilator;
  results.hospital = rpois(TOT_H + 1e-6) + results.ICU;

}
let params_log = ["betaI", "iota","beta_sd", "sigma", "kappa", "gammaI", "gammaH", "gammaC", "gammaV", "TF"];
let params_logit = ["rho", "theta", "dI", "dP", "dT", "dB","qP", "qI",  "qH", "qC", "qV"];
let statenamesFn = function() {
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

snippet.determineRW = function(run) {
  if(run === 1) {
    return ((time) => {
      let rw_size = 0.05;
      let t1 = 75;
      let betaI = time < t1 ? rw_size : 0
      let theta = time < t1 ? rw_size : 0
      let iota = time < t1 ? rw_size : 0
      let beta_sd = time < t1 ? rw_size : 0
      let dI = time < t1 ? 0 : rw_size
      let dP = time < t1 ? 0 : rw_size
      let dT = time < t1 ? 0 : rw_size
      let dB = time < t1 ? 0 : rw_size
      let qP = time < t1 ? rw_size : rw_size
      let qI = time < t1 ? rw_size : rw_size
      let qH = time < t1 ? rw_size : rw_size
      let qC = time < t1 ? rw_size : rw_size
      let qV = time < t1 ? rw_size : rw_size
      let sigma = time < t1 ? 0.5*rw_size : 0.5*rw_size
      let kappa = time < t1 ? 0.5*rw_size : 0.5*rw_size
      let gammaI = time < t1 ? 0.5*rw_size : 0.5*rw_size
      let gammaH = time < 93 ? 0 : rw_size
      let gammaC = time < 93 ? 0 : rw_size
      let gammaV = time < 93 ? 0 : rw_size
      let rho = time < t1 ? rw_size : rw_size
      let TF = time < 35 ? 0 : rw_size;
      return {
        betaI : betaI, theta : theta, iota : iota, beta_sd : beta_sd, dI : dI, dP : dP,  dT : dT,  dB : dB, qP : qP, qI : qI,  qH : qH,  qC : qC, qV : qV, sigma : sigma, kappa : kappa, gammaI : gammaI,  gammaH : gammaH, gammaC : gammaC, gammaV : gammaV,  rho : rho, TF: TF};
    }).toString()
  } else if (run === 2) {
    return ((time) => {
      let rw_size = 0.05;
      let t1 = 75;
      let betaI = time < t1? 0 : 0;
      let theta = time < t1? rw_size : 0;
      let iota = time < t1? rw_size : 0;
      let beta_sd = time < t1? rw_size : 0;
      let dI = time < t1? 0 : rw_size;
      let dP = time < t1? 0 : rw_size;
      let dT = time < t1? 0 : rw_size;
      let dB = time < t1? 0 : rw_size;
      let qP = time < t1? rw_size : rw_size;
      let qI = time < t1? rw_size : rw_size;
      let qH = time < t1? rw_size : rw_size;
      let qC = time < t1? rw_size : rw_size;
      let qV = time < t1? rw_size : rw_size;
      let sigma = time < t1? 0.5*rw_size : 0.5*rw_size;
      let kappa = time < t1? 0.5*rw_size : 0.5*rw_size;
      let gammaI = time < t1? 0.5*rw_size : 0.5*rw_size;
      let gammaH = time < 93? 0 : rw_size;
      let gammaC = time < 93? 0 : rw_size;
      let gammaV = time < 93? 0 : rw_size;
      let rho = time < t1? rw_size : rw_size;
      let TF = time < 35? 0 : rw_size;
      return {
        betaI : betaI, theta : theta, iota : iota, beta_sd : beta_sd, dI : dI, dP : dP,  dT : dT,  dB : dB, qP : qP, qI : qI,  qH : qH,  qC : qC, qV : qV, sigma : sigma, kappa : kappa, gammaI : gammaI,  gammaH : gammaH, gammaC : gammaC, gammaV : gammaV,  rho : rho, TF: TF};
    }).toString()
  } else if (run === 3) {
    return ((time) => {
      let rw_size = 0.05;
      let t1 = 75;
      let betaI = time < t1 ? rw_size : 0;
      let theta = time < t1 ? 0 : 0;
      let iota = time < t1 ? rw_size : 0;
      let beta_sd = time < t1 ? rw_size : 0;
      let dI = time < t1 ? 0 : rw_size;
      let dP = time < t1 ? 0 : rw_size;
      let dT = time < t1 ? 0 : rw_size;
      let dB = time < t1 ? 0 : rw_size;
      let qP = time < t1 ? rw_size : rw_size;
      let qI = time < t1 ? rw_size : rw_size;
      let qH = time < t1 ? rw_size : rw_size;
      let qC = time < t1 ? rw_size : rw_size;
      let qV = time < t1 ? rw_size : rw_size;
      let sigma = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let kappa = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let gammaI = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let gammaH = time < 93 ? 0 : rw_size;
      let gammaC = time < 93 ? 0 : rw_size;
      let gammaV = time < 93 ? 0 : rw_size;
      let rho = time < t1 ? rw_size : rw_size;
      let TF = time < 35 ? 0 : rw_size;
      return {
        betaI : betaI, theta : theta, iota : iota, beta_sd : beta_sd, dI : dI, dP : dP,  dT : dT,  dB : dB, qP : qP, qI : qI,  qH : qH,  qC : qC, qV : qV, sigma : sigma, kappa : kappa, gammaI : gammaI,  gammaH : gammaH, gammaC : gammaC, gammaV : gammaV,  rho : rho, TF: TF};
    }).toString()
  } else if (run === 4) {
    return ((time) => {
      let rw_size = 0.05;
      let t1 = 75;
      let betaI = time < t1 ? rw_size : 0;
      let theta = time < t1 ? rw_sd : 0;
      let iota = time < t1 ? 0 : 0;
      let beta_sd = time < t1 ? rw_size : 0;
      let dI = time < t1 ? 0 : rw_size;
      let dP = time < t1 ? 0 : rw_size;
      let dT = time < t1 ? 0 : rw_size;
      let dB = time < t1 ? 0 : rw_size;
      let qP = time < t1 ? rw_size : rw_size;
      let qI = time < t1 ? rw_size : rw_size;
      let qH = time < t1 ? rw_size : rw_size;
      let qC = time < t1 ? rw_size : rw_size;
      let qV = time < t1 ? rw_size : rw_size;
      let sigma = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let kappa = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let gammaI = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let gammaH = time < 93 ? 0 : rw_size;
      let gammaC = time < 93 ? 0 : rw_size;
      let gammaV = time < 93 ? 0 : rw_size;
      let rho = time < t1 ? rw_size : rw_size;
      let TF = time < 35 ? 0 : rw_size;
      return {
        betaI : betaI, theta : theta, iota : iota, beta_sd : beta_sd, dI : dI, dP : dP,  dT : dT,  dB : dB, qP : qP, qI : qI,  qH : qH,  qC : qC, qV : qV, sigma : sigma, kappa : kappa, gammaI : gammaI,  gammaH : gammaH, gammaC : gammaC, gammaV : gammaV,  rho : rho, TF: TF};
    }).toString()
  } else if (run === 5) {
    return ((time) => {
      let rw_size = 0.05;
      let t1 = 75;
      let betaI = time < t1 ? rw_size : 0;
      let theta = time < t1 ? rw_size : 0;
      let iota = time < t1 ? rw_size : 0;
      let beta_sd = time < t1 ? 0 : 0;
      let dI = time < t1 ? 0 : rw_size;
      let dP = time < t1 ? 0 : rw_size;
      let dT = time < t1 ? 0 : rw_size;
      let dB = time < t1 ? 0 : rw_size;
      let qP = time < t1 ? rw_size : rw_size;
      let qI = time < t1 ? rw_size : rw_size;
      let qH = time < t1 ? rw_size : rw_size;
      let qC = time < t1 ? rw_size : rw_size;
      let qV = time < t1 ? rw_size : rw_size;
      let sigma = time < t1 ? 0.5*rw_si : 0.5*rw_size;
      let kappa = time < t1 ? 0.5*rw_si : 0.5*rw_size;
      let gammaI = time < t1 ? 0.5*rw_si : 0.5*rw_size;
      let gammaH = time < 93 ? 0 : rw_size;
      let gammaC = time < 93 ? 0 : rw_size;
      let gammaV = time < 93 ? 0 : rw_size;
      let rho = time < t1 ? rw_size : rw_size;
      let TF = time < 35 ? 0 : rw_size;
      return {
        betaI : betaI, theta : theta, iota : iota, beta_sd : beta_sd, dI : dI, dP : dP,  dT : dT,  dB : dB, qP : qP, qI : qI,  qH : qH,  qC : qC, qV : qV, sigma : sigma, kappa : kappa, gammaI : gammaI,  gammaH : gammaH, gammaC : gammaC, gammaV : gammaV,  rho : rho, TF: TF};
    }).toString()
  } else if (run === 6) {
    return ((time) => {
      let rw_size = 0.05;
      let t1 = 75;
      let betaI = time < t1 ? rw_size : 0;
      let theta = time < t1 ? rw_size : 0;
      let iota = time < t1 ? rw_size : 0;
      let beta_sd = time < t1 ? rw_size : 0;
      let dI = time < t1 ? 0 : 0;
      let dP = time < t1 ? 0 : rw_size;
      let dT = time < t1 ? 0 : rw_size;
      let dB = time < t1 ? 0 : rw_size;
      let qP = time < t1 ? rw_size : rw_size;
      let qI = time < t1 ? rw_size : rw_size;
      let qH = time < t1 ? rw_size : rw_size;
      let qC = time < t1 ? rw_size : rw_size;
      let qV = time < t1 ? rw_size : rw_size;
      let sigma = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let kappa = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let gammaI = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let gammaH = time < 93 ? 0 : rw_size;
      let gammaC = time < 93 ? 0 : rw_size;
      let gammaV = time < 93 ? 0 : rw_size;
      let rho = time < t1 ? rw_size : rw_size;
      let TF = time < 35 ? 0 : rw_size;
      return {
        betaI : betaI, theta : theta, iota : iota, beta_sd : beta_sd, dI : dI, dP : dP,  dT : dT,  dB : dB, qP : qP, qI : qI,  qH : qH,  qC : qC, qV : qV, sigma : sigma, kappa : kappa, gammaI : gammaI,  gammaH : gammaH, gammaC : gammaC, gammaV : gammaV,  rho : rho, TF: TF};
    }).toString()
  } else if (run === 22) {
    return ((time) => {
      let rw_size = 0.05;
      let t1 = 75;
      let betaI = time < t1 ? rw_size : 0;
      let theta = time < t1 ? rw_size : 0;
      let iota = time < t1 ? rw_size : 0;
      let beta_sd = time < t1 ? rw_size : 0;
      let dI = time < t1 ? 0 : rw_size;
      let dP = time < t1 ? 0 : rw_size;
      let dT = time < t1 ? 0 : rw_size;
      let dB = time < t1 ? 0 : rw_size;
      let qP = time < t1 ? rw_size : rw_size;
      let qI = time < t1 ? rw_size : rw_size;
      let qH = time < t1 ? rw_size : rw_size;
      let qC = time < t1 ? rw_size : rw_size;
      let qV = time < t1 ? rw_size : rw_size;
      let sigma = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let kappa = time < t1 ? 0.5*rw_size : 0.5*rw_size;
      let gammaI = time < t1 ? 0 : 0;
      let gammaH = time < 93 ? 0 : rw_size;
      let gammaC = time < 93 ? 0 : rw_size;
      let gammaV = time < 93 ? 0 : rw_size;
      let rho = time < t1 ? rw_size : rw_size;
      let TF = time < 35 ? 0 : rw_size;
      return {
        betaI : betaI, theta : theta, iota : iota, beta_sd : beta_sd, dI : dI, dP : dP,  dT : dT,  dB : dB, qP : qP, qI : qI,  qH : qH,  qC : qC, qV : qV, sigma : sigma, kappa : kappa, gammaI : gammaI,  gammaH : gammaH, gammaC : gammaC, gammaV : gammaV,  rho : rho, TF: TF};
    }).toString()
  }
}


module.exports = snippet
