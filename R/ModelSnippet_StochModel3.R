library(pomp)
endTime <- "2020-07-16"
T0 <- 75
T1 <- 139
modeltype <- c(nstageE=3L, nstageI=3L, nstageP=3L, nstageH=3L, nstageC=3L, nstageV=3L)
pop <- 10e6
t1 <- 75
dObs <- Csnippet("
                 double *HOSP = &H1;
                 double *CARE = &C1;
                 double *VENT = &V1;
                 double lik_reports, lik_deaths, lik_hospital, lik_ICU, lik_ventilator;
                 double TOT_H, TOT_C, TOT_V;
                 int i;
                 if (R_FINITE(reports)) {
                 double reported_cases = casesH + casesIQ;
                 lik_reports = dpois(reports, reported_cases + 1e-6, 1);
                 } else {
                 lik_reports = 0;
                 }
                 if (R_FINITE(deaths)) {
                 double reported_deaths = deathsCV + deathsIIQ;
                 lik_deaths = dpois(deaths, reported_deaths + 1e-6, 1);
                 } else {
                 lik_deaths = 0;
                 }
                 for(i = 0, TOT_H = 0; i < nstageH; i++) {
                 TOT_H += HOSP[i];
                 }
                 for(i = 0, TOT_C = 0; i < nstageC; i++) {
                 TOT_C += CARE[i];
                 }
                 for(i = 0, TOT_V = 0; i < nstageV; i++) {
                 TOT_V += VENT[i];
                 }
                 if (R_FINITE(hospital) & R_FINITE(ICU)) {
                 lik_hospital = dpois(hospital - ICU, TOT_H + 1e-6, 1);
                 } else {
                 lik_hospital = 0;
                 }
                 if (R_FINITE(ICU) & R_FINITE(ventilator)) {
                 lik_ICU = dpois(ICU - ventilator, TOT_C + 1e-6, 1);
                 } else {
                 lik_ICU = 0;
                 }
                 if (R_FINITE(ventilator)) {
                 lik_ventilator = dpois(ventilator, TOT_V + 1e-6, 1);
                 } else {
                 lik_ventilator = 0;
                 }
                 lik = lik_reports + lik_deaths + lik_hospital + lik_ICU + lik_ventilator;
                 if (give_log == 0 ) {
                 lik = exp(lik);
                 }
                 ")

rObs <- Csnippet("
                 double *HOSP = &H1;
                 double *CARE = &C1;
                 double *VENT = &V1;
                 double TOT_H, TOT_C, TOT_V;
                 double reported_cases = casesH + casesIQ;
                 reports = rpois(reported_cases  + 1e-6);

                 double reported_deaths = deathsCV + deathsIIQ;
                 deaths = rpois(reported_deaths  + 1e-6);

                 int i;
                 for(i = 0, TOT_H = 0; i < nstageH; i++) {
                 TOT_H += HOSP[i];
                 }
                 for(i = 0, TOT_C = 0; i < nstageC; i++) {
                 TOT_C += CARE[i];
                 }
                 for(i = 0, TOT_V = 0; i < nstageV; i++) {
                 TOT_V += VENT[i];
                 }
                 ventilator = rpois(TOT_V + 1e-6);
                 ICU = rpois(TOT_C + 1e-6) + ventilator;
                 hospital = rpois(TOT_H + 1e-6) + ICU;
                 ")

rSim <- Csnippet("

                 double *SUSC = &S;
                 double *EXPD = &E1;
                 double *PRE = &P1;
                 double *INFD = &I1;
                 double *HOSP = &H1;
                 double *CARE = &C1;
                 double *VENT = &V1;
                 double *DEAD = &M;
                 double *RCVD = &R;

                 double *EXPDQ = &EQ1;
                 double *PREQ = &PQ1;
                 double *INFDQ = &IQ1;


                 // Different transmission rates

                 double dQdt;
                 double TOT_PRE, TOT_INFD;
                 int i;

                 for(i = 0, TOT_PRE = 0; i < nstageP; i++) {
                 TOT_PRE += PRE[i];
                 }
                 for(i = 0, TOT_INFD = 0; i < nstageI; i++) {
                 TOT_INFD += INFD[i];
                 }

                 double PD, lambdaI, lambdaP, lambdaPQ, lambda, lambdaQ;
                 if ( R_FINITE(tests) ) {
                 PD = rho * (tests / (tests + TF));
                 } else {
                 PD = rho * (15.0 / (15.0 + TF)) ;
                 }
                 lambdaI = betaI * TOT_INFD;
                 lambdaP = betaI * theta * TOT_PRE * (1-PD);
                 lambdaPQ = betaI * theta * TOT_PRE * PD;


                 if (t < T0) {
                 dQdt  = rgammawn(beta_sd, dt)/dt;
                 lambda = ( (lambdaI + lambdaP + iota) / pop ) * dQdt;
                 lambdaQ = ( (lambdaPQ) / pop ) * dQdt;
                 } else if (t < T0 + 7.0) {
                 double x = (t-T0) / 7.0;
                 double ss = 3*x*x - 2*x*x*x;
                 dQdt  = rgammawn(((1-ss) + ss*dB0)*beta_sd, dt)/dt;
                 lambda = ( ( ((1-ss) + ss*dI0) * lambdaI + ((1-ss) + ss*dP0) * lambdaP + ((1-ss) + ss*dT0)*iota ) / pop ) * dQdt;
                 lambdaQ = ( ( ((1-ss) + ss*dP0) * lambdaPQ  ) / pop ) * dQdt;
                 } else if (t < T1) {
                 dQdt  = rgammawn(dB0*beta_sd, dt)/dt;
                 lambda = ( ( dI0 * lambdaI + dP0 * lambdaP + dT0 * iota ) / pop ) * dQdt;
                 lambdaQ = ( (  dP0 * lambdaPQ  ) / pop ) * dQdt;
                 } else if (t < T1 +7.0) {
                 double x = (t-T1) / 7.0;
                 double ss = 3*x*x - 2*x*x*x;
                 dQdt  = rgammawn(((1-ss)*dB0 + ss*dB1)*beta_sd, dt)/dt;
                 lambda = ( ( ((1-ss)*dI0 + ss*dI1) * lambdaI + ((1-ss)*dP0 + ss*dP1) * lambdaP + ((1-ss)*dT0 + ss*dT1)*iota ) / pop ) * dQdt;
                 lambdaQ = ( ( ((1-ss)*dP0 + ss*dP1) * lambdaPQ  ) / pop ) * dQdt;
                 } else {
                 dQdt  = rgammawn(dB1*beta_sd, dt)/dt;
                 lambda = ( ( dI1 * lambdaI + dP1 * lambdaP + dT1 * iota ) / pop ) * dQdt;
                 lambdaQ = ( (  dP1 * lambdaPQ  ) / pop ) * dQdt;
                 }

                 // From class S
                 double transS[2];
                 double rateS[2];
                 rateS[0] = lambda ;
                 rateS[1] = lambdaQ;
                 reulermultinom(2,SUSC[0], &rateS[0], dt, &transS[0]);


                 // From class EQ
                 double transEQ[nstageE];
                 double rateEQ = nstageE * sigma;
                 for (i = 0; i < nstageE; i++) {
                 reulermultinom(1, EXPDQ[i], &rateEQ, dt, &transEQ[i]);
                 }

                 // From class PQ
                 double transPQ[nstageP+1];
                 double ratePQ = nstageP * kappa;
                 for (i = 0; i < nstageP-1; i++) {
                 reulermultinom(1, PREQ[i], &ratePQ, dt, &transPQ[i]);
                 }
                 double ratePQIQH[2];
                 ratePQIQH[0] = (1-qP) * nstageP * kappa;
                 ratePQIQH[1] = qP * nstageP * kappa;
                 reulermultinom(2, PREQ[nstageP-1], &ratePQIQH[0], dt, &transPQ[nstageP-1]);

                 // From class IQ
                 double transIQ[nstageI+1];
                 double rateIQ = nstageI * gammaI;
                 for (i = 0; i < nstageI-1; i++) {
                 reulermultinom(1, INFDQ[i], &rateIQ, dt, &transIQ[i]);
                 }

                 double rateIQRD[2];
                 rateIQRD[0] = (1-mI) * nstageI * gammaI;
                 rateIQRD[1] = mI * nstageI * gammaI;
                 reulermultinom(2, INFDQ[nstageI-1], &rateIQRD[0], dt, &transIQ[nstageI-1]);

                 // From class E
                 double transE[nstageE];
                 double rateE = nstageE * sigma;
                 for (i = 0; i < nstageE; i++) {
                 reulermultinom(1, EXPD[i], &rateE, dt, &transE[i]);
                 }

                 // From class P
                 double transP[nstageP+2];
                 double rateP = nstageP * kappa;
                 for (i = 0; i < nstageP-1; i++) {
                 reulermultinom(1, PRE[i], &rateP, dt, &transP[i]);
                 }

                 double ratePIHIQ[3];
                 ratePIHIQ[0] = (1.0-PD) * (1.0-qP) * nstageP * kappa;
                 ratePIHIQ[1] = qP * nstageP * kappa;
                 ratePIHIQ[2] = PD * (1.0-qP) * nstageP * kappa;
                 reulermultinom(3, PRE[nstageP-1], &ratePIHIQ[0], dt, &transP[nstageP-1]);

                 // From class I
                 double transI[nstageI+1];
                 double rateI = nstageI * gammaI;
                 for (i = 0; i < nstageI-1; i++) {
                 reulermultinom(1, INFD[i], &rateI, dt, &transI[i]);
                 }

                 double rateIRD[2];
                 rateIRD[0] = (1-mI) * nstageI * gammaI;
                 rateIRD[1] = mI * nstageI * gammaI;
                 reulermultinom(2, INFD[nstageI-1], &rateIRD[0], dt, &transI[nstageI-1]);


                 // From class H
                 double transH[nstageH+1];
                 double rateH = nstageH * gammaH;
                 for (i = 0; i < nstageH-1; i++) {
                 reulermultinom(1, HOSP[i], &rateH, dt, &transH[i]);
                 }
                 double rateHRC[2];
                 rateHRC[0] = (1-qH) * nstageH * gammaH;
                 rateHRC[1] = qH * nstageH * gammaH;
                 reulermultinom(2, HOSP[nstageH-1], &rateHRC[0], dt, &transH[nstageH-1]);


                 // From class C
                 double transC[nstageC+2];
                 double rateC = nstageC * gammaC;
                 for (i = 0; i < nstageC-1; i++) {
                 reulermultinom(1, CARE[i], &rateC, dt, &transC[i]);
                 }
                 double rateCRVM[3];
                 rateCRVM[0] = (1-mC) * (1-qC) * nstageC * gammaC;
                 rateCRVM[1] = qC * nstageC * gammaC;
                 rateCRVM[2] = mC * (1-qC) * nstageC * gammaC;
                 reulermultinom(3, CARE[nstageC-1], &rateCRVM[0], dt, &transC[nstageC-1]);

                 // From class V
                 double transV[nstageV+1];
                 double rateV = nstageV * gammaV;
                 for (i = 0; i < nstageV-1; i++) {
                 reulermultinom(1, VENT[i], &rateV, dt, &transV[i]);
                 }
                 double rateVRD[2];
                 rateVRD[0] = (1-mV) * nstageV * gammaV;
                 rateVRD[1] = mV * nstageV * gammaV;
                 reulermultinom(2, VENT[nstageV-1], &rateVRD[0], dt, &transV[nstageV-1]);



                 // Balance the equations
                 SUSC[0] -= transS[0];
                 EXPD[0] += transS[0];
                 for (i = 0; i < nstageE; i++) EXPD[i] -= transE[i];
                 for (i = 1; i < nstageE; i++) EXPD[i] += transE[i-1];
                 PRE[0] += transE[nstageE-1];
                 for (i = 0; i < nstageP; i++) PRE[i] -= transP[i];
                 for (i = 1; i < nstageP; i++) PRE[i] += transP[i-1];
                 INFD[0] += transP[nstageP-1];
                 HOSP[0] += transP[nstageP];
                 PRE[nstageP-1] -= transP[nstageP];
                 for (i = 0; i < nstageI; i++) INFD[i] -= transI[i];
                 for (i = 1; i < nstageI; i++) INFD[i] += transI[i-1];
                 for (i = 0; i < nstageH; i++) HOSP[i] -= transH[i];
                 for (i = 1; i < nstageH; i++) HOSP[i] += transH[i-1];
                 CARE[0] += transH[nstageH];
                 HOSP[nstageH-1] -= transH[nstageH];
                 INFD[nstageI-1] -= transI[nstageI];
                 for (i = 0; i < nstageC; i++) CARE[i] -= transC[i];
                 for (i = 1; i < nstageC; i++) CARE[i] += transC[i-1];
                 VENT[0] += transC[nstageC];
                 CARE[nstageC-1] -= transC[nstageC] + transC[nstageC+1];
                 for (i = 0; i < nstageV; i++) VENT[i] -= transV[i];
                 for (i = 1; i < nstageV; i++) VENT[i] += transV[i-1];
                 VENT[nstageV-1] -= transV[nstageV];
                 RCVD[0] += transI[nstageI-1] + transH[nstageH-1] + transC[nstageC-1] + transV[nstageV-1];
                 DEAD[0] += transI[nstageI] + transC[nstageC+1] + transV[nstageV];

                 SUSC[0] -= transS[1];
                 EXPDQ[0] += transS[1];
                 for (i = 0; i < nstageE; i++) EXPDQ[i] -= transEQ[i];
                 for (i = 1; i < nstageE; i++) EXPDQ[i] += transEQ[i-1];
                 PREQ[0] += transEQ[nstageE-1];
                 for (i = 0; i < nstageP; i++) PREQ[i] -= transPQ[i];
                 for (i = 1; i < nstageP; i++) PREQ[i] += transPQ[i-1];
                 PREQ[nstageP-1] -= transPQ[nstageP];
                 PRE[nstageP-1] -= transP[nstageP+1];
                 INFDQ[0] += transPQ[nstageP-1] + transP[nstageP+1];
                 HOSP[0] += transPQ[nstageP];
                 for (i = 0; i < nstageI; i++) INFDQ[i] -= transIQ[i];
                 for (i = 1; i < nstageI; i++) INFDQ[i] += transIQ[i-1];
                 INFDQ[nstageI-1] -= transIQ[nstageI];
                 RCVD[0] += transIQ[nstageI-1];
                 DEAD[0] += transIQ[nstageI];

                 casesI += transP[nstageP-1];
                 casesIQ += transPQ[nstageP-1] + transP[nstageP+1];
                 casesH += transPQ[nstageP] + transP[nstageP];
                 deathsIIQ += transI[nstageI] + transIQ[nstageI];
                 deathsCV += transC[nstageC+1] + transV[nstageV];
                 ")


skel <- Csnippet("
                 double dt = 0.1;
                 double *SUSC = &S;
                 double *EXPD = &E1;
                 double *PRE = &P1;
                 double *INFD = &I1;
                 double *HOSP = &H1;
                 double *CARE = &C1;
                 double *VENT = &V1;
                 double *DEAD = &M;
                 double *RCVD = &R;

                 double *DSUSC = &DS;
                 double *DEXPD = &DE1;
                 double *DPRE = &DP1;
                 double *DINFD = &DI1;
                 double *DHOSP = &DH1;
                 double *DCARE = &DC1;
                 double *DVENT = &DV1;
                 double *DDEAD = &DM;
                 double *DRCVD = &DR;


                 double *EXPDQ = &EQ1;
                 double *PREQ = &PQ1;
                 double *INFDQ = &IQ1;
                 double *DEXPDQ = &DEQ1;
                 double *DPREQ = &DPQ1;
                 double *DINFDQ = &DIQ1;


                 // Different transmission rates

                 double dQdt;
                 double TOT_INFD, TOT_PRE;
                 int i;

                 for(i = 0, TOT_PRE = 0; i < nstageP; i++) {
                 TOT_PRE += PRE[i];
                 }
                 for(i = 0, TOT_INFD = 0; i < nstageI; i++) {
                 TOT_INFD += INFD[i];
                 }


                 double PD, lambdaI, lambdaP, lambdaPQ, lambda, lambdaQ;
                 if ( R_FINITE(tests) ) {
                 PD = rho * (tests / (tests + TF));
                 } else {
                 PD = rho * (15.0 / (15.0 + TF)) ;
                 }
                 lambdaI = betaI * TOT_INFD;
                 lambdaP = betaI * theta * TOT_PRE * (1-PD);
                 lambdaPQ = betaI * theta * TOT_PRE * PD;

                 if (t < T0) {
                 dQdt  = rgammawn(beta_sd, dt)/dt;
                 lambda = ( (lambdaI + lambdaP + iota) / pop ) * dQdt;
                 lambdaQ = ( (lambdaPQ) / pop ) * dQdt;
                 } else if (t < T0 + 7.0) {
                 double x = (t-T0) / 7.0;
                 double ss = 3*x*x - 2*x*x*x;
                 dQdt  = rgammawn(((1-ss) + ss*dB0)*beta_sd, dt)/dt;
                 lambda = ( ( ((1-ss) + ss*dI0) * lambdaI + ((1-ss) + ss*dP0) * lambdaP + ((1-ss) + ss*dT0)*iota ) / pop ) * dQdt;
                 lambdaQ = ( ( ((1-ss) + ss*dP0) * lambdaPQ  ) / pop ) * dQdt;
                 } else if (t < T1) {
                 dQdt  = rgammawn(dB0*beta_sd, dt)/dt;
                 lambda = ( ( dI0 * lambdaI + dP0 * lambdaP + dT0 * iota ) / pop ) * dQdt;
                 lambdaQ = ( (  dP0 * lambdaPQ  ) / pop ) * dQdt;
                 } else if (t < T1 +7.0) {
                 double x = (t-T1) / 7.0;
                 double ss = 3*x*x - 2*x*x*x;
                 dQdt  = rgammawn(((1-ss)*dB0 + ss*dB1)*beta_sd, dt)/dt;
                 lambda = ( ( ((1-ss)*dI0 + ss*dI1) * lambdaI + ((1-ss)*dP0 + ss*dP1) * lambdaP + ((1-ss)*dT0 + ss*dT1)*iota ) / pop ) * dQdt;
                 lambdaQ = ( ( ((1-ss)*dP0 + ss*dP1) * lambdaPQ  ) / pop ) * dQdt;
                 } else {
                 dQdt  = rgammawn(dB1*beta_sd, dt)/dt;
                 lambda = ( ( dI1 * lambdaI + dP1 * lambdaP + dT1 * iota ) / pop ) * dQdt;
                 lambdaQ = ( (  dP1 * lambdaPQ  ) / pop ) * dQdt;
                 }

                 // From class S
                 double transS[2];
                 double rateS[2];
                 rateS[0] = lambda ;
                 rateS[1] = lambdaQ;
                 double rateS_tot = rateS[0] + rateS[1];
                 double transS_tot = (1.0 - exp(- rateS_tot * dt)) * SUSC[0];
                 transS[0] = rateS[0] / rateS_tot * transS_tot;
                 transS[1] = rateS[1] / rateS_tot * transS_tot;

                 // From class EQ
                 double transEQ[nstageE];
                 double rateEQ = nstageE * sigma;
                 for (i = 0; i < nstageE; i++) {
                 transEQ[i] =  (1.0 - exp(- rateEQ * dt))* EXPDQ[i];
                 }

                 // From class PQ
                 double transPQ[nstageP+1];
                 double ratePQ = nstageP * kappa;
                 for (i = 0; i < nstageP-1; i++) {
                 transPQ[i] =  (1.0 - exp(- ratePQ * dt))* PREQ[i];
                 }
                 double transPQIQH;
                 transPQIQH = (1.0 - exp(- ratePQ * dt))* PREQ[nstageP-1];
                 transPQ[nstageP-1] = (1-qP) * transPQIQH;
                 transPQ[nstageP] = qP * transPQIQH;

                 // From class IQ
                 double transIQ[nstageI+1];
                 double rateIQ = nstageI * gammaI;
                 for (i = 0; i < nstageI-1; i++) {
                 transIQ[i] =  (1.0 - exp(- rateIQ * dt))* INFDQ[i];
                 }
                 double transIQRD;
                 transIQRD = (1.0 - exp(- rateIQ * dt))* INFDQ[nstageI-1];
                 transIQ[nstageI-1] = (1-mI) * transIQRD;
                 transIQ[nstageI] = mI * transIQRD;


                 // From class E
                 double transE[nstageE];
                 double rateE = nstageE * sigma;
                 for (i = 0; i < nstageE; i++) {
                 transE[i] =  (1.0 - exp(- rateE * dt))* EXPD[i];
                 }

                 // From class P
                 double transP[nstageP+2];
                 double rateP = nstageP * kappa;
                 for (i = 0; i < nstageP-1; i++) {
                 transP[i] =  (1.0 - exp(- rateP * dt))* PRE[i];
                 }
                 double transPIHIQ;
                 transPIHIQ = (1.0 - exp(- rateP * dt))* PRE[nstageP-1];
                 transP[nstageP-1] = (1-PD) * (1-qP) * transPIHIQ;
                 transP[nstageP] = qP * transPIHIQ;
                 transP[nstageP+1] = PD * (1-qP) * transPIHIQ;

                 // From class I
                 double transI[nstageI+1];
                 double rateI = nstageI * gammaI;
                 for (i = 0; i < nstageI-1; i++) {
                 transI[i] =  (1.0 - exp(- rateI * dt))* INFD[i];
                 }

                 double transIRD;
                 transIRD = (1.0 - exp(- rateI * dt))* INFD[nstageI-1];
                 transI[nstageI-1] = (1-mI) * transIRD;
                 transI[nstageI] = mI * transIRD;


                 // From class H
                 double transH[nstageH+1];
                 double rateH = nstageH * gammaH;
                 for (i = 0; i < nstageH-1; i++) {
                 transH[i] =  (1.0 - exp(- rateH * dt))* HOSP[i];
                 }

                 double transHRC;
                 transHRC = (1.0 - exp(- rateH * dt))* HOSP[nstageH-1];
                 transH[nstageH-1] = (1-qH) * transHRC;
                 transH[nstageH] = qH * transHRC;


                 // From class C
                 double transC[nstageC+2];
                 double rateC = nstageC * gammaC;
                 for (i = 0; i < nstageC-1; i++) {
                 transC[i] =  (1.0 - exp(- rateC * dt))* CARE[i];
                 }
                 double transCRVM;
                 transCRVM = (1.0 - exp(- rateC * dt))* CARE[nstageC-1];
                 transC[nstageC-1] =  (1-mC) * (1-qC) * transCRVM;
                 transC[nstageC] = qC * transCRVM;
                 transC[nstageC+1] = mC * (1-qC) * transCRVM;

                 // From class V
                 double transV[nstageV+1];
                 double rateV = nstageV * gammaV;
                 for (i = 0; i < nstageV-1; i++) {
                 transV[i] =  (1.0 - exp(- rateV * dt))* VENT[i];
                 }
                 double transVRD;
                 transVRD = (1.0 - exp(- rateV * dt))* VENT[nstageV-1];
                 transV[nstageV-1] = (1-mV) * transVRD;
                 transV[nstageV] = mV * transVRD;

                 // Balance the equations
                 DSUSC[0] = SUSC[0];
                 for (i = 0; i < nstageE; i++) DEXPD[i] = EXPD[i];
                 for (i = 0; i < nstageP; i++) DPRE[i] = PRE[i];
                 for (i = 0; i < nstageI; i++) DINFD[i] = INFD[i];
                 for (i = 0; i < nstageH; i++) DHOSP[i] = HOSP[i];
                 for (i = 0; i < nstageC; i++) DCARE[i] = CARE[i];
                 for (i = 0; i < nstageV; i++) DVENT[i] = VENT[i];
                 DDEAD[0] = DEAD[0];
                 DRCVD[0] = RCVD[0];
                 DcasesI = casesI;
                 DcasesIQ = casesIQ;
                 DcasesH = casesH;
                 DdeathsIIQ = deathsIIQ;
                 DdeathsCV = deathsCV;

                 for (i = 0; i < nstageE; i++) DEXPDQ[i] = EXPDQ[i];
                 for (i = 0; i < nstageP; i++) DPREQ[i] = PREQ[i];
                 for (i = 0; i < nstageI; i++) DINFDQ[i] = INFDQ[i];



                 // Balance the equations
                 DSUSC[0] -= transS[0];
                 DEXPD[0] += transS[0];
                 for (i = 0; i < nstageE; i++) DEXPD[i] -= transE[i];
                 for (i = 1; i < nstageE; i++) DEXPD[i] += transE[i-1];
                 DPRE[0] += transE[nstageE-1];
                 for (i = 0; i < nstageP; i++) DPRE[i] -= transP[i];
                 for (i = 1; i < nstageP; i++) DPRE[i] += transP[i-1];
                 DINFD[0] += transP[nstageP-1];
                 DHOSP[0] += transP[nstageP];
                 DPRE[nstageP-1] -= transP[nstageP];
                 for (i = 0; i < nstageI; i++) DINFD[i] -= transI[i];
                 for (i = 1; i < nstageI; i++) DINFD[i] += transI[i-1];
                 for (i = 0; i < nstageH; i++) DHOSP[i] -= transH[i];
                 for (i = 1; i < nstageH; i++) DHOSP[i] += transH[i-1];
                 DCARE[0] += transH[nstageH];
                 DHOSP[nstageH-1] -= transH[nstageH];
                 DINFD[nstageI-1] -= transI[nstageI];
                 for (i = 0; i < nstageC; i++) DCARE[i] -= transC[i];
                 for (i = 1; i < nstageC; i++) DCARE[i] += transC[i-1];
                 DVENT[0] += transC[nstageC];
                 DCARE[nstageC-1] -= transC[nstageC] + transC[nstageC+1];
                 for (i = 0; i < nstageV; i++) DVENT[i] -= transV[i];
                 for (i = 1; i < nstageV; i++) DVENT[i] += transV[i-1];
                 DVENT[nstageV-1] -= transV[nstageV];
                 DRCVD[0] += transI[nstageI-1] + transH[nstageH-1] + transC[nstageC-1] + transV[nstageV-1];
                 DDEAD[0] += transI[nstageI] + transC[nstageC+1] + transV[nstageV];

                 DSUSC[0] -= transS[1];
                 DEXPDQ[0] += transS[1];
                 for (i = 0; i < nstageE; i++) DEXPDQ[i] -= transEQ[i];
                 for (i = 1; i < nstageE; i++) DEXPDQ[i] += transEQ[i-1];
                 DPREQ[0] += transEQ[nstageE-1];
                 for (i = 0; i < nstageP; i++) DPREQ[i] -= transPQ[i];
                 for (i = 1; i < nstageP; i++) DPREQ[i] += transPQ[i-1];
                 DPREQ[nstageP-1] -= transPQ[nstageP];
                 DPRE[nstageP-1] -= transP[nstageP+1];
                 DINFDQ[0] += transPQ[nstageP-1] + transP[nstageP+1];
                 DHOSP[0] += transPQ[nstageP];
                 for (i = 0; i < nstageI; i++) DINFDQ[i] -= transIQ[i];
                 for (i = 1; i < nstageI; i++) DINFDQ[i] += transIQ[i-1];
                 DINFDQ[nstageI-1] -= transIQ[nstageI];
                 DRCVD[0] += transIQ[nstageI-1];
                 DDEAD[0] += transIQ[nstageI];

                 DcasesI += transP[nstageP-1];
                 DcasesIQ += transPQ[nstageP-1] + transP[nstageP+1];
                 DcasesH += transPQ[nstageP] + transP[nstageP];
                 DdeathsIIQ += transI[nstageI] + transIQ[nstageI];
                 DdeathsCV += transC[nstageC+1] + transV[nstageV];
                 ")



rInit <- Csnippet("

                  double *SUSC = &S;
                  double *EXPD = &E1;
                  double *PRE = &P1;
                  double *INFD = &I1;
                  double *HOSP = &H1;
                  double *CARE = &C1;
                  double *VENT = &V1;
                  double *DEAD = &M;
                  double *RCVD = &R;

                  double *EXPDQ = &EQ1;
                  double *PREQ = &PQ1;
                  double *INFDQ = &IQ1;

                  int i;

                  double fS, fEQ, fPQ, fIQ, fE, fP, fI, fH, fC, fR, fM, fV;
                  fS = S0;
                  fEQ = EQ0/nstageE;
                  fPQ = PQ0/nstageP;
                  fIQ = IQ0/nstageI;
                  fE = E0/nstageE;
                  fP = P0/nstageP;
                  fI = I0/nstageI;
                  fH = H0/nstageH;
                  fC = C0/nstageC;
                  fV = V0/nstageV;
                  fM = M0;
                  fR = 1 - fS - nstageE*fE - nstageP*fP - nstageI*fI - nstageH*fH - nstageC*fC - nstageV*fV - nstageE*fEQ - nstageP*fPQ - nstageI*fIQ - fM;

                  SUSC[0] = nearbyint(pop*fS);
                  for (i = 0; i < nstageE; i++) EXPDQ[i] = nearbyint(pop*fEQ);
                  for (i = 0; i < nstageP; i++) PREQ[i] = nearbyint(pop*fPQ);
                  for (i = 0; i < nstageI; i++) INFDQ[i] = nearbyint(pop*fIQ);
                  for (i = 0; i < nstageE; i++) EXPD[i] = nearbyint(pop*fE);
                  for (i = 0; i < nstageP; i++) PRE[i] = nearbyint(pop*fP);
                  for (i = 0; i < nstageI; i++) INFD[i] = nearbyint(pop*fI);
                  for (i = 0; i < nstageH; i++) HOSP[i] = nearbyint(pop*fH);
                  for (i = 0; i < nstageC; i++) CARE[i] = nearbyint(pop*fC);
                  for (i = 0; i < nstageV; i++) VENT[i] = nearbyint(pop*fV);
                  DEAD[0] = nearbyint(pop*fM);
                  RCVD[0] = nearbyint(pop*fR);

                  casesI = 0.0;
                  casesIQ = 0.0;
                  casesH = 0.0;
                  deathsIIQ = 0.0;
                  deathsCV = 0.0;
                  ")

statenames <- function (nstageE=3L, nstageP=3L, nstageI=3L, nstageH=3L, nstageC=3L, nstageV=3L) {
  c("S",
    paste0("EQ",seq_len(nstageE)),
    paste0("PQ",seq_len(nstageP)),
    paste0("IQ",seq_len(nstageI)),
    paste0("E",seq_len(nstageE)),
    paste0("P",seq_len(nstageP)),
    paste0("I",seq_len(nstageI)),
    paste0("H",seq_len(nstageH)),
    paste0("C",seq_len(nstageC)),
    paste0("V",seq_len(nstageV)),
    "M", "R"
  )
}

icnames <- function () {
  c("S0","EQ0", "PQ0", "IQ0", "E0", "P0", "I0", "H0", "C0", "V0", "M0")
}

zeronames <- c("casesI", "casesIQ", "casesH", "deathsIIQ", "deathsCV")

params_log <- c("betaI", "iota","beta_sd",
                "sigma", "kappa", "gammaI", "gammaH", "gammaC", "gammaV",
                "TF")
params_logit <- c("rho", "theta",
                  "dI0", "dP0", "dT0", "dB0",
                  "dI1", "dP1", "dT1", "dB1",
                  "qP", "qH", "qC",
                  "mI", "mC", "mV")

params_mod <- c(params_log,params_logit)

params_ic <- icnames()
