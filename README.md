
# Looking Glass Epidemiological Fitting model function using DCP
## Introduction
This package includes three major functions from R package `pomp` to fit data to the model and provide the best set of parameters with the maximum liklihood.

**pomp**:
An R package for statistical inference on partially observed Markov processes.

**mif2**:
Maximum likelihood by iterated filtering. An iterated filtering algorithm for estimating the parameters of a partially-observed Markov process.
Running ***mif2*** causes the algorithm to perform a specified number of particle-filter iterations.
At each iteration, the particle filter is performed on a perturbed version of the model, in which the parameters to be estimated are subjected to random perturbations at each observation.
This extra variability effectively smooths the likelihood surface and combats particle depletion by introducing diversity into particle population.

**pfilter**:
A plain vanilla sequential Monte Carlo (particle filter) algorithm.

**trajMatch**:
In trajectory matching, one attempts to minimize the discrepancy between a POMP model's predictions and data under the assumption that the latent state process is deterministic and all discrepancies between model and data are due to measurement error. The measurement model likelihood (dmeasure), or rather its negative, is the natural measure of the discrepancy.


## Preparing input data
The first step is to provide `covar.csv` and `data.csv` in the samples folder.

`data.csv` includes columns of *time,	reports,	deaths,	hospital,	ICU,ventilator*. `covar.csv` includes columns *time* and *tests*. They are created using *Observed_data* and also CreateCovars.js and CreateDataset.js in the library folder.
`Observed_data` includes the following columns,

  **Date:** Reported date yyyy-mm-dd.

 **Cases:**  Number of new cases per date.
    

  **TotalTests:**  Total patients approved for testing as of reporting date.

  **Deaths:**  Number of death.

  **Hospital:**  Number of patients hospitalized with COVID-19.

  **ICU:**  Number of patients in ICU with COVID-19.

  **Ventilator:**  Number of patients in ICU on a ventilator with COVID-19.
    

## Installation & Requirements
Install the dependencies with `npm install`. 

## Running the application
Now you are able to start the model 

```
node epi-fit.js
```
By default, the model works on the [KDS Distributed Computer](https://portal.distributed.computer/) and uses the keystore located in `~/.dcp/default.keystore`. You can change these DCP-related defaults using the standard DCP-CLI arguments. Detailed usage is provided by the help command:
```
node epi-fit.js --help
```



the final set of parameters will calculated and be saved as  `modelParams.csv` in the results folder.  For a description of the parameters in modelParams, see epidemiological report. This is the list of `keys` in the object, the value is always a `number`.

```
betaI
iota
beta_sd
sigma	
kappa 
gammaI
gammaH
gammaC
gammaV
TF
rho
theta
dI0
dP0
dT0
dB0
dI1
dP1
dT1
dB1
dI2
dP2
dT2
dB2
qP
qH
qC
mI
mC
mV
S0
EQ0
PQ0
IQ0
E0
P0
I0
H0
C0
V0
M0
```
`epi-fit.js`  also returns one more CSV file (savedStates.csv) that will be used later as an input in the forward model.

## Forward Function in R

In the R-forward-function folder `forwardFunction.R` returns `predicion.csv` and saves it in the simulations folder. It uses data.csv, covar.csv and also savedStates.csv.

```
forward_model(nsim, mainDir, path_to_params, endTime, predTime)
```
**nsim:**  Number of simulations.

**mainDir** Main working directory.

**path_to_params** Path to the CSV file includes modified parameters from RT.

**endTime** End time in data and start time of the prediction (yyyy-mm-dd).

**predTime** End time of prediction (yyyy-mm-dd).

