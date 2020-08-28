
source("CreateModel.R")
source("CreateCovars.R")
source("CreateDataset.R")


data <- create_dataset(endTime=endTime)
covars <- create_covars(endTime=endTime)
out <- create_pomp_model(data = data, covars=covars, t0 = 0, dt=0.1) 

model <- out$model
rm(out)
