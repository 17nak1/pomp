rm(list=ls())
mainDir <- "~/Git/pomp/R/stochmodel3"
# mainDir <- getwd()
setwd(mainDir)

library(reshape2)
library(pomp)
library(ggplot2)

# loc <- file.path("/home/fmagpant/pomp_devel")
# library("digest",lib.loc=loc)
# library("mvtnorm",lib.loc=loc)
# library("deSolve",lib.loc=loc)
# library("coda",lib.loc=loc)
# library("subplex",lib.loc=loc)
# library("nloptr",lib.loc=loc)
# library("pomp",lib.loc=loc)

modelnumber <- 3
modelname <- paste0("StochModel",modelnumber)
source(paste0("ModelSnippet_",modelname,".R"))

