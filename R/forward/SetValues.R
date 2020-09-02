rm(list=ls())
mainDir <- setwd("~/Git/pomp/R/forward")
setwd(mainDir)

library(reshape2)
library(pomp)
library(ggplot2)

modelnumber <- 3
modelname <- paste0("StochModel",modelnumber)
source(paste0("ModelSnippet_",modelname,".R"))

