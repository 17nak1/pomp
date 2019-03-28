rm(list=ls())

library(pomp)
library(plyr)
library(reshape2)
library(magrittr)
mainDir <- "~/Documents/R/StochasticSEIR"


# loc <- "/home/fmagpant/pomp_devel"
# library("pomp",lib.loc=loc)
# library("plyr",lib.loc=loc)
# library("reshape2",lib.loc=loc)
# library("magrittr",lib.loc=loc)
# mainDir <- getwd()

setwd(mainDir)
source("ModelSnippet.R")