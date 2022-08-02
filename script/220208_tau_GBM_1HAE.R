#SR_TAU_CELL
.libPaths("/data/user/jfisher7/.conda/envs/SR_TAU_CELL/lib/R/library")

#make sure the environment is clean
rm(list=ls())


args <- commandArgs(trailingOnly = TRUE)

num1<- args[1]
num2<- args[2]

library(MASS)
library(parallel)
library(signatureSearch)
library(tidyverse)
library(HDF5Array)
library(SummarizedExperiment)
library(BiocParallel)
library(ExperimentHub)


file <- paste0("/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/result_GBM_", 1:69)
tmpDF <- read.delim(file[1], row.names=1, check.names=FALSE)
tmpDF <- round(tmpDF, 2)
for(i in 2:69){
    data<- read.delim(file[i], row.names=1, check.names=FALSE)
    data <- round(data, 2)
    tmpDF<- cbind(tmpDF, data)
}


saveRDS(tmpDF, file="/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/GI1_taudf.rds")

file <- paste0("/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/result_1HAE_", 1:32)
tmpDF <- read.delim(file[1], row.names=1, check.names=FALSE)
tmpDF <- round(tmpDF, 2)
for(i in 2:32){
    data<- read.delim(file[i], row.names=1, check.names=FALSE)
    data <- round(data, 2)
    tmpDF<- cbind(tmpDF, data)
}
saveRDS(tmpDF, file="/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/1HAE_taudf.rds")



info<- sessionInfo()
print(info)
