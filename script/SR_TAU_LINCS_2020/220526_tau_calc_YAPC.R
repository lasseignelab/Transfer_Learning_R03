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

query_list <- readRDS("/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/query_list.rds")
se <- readRDS("/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/se.rds")
db_path<- "/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/lincs_2020.h5"


query_list_v2 <- query_list[grep("YAPC", names(query_list))]
print("change query list to gbm")
print(length(query_list_v2))
numCores<- detectCores()
print(numCores)
###  Query signatures in the query_list against lincs database
###  To save time, the processing is paralleled with BiocParallel to run on 
###  CPU cores of a computer cluster with a scheduler (e.g. Slurm). 
#### Define submission function
f <- function(x, se, query_list_v2, query_list, dest_dir) {
  require(signatureSearch)
  chunkno <- x 
  sz <- 10 # small enough to use short queue 
  qlist <- split(query_list_v2, ceiling(seq_along(names(query_list_v2))/sz))
  #print(qlist)
  myMA <- matrix(NA, length(query_list), sz, 
                 dimnames=list(names(query_list), seq_len(sz)))
  qlistone <- qlist[[chunkno]] 
  #print(qlist[[1]])
  #print(qlistone[[1]])
  for(i in seq_along(qlistone)){
  print(i)
    qsig <- suppressMessages(qSig(query=list(upset=qlistone[[i]]$upset, 
                                            downset=qlistone[[i]]$downset), 
                 gess_method = "LINCS", refdb = db_path))
    lincs <- gess_lincs(qsig, sortby=NA, workers= 4)
    print("lincs done")
    resultDF <- result(lincs)
    ncs <- resultDF$NCS
    mynames <- paste(resultDF$pert, resultDF$cell, resultDF$type, sep="__")
    names(ncs) <- mynames
    myMA[,i] <- ncs[rownames(myMA)]
    colnames(myMA)[i] <- names(qlistone[i])
  }
  myMA <- myMA[, seq_along(qlistone)] 
  ## Only relevant for last entry that may not have as many columns as sz
  if(! dir.exists(dest_dir)) dir.create(dest_dir)
  write.table(as.data.frame(myMA), 
              file=paste0(dest_dir, "/result_YAPC_", as.character(x)), 
              col.names=NA, quote=FALSE, sep="\t")
}

print("start run")
#run through function 
f(num1, se=se,  query_list_v2=query_list_v2, query_list=query_list, "/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020")
#mclapply(num1:num2, f, se, query_list, "/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/")
print("Done run")

info<- sessionInfo()
print(info)