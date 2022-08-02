

.libPaths("/data/user/jfisher7/.conda/envs/DR_Methods/lib/R/library")

#make sure the environment is clean
rm(list=ls())


library(MASS)
library(parallel)
library(signatureSearch)
library(tidyverse)


meta42 <- readr::read_tsv("/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/siginfo_beta.txt")
dose <- "10 uM"
## filter rows by 'pert_type' as compound, 10uM concentration, and 24h treatment time
#meta42_filter <- sig_filter(meta42, pert_type="trt_cp", dose=dose, time="24 h") 


sig_filter_JLF<- function (meta, pert_type = "trt_cp", dose, time = "24 h") {
    meta %<>% dplyr::filter(pert_type == pert_type & pert_idose == 
        dose & pert_itime == time)
    meta %<>% bind_cols(alt_id = paste(meta$cmap_name, meta$cell_iname, 
        sep = "_")) %>% bind_cols(pert_cell_factor = paste(meta$cmap_name, 
        meta$cell_iname, meta$pert_type, sep = "__")) %>% distinct(alt_id, 
        .keep_all = TRUE)
    return(meta)
}

library(magrittr)
library(dplyr)
meta42_filter <- sig_filter_JLF(meta42, pert_type="trt_cp", dose=dose, time="24 h") 


#########################
## make taurefList.rds ##
#########################

# It uses all signatures in the reference database (such as LINCS) to query 
# against itself as Qref to compute tau score of `gess_lincs` 
# method in `signatureSearch` package. Tau score compares observed enrichment 
# score to all others in Qref. It represents the percentage of reference queries
# with a lower |NCS| than |NCSq,r|, adjusted to retain the sign of NCSq,r. 
# NCSq,r is the normalized connectivity score for signature r relative to 
# query q. A tau of 90 indicates that only 10 percent of reference perturbations
# showed stronger connectivity to the query. For more details, please refer to 
# Subramanian et al., 2017, Cell, A Next Generation Connectivity Map: L1000 
# Platform and the First 1,000,000 Profiles

print("create se file")
## Create Query Reference DB for Tau Score Computation of `gess_lincs` method
### Load `lincs` database created above
library(HDF5Array); library(SummarizedExperiment)
#se = loadHDF5SummarizedExperiment("/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020")
se<- SummarizedExperiment(HDF5Array("/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/lincs_2020.h5",name= "assay"))


#missing colData
#se@colData<- as.matrix(meta42_filter)
score_mat <- assay(se)
se <- SummarizedExperiment(assays=list(assay=score_mat), colData=meta42_filter)
rownames(se) <- HDF5Array("/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/lincs_2020.h5", name="rownames")
colnames(se) <- HDF5Array("/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/lincs_2020.h5", name="colnames")
score_mat <- assay(se)
### Create query list for all signatures in se
query_list <- lapply(colnames(score_mat), function(x) {
  vec <- as.matrix(score_mat[,x])
  names(vec) <- rownames(score_mat)
  sigvec = sort(vec, decreasing = TRUE)
  list(upset=utils::head(names(sigvec), 150), 
       downset=utils::tail(names(sigvec), 150))
})
names(query_list) = colData(se)$pert_cell_factor
print("done with query object")

saveRDS(query_list, "/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/query_list.rds")
saveRDS(se, "/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/data/LINCS_2020/se.rds")
info<- SessionInfo()
print(info)


