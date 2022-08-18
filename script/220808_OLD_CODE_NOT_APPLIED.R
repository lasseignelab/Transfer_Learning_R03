#RANDOM CODE JUST STORE BUT NOT USED IN PROJECT CURRENTLY
test the go term stuff

go_gbm_up_sim

```{r}
go_sim_median<- function(go_sim_matrix, path_list,  title , method){
  cos_GI1_cand_tfl <- go_sim_matrix[grep(paste(path_list,collapse="|"), colnames(go_sim_matrix)), grep(paste(path_list,collapse="|"), colnames(go_sim_matrix))]
  
  cos_GI1_cand_tfl[cos_GI1_cand_tfl == 1] <- NA
  
  cos_median <- rowMedians(cos_GI1_cand_tfl, na.rm=TRUE)
  meth<- rep(method, nrow(cos_GI1_cand_tfl))
  
  #cos_median <- median(cos_GI1_cand_tfl, na.rm=TRUE)
  #meth<- rep(method, length(cos_GI1_cand_tfl))
  
  
  
  cos_median<- data.frame(drug = rownames(cos_GI1_cand_tfl), cos_sim_med= cos_median, method= meth)
  
  cos_GI1_cand_tfl<- as.data.frame(cos_GI1_cand_tfl)
  cos_GI1_cand_tfl$names <-  rownames(cos_GI1_cand_tfl)
  cos_GI1_cand_tfl_long <- cos_GI1_cand_tfl %>% pivot_longer(!names, names_to = "other")
  #p <- ggplot(cos_GI1_cand_tfl_long, aes(x=names, y=value)) + 
  #geom_boxplot()+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  #labs(title=title,x="pathways", y = "Tanimoto Coefficient")
  #print(p)
  
  return(cos_median)
}
```


```{r}
go_gbm_med_tfl <- go_sim_median(go_gbm_up_sim, TFL_bp_up, "gbm_up_paths", "Transfer Learning")
```

```{r}
go_gbm_med_limma <- go_sim_median(go_gbm_up_sim, limma_bp_up, "gbm_up_paths", "limma")
```

```{r}
go_gbm_med_deseq2 <- go_sim_median(go_gbm_up_sim, deseq2_bp_up, "gbm_up_paths", "deseq2")
```


```{r}
go_gbm_med <- rbind(go_gbm_med_tfl, go_gbm_med_limma, go_gbm_med_deseq2)
```

```{r}
ggplot(go_gbm_med, aes(x=method, y=cos_sim_med, fill=method)) + 
  geom_violin()+ geom_point() +labs(title="GOSim (Wang) of GBM Input Enriched Pathways",x="Method", y = "Median GOSim (Wang)", color= "Method") + scale_colour_manual(values =  c("#440154FF", "#31688EFF", "#FDE725FF"), aesthetics = c("colour", "fill"))
```
```{r}
go_sim_median_others <- function(go_sim_matrix, drug_list_testing_method, drug_list_method1, drug_list_method2,  title , method){
  drug_list_others<- c(drug_list_method1, drug_list_method2)
  cos_GI1_cand_tfl <- go_sim_matrix[grep(paste(drug_list_testing_method,collapse="|"), colnames(go_sim_matrix)), grep(paste(drug_list_others,collapse="|"), colnames(go_sim_matrix))]
  
  cos_GI1_cand_tfl[cos_GI1_cand_tfl == 1] <- NA
  cos_median <- rowMedians(cos_GI1_cand_tfl, na.rm=TRUE)
  meth<- rep(method, nrow(cos_GI1_cand_tfl))
  
  cos_median<- data.frame(drug = rownames(cos_GI1_cand_tfl), cos_sim_med= cos_median, method= meth)
  
  cos_GI1_cand_tfl<- as.data.frame(cos_GI1_cand_tfl)
  cos_GI1_cand_tfl$names <-  rownames(cos_GI1_cand_tfl)
  cos_GI1_cand_tfl_long <- cos_GI1_cand_tfl %>% pivot_longer(!names, names_to = "other")
  #p <- ggplot(cos_GI1_cand_tfl_long, aes(x=names, y=value)) + 
  #geom_boxplot()+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  #labs(title=title,x="Drug Candidates", y = "Tanimoto Coefficient of drugs of other methods")
  #print(p)
  
  return(cos_median)
}
```

```{r}
go_gbm_med_tfl_others <- go_sim_median_others(go_gbm_up_sim, TFL_bp_up, limma_bp_up, deseq2_bp_up, "gbm_up_paths", "Transfer Learning_other")
```

```{r}
go_gbm_med_limma_others <- go_sim_median_others(go_gbm_up_sim, limma_bp_up, TFL_bp_up, deseq2_bp_up,"gbm_up_paths", "limma_other")
```

```{r}
go_gbm_med_deseq2_others <- go_sim_median_others(go_gbm_up_sim, deseq2_bp_up,TFL_bp_up, limma_bp_up, "gbm_up_paths", "deseq2_other")
```


```{r}
go_gbm_med_others <- rbind(go_gbm_med_tfl_others, go_gbm_med_limma_others, go_gbm_med_deseq2_others)
```

```{r}
go_sim_all_gbm<- rbind(go_gbm_med_others, go_gbm_med)
```

```{r}
ggplot(go_sim_all_gbm, aes(x=method, y=cos_sim_med, fill=method)) + 
  geom_violin()+ geom_point() +labs(x="Method", y = "Median GOSim (Wang Method)", color= "Method") + scale_fill_viridis(discrete=TRUE)
```


```{r}
compare_methods_go_sim<- function(limma_pathways, deseq2_pathways, TFL_pathways, list_type){
  
  limma_bp<- limma_pathways$term_id[ limma_pathways$source == "GO:BP" & limma_pathways$set == list_type]
  deseq2_bp <- deseq2_pathways$term_id[ deseq2_pathways$source == "GO:BP" & deseq2_pathways$set == list_type]
  #print(deseq2_bp)
  TFL_bp <- TFL_pathways$term_id[ TFL_pathways$source == "GO:BP" & TFL_pathways$set == list_type]
  
  go1 <-  unique(c(TFL_bp, deseq2_bp, limma_bp ))
  
  
  go_gbm_up_sim <- mgoSim(go1, go1, semData=hsGO, measure="Wang", combine=NULL)
  
  # tfl_drugs <- grepl(paste(TFL_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  # deseq2_drugs<- grepl(paste(deseq2_bp,collapse="|"), colnames(go_gbm_up_sim))
  # limma_drugs <- grepl(paste(limma_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  # 
  # if(length(deseq2_bp) ==0){ deseq2_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  # if(length(TFL_bp) ==0){ tfl_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  # if(length(limma_bp) ==0){ limma_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  
  
  
  #within
  go_gbm_med_tfl <- go_sim_median(go_gbm_up_sim, TFL_bp, "gbm_up_paths", "Transfer Learning")
  
  go_gbm_med_limma <- go_sim_median(go_gbm_up_sim, limma_bp, "gbm_up_paths", "limma")
  
  go_gbm_med_deseq2 <- go_sim_median(go_gbm_up_sim, deseq2_bp, "gbm_up_paths", "deseq2")
  #print(go_gbm_med_deseq2)
  go_gbm_med <- rbind(go_gbm_med_tfl, go_gbm_med_limma, go_gbm_med_deseq2)
  
  #ggplot(go_gbm_med, aes(x=method, y=cos_sim_med, fill=method)) + 
  #  geom_violin()+ geom_point() +labs(title="GOSim (Wang) of GBM Input Enriched Pathways",x="Method", y = "Median GOSim (Wang)", color= "Method") + scale_colour_manual(values =  c("#440154FF", "#31688EFF", "#FDE725FF"), aesthetics = c("colour", "fill"))
  
  #outside
  go_gbm_med_tfl_others <- go_sim_median_others(go_gbm_up_sim, TFL_bp, limma_bp, deseq2_bp, "gbm_up_paths", "Transfer Learning_other")
  
  go_gbm_med_limma_others <- go_sim_median_others(go_gbm_up_sim, limma_bp, TFL_bp, deseq2_bp,"gbm_up_paths", "limma_other")
  
  go_gbm_med_deseq2_others <- go_sim_median_others(go_gbm_up_sim, deseq2_bp,TFL_bp, limma_bp, "gbm_up_paths", "deseq2_other")
  
  go_gbm_med_others <- rbind(go_gbm_med_tfl_others, go_gbm_med_limma_others, go_gbm_med_deseq2_others)
  
  go_sim_all_gbm<- rbind(go_gbm_med_others, go_gbm_med)
  
  if(length(deseq2_bp) == 0){go_sim_all_gbm<- go_sim_all_gbm[!go_sim_all_gbm$method %in% c("deseq2", "deseq2_other"),]}
  if(length(limma_bp) == 0){go_sim_all_gbm<- go_sim_all_gbm[!go_sim_all_gbm$method %in% c("limma", "limma_other"),]}
  if(length(TFL_bp) == 0){go_sim_all_gbm<- go_sim_all_gbm[!go_sim_all_gbm$method %in% c("Transfer Learning", "Transfer Learning_other"),]}
  
  
  ggplot(go_sim_all_gbm, aes(x=method, y=cos_sim_med, fill=method)) + 
    geom_violin()+ geom_point() +labs(x="Method", y = "Median GOSim (Wang Method)", color= "Method") + scale_fill_viridis(discrete=TRUE)
  
}
```
# liver 
```{r}
TFL_pathways <- read_delim("~/output/liver_cancer/TFL_pathways.csv", 
                           delim = ";", escape_double = FALSE, trim_ws = TRUE)

deseq2_pathways1 <-read_delim("~/output/liver_cancer/DESEq2_pathways.csv", 
                              delim = ";", escape_double = FALSE, trim_ws = TRUE)

limma_pathways <- read_delim("~/output/liver_cancer/limma_pathways.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
```

```{r}
compare_methods_go_sim(limma_pathways, deseq2_pathways1, TFL_pathways, "Down")
```

```{r}
compare_methods_go_sim(limma_pathways, deseq2_pathways1, TFL_pathways, "Up")
```

# pancreatic 

```{r}
TFL_pathways <- read_delim("~/output/pancreas_cancer/TFL_pathways.csv", 
                           delim = ";", escape_double = FALSE, trim_ws = TRUE)

deseq2_pathways1 <-read_delim("~/output/pancreas_cancer/DESEq2_pathways.csv", 
                              delim = ";", escape_double = FALSE, trim_ws = TRUE)

limma_pathways <- read_delim("~/output/pancreas_cancer/limma_pathways.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
```

```{r}
compare_methods_go_sim(limma_pathways, deseq2_pathways1, TFL_pathways, "Down")
```

```{r}
compare_methods_go_sim(limma_pathways= limma_pathways,deseq2_pathways1, TFL_pathways,  list_type = "Up")
```

GBM
```{r}
TFL_pathways <- read_delim("~/output/gbm_plots/TFL_pathways.csv", 
                           delim = ";", escape_double = FALSE, trim_ws = TRUE)

deseq2_pathways <- read_delim("~/output/gbm_plots/DESEq2_pathways.csv", 
                              delim = ";", escape_double = FALSE, trim_ws = TRUE)

limma_pathways <- read_delim("~/output/gbm_plots/limma_pathways.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
```

```{r}
compare_methods_go_sim(limma_pathways, deseq2_pathways, TFL_pathways, "Down")
```

```{r}
compare_methods_go_sim(limma_pathways= limma_pathways,deseq2_pathways, TFL_pathways,  list_type = "Up")
```
#lungs
```{r}
TFL_pathways <- read_delim("~/output/lung_cancer/TFL_pathways.csv", 
                           delim = ";", escape_double = FALSE, trim_ws = TRUE)

deseq2_pathways1 <-read_delim("~/output/lung_cancer/DESEq2_pathways.csv", 
                              delim = ";", escape_double = FALSE, trim_ws = TRUE)

limma_pathways <- read_delim("~/output/lung_cancer/limma_pathways.csv", 
                             delim = ";", escape_double = FALSE, trim_ws = TRUE)
```

```{r}
compare_methods_go_sim(limma_pathways, deseq2_pathways1, TFL_pathways, "Down")
```

```{r}
compare_methods_go_sim_quick_fix <- function(limma_pathways, deseq2_pathways, TFL_pathways, list_type){
  
  limma_bp<- limma_pathways$term_id[ limma_pathways$source == "GO:BP" & limma_pathways$set == list_type]
  deseq2_bp <- deseq2_pathways$term_id[ deseq2_pathways$source == "GO:BP" & deseq2_pathways$set == list_type]
  #print(deseq2_bp)
  TFL_bp <- TFL_pathways$term_id[ TFL_pathways$source == "GO:BP" & TFL_pathways$set == list_type]
  
  go1 <-  unique(c(TFL_bp, deseq2_bp, limma_bp ))
  
  
  go_gbm_up_sim <- mgoSim(go1, go1, semData=hsGO, measure="Wang", combine=NULL)
  
  # tfl_drugs <- grepl(paste(TFL_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  # deseq2_drugs<- grepl(paste(deseq2_bp,collapse="|"), colnames(go_gbm_up_sim))
  # limma_drugs <- grepl(paste(limma_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  # 
  # if(length(deseq2_bp) ==0){ deseq2_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  # if(length(TFL_bp) ==0){ tfl_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  # if(length(limma_bp) ==0){ limma_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  
  
  
  #within
  go_gbm_med_tfl <- go_sim_median(go_gbm_up_sim, TFL_bp, "gbm_up_paths", "Transfer Learning")
  
  #go_gbm_med_limma <- go_sim_median(go_gbm_up_sim, limma_bp, "gbm_up_paths", "limma")
  
  go_gbm_med_deseq2 <- go_sim_median(go_gbm_up_sim, deseq2_bp, "gbm_up_paths", "deseq2")
  #print(go_gbm_med_deseq2)
  go_gbm_med <- rbind(go_gbm_med_tfl, go_gbm_med_deseq2)
  
  #ggplot(go_gbm_med, aes(x=method, y=cos_sim_med, fill=method)) + 
  #  geom_violin()+ geom_point() +labs(title="GOSim (Wang) of GBM Input Enriched Pathways",x="Method", y = "Median GOSim (Wang)", color= "Method") + scale_colour_manual(values =  c("#440154FF", "#31688EFF", "#FDE725FF"), aesthetics = c("colour", "fill"))
  
  
  print(TFL_bp)
  print(limma_bp)
  print(deseq2_bp)
  #outside
  go_gbm_med_tfl_others <- go_sim_median_others(go_gbm_up_sim, TFL_bp, limma_bp, deseq2_bp, "gbm_up_paths", "Transfer Learning_other")
  
  #go_gbm_med_limma_others <- go_sim_median_others(go_gbm_up_sim, limma_bp, TFL_bp, deseq2_bp,"gbm_up_paths", "limma_other")
  
  go_gbm_med_deseq2_others <- go_sim_median_others(go_gbm_up_sim, deseq2_bp,TFL_bp, limma_bp, "gbm_up_paths", "deseq2_other")
  
  go_gbm_med_others <- rbind(go_gbm_med_tfl_others,  go_gbm_med_deseq2_others)
  
  go_sim_all_gbm<- rbind(go_gbm_med_others, go_gbm_med)
  
  if(length(deseq2_bp) == 0){go_sim_all_gbm<- go_sim_all_gbm[!go_sim_all_gbm$method %in% c("deseq2", "deseq2_other"),]}
  # if(length(limma_bp) == 0){go_sim_all_gbm<- go_sim_all_gbm[!go_sim_all_gbm$method %in% c("limma", "limma_other"),]}
  if(length(TFL_bp) == 0){go_sim_all_gbm<- go_sim_all_gbm[!go_sim_all_gbm$method %in% c("Transfer Learning", "Transfer Learning_other"),]}
  
  
  ggplot(go_sim_all_gbm, aes(x=method, y=cos_sim_med, fill=method)) + 
    geom_violin()+ geom_point() +labs(x="Method", y = "Median GOSim (Wang Method)", color= "Method") + scale_fill_viridis(discrete=TRUE)
  
}
```

```{r}
compare_methods_go_sim_quick_fix(limma_pathways= limma_pathways,deseq2_pathways1, TFL_pathways,  list_type = "Up")
```

```{r}
library(readr)
test <- read_csv("~/output/liver_cancer/deseq2_res/220602_SR_LINCS_LIHC_DESEQ2_RES_tau.csv")
```

```{r}

```


```{r}
test_fda <- test[test$fda_approved,]
lincs_commpound_info[lincs_commpound_info$cmap_name %in% test_fda$pert,]
```

```{r}
summ_df <- as.data.frame(matrix(nrow=1, ncol=3))
colnames(summ_df)<- c("Var1", "Var2", "Var3")
for (i in 1:nrow(test_fda)){
  #get the target from the data results data frame
  targets1 <- as.vector(unlist(strsplit(test_fda$t_gn_sym[i], "; ")))
  
  targets2 <- lincs_commpound_info$target[lincs_commpound_info$cmap_name == test_fda$pert[i]]
  targets2<- ifelse(targets2 == "", NA, targets2)
  targets<- unique(c(targets1, targets2))
  
  moa <- unique(lincs_commpound_info$moa[lincs_commpound_info$cmap_name == test_fda$pert[i]])
  moa<- ifelse(moa == "", NA, moa)
  
  df<- expand.grid(test_fda$pert[i],targets, moa )
  
  summ_df<- rbind(summ_df, df)
  
}
summ_df<- summ_df[-1,]
colnames(summ_df)<- c( "cmap_name", "target", "moa")
```

```{r}
summ_df<- summ_df[!duplicated(summ_df), ]
```

```{r}
TEST_LODES<- to_lodes_form(as.data.frame(summ_df),
                           axes = 1:3,
                           id = "INDEX")
```

```{r}
TEST_LODES$stratum = str_wrap(TEST_LODES$stratum, width = 5)
TEST_LODES$x <- factor(TEST_LODES$x, levels= c("moa", "cmap_name", "target"))
```

```{r}
TEST_LODES$drug <- rep(NA, nrow(summ_df))
drugs<- unique(summ_df$cmap_name)
for (i in 1:length(drugs)){
  ind <- TEST_LODES$INDEX[TEST_LODES$stratum == drugs[i]]
  TEST_LODES$drug<- ifelse(TEST_LODES$INDEX %in% ind, drugs[i], TEST_LODES$drug)
}
```

```{r}
#not need for this analysis but many for future.
TEST_LODES$drug<- factor(TEST_LODES$drug, levels= test_fda$pert)
```

```{r}
#moa 
moa_df <- TEST_LODES[TEST_LODES$x == "moa", 2:4]
moa_df<- moa_df[!duplicated(moa_df), ]
counts<- as.data.frame(table(moa_df$stratum))

#counts<- as.data.frame(table(TEST_LODES$stratum[TEST_LODES$x == "moa"]))
moa_order <- counts$Var1[ order(-counts$Freq)]
counts_moa<- counts
#targets
target_df <- TEST_LODES[TEST_LODES$x == "target", 2:4]
target_df<- target_df[!duplicated(target_df), ]
counts<- as.data.frame(table(target_df$stratum))
#counts
target_order <- counts$Var1[ order(-counts$Freq)]

#factor the stratum column
TEST_LODES$stratum<- factor(TEST_LODES$stratum, levels= unique(c(test_fda$pert, as.character(moa_order), as.character(target_order))))
```

```{r}
TEST_LODES_com <- TEST_LODES[ complete.cases(TEST_LODES),]
#TEST_LODES_com$index2 <- as.numeric(TEST_LODES_com$drug)
#TEST_LODES_com_v2 <- TEST_LODES_com[!TEST_LODES_com$x == "cmap_name",]

ggplot(TEST_LODES_com,
       aes(x = x, stratum = stratum, alluvium = INDEX, label= stratum))+ geom_flow( aes(fill= drug), stat = "alluvium") +
  geom_stratum(aes(fill= drug)) +geom_label(stat = "stratum", fill="white") + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.position="none") + scale_y_continuous(breaks=NULL) + scale_x_discrete(labels=c("moa" = "Mechanism of Action", "cmap_name" = "Drug", "target" = "Target")) +
  scale_fill_manual(
    values=viridis(n=length(drugs)),
    breaks=drugs,
    labels=drugs,
    na.value = NA
  )
```

```{r}
counts_moa_v2<- counts_moa[counts_moa$Freq >0,]
ggplot(counts_moa_v2, aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Mechanism of Action") + xlab("Number of Drugs") 
```


```{r}
ggplot(counts, aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Targets") + xlab("Number of Drugs") 
```


```{r}
summ_df_moa<- summ_df[,c(1,3)]
summ_df_moa<- summ_df_moa[!duplicated(summ_df_moa), ]
TEST_LODES_moa<- to_lodes_form(as.data.frame(summ_df_moa),
                               axes = 1:2,
                               id = "INDEX")
```

```{r}
TEST_LODES_moa$stratum = str_wrap(TEST_LODES_moa$stratum, width = 5)
```

```{r}
TEST_LODES_moa$drug <- rep(NA, nrow(summ_df_moa))
drugs<- unique(summ_df_moa$cmap_name)
for (i in 1:length(drugs)){
  ind <- TEST_LODES_moa$INDEX[TEST_LODES_moa$stratum == drugs[i]]
  TEST_LODES_moa$drug<- ifelse(TEST_LODES_moa$INDEX %in% ind, drugs[i], TEST_LODES_moa$drug)
}
```

```{r}
#not need for this analysis but many for future.
TEST_LODES_moa$drug<- factor(TEST_LODES_moa$drug, levels= test_fda$pert)
```

```{r}
#moa 
moa_df <- TEST_LODES_moa[TEST_LODES_moa$x == "moa", 2:4]
moa_df<- moa_df[!duplicated(moa_df), ]
counts<- as.data.frame(table(moa_df$stratum))

#counts<- as.data.frame(table(TEST_LODES_moa$stratum[TEST_LODES_moa$x == "moa"]))
moa_order <- counts$Var1[ order(-counts$Freq)]
counts_moa<- counts

#factor the stratum column
TEST_LODES_moa$stratum<- factor(TEST_LODES_moa$stratum, levels= unique(c(test_fda$pert, as.character(moa_order))))
```

```{r}
TEST_LODES_moa_com <- TEST_LODES_moa[ complete.cases(TEST_LODES_moa),]
#TEST_LODES_com$index2 <- as.numeric(TEST_LODES_com$drug)
#TEST_LODES_com_v2 <- TEST_LODES_com[!TEST_LODES_com$x == "cmap_name",]

ggplot(TEST_LODES_moa_com,
       aes(x = x, stratum = stratum, alluvium = INDEX, label= stratum))+ geom_flow( aes(fill= drug), stat = "alluvium", color= "black") +
  geom_stratum(aes(fill= drug)) +geom_label(stat = "stratum", fill="white") + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.position="none") + scale_y_continuous(breaks=NULL) + scale_x_discrete(labels=c("cmap_name" = "Drug", "moa" = "Mechanism of Action")) +
  scale_fill_manual(
    values=viridis(n=length(drugs)),
    breaks=drugs,
    labels=drugs,
    na.value = NA
  )
```

```{r}
summ_df_tar<- summ_df[,c(1,2)]
summ_df_tar<- summ_df_tar[!duplicated(summ_df_tar), ]
TEST_LODES_tar<- to_lodes_form(as.data.frame(summ_df_tar),
                               axes = 1:2,
                               id = "INDEX")
```

```{r}
TEST_LODES_tar$stratum = str_wrap(TEST_LODES_tar$stratum, width = 5)
```

```{r}
TEST_LODES_tar$drug <- rep(NA, nrow(summ_df_tar))
drugs<- unique(summ_df_tar$cmap_name)
for (i in 1:length(drugs)){
  ind <- TEST_LODES_tar$INDEX[TEST_LODES_tar$stratum == drugs[i]]
  TEST_LODES_tar$drug<- ifelse(TEST_LODES_tar$INDEX %in% ind, drugs[i], TEST_LODES_tar$drug)
}
```

```{r}
#not need for this analysis but many for future.
TEST_LODES_tar$drug<- factor(TEST_LODES_tar$drug, levels= test_fda$pert)
```

```{r}
#targets
target_df <- TEST_LODES[TEST_LODES$x == "target", 2:4]
target_df<- target_df[!duplicated(target_df), ]
counts<- as.data.frame(table(target_df$stratum))
#counts
target_order <- counts$Var1[ order(-counts$Freq)]


#factor the stratum column
TEST_LODES_tar$stratum<- factor(TEST_LODES_tar$stratum, levels= unique(c(test_fda$pert, as.character(target_order))))
```

```{r}
TEST_LODES_tar_com <- TEST_LODES_tar[ complete.cases(TEST_LODES_tar),]
#TEST_LODES_com$index2 <- as.numeric(TEST_LODES_com$drug)
#TEST_LODES_com_v2 <- TEST_LODES_com[!TEST_LODES_com$x == "cmap_name",]

ggplot(TEST_LODES_tar_com,
       aes(x = x, stratum = stratum, alluvium = INDEX, label= stratum))+ geom_flow( aes(fill= drug), stat = "alluvium", color= "black") +
  geom_stratum(aes(fill= drug)) +geom_label(stat = "stratum", fill="white") + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.position="none") + scale_y_continuous(breaks=NULL) + scale_x_discrete(labels=c("cmap_name" = "Drug", "moa" = "Mechanism of Action")) +
  scale_fill_manual(
    values=viridis(n=length(drugs)),
    breaks=drugs,
    labels=drugs,
    na.value = NA
  )
```

for each of the drug from the same method what is the average cosine siamilirty
```{r}
cos_GI1_cand_tfl <- cos_GI1[grep(paste(lincs_tfl_gbm_res_GI1_fda_approved$pert,collapse="|"), colnames(cos_GI1)), grep(paste(lincs_tfl_gbm_res_GI1_fda_approved$pert,collapse="|"), colnames(cos_GI1))]
```

```{r}
cos_GI1_cand_tfl[cos_GI1_cand_tfl == 1] <- NA
```

```{r}
rowMedians(cos_GI1_cand_tfl, na.rm=TRUE)
boxplot(cos_GI1_cand_tfl)
```
```{r}
cos_GI1_cand_tfl<- as.data.frame(cos_GI1_cand_tfl)
cos_GI1_cand_tfl$names <-  str_split(rownames(cos_GI1_cand_tfl), "__GI1", simplify = TRUE)[,1]
cos_GI1_cand_tfl_long <- cos_GI1_cand_tfl %>% pivot_longer(!names, names_to = "other")
```

```{r}
p <- ggplot(cos_GI1_cand_tfl_long, aes(x=names, y=value)) + 
  geom_boxplot()+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  labs(title="Cosine Similarity of Transfer Learning GBM candidates",x="Drug Candidates", y = "Cosine Similarity")
p
```

```{r}
plot_cosine_similarity_boxplot<- function(cos_sim_matrix, drug_list, pattern_cell= "__GI1",  title , method){
  cos_GI1_cand_tfl <- cos_sim_matrix[grep(paste(drug_list,collapse="|"), colnames(cos_sim_matrix)), grep(paste(drug_list,collapse="|"), colnames(cos_sim_matrix))]
  
  cos_GI1_cand_tfl[cos_GI1_cand_tfl == 1] <- NA
  cos_median <- rowMedians(cos_GI1_cand_tfl, na.rm=TRUE)
  meth<- rep(method, nrow(cos_GI1_cand_tfl))
  
  cos_median<- data.frame(drug = str_split(rownames(cos_GI1_cand_tfl), pattern_cell, simplify = TRUE)[,1], cos_sim_med= cos_median, method= meth)
  
  cos_GI1_cand_tfl<- as.data.frame(cos_GI1_cand_tfl)
  cos_GI1_cand_tfl$names <-  str_split(rownames(cos_GI1_cand_tfl), pattern_cell, simplify = TRUE)[,1]
  cos_GI1_cand_tfl_long <- cos_GI1_cand_tfl %>% pivot_longer(!names, names_to = "other")
  p <- ggplot(cos_GI1_cand_tfl_long, aes(x=names, y=value)) + 
    geom_boxplot()+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
    labs(title=title,x="Drug Candidates", y = "Cosine Similarity")
  print(p)
  
  return(cos_median)
}
```

```{r}
cos_median_tfl<- plot_cosine_similarity_boxplot(cos_GI1,lincs_tfl_gbm_res_GI1_fda_approved$pert, "__GI1", "Cosine Similarity of Transfer Learning GBM candidates" , "Transfer Learning")
```
```{r}
cos_median_limma<- plot_cosine_similarity_boxplot(cos_GI1,lincs_limma_gbm_res_GI1_fda_approved$pert, "__GI1", "Cosine Similarity of limma GBM candidates", "limma" )
```

```{r}
cos_median_deseq2<- plot_cosine_similarity_boxplot(cos_GI1,lincs_deseq_gbm_res_GI1_fda_approved$pert, "__GI1", "Cosine Similarity of DESeq2 GBM candidates", "DESeq2" )
```

```{r}
cos_median_gbm <- rbind(cos_median_deseq2, cos_median_limma, cos_median_tfl)
```

```{r}
ggplot(cos_median_gbm, aes(x=method, y=cos_sim_med, fill=method)) + 
  geom_violin()+ geom_point() +labs(title="Median Cosine Similarity of GBM candidates",x="Method", y = "Median Cosine Similarity", color= "Method") + scale_colour_manual(values =  c("#440154FF", "#31688EFF", "#FDE725FF"), aesthetics = c("colour", "fill"))
```
```{r}
kruskal.test(cos_sim_med ~ method, data = cos_median_gbm)
```

```{r}
pairwise.wilcox.test(cos_median_gbm$cos_sim_med, cos_median_gbm$method,
                     p.adjust.method = "BH")
```


```{r}
cosine_similarity_others <- function(cos_sim_matrix, drug_list, pattern_cell= "__GI1",  title , method){
  cos_median<- c()
  for (i in 1:length(drug_list)){
    drug<- drug_list[i]
    not_list <- drug_list[-i]
    #print(not_list)
    cos_GI1_cand_tfl <- cos_sim_matrix[!grepl(paste(not_list ,collapse="|"), colnames(cos_sim_matrix)), grepl(drug, colnames(cos_sim_matrix))]
    
    cos_GI1_cand_tfl[cos_GI1_cand_tfl == 1] <- NA
    cos_median[i] <- median(cos_GI1_cand_tfl, na.rm=TRUE)
    #print(cos_median[i])
  }
  
  
  meth<- rep(method, length(drug_list))
  #make a dataframe
  cos_median<- data.frame(drug =drug_list, cos_sim_med= cos_median, method= meth)
  
  
  return(cos_median)
}
```


```{r}
cos_median_tfl_v2<- cosine_similarity_others(cos_GI1,lincs_tfl_gbm_res_GI1_fda_approved$pert, "__GI1", "Cosine Similarity of Transfer Learning GBM candidates" , "Transfer Learning")
```



```{r}
cos_median_limma_v2 <- cosine_similarity_others(cos_GI1,lincs_limma_gbm_res_GI1_fda_approved$pert, "__GI1", "Cosine Similarity of limma GBM candidates", "limma" )
```

```{r}
cos_median_deseq2_v2<- cosine_similarity_others(cos_GI1,lincs_deseq_gbm_res_GI1_fda_approved$pert, "__GI1", "Cosine Similarity of DESeq2 GBM candidates", "DESeq2" )
```

```{r}
cos_median_gbm_v2 <- rbind(cos_median_deseq2_v2, cos_median_limma_v2, cos_median_tfl_v2)
```

```{r}
ggplot(cos_median_gbm_v2, aes(x=method, y=cos_sim_med, fill=method)) + 
  geom_violin()+ geom_point() +labs(title="Median Cosine Similarity of GBM candidates",x="Method", y = "Median Cosine Similarity", color= "Method") + scale_colour_manual(values =  c("#440154FF", "#31688EFF", "#FDE725FF"), aesthetics = c("colour", "fill"))
```
```{r}
kruskal.test(cos_sim_med ~ method, data = cos_median_gbm_v2)
```

```{r}
pairwise.wilcox.test(cos_median_gbm_v2$cos_sim_med, cos_median_gbm_v2$method,
                     p.adjust.method = "BH")
```


```{r}
cos_median_gbm_v3<- cos_median_gbm_v2
```
```{r}
cos_median_gbm_v3$diff_cos <- cos_median_gbm$cos_sim_med - cos_median_gbm_v2$cos_sim_med
```

```{r}
ggplot(cos_median_gbm_v3, aes(x=method, y=diff_cos, fill=method)) + 
  geom_violin()+ geom_point() +labs(title="Difference of Median Cosine Similarity of GBM candidates",x="Method", y = "Difference of Median Cosine Similarity", color= "Method") + scale_colour_manual(values =  c("#440154FF", "#31688EFF", "#FDE725FF"), aesthetics = c("colour", "fill"))
```

```{r}
kruskal.test(diff_cos ~ method, data = cos_median_gbm_v3)
```

```{r}
pairwise.wilcox.test(cos_median_gbm_v3$diff_cos, cos_median_gbm_v3$method,
                     p.adjust.method = "BH")
```


```{r}
cos_GI1_table<-  as.data.frame.table(cos_GI1)
```

```{r}
ggplot(cos_GI1_table, aes(x= Freq)) +geom_histogram()
```


```{r}
fda_info <- read.table("~/data/fda_product_info_df.csv" )
fda_approve <- fda_info[fda_info$market_status_v2 %in% c("Over-the-Counter", "Prescription" ),]
```

```{r}
FDA_APPROVAL_CHECK<- function(drug_list){
  return_list <- c()
  for (i in 1:length(drug_list)){
    test <- grep( drug_list[i], fda_approve$ActiveIngredient, ignore.case= TRUE)
    if (length(test)>0 ){
      return_list[i] <- TRUE
    }else{
      return_list[i] <- FALSE
    }
    if (nchar(drug_list[i])< 3){
      return_list[i] <- FALSE
    }
  }
  names(return_list)<- drug_list
  return(return_list)
}
```



```{r}
cosine_sim_random_test<- function(med_cos_sim, num_drug, cos_sim_matrix, seed= 101){
  
  gbm_lincs_drugs <- as.vector(str_split(colnames(cos_sim_matrix),  "__", simplify = TRUE)[,1])
  
  gbm_lincs_drugs<- gbm_lincs_drugs[FDA_APPROVAL_CHECK(gbm_lincs_drugs)]
  number_sig_drugs<- num_drug
  # do that again 10,000 times 
  set.seed(seed)
  fraction<-c()
  for (i in 1:10000){
    
    test <- gbm_lincs_drugs[ sample(x = 1:length(gbm_lincs_drugs),size = number_sig_drugs ,replace = FALSE)]
    
    cos_GI1_cand_tfl <- cos_sim_matrix[grep(paste(test,collapse="|"), colnames(cos_sim_matrix)), grep(paste(test,collapse="|"), colnames(cos_sim_matrix))]
    
    cos_GI1_cand_tfl[cos_GI1_cand_tfl == 1] <- NA
    cos_median <- rowMedians(cos_GI1_cand_tfl, na.rm=TRUE)
    
    
    fraction[i] <- median( cos_median )
    
  }
  
  plot_data<- as.data.frame(cbind(1:10000, fraction))
  
  plot <- ggplot(plot_data, aes(x=fraction)) + geom_histogram() + geom_vline(xintercept = med_cos_sim)+ labs( x = "median cosine similarity of all drugs to the other candidates") + annotate("text", x=(med_cos_sim + 0.1) , y=800, label= med_cos_sim)
  
  print(plot)
  
  return(wilcox.test(fraction, mu= med_cos_sim , alternative="two.sided" ))
  #return(fraction_adj)
  
}
```

```{r}

cos_tfl_test<- cosine_sim_random_test(median(cos_median_tfl$cos_sim_med), length(cos_median_tfl), cos_GI1)
```

```{r}
cos_limma_test<- cosine_sim_random_test(median(cos_median_limma$cos_sim_med), length(cos_median_limma), cos_GI1)
```

```{r}
cos_deseq2_test<- cosine_sim_random_test(median(cos_median_deseq2$cos_sim_med), length(cos_median_deseq2), cos_GI1)
```

```{r}
remove_list <- ls()[grep("GI1|cos", ls())]
```

```{r}
rm( list=remove_list)
gc()
```


heatmap for lung 
```{r}
se_assay_A549<- se_assay[, grep( "A549", colnames(se_assay)) ]
```

```{r}
cos_A549  <- cosine(as.matrix(se_assay_A549))
#<- test
```

```{r}
saveRDS(cos_A549, "~/output/lung_cancer/220616_cosine_sim_LINCS_A549.rds")
```


```{r eval=FALSE}
#comment.   
lincs_limma_gbm_res_A549_fda_approved<- read_csv("~/output/lung_cancer/limma_res/220602_SR_LINCS_LUAD_limma_RES_tau.csv")
lincs_limma_gbm_res_A549_fda_approved<- lincs_limma_gbm_res_A549_fda_approved[lincs_limma_gbm_res_A549_fda_approved$fda_approved == TRUE, ]
lincs_deseq_gbm_res_A549_fda_approved<- read_csv("~/output/lung_cancer/deseq2_res/220602_SR_LINCS_LUAD_DESEQ2_RES_tau.csv")
lincs_deseq_gbm_res_A549_fda_approved<- lincs_deseq_gbm_res_A549_fda_approved[lincs_deseq_gbm_res_A549_fda_approved$fda_approved == TRUE, ]
lincs_tfl_gbm_res_A549_fda_approved<- read_csv("~/output/lung_cancer/TFL_res/220602_SR_LINCS_LUAD_TFL_RES_tau.csv")
lincs_tfl_gbm_res_A549_fda_approved<- lincs_tfl_gbm_res_A549_fda_approved[lincs_tfl_gbm_res_A549_fda_approved$fda_approved==TRUE,]
```

```{r}
all_drugs<- c(lincs_limma_gbm_res_A549_fda_approved$pert, lincs_deseq_gbm_res_A549_fda_approved$pert, lincs_tfl_gbm_res_A549_fda_approved$pert)
```

```{r}
#grep(all_drugs, colnames(cos_A549))
assay_A549_cand<- se_assay_A549[, grep(paste(all_drugs,collapse="|"), colnames(se_assay_A549))]
```



```{r}
cos_A549_cand <- cosine(as.matrix(assay_A549_cand))
```



```{r}

tfl_drugs <- grepl(paste(lincs_tfl_gbm_res_A549_fda_approved$pert,collapse="|"), colnames(cos_A549_cand)) 
deseq2_drugs<- grepl(paste(lincs_deseq_gbm_res_A549_fda_approved$pert,collapse="|"), colnames(cos_A549_cand)) 
limma_drugs <- grepl(paste(lincs_limma_gbm_res_A549_fda_approved$pert,collapse="|"), colnames(cos_A549_cand)) 

row_ha = HeatmapAnnotation(Transfer_Learning=tfl_drugs,DESeq2= deseq2_drugs, limma= limma_drugs , col = list(Transfer_Learning = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),DESeq2 = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),limma = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF") ))
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "black", "yellow"))
rownames(cos_A549_cand) <- str_split(rownames(cos_A549_cand), "__A549", simplify = TRUE)[,1]
Heatmap(cos_A549_cand, nam= "Cosine Similarity of LINCS Profiles",  col = col_fun, show_column_names = FALSE, show_row_names =FALSE, top_annotation = row_ha, 
        clustering_distance_rows= "euclidean",
        clustering_distance_columns=  "euclidean",
        clustering_method_rows = "ward.D2" ,
        clustering_method_columns="ward.D2")
```

for each of the drug from the same method what is the average cosine siamilirty
```{r}
cos_A549_cand_tfl <- cos_A549[grep(paste(lincs_tfl_gbm_res_A549_fda_approved$pert,collapse="|"), colnames(cos_A549)), grep(paste(lincs_tfl_gbm_res_A549_fda_approved$pert,collapse="|"), colnames(cos_A549))]
```

```{r}
cos_A549_cand_tfl[cos_A549_cand_tfl == 1] <- NA
```

```{r}
rowMedians(cos_A549_cand_tfl, na.rm=TRUE)
boxplot(cos_A549_cand_tfl)
```
```{r}
cos_A549_cand_tfl<- as.data.frame(cos_A549_cand_tfl)
cos_A549_cand_tfl$names <-  str_split(rownames(cos_A549_cand_tfl), "__A549", simplify = TRUE)[,1]
cos_A549_cand_tfl_long <- cos_A549_cand_tfl %>% pivot_longer(!names, names_to = "other")
```

```{r}
p <- ggplot(cos_A549_cand_tfl_long, aes(x=names, y=value)) + 
  geom_boxplot()+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
  labs(title="Cosine Similarity of Transfer Learning GBM candidates",x="Drug Candidates", y = "Cosine Similarity")
p
```

```{r}
plot_cosine_similarity_boxplot<- function(cos_sim_matrix, drug_list, pattern_cell= "__A549",  title , method){
  cos_A549_cand_tfl <- cos_sim_matrix[grep(paste(drug_list,collapse="|"), colnames(cos_sim_matrix)), grep(paste(drug_list,collapse="|"), colnames(cos_sim_matrix))]
  
  cos_A549_cand_tfl[cos_A549_cand_tfl == 1] <- NA
  cos_median <- rowMedians(cos_A549_cand_tfl, na.rm=TRUE)
  meth<- rep(method, nrow(cos_A549_cand_tfl))
  
  cos_median<- data.frame(drug = str_split(rownames(cos_A549_cand_tfl), pattern_cell, simplify = TRUE)[,1], cos_sim_med= cos_median, method= meth)
  
  cos_A549_cand_tfl<- as.data.frame(cos_A549_cand_tfl)
  cos_A549_cand_tfl$names <-  str_split(rownames(cos_A549_cand_tfl), pattern_cell, simplify = TRUE)[,1]
  cos_A549_cand_tfl_long <- cos_A549_cand_tfl %>% pivot_longer(!names, names_to = "other")
  p <- ggplot(cos_A549_cand_tfl_long, aes(x=names, y=value)) + 
    geom_boxplot()+ geom_dotplot(binaxis='y', stackdir='center', dotsize=0.5) +
    labs(title=title,x="Drug Candidates", y = "Cosine Similarity")
  print(p)
  
  return(cos_median)
}
```

```{r}
cos_median_tfl<- plot_cosine_similarity_boxplot(cos_A549,lincs_tfl_gbm_res_A549_fda_approved$pert, "__A549", "Cosine Similarity of Transfer Learning GBM candidates" , "Transfer Learning")
```
```{r}
cos_median_limma<- plot_cosine_similarity_boxplot(cos_A549,lincs_limma_gbm_res_A549_fda_approved$pert, "__A549", "Cosine Similarity of limma GBM candidates", "limma" )
```

```{r}
cos_median_deseq2<- plot_cosine_similarity_boxplot(cos_A549,lincs_deseq_gbm_res_A549_fda_approved$pert, "__A549", "Cosine Similarity of DESeq2 GBM candidates", "DESeq2" )
```

```{r}
cos_median_gbm <- rbind(cos_median_deseq2, cos_median_limma, cos_median_tfl)
```

```{r}
ggplot(cos_median_gbm, aes(x=method, y=cos_sim_med, fill=method)) + 
  geom_violin()+ geom_point() +labs(title="Median Cosine Similarity of GBM candidates",x="Method", y = "Median Cosine Similarity", color= "Method") + scale_colour_manual(values =  c("#440154FF", "#31688EFF", "#FDE725FF"), aesthetics = c("colour", "fill"))
```
```{r}
kruskal.test(cos_sim_med ~ method, data = cos_median_gbm)
```

```{r}
pairwise.wilcox.test(cos_median_gbm$cos_sim_med, cos_median_gbm$method,
                     p.adjust.method = "BH")
```


```{r}
cosine_sim_random_test<- function(med_cos_sim, num_drug, cos_sim_matrix, seed= 101){
  
  gbm_lincs_drugs <- as.vector(str_split(colnames(cos_sim_matrix),  "__", simplify = TRUE)[,1])
  
  gbm_lincs_drugs<- gbm_lincs_drugs[FDA_APPROVAL_CHECK(gbm_lincs_drugs)]
  number_sig_drugs<- num_drug
  # do that again 10,000 times 
  set.seed(seed)
  fraction<-c()
  for (i in 1:10000){
    
    test <- gbm_lincs_drugs[ sample(x = 1:length(gbm_lincs_drugs),size = number_sig_drugs ,replace = FALSE)]
    
    cos_A549_cand_tfl <- cos_sim_matrix[grep(paste(test,collapse="|"), colnames(cos_sim_matrix)), grep(paste(test,collapse="|"), colnames(cos_sim_matrix))]
    
    cos_A549_cand_tfl[cos_A549_cand_tfl == 1] <- NA
    cos_median <- rowMedians(cos_A549_cand_tfl, na.rm=TRUE)
    
    
    fraction[i] <- median( cos_median )
    
  }
  
  plot_data<- as.data.frame(cbind(1:10000, fraction))
  
  plot <- ggplot(plot_data, aes(x=fraction)) + geom_histogram() + geom_vline(xintercept = med_cos_sim)+ labs( x = "median cosine similarity of all drugs to the other candidates") + annotate("text", x=(med_cos_sim + 0.1) , y=800, label= med_cos_sim)
  
  print(plot)
  
  return(wilcox.test(fraction, mu= med_cos_sim , alternative="two.sided" ))
  #return(fraction_adj)
  
}
```

```{r}

cos_tfl_test<- cosine_sim_random_test(median(cos_median_tfl$cos_sim_med), length(cos_median_tfl), cos_A549)
```

```{r}
cos_limma_test<- cosine_sim_random_test(median(cos_median_limma$cos_sim_med), length(cos_median_limma), cos_A549)
```

```{r}
cos_deseq2_test<- cosine_sim_random_test(median(cos_median_deseq2$cos_sim_med), length(cos_median_deseq2), cos_A549)
```



heatmap for liver 

```{r}
lincs_deseq_gbm_res_HEPG2_fda_approved <- read_csv("~/output/liver_cancer/deseq2_res/220602_SR_LINCS_LIHC_DESEQ2_RES_tau.csv")
lincs_deseq_gbm_res_HEPG2_fda_approved <- lincs_deseq_gbm_res_HEPG2_fda_approved [lincs_deseq_gbm_res_HEPG2_fda_approved$fda_approved== TRUE,]

lincs_limma_gbm_res_HEPG2_fda_approved <- read_csv("~/output/liver_cancer/limma_res/220602_SR_LINCS_LIHC_limma_RES_tau.csv")

lincs_limma_gbm_res_HEPG2_fda_approved<- lincs_limma_gbm_res_HEPG2_fda_approved[ lincs_limma_gbm_res_HEPG2_fda_approved$fda_approved== TRUE,]

lincs_tfl_gbm_res_HEPG2_fda_approved <- read_csv("~/output/liver_cancer/TFL_res/220602_SR_LINCS_LIHC_TFL_RES_tau.csv")
lincs_tfl_gbm_res_HEPG2_fda_approved<- lincs_tfl_gbm_res_HEPG2_fda_approved[ lincs_tfl_gbm_res_HEPG2_fda_approved$fda_approved == TRUE, ]
```



```{r}
all_drugs<- c(lincs_limma_gbm_res_HEPG2_fda_approved$pert, lincs_deseq_gbm_res_HEPG2_fda_approved$pert, lincs_tfl_gbm_res_HEPG2_fda_approved$pert)
```

```{r}
se_assay_HEPG2<- se_assay[ , grep( "HEPG2", colnames(se_assay)) ]
```

```{r}
cos_HEPG2<- cosine(as.matrix(se_assay_HEPG2))
```

```{r}
saveRDS(cos_HEPG2, "~/output/liver_cancer/220616_cosine_sim_LINCS_HEPG2.rds")
```


```{r}
#grep(all_drugs, colnames(cos_HEPG2))
assay_HEPG2_cand<- se_assay_HEPG2[, grep(paste(all_drugs,collapse="|"), colnames(se_assay_HEPG2))]
```



```{r}
cos_HEPG2_cand <- cosine(as.matrix(assay_HEPG2_cand))
```



```{r}

tfl_drugs <- grepl(paste(lincs_tfl_gbm_res_HEPG2_fda_approved$pert,collapse="|"), colnames(cos_HEPG2_cand)) 
deseq2_drugs<- grepl(paste(lincs_deseq_gbm_res_HEPG2_fda_approved$pert,collapse="|"), colnames(cos_HEPG2_cand)) 
limma_drugs <- grepl(paste(lincs_limma_gbm_res_HEPG2_fda_approved$pert,collapse="|"), colnames(cos_HEPG2_cand)) 

row_ha = HeatmapAnnotation(Transfer_Learning=tfl_drugs,DESeq2= deseq2_drugs, limma= limma_drugs , col = list(Transfer_Learning = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),DESeq2 = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),limma = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF") ))
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "black", "yellow"))
rownames(cos_HEPG2_cand) <- str_split(rownames(cos_HEPG2_cand), "__HEPG2", simplify = TRUE)[,1]
Heatmap(cos_HEPG2_cand, nam= "Cosine Similarity of LINCS Profiles",  col = col_fun, show_column_names = FALSE,  top_annotation = row_ha, 
        clustering_distance_rows= "euclidean",
        clustering_distance_columns=  "euclidean",
        clustering_method_rows = "ward.D2" ,
        clustering_method_columns="ward.D2")
```
```{r}
cos_median_tfl<- plot_cosine_similarity_boxplot(cos_HEPG2_cand,lincs_tfl_gbm_res_HEPG2_fda_approved$pert, "__HEPG2", "Cosine Similarity of Transfer Learning GBM candidates" , "Transfer Learning")
```
```{r}
cos_median_limma<- plot_cosine_similarity_boxplot(cos_HEPG2_cand,lincs_limma_gbm_res_HEPG2_fda_approved$pert, "__HEPG2", "Cosine Similarity of limma GBM candidates", "limma" )
```

```{r}
cos_median_deseq2<- plot_cosine_similarity_boxplot(cos_HEPG2_cand,lincs_deseq_gbm_res_HEPG2_fda_approved$pert, "__HEPG2", "Cosine Similarity of DESeq2 GBM candidates", "DESeq2" )
```

```{r}
cos_median_gbm <- rbind(cos_median_deseq2, cos_median_limma, cos_median_tfl)
```

```{r}

cos_tfl_test<- cosine_sim_random_test(median(cos_median_tfl$cos_sim_med), length(cos_median_tfl), cos_HEPG2)
```

```{r}
cos_limma_test<- cosine_sim_random_test(median(cos_median_limma$cos_sim_med), length(cos_median_limma), cos_HEPG2)
```

```{r}
cos_deseq2_test<- cosine_sim_random_test(median(cos_median_deseq2$cos_sim_med), length(cos_median_deseq2), cos_HEPG2)
```


```{r}
ggplot(cos_median_gbm, aes(x=method, y=cos_sim_med, fill=method)) + 
  geom_violin()+ geom_point() +labs(title="Median Cosine Similarity of GBM candidates",x="Method", y = "Median Cosine Similarity", color= "Method") + scale_colour_manual(values =  c("#440154FF", "#31688EFF", "#FDE725FF"), aesthetics = c("colour", "fill"))
```

```{r}
kruskal.test(cos_sim_med ~ method, data = cos_median_gbm)
```

```{r}
pairwise.wilcox.test(cos_median_gbm$cos_sim_med, cos_median_gbm$method,
                     p.adjust.method = "BH")
```


heatmap for pancreatic 


```{r}

lincs_deseq_gbm_res_YAPC_fda_approved <- read_csv("~/output/pancreas_cancer/deseq2_res/220602_SR_LINCS_PAAD_DESEQ2_RES_tau.csv")
lincs_deseq_gbm_res_YAPC_fda_approved<- lincs_deseq_gbm_res_YAPC_fda_approved[lincs_deseq_gbm_res_YAPC_fda_approved$fda_approved == TRUE,]

lincs_limma_gbm_res_YAPC_fda_approved <- read_csv("~/output/pancreas_cancer/limma_res/220602_SR_LINCS_PAAD_limma_RES_tau.csv")
lincs_limma_gbm_res_YAPC_fda_approved<- lincs_limma_gbm_res_YAPC_fda_approved[lincs_limma_gbm_res_YAPC_fda_approved$fda_approved ==TRUE, ]

lincs_tfl_gbm_res_YAPC_fda_approved <- read_csv("~/output/pancreas_cancer/TFL_res/220602_SR_LINCS_PAAD_TFL_RES_tau.csv")
lincs_tfl_gbm_res_YAPC_fda_approved<- lincs_tfl_gbm_res_YAPC_fda_approved[ lincs_tfl_gbm_res_YAPC_fda_approved$fda_approved ==TRUE,]
```



```{r}
all_drugs<- c(lincs_limma_gbm_res_YAPC_fda_approved$pert, lincs_deseq_gbm_res_YAPC_fda_approved$pert, lincs_tfl_gbm_res_YAPC_fda_approved$pert)
```

```{r}
se_assay_YAPC<- se_assay[ , grep( "YAPC", colnames(se_assay)) ]
```

```{r}
cos_YAPC<- cosine(as.matrix(se_assay_YAPC))
saveRDS(cos_YAPC, "~/output/pancreas_cancer/220616_cosine_sim_LINCS_YAPC.rds")
```


```{r}
#grep(all_drugs, colnames(cos_YAPC))
assay_YAPC_cand<- se_assay_YAPC[, grep(paste(all_drugs,collapse="|"), colnames(se_assay_YAPC))]
```



```{r}
cos_YAPC_cand <- cosine(as.matrix(assay_YAPC_cand))
```



```{r}

tfl_drugs <- grepl(paste(lincs_tfl_gbm_res_YAPC_fda_approved$pert,collapse="|"), colnames(cos_YAPC_cand)) 
deseq2_drugs<- grepl(paste(lincs_deseq_gbm_res_YAPC_fda_approved$pert,collapse="|"), colnames(cos_YAPC_cand)) 
limma_drugs <- grepl(paste(lincs_limma_gbm_res_YAPC_fda_approved$pert,collapse="|"), colnames(cos_YAPC_cand)) 

row_ha = HeatmapAnnotation(Transfer_Learning=tfl_drugs,DESeq2= deseq2_drugs, limma= limma_drugs , col = list(Transfer_Learning = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),DESeq2 = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),limma = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF") ))
col_fun = colorRamp2(c(-1, 0, 1), c("blue", "black", "yellow"))
rownames(cos_YAPC_cand) <- str_split(rownames(cos_YAPC_cand), "__YAPC", simplify = TRUE)[,1]
Heatmap(cos_YAPC_cand, nam= "Cosine Similarity of LINCS Profiles",  col = col_fun, show_column_names = FALSE,show_row_names= FALSE, top_annotation = row_ha, 
        clustering_distance_rows= "euclidean",
        clustering_distance_columns=  "euclidean",
        clustering_method_rows = "ward.D2" ,
        clustering_method_columns="ward.D2")
```

```{r}
cos_median_tfl<- plot_cosine_similarity_boxplot(cos_YAPC_cand,lincs_tfl_gbm_res_YAPC_fda_approved$pert, "__YAPC", "Cosine Similarity of Transfer Learning GBM candidates" , "Transfer Learning")
```
```{r}
cos_median_limma<- plot_cosine_similarity_boxplot(cos_YAPC_cand,lincs_limma_gbm_res_YAPC_fda_approved$pert, "__YAPC", "Cosine Similarity of limma GBM candidates", "limma" )
```

```{r}
cos_median_deseq2<- plot_cosine_similarity_boxplot(cos_YAPC_cand,lincs_deseq_gbm_res_YAPC_fda_approved$pert, "__YAPC", "Cosine Similarity of DESeq2 GBM candidates", "DESeq2" )
```
```{r}

cos_tfl_test<- cosine_sim_random_test(median(cos_median_tfl$cos_sim_med), length(cos_median_tfl), cos_YAPC)
```

```{r}
cos_limma_test<- cosine_sim_random_test(median(cos_median_limma$cos_sim_med), length(cos_median_limma), cos_YAPC)
```

```{r}
cos_deseq2_test<- cosine_sim_random_test(median(cos_median_deseq2$cos_sim_med), length(cos_median_deseq2), cos_YAPC)
```

```{r}
cos_median_gbm <- rbind(cos_median_deseq2, cos_median_limma, cos_median_tfl)
```

```{r}
ggplot(cos_median_gbm, aes(x=method, y=cos_sim_med, fill=method)) + 
  geom_violin()+ geom_point() +labs(title="Median Cosine Similarity of GBM candidates",x="Method", y = "Median Cosine Similarity", color= "Method") + scale_colour_manual(values =  c("#440154FF", "#31688EFF", "#FDE725FF"), aesthetics = c("colour", "fill"))
```

```{r}
kruskal.test(cos_sim_med ~ method, data = cos_median_gbm)
```

```{r}
pairwise.wilcox.test(cos_median_gbm$cos_sim_med, cos_median_gbm$method,
                     p.adjust.method = "BH")
```


```{r}
tiagabine_1_test_v2 <- tiagabine_1_test[order(tiagabine_1_test$logfoldchange),]
#tiagabine_1_test_v2$cell_line <- as.factor(tiagabine_1_test_v2$cell_line)
tiagabine_1_test_v2$cell_line <- factor(tiagabine_1_test_v2$cell_line, levels= as.character(tiagabine_1_test_v2$cell_line))

ggplot(tiagabine_1_test_v2, aes(x=logfoldchange, y= cell_line))+ geom_point(stat = "identity") +theme_minimal()
```

```{r}
ggplot(tiagabine_1_test_v2, aes(x=logfoldchange)) + geom_density(alpha=.2, fill="#440154FF") +theme_minimal()
```


```{r}
p2 <- ggplot(tiagabine_1_test_v2, aes(x=logfoldchange, y= cell_line))+ geom_point(stat = "identity") +theme_minimal() + xlab("PRISM Primary Screen log(fold change)") + ylab("GBM Cell Lines") + geom_vline(xintercept = 0.3)
p1 <- ggplot(tiagabine_1_test_v2, aes(x=logfoldchange)) + geom_density(alpha=.2, fill="#440154FF") +theme_minimal() +  theme(
  axis.title.x = element_blank(),
)+ geom_vline(xintercept = 0.3)
plot_grid(p1, p2, ncol = 1, align = "v")
```

```{r}
tiagabine_secondary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/tiagabine AUC (PRISM Repurposing Secondary Screen)  tiagabine log2 fold change Drug sensitivity.csv")
colnames(tiagabine_secondary)[2]<- "logfoldchange"
colnames(tiagabine_secondary)[3]<- "AUC"
colnames(tiagabine_secondary)[5]<- "cell_line"
```

```{r}
ggplot(tiagabine_secondary, aes(x=logfoldchange, y= AUC))+geom_point()+ geom_text(
  label=tiagabine_secondary$cell_line, nudge_x = 0.005, nudge_y = 0.005, 
  check_overlap = T) + xlab("PRISM Primary Screen log(fold change)") + ylab("PRISM Secondary Screen AUC") + geom_vline(xintercept = 0.3) +geom_hline(yintercept = 0.9) +theme_minimal()
```

```{r}
vemurafenib_primary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/vemurafenib (BRDBRD-K56343971-001-10-6) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
colnames(vemurafenib_primary)[2]<- "logfoldchange"
colnames(vemurafenib_primary)[4]<- "cell_line"
```

```{r}
vemurafenib_primary_v2 <- vemurafenib_primary[order(vemurafenib_primary$logfoldchange),]
#tiagabine_1_test_v2$cell_line <- as.factor(tiagabine_1_test_v2$cell_line)
vemurafenib_primary_v2$cell_line <- factor(vemurafenib_primary_v2$cell_line, levels= as.character(vemurafenib_primary_v2$cell_line))

```

```{r}
p2 <- ggplot(vemurafenib_primary_v2, aes(x=logfoldchange, y= cell_line))+ geom_point(stat = "identity") +theme_minimal() + xlab("PRISM Primary Screen log(fold change)") + ylab("GBM Cell Lines") + geom_vline(xintercept = 0.3)
p1 <- ggplot(vemurafenib_primary_v2, aes(x=logfoldchange)) + geom_density(alpha=.2, fill="#440154FF") +theme_minimal() +  theme(
  axis.title.x = element_blank(),
)+ geom_vline(xintercept = 0.3)
plot_grid(p1, p2, ncol = 1, align = "v")
```

```{r}
vemurafenib_secondary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/vemurafenib (BRDBRD-K56343971-001-14-8) AUC Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4 vs vemurafenib (BRDBRD-K56343971-001-10-6) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
```

```{r}
colnames(vemurafenib_secondary)[2]<- "logfoldchange"
colnames(vemurafenib_secondary)[3]<- "AUC"
colnames(vemurafenib_secondary)[5]<- "cell_line"
```

```{r}
ggplot(vemurafenib_secondary, aes(x=logfoldchange, y= AUC))+geom_point()+ geom_text(
  label=vemurafenib_secondary$cell_line, nudge_x = 0.005, nudge_y = 0.005, 
  check_overlap = T) + xlab("PRISM Primary Screen log(fold change)") + ylab("PRISM Secondary Screen AUC") + geom_vline(xintercept = 0.3) +geom_hline(yintercept = 0.9) +theme_minimal()
```


```{r}
saxagliptin_primary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/saxagliptin (BRDBRD-A81513827-001-03-6) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
```

```{r}
colnames(saxagliptin_primary)[2]<- "logfoldchange"
colnames(saxagliptin_primary)[4]<- "cell_line"
```

```{r}
saxagliptin_primary_v2 <- saxagliptin_primary[order(saxagliptin_primary$logfoldchange),]
#tiagabine_1_test_v2$cell_line <- as.factor(tiagabine_1_test_v2$cell_line)
saxagliptin_primary_v2$cell_line <- factor(saxagliptin_primary_v2$cell_line, levels= as.character(saxagliptin_primary_v2$cell_line))
```

```{r}
p2 <- ggplot(saxagliptin_primary_v2, aes(x=logfoldchange, y= cell_line))+ geom_point(stat = "identity") +theme_minimal() + xlab("PRISM Primary Screen log(fold change)") + ylab("GBM Cell Lines") + geom_vline(xintercept = 0.3)
p1 <- ggplot(saxagliptin_primary_v2, aes(x=logfoldchange)) + geom_density(alpha=.2, fill="#440154FF") +theme_minimal() +  theme(
  axis.title.x = element_blank(),
)+ geom_vline(xintercept = 0.3)
plot_grid(p1, p2, ncol = 1, align = "v")
```

```{r}
saxagliptin_secondary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/saxagliptin (BRDBRD-A81513827-001-03-6) AUC Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4 vs saxagliptin (BRDBRD-A81513827-001-03-6) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
```


```{r}
colnames(saxagliptin_secondary)[2]<- "logfoldchange"
colnames(saxagliptin_secondary)[3]<- "AUC"
colnames(saxagliptin_secondary)[5]<- "cell_line"
```

```{r}
ggplot(saxagliptin_secondary, aes(x=logfoldchange, y= AUC))+geom_point()+ geom_text(
  label=saxagliptin_secondary$cell_line, nudge_x = 0.005, nudge_y = 0.005, 
  check_overlap = T) + xlab("PRISM Primary Screen log(fold change)") + ylab("PRISM Secondary Screen AUC") + geom_vline(xintercept = 0.3) +geom_hline(yintercept = 0.9) +theme_minimal()
```
```{r}
spironolactone_primary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/spironolactone (BRDBRD-K90027355-001-13-3) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
```


```{r}
colnames(spironolactone_primary)[2]<- "logfoldchange"
colnames(spironolactone_primary)[4]<- "cell_line"
```

```{r}
spironolactone_primary_v2 <- spironolactone_primary[order(spironolactone_primary$logfoldchange),]
#tiagabine_1_test_v2$cell_line <- as.factor(tiagabine_1_test_v2$cell_line)
spironolactone_primary_v2$cell_line <- factor(spironolactone_primary_v2$cell_line, levels= as.character(spironolactone_primary_v2$cell_line))
```

```{r}
p2 <- ggplot(spironolactone_primary_v2, aes(x=logfoldchange, y= cell_line))+ geom_point(stat = "identity") +theme_minimal() + xlab("PRISM Primary Screen log(fold change)") + ylab("GBM Cell Lines") + geom_vline(xintercept = 0.3)
p1 <- ggplot(spironolactone_primary_v2, aes(x=logfoldchange)) + geom_density(alpha=.2, fill="#440154FF") +theme_minimal() +  theme(
  axis.title.x = element_blank(),
)+ geom_vline(xintercept = 0.3)
plot_grid(p1, p2, ncol = 1, align = "v")
```

```{r}
spironolactone_secondary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/spironolactone (BRDBRD-K90027355-001-13-3) AUC Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4 vs spironolactone (BRDBRD-K90027355-001-13-3) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
```


```{r}
colnames(spironolactone_secondary)[2]<- "logfoldchange"
colnames(spironolactone_secondary)[3]<- "AUC"
colnames(spironolactone_secondary)[5]<- "cell_line"
```

```{r}
ggplot(spironolactone_secondary, aes(x=logfoldchange, y= AUC))+geom_point()+ geom_text(
  label=spironolactone_secondary$cell_line, nudge_x = 0.005, nudge_y = 0.005, 
  check_overlap = T) + xlab("PRISM Primary Screen log(fold change)") + ylab("PRISM Secondary Screen AUC") + geom_vline(xintercept = 0.3) +geom_hline(yintercept = 0.9) +theme_minimal()
```


```{r}
icosapent_primary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/icosapent (BRDBRD-K47192521-001-05-4) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
```


```{r}
colnames(icosapent_primary )[2]<- "logfoldchange"
colnames(icosapent_primary )[4]<- "cell_line"
```

```{r}
icosapent_primary_v2 <- icosapent_primary [order(icosapent_primary$logfoldchange),]
#tiagabine_1_test_v2$cell_line <- as.factor(tiagabine_1_test_v2$cell_line)
icosapent_primary_v2$cell_line <- factor(icosapent_primary_v2$cell_line, levels= as.character(icosapent_primary_v2$cell_line))
```

```{r}
p2 <- ggplot(icosapent_primary_v2, aes(x=logfoldchange, y= cell_line))+ geom_point(stat = "identity") +theme_minimal() + xlab("PRISM Primary Screen log(fold change)") + ylab("GBM Cell Lines") + geom_vline(xintercept = 0.3)
p1 <- ggplot(icosapent_primary_v2, aes(x=logfoldchange)) + geom_density(alpha=.2, fill="#440154FF") +theme_minimal() +  theme(
  axis.title.x = element_blank(),
)+ geom_vline(xintercept = 0.3)
plot_grid(p1, p2, ncol = 1, align = "v")
```
```{r}
pamidronate_primary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/pamidronate (BRDBRD-K58513245-304-02-7) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
```


```{r}
colnames(pamidronate_primary )[2]<- "logfoldchange"
colnames(pamidronate_primary )[4]<- "cell_line"
```

```{r}
pamidronate_primary_v2 <- pamidronate_primary[order(pamidronate_primary$logfoldchange),]
#tiagabine_1_test_v2$cell_line <- as.factor(tiagabine_1_test_v2$cell_line)
pamidronate_primary_v2$cell_line <- factor(pamidronate_primary_v2$cell_line, levels= as.character(pamidronate_primary_v2$cell_line))
```

```{r}
p2 <- ggplot(pamidronate_primary_v2, aes(x=logfoldchange, y= cell_line))+ geom_point(stat = "identity") +theme_minimal() + xlab("PRISM Primary Screen log(fold change)") + ylab("GBM Cell Lines") + geom_vline(xintercept = 0.3)
p1 <- ggplot(pamidronate_primary_v2, aes(x=logfoldchange)) + geom_density(alpha=.2, fill="#440154FF") +theme_minimal() +  theme(
  axis.title.x = element_blank(),
)+ geom_vline(xintercept = 0.3)
plot_grid(p1, p2, ncol = 1, align = "v")
```

also need to plot tmz
```{r}
temozolomide_primary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/temozolomide (BRDBRD-K32107296-001-16-9) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
```


```{r}
colnames(temozolomide_primary )[2]<- "logfoldchange"
colnames(temozolomide_primary)[4]<- "cell_line"
```

```{r}
temozolomide_primary_v2 <- temozolomide_primary[order(temozolomide_primary$logfoldchange),]
#tiagabine_1_test_v2$cell_line <- as.factor(tiagabine_1_test_v2$cell_line)
temozolomide_primary_v2$cell_line <- factor(temozolomide_primary_v2$cell_line, levels= as.character(temozolomide_primary_v2$cell_line))
```

```{r}
p2 <- ggplot(temozolomide_primary_v2, aes(x=logfoldchange, y= cell_line))+ geom_point(stat = "identity") +theme_minimal() + xlab("PRISM Primary Screen log(fold change)") + ylab("GBM Cell Lines") + geom_vline(xintercept = 0.3)
p1 <- ggplot(temozolomide_primary_v2, aes(x=logfoldchange)) + geom_density(alpha=.2, fill="#440154FF") +theme_minimal() +  theme(
  axis.title.x = element_blank(),
)+ geom_vline(xintercept = 0.3)
plot_grid(p1, p2, ncol = 1, align = "v")
```

```{r}
temozolomide_secondary <- read_csv("~/output/TF_L_GBM/220503_PRISM_DEPMAP_candidate_data/temozolomide (BRDBRD-K32107296-001-16-9) AUC Drug sensitivity AUC (PRISM Repurposing Secondary Screen) 19Q4 vs temozolomide (BRDBRD-K32107296-001-16-9) log2 fold change Drug sensitivity (PRISM Repurposing Primary Screen) 19Q4.csv")
```

```{r}
colnames(temozolomide_secondary)[2]<- "logfoldchange"
colnames(temozolomide_secondary)[3]<- "AUC"
colnames(temozolomide_secondary)[5]<- "cell_line"
```

```{r}
ggplot(temozolomide_secondary, aes(x=logfoldchange, y= AUC))+geom_point()+ geom_text(
  label=temozolomide_secondary$cell_line, nudge_x = 0.00008, nudge_y = 0.00008, 
  check_overlap = T) + xlab("PRISM Primary Screen log(fold change)") + ylab("PRISM Secondary Screen AUC") + geom_vline(xintercept = 0.3) +theme_minimal()
```



