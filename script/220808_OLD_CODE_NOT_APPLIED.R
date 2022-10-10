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


PRISM primary result differences between the three methods

```{r}
library(readr)
primary_screen <- read_csv("~/data/DEPMAP_PRISM_220228/primary-screen-replicate-collapsed-logfold-change.csv")
sample_info <- read_csv("~/data/DEPMAP_PRISM_220228/sample_info.csv")
primary_screen_treatment_info <- read_csv("~/data/DEPMAP_PRISM_220228/primary-screen-replicate-collapsed-treatment-info.csv")
```
different order

```{r}
names_vector<- c()
for ( i in 2:ncol(primary_screen)){
  names_vector[i] <- primary_screen_treatment_info$name[ grep(colnames(primary_screen)[i], primary_screen_treatment_info$column_name)]
}
```

```{r}
good_id <- c()
for (i in 1:nrow(primary_screen)){
  good_id[i]<- nchar(primary_screen[i,1]) == 10 
}

```

```{r}
table(good_id)
```
```{r}
primary_screen_v2 <- primary_screen[good_id,]
```


```{r}
cell_vector<- c()
for ( i in 1:nrow(primary_screen_v2)){
  cell_vector[i] <-sample_info$stripped_cell_line_name[ grep(primary_screen_v2[i,1], sample_info$DepMap_ID)]
}
```

```{r}
primary_screen_v3<-  as.matrix(primary_screen_v2[,-1])
rownames(primary_screen_v3)<- cell_vector
```


```{r}
colnames(primary_screen_v3) <- names_vector[-1]
```

```{r}
test2<- as.data.frame.table(primary_screen_v3, responseName = "log_fold_change")
```


filter by gbm cell lines

```{r}
#table(sample_info$Subtype)
gbm_sample_info <- sample_info[sample_info$Subtype == "Glioblastoma",]
```


```{r}
primary_gbm_res <- test2[test2$Var1 %in% gbm_sample_info$stripped_cell_line_name,]
```


deseq2 only 
```{r}
primary_gbm_deseq2<- primary_gbm_res[primary_gbm_res$Var2 %in% c("suprofen", "trandolapril", "felbamate", "pralidoxime", "abiraterone", "ethotoin", "floxuridine", "vardenafil", "moxifloxacin", "nimodipine", "diltiazem"),]
```

```{r}
unique(primary_gbm_deseq2$Var2)
```

```{r}
primary_gbm_deseq2$method <- rep("DESeq2", nrow(primary_gbm_deseq2))
```

limma- only 
```{r}
primary_gbm_limma<- primary_gbm_res[primary_gbm_res$Var2 %in% c("lonafarnib", "dabrafenib", "apremilast", "rivastigmine", "crizotinib", "ixazomib", "imatinib", "sulfasalazine", "maraviroc", "amiodarone"),]
```

```{r}
unique(primary_gbm_limma$Var2)
```

```{r}
primary_gbm_limma$method <- rep("limma", nrow(primary_gbm_limma))
```

tf-only 
```{r}
primary_gbm_tfl<- primary_gbm_res[primary_gbm_res$Var2 %in% c("fludarabine", "tiagabine", "vemurafenib", "dasatinib", "nifedipine", "spironolactone", "lapatinib", "pentoxifylline", "erlotinib"),]
```
nifedipine not in prism
```{r}
unique(primary_gbm_tfl$Var2)
```

```{r}
primary_gbm_tfl$method <- rep("Transfer Learning", nrow(primary_gbm_tfl))
```

deseq2/tf 
```{r}
primary_gbm_deseq2_tfl<- primary_gbm_res[primary_gbm_res$Var2 %in% c("diflunisal", "rucaparib"),]
```
nifedipine not in prism
```{r}
unique(primary_gbm_deseq2_tfl$Var2)
```

```{r}
primary_gbm_deseq2_tfl$method <- rep("Transfer Learning and DESeq2", nrow(primary_gbm_deseq2_tfl))
```


limma/tf
```{r}
primary_gbm_limma_tfl<- primary_gbm_res[primary_gbm_res$Var2 %in% c("vorinostat", "thioridazine", "saxagliptin", "icosapent"),]
```
nifedipine not in prism
```{r}
unique(primary_gbm_limma_tfl$Var2)
```

```{r}
primary_gbm_limma_tfl$method <- rep("Transfer Learning and limma", nrow(primary_gbm_limma_tfl))
```

all 
```{r}
primary_gbm_all<- primary_gbm_res[primary_gbm_res$Var2 %in% c("pamidronate"),]
```
nifedipine not in prism
```{r}
unique(primary_gbm_all$Var2)
```

```{r}
primary_gbm_all$method <- rep("ALL", nrow(primary_gbm_all))
```

TMZ
```{r}
primary_gbm_tmz<- primary_gbm_res[primary_gbm_res$Var2 %in% c("temozolomide"),]
```
nifedipine not in prism
```{r}
unique(primary_gbm_tmz$Var2)
```

```{r}
primary_gbm_tmz$method <- rep("TMZ Treatment", nrow(primary_gbm_tmz))
```

```{r}
primary_gbm_candidates <- rbind(primary_gbm_deseq2, primary_gbm_limma, primary_gbm_tfl, primary_gbm_deseq2_tfl, primary_gbm_limma_tfl, primary_gbm_all, primary_gbm_tmz)
```

```{r}
library(ggplot2)
library(viridis)
```


```{r}
primary_gbm_candidates$Var2 <- factor(primary_gbm_candidates$Var2, levels=unique(primary_gbm_candidates$Var2))
ggplot(primary_gbm_candidates, aes(x=log_fold_change, y= Var2))+ geom_violin(aes(fill=method)) +geom_vline(xintercept = -1) +scale_fill_viridis_d()+ theme_minimal()
```

```{r}
primary_gbm_candidates_v2 <- rbind(primary_gbm_deseq2, primary_gbm_limma, primary_gbm_tfl,  primary_gbm_tmz)
```
```{r}
primary_gbm_candidates_v2$Var2 <- factor(primary_gbm_candidates_v2$Var2, levels=unique(primary_gbm_candidates_v2$Var2))
ggplot(primary_gbm_candidates_v2, aes(x=log_fold_change, y= Var2))+ geom_violin(aes(fill=method)) +theme_minimal()
```

look to see how many gbm cell lines are senstive
```{r}
primary_gbm_candidates$Sensitive <- primary_gbm_candidates$log_fold_change < 0.3
```

```{r}
table(primary_gbm_candidates$Sensitive)
```
```{r}
drug_fraction <- c()
drug_list <- as.vector(unique(primary_gbm_candidates$Var2))
for (i in 1:length(drug_list)){
  cells <- primary_gbm_candidates$Sensitive[primary_gbm_candidates$Var2 %in% drug_list[i]]
  cells<- cells[!is.na(cells)]
  #x<- table(cells)
  drug_fraction[i] <- length(cells[cells== TRUE])/ length(cells)
}
```

```{r}
hist(drug_fraction)
```

```{r}
drug_senstive_precentage <- data.frame(drug_fraction, drug_list)
```

```{r}
ggplot(drug_senstive_precentage, aes( drug_fraction, drug_list ))+ geom_point()
```

```{r}
deseq2<- c("suprofen", "trandolapril", "felbamate", "pralidoxime", "abiraterone", "ethotoin", "floxuridine", "vardenafil", "moxifloxacin", "nimodipine", "diltiazem")
limma<- c("lonafarnib", "dabrafenib", "apremilast", "rivastigmine", "crizotinib", "ixazomib", "imatinib", "sulfasalazine", "maraviroc", "amiodarone")
transfer_learning<- c("fludarabine", "tiagabine", "vemurafenib", "dasatinib", "nifedipine", "spironolactone", "lapatinib", "pentoxifylline", "erlotinib")
limma_transfer_learning<- c("vorinostat", "thioridazine", "saxagliptin", "icosapent")
deseq2_transfer_learning <- c("diflunisal", "rucaparib")

```

```{r}
drug_senstive_precentage$method<- ifelse(drug_senstive_precentage$drug_list %in% deseq2, "DESeq2", ifelse(
  drug_senstive_precentage$drug_list %in% limma, "limma", ifelse(
    drug_senstive_precentage$drug_list %in% transfer_learning, "Transfer Learning", ifelse(
      drug_senstive_precentage$drug_list %in% limma_transfer_learning, "limma and Transfer Learning", ifelse(
        drug_senstive_precentage$drug_list %in% deseq2_transfer_learning, "DESeq2 and Transfer Learning", ifelse(
          drug_senstive_precentage$drug_list == "pamidronate", "ALL", "temozolomide"
        )
      )
    )
  )
) )
```


```{r}
drug_senstive_precentage$drug_list<- factor(drug_senstive_precentage$drug_list,levels= drug_senstive_precentage$drug_list[order(drug_senstive_precentage$drug_fraction)])
ggplot(drug_senstive_precentage, aes( drug_fraction, drug_list, fill=method ))+ geom_bar(stat="identity")
```

```{r}
drug_senstive_precentage_v2 <- drug_senstive_precentage[! drug_senstive_precentage$drug_list == "temozolomide",]
ggplot(drug_senstive_precentage_v2, aes( drug_fraction, drug_list, fill=method ))+  geom_bar(stat="identity", color = "black")+ scale_fill_viridis_d()+ geom_vline(xintercept = 0.75)+  xlab("PRISM Primary Screen log(fold change)") + ylab("Drug Candidates (All Methods)")
```

```{r}
primary_gbm_res$Sensitive <- primary_gbm_res$log_fold_change < 0.3
drug_fraction <- c()
drug_list <- as.vector(unique(primary_gbm_res$Var2))
for (i in 1:length(drug_list)){
  cells <- primary_gbm_res$Sensitive[primary_gbm_res$Var2 %in% drug_list[i]]
  cells<- cells[!is.na(cells)]
  if( length(cells) < 12){
    print(drug_list[i])
  }
  #x<- table(cells)
  drug_fraction[i] <- length(cells[cells== TRUE])/ length(cells)
}
drug_senstive_precentage_all_drugs <- data.frame(drug_fraction, drug_list)
```

```{r}
table(drug_senstive_precentage_all_drugs$drug_fraction >= 0.75)
```

```{r}
drug_senstive_precentage_all_drugs$tmz <-drug_senstive_precentage_all_drugs$drug_fraction >= 0.75
```

```{r}
drug_senstive_precentage_all_drugs<- drug_senstive_precentage_all_drugs[!drug_senstive_precentage_all_drugs$drug_list == "temozolomide",]
```



```{r}
drug_senstive_precentage_all_drugs_tfl<- drug_senstive_precentage_all_drugs[!drug_senstive_precentage_all_drugs$drug_list %in% c("fludarabine", "tiagabine", "vemurafenib", "dasatinib", "nifedipine", "spironolactone", "lapatinib", "pentoxifylline", "erlotinib", "vorinostat", "thioridazine", "saxagliptin", "icosapent","diflunisal", "rucaparib", "pamidronate"),]
number_of_drugs<- c()
for ( i in 1:10000){
  test3 <- table(drug_senstive_precentage_all_drugs_tfl$tmz[sample(1:nrow(drug_senstive_precentage_all_drugs_tfl), 15)])
  number_of_drugs[i]<- test3[names(test3) == TRUE]
}
```

```{r}
hist(number_of_drugs)
```

```{r}
t.test(mu=11, number_of_drugs)
```

```{r}

drug_senstive_precentage_all_drugs_limma <- drug_senstive_precentage_all_drugs[drug_senstive_precentage_all_drugs$drug_list %in%c("lonafarnib", "dabrafenib", "apremilast", "rivastigmine", "crizotinib", "ixazomib", "imatinib", "sulfasalazine", "maraviroc", "amiodarone", "vorinostat", "thioridazine", "saxagliptin", "icosapent", "pamidronate"),]

drug_senstive_precentage_all_drugs_deseq2 <- drug_senstive_precentage_all_drugs[drug_senstive_precentage_all_drugs$drug_list %in%c("suprofen", "trandolapril", "felbamate", "pralidoxime", "abiraterone", "ethotoin", "floxuridine", "vardenafil", "moxifloxacin", "nimodipine", "diltiazem", "diflunisal", "rucaparib", "pamidronate"),]

```


```{r}
number_of_drugs_v2 <- c()
for ( i in 1:10000){
  test3 <- table(drug_senstive_precentage_all_drugs_tfl$tmz[sample(1:nrow(drug_senstive_precentage_all_drugs_tfl), 14)])
  number_of_drugs_v2[i]<- test3[names(test3) == TRUE]
}
```

```{r}
hist(number_of_drugs_v2)
```

```{r}
t.test(mu=9, number_of_drugs_v2, alternative = "less")
```

```{r}
hist(drug_senstive_precentage_all_drugs$drug_fraction)
```



```{r}
drug_senstive_precentage$bi_method<- grepl("Transfer Learning",drug_senstive_precentage$method)
```

```{r}
ggplot(drug_senstive_precentage, aes( drug_fraction, drug_list, fill=bi_method ))+ geom_bar(stat="identity") +geom_vline(xintercept = 0.75)
```

```{r}
drug_senstive_precentage$bi_method<- grepl("DES",drug_senstive_precentage$method)
```

```{r}
ggplot(drug_senstive_precentage, aes( drug_fraction, drug_list, fill=bi_method ))+ geom_bar(stat="identity")
```

```{r}
drug_senstive_precentage$bi_method<- grepl("limma",drug_senstive_precentage$method)
```

```{r}
ggplot(drug_senstive_precentage, aes( drug_fraction, drug_list, fill=bi_method ))+ geom_bar(stat="identity") 
```

```{r}
drug_senstive_precentage$drug_list
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

#tfl_Prism
```{r}
primary_paad_res_sub <-primary_paad_res[ primary_paad_res$Var2 %in% lincs_deseq_gbm_res_GI1C_fda_approved$pert,]
```


```{r eval=FALSE}
write.csv(primary_paad_res_sub , file= "~/output/TF_L_GBM/220602_prism_tfl_candidates.csv" )
```

#deseq2_Prism

```{r}
primary_paad_res_sub <-primary_paad_res[ primary_paad_res$Var2 %in% lincs_deseq_gbm_res_GI1_fda_approved$pert,]
```


```{r eval=FALSE}
write.csv(primary_paad_res_sub , file= "~/output/deseq2_gbm/220602_prism_deseq2_candidates.csv" )
```

#limma_Prism

```{r}
primary_paad_res_sub <-primary_paad_res[ primary_paad_res$Var2 %in% lincs_limma_gbm_res_GI1_fda_approved$pert,]
```


```{r eval=FALSE}
write.csv(primary_paad_res_sub , file= "~/output/limma_gbm/220602_prism_limma_candidates.csv" )
```

```{r}
dim(drug_senstive_precentage_all_drugs)
```



```{r}
fraction_cell_line <- drug_senstive_precentage_all_drugs[1:30,]
```


```{r}
fraction_cell_line
```

# fraction first

data wrangle the cancers
```{r}
#Glioblastoma
##deseq2
drug_senstive_precentage_gbm_deseq2_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% lincs_deseq_gbm_res_GI1_fda_approved$pert,]
drug_senstive_precentage_gbm_deseq2_drugs$method <- rep("DESeq2", nrow(drug_senstive_precentage_gbm_deseq2_drugs))
##limma
drug_senstive_precentage_gbm_limma_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% lincs_limma_gbm_res_GI1_fda_approved$pert,]
drug_senstive_precentage_gbm_limma_drugs$method <- rep("limma", nrow(drug_senstive_precentage_gbm_limma_drugs))
##tfl
drug_senstive_precentage_gbm_tfl_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% lincs_deseq_gbm_res_GI1C_fda_approved$pert,]
drug_senstive_precentage_gbm_tfl_drugs$method <- rep("Transfer Learning", nrow(drug_senstive_precentage_gbm_tfl_drugs))

drug_senstive_precentage_gbm_methods<- rbind(drug_senstive_precentage_gbm_tfl_drugs, drug_senstive_precentage_gbm_limma_drugs, drug_senstive_precentage_gbm_deseq2_drugs)
drug_senstive_precentage_gbm_methods$cancer <- rep("Glioblastima", nrow(drug_senstive_precentage_gbm_methods))
```

```{r}
#Liver hepatocellular carcinoma
liver_deseq_tau <- read_csv("~/output/liver_cancer/deseq2_res/220602_SR_LINCS_LIHC_DESEQ2_RES_tau.csv")
liver_deseq_tau_fda_approved <- liver_deseq_tau[liver_deseq_tau$fda_approved == TRUE, ]

liver_limma_tau <- read_csv("~/output/liver_cancer/limma_res/220602_SR_LINCS_LIHC_limma_RES_tau.csv")
liver_limma_tau_fda_approved <- liver_limma_tau[liver_limma_tau$fda_approved == TRUE,]

liver_limma_tau <- read_csv("~/output/liver_cancer/limma_res/220602_SR_LINCS_LIHC_limma_RES_tau.csv")
liver_limma_tau_fda_approved <- liver_limma_tau[liver_limma_tau$fda_approved == TRUE,]

liver_tfl_tau <- read_csv("~/output/liver_cancer/TFL_res/220602_SR_LINCS_LIHC_TFL_RES_tau.csv")
liver_tfl_tau_fda_approved <- liver_tfl_tau[liver_tfl_tau$fda_approved == TRUE,]

##deseq2
drug_senstive_precentage_lihc_deseq2_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% liver_deseq_tau_fda_approved$pert,]
drug_senstive_precentage_lihc_deseq2_drugs$method <- rep("DESeq2", nrow(drug_senstive_precentage_lihc_deseq2_drugs))
##limma
drug_senstive_precentage_lihc_limma_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% liver_limma_tau_fda_approved$pert,]
drug_senstive_precentage_lihc_limma_drugs$method <- rep("limma", nrow(drug_senstive_precentage_lihc_limma_drugs))
##tfl
drug_senstive_precentage_lihc_tfl_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% liver_tfl_tau_fda_approved $pert,]
drug_senstive_precentage_lihc_tfl_drugs$method <- rep("Transfer Learning", nrow(drug_senstive_precentage_lihc_tfl_drugs))

drug_senstive_precentage_lihc_methods<- rbind(drug_senstive_precentage_lihc_tfl_drugs, drug_senstive_precentage_lihc_limma_drugs, drug_senstive_precentage_lihc_deseq2_drugs)
drug_senstive_precentage_lihc_methods$cancer <- rep("Liver hepatocellular carcinoma", nrow(drug_senstive_precentage_lihc_methods))
```

```{r}
#lung adenocaracionoma
lung_deseq_tau <- read_csv("~/output/lung_cancer/deseq2_res/220602_SR_LINCS_LUAD_DESEQ2_RES_tau.csv")
lung_deseq_tau_fda_approved <- lung_deseq_tau[lung_deseq_tau$fda_approved == TRUE, ]

lung_limma_tau <- read_csv("~/output/lung_cancer/limma_res/220602_SR_LINCS_LUAD_limma_RES_tau.csv")
lung_limma_tau_fda_approved <- lung_limma_tau[lung_limma_tau$fda_approved == TRUE,]

lung_limma_tau <- read_csv("~/output/lung_cancer/limma_res/220602_SR_LINCS_LUAD_limma_RES_tau.csv")
lung_limma_tau_fda_approved <- lung_limma_tau[lung_limma_tau$fda_approved == TRUE,]

lung_tfl_tau <- read_csv("~/output/lung_cancer/TFL_res/220602_SR_LINCS_LUAD_TFL_RES_tau.csv")
lung_tfl_tau_fda_approved <- lung_tfl_tau[lung_tfl_tau$fda_approved == TRUE,]

##deseq2
drug_senstive_precentage_luad_deseq2_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% lung_deseq_tau_fda_approved$pert,]
drug_senstive_precentage_luad_deseq2_drugs$method <- rep("DESeq2", nrow(drug_senstive_precentage_luad_deseq2_drugs))
##limma
drug_senstive_precentage_luad_limma_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% lung_limma_tau_fda_approved$pert,]
drug_senstive_precentage_luad_limma_drugs$method <- rep("limma", nrow(drug_senstive_precentage_luad_limma_drugs))
##tfl
drug_senstive_precentage_luad_tfl_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% lung_tfl_tau_fda_approved $pert,]
drug_senstive_precentage_luad_tfl_drugs$method <- rep("Transfer Learning", nrow(drug_senstive_precentage_luad_tfl_drugs))

drug_senstive_precentage_luad_methods<- rbind(drug_senstive_precentage_luad_tfl_drugs, drug_senstive_precentage_luad_limma_drugs, drug_senstive_precentage_luad_deseq2_drugs)
drug_senstive_precentage_luad_methods$cancer <- rep("Lung adenocaracionoma", nrow(drug_senstive_precentage_luad_methods))



```


```{r}
#pancreatic adenicarcinoma
pancreas_deseq_tau <- read_csv("~/output/pancreas_cancer/deseq2_res/220602_SR_LINCS_PAAD_DESEQ2_RES_tau.csv")
pancreas_deseq_tau_fda_approved <- pancreas_deseq_tau[pancreas_deseq_tau$fda_approved == TRUE, ]

pancreas_limma_tau <- read_csv("~/output/pancreas_cancer/limma_res/220602_SR_LINCS_PAAD_limma_RES_tau.csv")
pancreas_limma_tau_fda_approved <- pancreas_limma_tau[pancreas_limma_tau$fda_approved == TRUE,]

pancreas_limma_tau <- read_csv("~/output/pancreas_cancer/limma_res/220602_SR_LINCS_PAAD_limma_RES_tau.csv")
pancreas_limma_tau_fda_approved <- pancreas_limma_tau[pancreas_limma_tau$fda_approved == TRUE,]

pancreas_tfl_tau <- read_csv("~/output/pancreas_cancer/TFL_res/220602_SR_LINCS_PAAD_TFL_RES_tau.csv")
pancreas_tfl_tau_fda_approved <- pancreas_tfl_tau[pancreas_tfl_tau$fda_approved == TRUE,]

##deseq2
drug_senstive_precentage_paad_deseq2_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% pancreas_deseq_tau_fda_approved$pert,]
drug_senstive_precentage_paad_deseq2_drugs$method <- rep("DESeq2", nrow(drug_senstive_precentage_paad_deseq2_drugs))
##limma
drug_senstive_precentage_paad_limma_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% pancreas_limma_tau_fda_approved$pert,]
drug_senstive_precentage_paad_limma_drugs$method <- rep("limma", nrow(drug_senstive_precentage_paad_limma_drugs))
##tfl
drug_senstive_precentage_paad_tfl_drugs<-drug_senstive_precentage_all_drugs[ drug_senstive_precentage_all_drugs$drug_list %in% pancreas_tfl_tau_fda_approved $pert,]
drug_senstive_precentage_paad_tfl_drugs$method <- rep("Transfer Learning", nrow(drug_senstive_precentage_paad_tfl_drugs))

drug_senstive_precentage_paad_methods<- rbind(drug_senstive_precentage_paad_tfl_drugs, drug_senstive_precentage_paad_limma_drugs, drug_senstive_precentage_paad_deseq2_drugs)
drug_senstive_precentage_paad_methods$cancer <- rep("Pancreatic adenicarcinoma", nrow(drug_senstive_precentage_paad_methods))

```

```{r}
drug_senstive_precentage_all_methods_cancers <- rbind(drug_senstive_precentage_paad_methods, drug_senstive_precentage_gbm_methods, drug_senstive_precentage_lihc_methods, drug_senstive_precentage_luad_methods)
```

```{r}
t <- ggplot(drug_senstive_precentage_all_methods_cancers,aes(x=method, y=drug_fraction ,fill=method))+
  #geom_point()+
  geom_violin() +
  geom_boxplot(width = 0.1, fill = "grey", color = "black") +
  facet_wrap(~cancer,nrow=1) + scale_fill_manual(values= c("#440154FF" ,"#21908CFF" ,"#FDE725FF")) + ylab("Fraction of Sensitive Cell Lines")
#t + geom_text(aes(label=labels_v2),vjust=- 15,  size= 10) +ylim(0, 90)
t
```

```{r}
primary_paad_res_v2 <- primary_paad_res[!is.na(primary_paad_res$log_fold_change),]


drug_list <- as.vector(unique(primary_paad_res_v2$Var2))
drug_median<- c()
for (i in 1:length(drug_list)){
  cells <-primary_paad_res_v2$log_fold_change[primary_paad_res_v2$Var2 %in% drug_list[i]]
  drug_median[i] <- median(cells)
}
drug_prism_median_all_drugs <- data.frame( drug_list, drug_median)
```


```{r}
#Glioblastoma
##deseq2
drug_senstive_precentage_gbm_deseq2_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% lincs_deseq_gbm_res_GI1_fda_approved$pert,]
drug_senstive_precentage_gbm_deseq2_drugs$method <- rep("DESeq2", nrow(drug_senstive_precentage_gbm_deseq2_drugs))
##limma
drug_senstive_precentage_gbm_limma_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% lincs_limma_gbm_res_GI1_fda_approved$pert,]
drug_senstive_precentage_gbm_limma_drugs$method <- rep("limma", nrow(drug_senstive_precentage_gbm_limma_drugs))
##tfl
drug_senstive_precentage_gbm_tfl_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% lincs_deseq_gbm_res_GI1C_fda_approved$pert,]
drug_senstive_precentage_gbm_tfl_drugs$method <- rep("Transfer Learning", nrow(drug_senstive_precentage_gbm_tfl_drugs))

drug_senstive_precentage_gbm_methods<- rbind(drug_senstive_precentage_gbm_tfl_drugs, drug_senstive_precentage_gbm_limma_drugs, drug_senstive_precentage_gbm_deseq2_drugs)
drug_senstive_precentage_gbm_methods$cancer <- rep("Glioblastima", nrow(drug_senstive_precentage_gbm_methods))
```

```{r}
#Liver hepatocellular carcinoma

##deseq2
drug_senstive_precentage_lihc_deseq2_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% liver_deseq_tau_fda_approved$pert,]
drug_senstive_precentage_lihc_deseq2_drugs$method <- rep("DESeq2", nrow(drug_senstive_precentage_lihc_deseq2_drugs))
##limma
drug_senstive_precentage_lihc_limma_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% liver_limma_tau_fda_approved$pert,]
drug_senstive_precentage_lihc_limma_drugs$method <- rep("limma", nrow(drug_senstive_precentage_lihc_limma_drugs))
##tfl
drug_senstive_precentage_lihc_tfl_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% liver_tfl_tau_fda_approved $pert,]
drug_senstive_precentage_lihc_tfl_drugs$method <- rep("Transfer Learning", nrow(drug_senstive_precentage_lihc_tfl_drugs))

drug_senstive_precentage_lihc_methods<- rbind(drug_senstive_precentage_lihc_tfl_drugs, drug_senstive_precentage_lihc_limma_drugs, drug_senstive_precentage_lihc_deseq2_drugs)
drug_senstive_precentage_lihc_methods$cancer <- rep("Liver hepatocellular carcinoma", nrow(drug_senstive_precentage_lihc_methods))
```

```{r}
#lung adenocaracionoma

##deseq2
drug_senstive_precentage_luad_deseq2_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% lung_deseq_tau_fda_approved$pert,]
drug_senstive_precentage_luad_deseq2_drugs$method <- rep("DESeq2", nrow(drug_senstive_precentage_luad_deseq2_drugs))
##limma
drug_senstive_precentage_luad_limma_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% lung_limma_tau_fda_approved$pert,]
drug_senstive_precentage_luad_limma_drugs$method <- rep("limma", nrow(drug_senstive_precentage_luad_limma_drugs))
##tfl
drug_senstive_precentage_luad_tfl_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% lung_tfl_tau_fda_approved $pert,]
drug_senstive_precentage_luad_tfl_drugs$method <- rep("Transfer Learning", nrow(drug_senstive_precentage_luad_tfl_drugs))

drug_senstive_precentage_luad_methods<- rbind(drug_senstive_precentage_luad_tfl_drugs, drug_senstive_precentage_luad_limma_drugs, drug_senstive_precentage_luad_deseq2_drugs)
drug_senstive_precentage_luad_methods$cancer <- rep("Lung adenocaracionoma", nrow(drug_senstive_precentage_luad_methods))



```


```{r}
#pancreatic adenicarcinoma

##deseq2
drug_senstive_precentage_paad_deseq2_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% pancreas_deseq_tau_fda_approved$pert,]
drug_senstive_precentage_paad_deseq2_drugs$method <- rep("DESeq2", nrow(drug_senstive_precentage_paad_deseq2_drugs))
##limma
drug_senstive_precentage_paad_limma_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% pancreas_limma_tau_fda_approved$pert,]
drug_senstive_precentage_paad_limma_drugs$method <- rep("limma", nrow(drug_senstive_precentage_paad_limma_drugs))
##tfl
drug_senstive_precentage_paad_tfl_drugs<-drug_prism_median_all_drugs[ drug_prism_median_all_drugs$drug_list %in% pancreas_tfl_tau_fda_approved $pert,]
drug_senstive_precentage_paad_tfl_drugs$method <- rep("Transfer Learning", nrow(drug_senstive_precentage_paad_tfl_drugs))

drug_senstive_precentage_paad_methods<- rbind(drug_senstive_precentage_paad_tfl_drugs, drug_senstive_precentage_paad_limma_drugs, drug_senstive_precentage_paad_deseq2_drugs)
drug_senstive_precentage_paad_methods$cancer <- rep("Pancreatic adenicarcinoma", nrow(drug_senstive_precentage_paad_methods))

```

```{r}
drug_senstive_precentage_all_methods_cancers <- rbind(drug_senstive_precentage_paad_methods, drug_senstive_precentage_gbm_methods, drug_senstive_precentage_lihc_methods, drug_senstive_precentage_luad_methods)
```



```{r}
t <- ggplot(drug_senstive_precentage_all_methods_cancers,aes(x=method, y=drug_median ,fill=method))+
  #geom_point()+
  geom_violin() + geom_hline(yintercept=0.3, linetype="dashed", color = "red") +
  geom_boxplot(width = 0.1, fill = "grey", color = "black") +
  facet_wrap(~cancer,nrow=1) + scale_fill_manual(values= c("#440154FF" ,"#21908CFF" ,"#FDE725FF")) + ylab("Median log2fold change (PRISM)")
#t + geom_text(aes(label=labels_v2),vjust=- 15,  size= 10) +ylim(0, 90)
t
```


add the limma gene set list 

```{r}
GBM_GTEX_gene_limma_res <- readRDS("~/output/limma_gbm/220421_GBM_GTEX_gene_limma_res.rds")
```


```{r}
GBM_gene_limma_sig<- GBM_GTEX_gene_limma_res$limma[abs(GBM_GTEX_gene_limma_res$limma$logFC) > 2 & GBM_GTEX_gene_limma_res$limma$adj.P.Val < 0.05,]
```
filter to genes only in latent varaibles 

GBM_gene_limma_sig
```{r}
limma_genes_v2 <- GBM_gene_limma_sig$LV[GBM_gene_limma_sig$LV %in% rownames(Z_mtx)]
```


compare gene list across different methods

```{r}
library(VennDiagram)
```
```{r}
# Prepare a palette of 3 colors with R colorbrewer:
myCol <- c("#440154FF" , "#31688EFF" ,"#35B779FF")
```

```{r}
venn.diagram(
  x = list(limma_genes_v2, deseq2_genes_v2, topgenes),
  category.names = c("limma" , "DESeq2" , "Transfer Learning"),
  filename = '~/output/gene_gbm_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.3,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
```


pathway analysis for each unique gene set
```{r}
un_tf_genes <- setdiff(topgenes, c(limma_genes_v2, deseq2_genes_v2))
```

```{r}
deseq_results<- deseq_results[complete.cases(deseq_results),]
up_genes<- deseq_results$Symbol[deseq_results$log2FoldChange >0]
down_genes<- deseq_results$Symbol[deseq_results$log2FoldChange <0]
```


```{r}
pathway_geneset_analysis<- function(geneset , method ){
  #get up down gene set 
  up <- geneset[geneset %in% up_genes]
  #print(up)
  down <- geneset[geneset %in% down_genes]
  downset_pathway_results <- gost(query = down, 
                                  organism = "hsapiens", ordered_query = TRUE, 
                                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                                  user_threshold = 0.05, correction_method = "g_SCS", 
                                  domain_scope = "annotated", custom_bg = rownames(Z_mtx), 
                                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
  downset_pathway_results<- downset_pathway_results$result
  
  upset_pathway_results <- gost(query = up, 
                                organism = "hsapiens", ordered_query = TRUE, 
                                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                measure_underrepresentation = FALSE, evcodes = FALSE, 
                                user_threshold = 0.05, correction_method = "g_SCS", 
                                domain_scope = "annotated", custom_bg = rownames(Z_mtx), 
                                numeric_ns = "", sources = NULL, as_short_link = FALSE)
  upset_pathway_results<- upset_pathway_results$result
  
  pathway_results<- rbind(upset_pathway_results, downset_pathway_results)
  pathway_results$set<- c(rep("Up", nrow(upset_pathway_results)), rep("Down", nrow(downset_pathway_results)))
  file_name<- paste0("~/output/gbm_",  method, "_gene_pathways.csv")
  write.csv2(as.data.frame(pathway_results[,-14]), file_name )
  
  file_name<- paste0("~/output/gbm_", method, "_gene_pathways.png")
  if(nrow(pathway_results)<50){
    ggplot(pathway_results, aes(x=set, y=term_name, fill=p_value))+geom_tile()+scale_fill_viridis(direction=-1)+theme_classic()+ labs(title=name,x="Drug Gene Sets", y = "g:Profiler Gene Sets", fill= "p-value")
    
  }else{
    pathway_results<- pathway_results[pathway_results$source %in% c("KEGG","REAC","GO:MF" ),]
    ggplot(pathway_results, aes(x=set, y=term_name, fill=p_value)) + geom_tile()+ scale_fill_viridis(direction=-1)+ theme_classic() + labs(title=method,x="Drug Gene Sets", y = "g:Profiler Gene Sets", fill= "p-value")
  }
  ggsave(file_name, width = 12, height = 14)
}
```

note due to the overlapping of genes between the different methods. A direct statistical test can not be done due to the independence assumption. 

pathways analysi of genes from limma and DESEQ2
```{r}
library(gprofiler2)
#no signifcant pathways
#pathway_geneset_analysis(un_tf_genes, "transfer_learning")
```

```{r}
pathway_geneset_analysis(limma_genes_v2, "limma")
```

```{r}
pathway_geneset_analysis(deseq2_genes_v2, "DESeq2")
```


To compare the different centrality metrics for the GBM output genes from limma, deseq2, and the transfer learning approach. For each method, I randomly selected the same number of genes and calculated the median of each of the three log2 transformed centrality metrics. This was done 100,000 times. A one-way t-test was conducted to determine if the genes determined if one of the methods had a higher median of the centrality metrics than by chance. We used this approach over a comparison of one method vs another via Wilcox test and ANOVA due to the assumption of independence because of the overlap of genes between methods. 


```{r}
network_metric_anaylsis<- function(geneset){
  #calculate the median for each metric for method
  degree_med <- median(string_ppi_details$degree_log[ rownames(string_ppi_details) %in% geneset ])
  print(degree_med)
  betweenness_med <- median(string_ppi_details$log_betweness[ rownames(string_ppi_details) %in% geneset ])
  print(betweenness_med)
  eign_cent_med <- median(string_ppi_details$eign_cent_log[ rownames(string_ppi_details) %in% geneset ])
  print(eign_cent_med)
  
  random_degree_med <- c()
  random_betweenness_med<- c()
  random_eign_cent_med<- c() 
  
  for( i in 1:100000){
    sub_genes <- rownames(Z_mtx)[sample(1:nrow(Z_mtx), length(geneset))]
    sub_ppi <- string_ppi_details[rownames(string_ppi_details) %in% sub_genes,]
    random_degree_med[i] <- median(sub_ppi$degree_log)
    random_betweenness_med[i] <- median(sub_ppi$log_betweness)
    random_eign_cent_med[i] <- median(sub_ppi$eign_cent_log)
  }
  #print(hist(random_degree_med))
  #print("wilcox degree")
  #print(random_degree_med[1:5])
  print(t.test(random_degree_med,mu= degree_med, alternative = "less" ))
  
  #print(hist(random_betweenness_med))
  # print("wilcox betweenness")
  print(t.test(random_betweenness_med,mu=betweenness_med, alternative = "less" ))
  
  #print(hist(random_eign_cent_med))
  #print("wilcox betweenness")
  print(t.test(random_eign_cent_med, mu=eign_cent_med, alternative ="less"))
  
  #degree
  
  #Print plot   
  
  #betweeness 
  
  #Print plot
  
  #eignvector 
  #print plot 
  
  
}
```

print("kruskal and pairwise.wilcox degree")
print(kruskal.test(degree ~ method, data = ppi_details_sig_3_methods_filtered))
print(pairwise.wilcox.test(ppi_details_sig_3_methods_filtered$degree, ppi_details_sig_3_methods_filtered$method, p.adjust.method = "bonf"))

print("kruskal and pairwise.wilcox betweenness")
print(kruskal.test(betweenness ~ method, data = ppi_details_sig_3_methods_filtered))
print(pairwise.wilcox.test(ppi_details_sig_3_methods_filtered$betweenness, ppi_details_sig_3_methods_filtered$method, p.adjust.method = "bonf"))

print("kruskal and pairwise.wilcox betweenness")
print(kruskal.test(eign_cent ~ method, data = ppi_details_sig_3_methods_filtered))
print(pairwise.wilcox.test(ppi_details_sig_3_methods_filtered$eign_cent, ppi_details_sig_3_methods_filtered$method, p.adjust.method = "bonf"))


#plot each metric 
#ppi_details_sig_3_methods$degree_log <- log(ppi_details_sig_3_methods$degree)
p1 <- ggplot(ppi_details_sig_3_methods, aes(x=method , y=degree_log, fill=method )) + 
  geom_violin()+ geom_boxplot(width = 0.1, fill = "grey", color = "black")+labs(title="PPI Degree",x="Genes from Disease Signature for Signature Reversion", y = "PPI log(Degree)", color= "Method") + scale_colour_manual(values = c("#440154FF" ,"#21908CFF" ,"#FDE725FF"), aesthetics = c("colour", "fill"))+ theme_minimal()


#ppi_details_sig_3_methods$log_betweness <- log(ppi_details_sig_3_methods$betweenness)
p2<- ggplot(ppi_details_sig_3_methods, aes(x=method, y=log_betweness , fill=method)) + 
  geom_violin()+ geom_boxplot(width = 0.1, fill = "grey", color = "black")+labs(title="PPI Betweeness",x="Genes from Disease Signature for Signature Reversion", y = "PPI log(Betweeness)", color= "Method") + scale_colour_manual(values =  c("#440154FF" ,"#21908CFF" ,"#FDE725FF"), aesthetics = c("colour", "fill"))+ theme_minimal()


#ppi_details_sig_3_methods$eign_cent_log <- log(ppi_details_sig_3_methods$eign_cent)
p3<- ggplot(ppi_details_sig_3_methods, aes(x=method, y=eign_cent_log, fill=method)) + 
  geom_violin()+ geom_boxplot(width = 0.1, fill = "grey", color = "black")  +labs(title="PPI Eigenvector Centrality Scores",x="Genes from Disease Signature for Signature Reversion", y = "PPI log(Eigenvector Centrality Scores)", color= "Method") + scale_colour_manual(values =  c("#440154FF" ,"#21908CFF" ,"#FDE725FF"), aesthetics = c("colour", "fill"))+ theme_minimal()

#print(p1)
#print(p2)
#print(p3)


#return(ppi_details_sig_3_methods)
print("wilcox limma")
#print(limma_genes_v2$target)
network_metric_anaylsis(limma_genes_v2$target)
print("wilcox deseq2")
network_metric_anaylsis(deseq2_genes_v2$target)
print("wilcox transfer learning")
network_metric_anaylsis(topgenes$target)


```{r}
transfer_gbm_res_GI1_tau_fda_approved <- read_csv("~/Documents/Transfer_Learning_R03/output/TF_L_GBM/220207_lincs_transfer_gbm_res_GI1_tau_fda_approved.csv")
```


```{r}
#some deprecated or multiple ids 
tfl_gbm_ids <- unique(lincs_commpound_info$pert_id[grepl(paste(transfer_gbm_res_GI1_tau_fda_approved$pert,collapse="|"), lincs_commpound_info$cmap_name) ])
```


```{r}
qres <- queryAnnotDB(tfl_gbm_ids, annot=c("lincsAnnot"))
qres
```

```{r}
#if there are multiples keep the toouchstone case
qres[qres$is_touchstone == 1,]
```



```{r}
IDs <- lincs_sdfset@ID[grepl(paste(transfer_gbm_res_GI1_tau_fda_approved$pert,collapse="|"), lincs_sdfset@ID)]
```
two of the candidates sxagliptin and pamidronte were not in the dataset 


```{r}
tfl_gbm_sdfset <- lincs_sdfset[IDs]
```

```{r}

```


```{r}
d <- sapply(cid(tfl_gbm_sdfset), function(x) fmcsBatch(tfl_gbm_sdfset[x], tfl_gbm_sdfset, au=0, bu=0, matching.mode = "aromatic")[,"Overlap_Coefficient"])
```

```{r}
hc <- hclust(as.dist(1-d), method="complete")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE) 
```


```{r}
drug_list <- unique(c(lincs_tfl_gbm_res_GI1_fda_approved$pert,lincs_limma_gbm_res_GI1_fda_approved$pert,  lincs_deseq2_gbm_res_GI1_fda_approved$pert))
```


```{r}
IDs <- lincs_sdfset@ID[grepl(paste(drug_list,collapse="|"), lincs_sdfset@ID)]
```

```{r}
drug_gbm_sdfset <- lincs_sdfset[IDs]
```

```{r}
drug_gbm_res <- sapply(cid(drug_gbm_sdfset ), function(x) fmcsBatch(drug_gbm_sdfset[x], drug_gbm_sdfset , au=0, bu=0, matching.mode = "aromatic")[,"Tanimoto_Coefficient"])
```
```{r}
hc <- hclust(as.dist(1-drug_gbm_res ), method="complete")
plot(as.dendrogram(hc), edgePar=list(col=4, lwd=2), horiz=TRUE) 
```

```{r}
library(dendextend)
dendro <- as.dendrogram(hc)
tfl_drugs2 <-  IDs %in% lincs_tfl_gbm_res_GI1_fda_approved$pert
limma_drugs2 <- IDs %in% lincs_limma_gbm_res_GI1_fda_approved$pert
deseq2_drugs2<- IDs %in% lincs_deseq2_gbm_res_GI1_fda_approved$pert

test<- ifelse(tfl_drugs2== TRUE, "#440154FF", "#228C8DFF")
test2<- ifelse(limma_drugs2== TRUE, "#440154FF", "#228C8DFF")
test3<- ifelse(deseq2_drugs2== TRUE, "#440154FF", "#228C8DFF")
#test<- ifelse(tfl_drugs== TRUE, "#440154FF", "#228C8DFF")
#test2<- ifelse(limma_drugs== TRUE, "#440154FF", "#228C8DFF")
#test3<- ifelse(deseq2_drugs== TRUE, "#440154FF", "#228C8DFF")

par(oma=c(0.5,0.5,0.5,20), mar = c(4, 0.01, 0.01, 12))
#dend2 <-prune( dend, med_bottom, reindex_dend=FALSE)
dendro %>%
  set("labels_col",value =c( "#117733", "#661100",  "#0072B2", "#D55E00", "#AA4499"), k=5) %>%
  set("branches_k_color",value =c("#117733", "#661100",  "#0072B2", "#D55E00", "#AA4499"),  k = 5) %>%
  set("leaves_pch", 19)  %>% 
  set("nodes_cex", 0.6) %>% 
  set("branches_lwd", 3) %>% 
  set("labels_cex", 0.5) %>%
  plot( horiz=TRUE, axes=FALSE)
#addjust eh y shift
colored_bars(cbind(test, test2, test3), dendro ,rowLabels = c("Transfer", "limma", "DESeq2"), y_shift = 4, horiz = TRUE)
#note have to save manually via export function as a pdf
```


```{r}
library(ComplexHeatmap)
library(circlize)
```


```{r}
row_ha = HeatmapAnnotation(Transfer_Learning=tfl_drugs2,DESeq2= deseq2_drugs2, limma= limma_drugs2 , col = list(Transfer_Learning = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),DESeq2 = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),limma = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF") ))
col_fun = colorRamp2(c(0,  1), c( "black", "yellow"))
#ha = rowAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 97:100), labels = TFL_bp_up[1:10]))

Heatmap(drug_gbm_res, nam= "Overlap_Coefficient drug structure", col = col_fun, show_column_names = FALSE,  top_annotation = row_ha,  
        clustering_distance_rows= "euclidean",
        clustering_distance_columns=  "euclidean",
        clustering_method_rows = "ward.D2" ,
        clustering_method_columns="ward.D2")
```

