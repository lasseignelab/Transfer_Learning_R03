#functions for cancer signature reversion paper 
# Jennifer Fisher
# jfisher7@uab.edu 

# the following function was from the MultiPLIER paper and github. I just added a few print statements in function for updates during run
PLIERNewData_JLF <- function(exprs.mat, seed = 12345) {
  # A wrapper function for applying PLIER to a data set. We use the following
  # genesets that come with PLIER: bloodCellMarkersIRISDMAP, svmMarkers, 
  # and canonicalPathways. We set the k parameter for the PLIER model by
  # identifying the number of "significant PCs" with PLIER::num.pc and then 
  # using sig PCs * 0.3. This is consistent with recommendations from the 
  # PLIER authors.
  # 
  # Args:
  #   exprs.mat: a gene expression matrix, rows are genes, columns are samples
  #   seed: an integer to be supplied to set.seed() for reproducibility 
  #         purposes, default is 12345
  #         
  # Returns:
  #   plier.res: output from PLIER::PLIER()
  #
  require(PLIER)
  
  set.seed(seed)
  
  # load PLIER pathway and cell type data
  data(bloodCellMarkersIRISDMAP)
  data(svmMarkers)
  data(canonicalPathways)
  
  print("combine the pathway data from PLIER")
  all.paths <- PLIER::combinePaths(bloodCellMarkersIRISDMAP, svmMarkers, 
                                   canonicalPathways)
  
  # what genes are common to the pathway data and the expression matrix
  cm.genes <- PLIER::commonRows(all.paths, exprs.mat)
  
  print(" row normalize")
  exprs.norm <- PLIER::rowNorm(exprs.mat)
  
  # what should we set the minimum k parameter to in PLIER? estimate the number 
  # of PC for the SVD decomposition 
  print("num.pc")
  mtx <- exprs.norm[cm.genes, ]
  mtx<- mtx[complete.cases(mtx),]
  print(class(mtx))
  print(dim(mtx))
  set.k <- PLIER::num.pc(mtx)
  
  # PLIER main function + return results
  print("PLIER main function")
  plier.res <- PLIER::PLIER(mtx, all.paths[cm.genes, ], scale= F,
                            k = round((set.k + set.k * 0.3), 0), trace = TRUE)
  
  return(plier.res)
  
}

#function to make raw counts to TPM 
Counts_to_tpm <- function(counts, featureLength) {
  #counts the data frame
  #featureLength is a vector with the gene lengths from Recount3
  
  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  # Compute effective lengths of features in each library.
  effLen <- featureLength
  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    rate = log(counts[,i]) - log(effLen)
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(counts)
  rownames(tpm) <- rownames(counts)
  return(tpm)
}


#deseq2 function
deseq2_function <- function(dds, file_name, condition= "status", B= "Primary Tumor", A= "Solid Tissue Normal"){
  #dds - deseq2 object from the DESeqDataSetFromMatrix function. This will contain the counts and metadata
  #file name - ex. "/data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/output/deseq2_gbm/Deseq2_gbm_normal_gtx_res.csv"
  #condition is the condition tested for the results 
  #B and A are the different groups to compare to
  #B-TUMOR 
  #A-CONTROL
  dds$status <- factor(dds$status, levels = c("Solid Tissue Normal","Primary Tumor"))
  #run deseq2 model
  dds <- DESeq(dds)
  print(resultsNames(dds))
  #get results and save them and return it 
  #dds_res <- results(dds, c(condition,B,A))
  #res_df <- as.data.frame(dds_res)
  res.ape <- as.data.frame(lfcShrink(dds=dds, coef=3, type="apeglm"))
  write.table(res.ape, file=file_name, sep=",", row.names =TRUE)
  return(res.ape)
}


#Transfer Learning functions 

#taken from the multiplier paper but added the blocking variable for GTEx vs TCGA databases
TestLVDifferencesJLF <- function(b.matrix, phenotype, blocking,
                                 use.bonferroni = FALSE) {
  # This function tests for differential expression of PLIER-derived latent
  # variables between groups specified by the phenotype argument
  # using limma. Examples include different disease groups or disease v. 
  # control. It takes the B matrix from a PLIER model and a factor vector with 
  # the group labels. It will reorder the phenotype vector if control.phenotype
  # is specified. The resulting p-values are Benjamini-Hochberg corrected by 
  # default or Bonferroni corrected if use.bonferroni = TRUE
  # 
  # Args:
  #   b.matrix: a B matrix (latent variables are rows, samples are columns) 
  #             from a PLIER model. Can be the B element of the list 
  #             returned by PLIER::PLIER or the output of GetNewDataB
  #   phenotype: a named factor vector that contains group labels to be used
  #              for contrasts
  #   blocking: a names facotor vector that contains group labels to be used for blocking
  #   use.bonferroni: logical - should bonferroni correction be used? if FALSE
  #                   (default), will use "BH"
  #   
  # Returns:
  #   A limma::topTable, where the first column is the latent variable name and
  #   all pathways are returned (without filtering or sorting)
  
  ## error-handling ##
  
  if (is.null(names(phenotype))) {
    stop("phenotype should be a named factor vector -- the names ensure that
         the vector is correctly ordered")
  }
  if (is.null(names(blocking))) {
    stop("blocking should be a named factor vector -- the names ensure that
         the vector is correctly ordered")
  }  
  # no names should be "missing"
  check.names <- all(colnames(b.matrix) %in% names(phenotype)) & 
    all(names(phenotype) %in% colnames(b.matrix)) 
  
  if (!check.names) {
    stop("Some sample(s) is missing from colnames(b.matrix) or 
         names(phenotype)")
  }
  
  # the phenotype labels should be in the same order as the b.matrix samples
  ordered.phenotype <- as.factor(phenotype[colnames(b.matrix)])
  
  # the blocking labels should be in the same order as the b.matrix samples
  ordered.blocking <- as.factor(blocking[colnames(b.matrix)])
  
  # get contrasts (all pairwise)  
  num.levels <- length(levels(ordered.phenotype))
  contrast.vector <- c()
  for (lvl in levels(ordered.phenotype)[1:(num.levels - 1)]) {
    for (lvl.2 in levels(ordered.phenotype)[2:num.levels]) {
      if (lvl != lvl.2) {
        contrast.vector <- append(contrast.vector, paste(lvl, lvl.2, sep = "-"))
      }
    }
  }
  contrast.vector <- paste(contrast.vector, collapse = ", ")
  
  # prep design matrix
  design <- model.matrix(~ 0+ ordered.phenotype+ ordered.blocking )
  colnames(design)[seq_len(nlevels(ordered.phenotype))] <- levels(ordered.phenotype)
  
  
  
  
  #colnames(design) <- levels(ordered.phenotype)
  
  # fit linear model
  fit <- limma::lmFit(b.matrix, design)
  contrast.matrix <-
    eval(parse(
      text = paste0(
        'limma::makeContrasts(', contrast.vector, ', levels = design)')))
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2)
  
  # extract results as a data.frame
  if (use.bonferroni) {  # if specified, use Bonferroni correction
    limma.result <- limma::topTable(fit2, number = nrow(b.matrix),
                                    adjust.method = "bonferroni", 
                                    sort.by = "none")
  } else {  # calculate FDR
    limma.result <- limma::topTable(fit2, number = nrow(b.matrix),
                                    adjust.method = "BH", sort.by = "none")
  }
  
  # want feature (latent variable) names as column
  limma.result <- tibble::rownames_to_column(limma.result, var = "LV")
  
  return(limma.result)
  
}
#taken from the multiplier paper but added the blocking variable for GTEx vs TCGA databases
LVTestWrapperJLF <- function(b.matrix,
                             sample.info.df,
                             phenotype.col,
                             blocking,
                             file.lead,
                             plot.dir = "plots",
                             results.dir = "results",
                             use.bonferroni = FALSE,
                             significant.only = FALSE,
                             sig.threshold = 0.05) {
  # A wrapper function for TestLVDifferences and BoxplotDiffLV; does the 
  # reshaping required for plotting. Produces the following files: 1) tsv of the
  # differential expression results 2) a long form of the B matrix joined with
  # the sample information (sample.info.df) and 3) a PDF of boxplots
  # 
  # Args:
  #   b.matrix: a B matrix (latent variables are rows, samples are columns) 
  #             from a PLIER model. Can be the B element of the list 
  #             returned by PLIER::PLIER or the output of GetNewDataB
  #   sample.info.df: a long form data.frame that contains sample information,
  #                   sample names that match the B matrix sample identifiers
  #                   must be in a "Sample" column & it also must contain
  #                   the factor to group by for testing differential expression
  #                   (DE) and plotting
  #   phenotype.col: the column name of the column in sample.info.df to be used
  #                  for DE and plotting; character
  #   file.lead: string that designates the "beginning" of filenames
  #   plot.dir: plot directory where the boxplots PDF should be saved
  #   results.dir: results directory for DE results and reshaped B data.frame
  #   use.bonferroni: logical - should bonferroni correction be used for DE? 
  #                   if FALSE (default), will use "BH"
  #   significant.only: logical - should only differentially expressed LVs be
  #                     plotted? if FALSE (default), all will be plotted
  #   sig.threshold: the adj. P cutoff to be used if only plotting significant 
  #                  results; default is 0.05
  #                  
  # Returns:
  #   A list with the following elements
  #       limma: the results from TestLVDifferences 
  #       b.df: sample.b.df, prior to any filtering for plotting (if applicable)
  # 
  #   the following files are written by this function (see above):
  #     1) <results.dir>/<file.lead>_LV_limma_results.tsv
  #     2) <results.dir>/<file.lead>_B_long_sample_info.tsv
  #     3) <plot.dir>/<file.lead>_LV_boxplots.pdf
  
  `%>%` <- dplyr::`%>%`
  
  # error-handling
  # we need to join by "Sample" column for boxplots
  if (!("Sample" %in% colnames(sample.info.df))) {
    stop("'Sample' must be a column in sample.info.df")
  }
  # phenotype.col needs to be in column names
  if (!(phenotype.col %in% colnames(sample.info.df))) {
    stop("phenotype.col must be a column name in sample.info.df")
  }
  
  # initialize list to hold results to be returned
  return.list <- list()
  
  #### Differential Expression ####
  # get the named vector to use as the phenotype for testing differential 
  # expression
  phenotype.vector <- as.factor(make.names(sample.info.df[[phenotype.col]]))
  names(phenotype.vector) <- sample.info.df$Sample
  # test itself
  limma.df <- TestLVDifferencesJLF(b.matrix = b.matrix,
                                   phenotype = phenotype.vector,
                                   blocking= blocking,
                                   use.bonferroni = use.bonferroni)
  # write to file
  dlve.file <- file.path(results.dir, 
                         paste0(file.lead, "_LV_limma_results.tsv"))
  readr::write_tsv(limma.df, path = dlve.file)
  # add to list to be returned
  return.list[["limma"]] <- limma.df
  
  #### Reshape & join with sample information ####
  b.df <- reshape2::melt(b.matrix)
  colnames(b.df) <- c("LV", "Sample", "Value")
  sample.b.df <- dplyr::inner_join(b.df, sample.info.df, by = "Sample")
  long.file <- file.path(results.dir, 
                         paste0(file.lead, "_B_long_sample_info.tsv"))
  readr::write_tsv(sample.b.df, long.file)
  # add to list to be returned
  return.list[["b.df"]] <- sample.b.df
  
  #### Plotting ####
  plot.file <- file.path(plot.dir,
                         paste0(file.lead, "_LV_boxplots.pdf"))
  
  # if we only want significant LVs plotted, filter sample.b.df using the
  # adj.P.Val cutoff sig.threshold
  if (significant.only) {
    sig.lvs <- limma.df$LV[which(limma.df$adj.P.Val < sig.threshold)]
    sample.b.df <- sample.b.df %>%
      dplyr::filter(LV %in% sig.lvs)
  }
  
  BoxplotDiffLV(tidy.b.df = sample.b.df, 
                phenotype.column = phenotype.col,
                pdf.path = plot.file)
  
  
  return(return.list)
}

#limma wrapper- smiliar to Transfer Learning wrapper, but adjusted to not plot the boxplots because it genes in this case it is just a ton plots
LVTestWrapperJLF_v2 <- function(b.matrix,
                                sample.info.df,
                                phenotype.col,
                                blocking,
                                file.lead,
                                plot.dir = "plots",
                                results.dir = "results",
                                use.bonferroni = FALSE,
                                significant.only = FALSE,
                                sig.threshold = 0.05) {
  # A wrapper function for TestLVDifferences and BoxplotDiffLV; does the 
  # reshaping required for plotting. Produces the following files: 1) tsv of the
  # differential expression results 2) a long form of the B matrix joined with
  # the sample information (sample.info.df) and 3) a PDF of boxplots
  # 
  # Args:
  #   b.matrix: a B matrix (latent variables are rows, samples are columns) 
  #             from a PLIER model. Can be the B element of the list 
  #             returned by PLIER::PLIER or the output of GetNewDataB
  #   sample.info.df: a long form data.frame that contains sample information,
  #                   sample names that match the B matrix sample identifiers
  #                   must be in a "Sample" column & it also must contain
  #                   the factor to group by for testing differential expression
  #                   (DE) and plotting
  #   phenotype.col: the column name of the column in sample.info.df to be used
  #                  for DE and plotting; character
  #   file.lead: string that designates the "beginning" of filenames
  #   plot.dir: plot directory where the boxplots PDF should be saved
  #   results.dir: results directory for DE results and reshaped B data.frame
  #   use.bonferroni: logical - should bonferroni correction be used for DE? 
  #                   if FALSE (default), will use "BH"
  #   significant.only: logical - should only differentially expressed LVs be
  #                     plotted? if FALSE (default), all will be plotted
  #   sig.threshold: the adj. P cutoff to be used if only plotting significant 
  #                  results; default is 0.05
  #                  
  # Returns:
  #   A list with the following elements
  #       limma: the results from TestLVDifferences 
  #       b.df: sample.b.df, prior to any filtering for plotting (if applicable)
  # 
  #   the following files are written by this function (see above):
  #     1) <results.dir>/<file.lead>_LV_limma_results.tsv
  #     2) <results.dir>/<file.lead>_B_long_sample_info.tsv
  #     3) <plot.dir>/<file.lead>_LV_boxplots.pdf
  
  `%>%` <- dplyr::`%>%`
  
  # error-handling
  # we need to join by "Sample" column for boxplots
  if (!("Sample" %in% colnames(sample.info.df))) {
    stop("'Sample' must be a column in sample.info.df")
  }
  # phenotype.col needs to be in column names
  if (!(phenotype.col %in% colnames(sample.info.df))) {
    stop("phenotype.col must be a column name in sample.info.df")
  }
  
  # initialize list to hold results to be returned
  return.list <- list()
  
  #### Differential Expression ####
  # get the named vector to use as the phenotype for testing differential 
  # expression
  phenotype.vector <- as.factor(make.names(sample.info.df[[phenotype.col]]))
  names(phenotype.vector) <- sample.info.df$Sample
  # test itself
  limma.df <- TestLVDifferencesJLF(b.matrix = b.matrix,
                                   phenotype = phenotype.vector,
                                   blocking= blocking,
                                   use.bonferroni = use.bonferroni)
  # write to file
  dlve.file <- file.path(results.dir, 
                         paste0(file.lead, "_limma_results.tsv"))
  readr::write_tsv(limma.df, path = dlve.file)
  # add to list to be returned
  return.list[["limma"]] <- limma.df
  
  #### Reshape & join with sample information ####
  b.df <- reshape2::melt(b.matrix)
  #print(b.df)
  colnames(b.df) <- c("LV", "Sample", "Value")
  sample.b.df <- dplyr::inner_join(b.df, sample.info.df, by = "Sample")
  #long.file <- file.path(results.dir, 
  #paste0(file.lead, "_B_long_sample_info.tsv"))
  #readr::write_tsv(sample.b.df, long.file)
  # add to list to be returned
  return.list[["b.df"]] <- sample.b.df
  
  #### Plotting ####
  plot.file <- file.path(plot.dir,
                         paste0(file.lead, "_boxplots.pdf"))
  
  # if we only want significant LVs plotted, filter sample.b.df using the
  # adj.P.Val cutoff sig.threshold
  if (significant.only) {
    sig.lvs <- limma.df$LV[which(limma.df$adj.P.Val < sig.threshold)]
    sample.b.df <- sample.b.df %>%
      dplyr::filter(LV %in% sig.lvs)
  }
  
  # BoxplotDiffLV(tidy.b.df = sample.b.df, 
  #  phenotype.column = phenotype.col,
  #   pdf.path = plot.file)
  
  
  return(return.list)
}

gene_reordering_testing <- function(seed_num, tpm_table, test_vector, groups){
  #input 
  #seed_num- the seed number 
  #tpm_table- the gene expression tpm normalized table
  #test_vector- the vector that matches the colnames in in the tpm table that say what group the samples belong in 
  #groups- the groups in the comparision 
  
  #output
  #noise_number- the number of latent variabels are different
  
  limma_list <- vector(mode = "list", length = 10)
  #print(test)
  set.seed(seed_num)
  
  for(i in 1:10){
    #print(groups[[i]])
    tpm_table_sub <- tpm_table[sample(1:nrow(tpm_table), size= groups[[i]]) ,]
    tpm_table_left<- tpm_table[! rownames(tpm_table) %in% rownames(tpm_table_sub ) ,]
    test_index <- sample(1:groups[[i]], size= groups[[i]])
    reorder_genes <-  rownames(tpm_table_sub)[test_index ]
    rownames(tpm_table_sub) <- reorder_genes
    tpm_table_v3 <- rbind(tpm_table_sub, tpm_table_left)
    test <- GetOrderedRowNorm(as.matrix(tpm_table_v3),recount3_plier)
    zeros<- rownames(test)[!complete.cases(test)]
    tpm_table_v4 <- tpm_table_v3[!rownames(tpm_table_v3) %in% zeros,]
    b.matrix <- GetNewDataB(as.matrix(tpm_table_v4), recount3_plier)
    test_res <- TestLVDifferences(b.matrix, test_vector ,use.bonferroni = TRUE)
    limma_list[[i]] <- test_res$LV[test_res$adj.P.Val < 0.05 & abs(test_res$logFC)> 0.05]
    #print("Done")
    
  }
  noise_number<- c()
  for(j in 1:10){
    noise_number[j] <- length(limma_list[[j]])
  }
  return(noise_number)    
}



#LINCS match the ids to gene ids
match_to_LINCS_genes <- function(gene_list){
  #gene_list- gene symbols to matching 
  
  #upset is the mapping to the lincs data
  gene_info_LINCS<- read.table("~/data/LINCS_210914_LEVEL3/geneinfo_beta.txt", sep="\t", header=TRUE)
  upset_v2<- c()
  for (i in 1:length(gene_list)){
    possibleError <- tryCatch(
      gene_info_LINCS$gene_id[gene_list[i] == gene_info_LINCS$gene_symbol],
      error=function(e) e
    )
    if(length(possibleError)> 0){
      #REAL WORK
      upset_v2[i]<- gene_info_LINCS$gene_id[gene_list[i] == gene_info_LINCS$gene_symbol]
      }else{
        upset_v2[i]<- NA
    }
  
  }
  names(upset_v2)<- gene_list
  return(upset_v2)
}

#calc the tau scores for the results
#au values were calculated based on the calculation conducted by the SignatureSearch package for a later version of LINCS(Duan et al. 2020). 
#A Tau value of -100 means that the negative connectivity score (ie., a more inverse signature) from the comparison between the perturbation 
#and disease signature was more negative than that same perturbation signature to other perturbation signatures in the LINCS 2020 database.
#A Tau value of 100 indicates that the positive connectivity score (ie., a more matching signature) from the comparison between the 
#perturbation and disease signature was more positive than that same perturbation signature to other perturbation signatures in the 
#LINCS 2020 database. Tau values close to 0 suggest that the overlap between disease signatures and perturbation signatures is not a unique 
#overlap for the perturbation. Overall, the Tau value helps identify drugs uniquely inverse to the disease signature to avoid off-target 
#and adverse events (Subramanian et al., 2017).

tau_scores<- function(results_df, tau_score_df){
  #input
  #results_df- this is results dataframe from LINCS that has been filtered to only include the cell of interest
  #tau_score_df - this is the dataframe of the connectivity score of the cell line preturbations profiles compared to other profiles
  
  #output
  #tau_list- is a vector of tau values for each of the results in the results_df
  
  results_df$id <- paste(results_df$pert, results_df$cell, results_df$type, sep= "__")
  #print(results_df$id)
  tau_list<- c()
  for (i in 1:nrow(results_df)){
    #print( results_df$id[i])
    #print(colnames(tau_score_df) == results_df$id[i])
    index <- grepl(results_df$pert[i], colnames(tau_score_df))
    scores <- tau_score_df[,colnames(tau_score_df) == results_df$id[i]]
    #print(scores)
    tmpsum <- sum( abs(scores) < abs(results_df$NCS[i]))
    secd<- (tmpsum/length(scores)) *100
    tau<- secd* sign( results_df$NCS[i])
    #print(tau)
    tau_list[i]<- tau
  }
  return(tau_list)
}


#This function is check the FDA-approval of candidates. The approved list comes from Drugs@FDA
FDA_APPROVAL_CHECK<- function(drug_list){
  #input 
  #drug_list- a vector of drugs 
  
  #output
  #return_list - a vector of true or false of drug list being FDA-approved 
  fda_info <- read.table("~/data/fda_product_info_df.csv" )
  fda_approve <- fda_info[fda_info$market_status_v2 %in% c("Over-the-Counter", "Prescription" ),]
  return_list <- c()
  for (i in 1:length(drug_list)){
    drug <- str_replace_all(drug_list[i], "[^[:alnum:]]", " ")
    drug <- str_squish(drug)
    test <- grep( drug, fda_approve$ActiveIngredient, ignore.case= TRUE)
    if (length(test)>0 ){
      return_list[i] <- TRUE
    }else{
      return_list[i] <- FALSE
    }
    #more than 5 character are needed for matching
    if (nchar(drug_list[i])< 6){
      return_list[i] <- FALSE
    }
  }
  names(return_list)<- drug_list
  return(return_list)
}

#this function detemines if a drug has been/is in a clinical trial. this data was download from clinicaltrials.gov
clinical_trial_check<- function(drug_list, clinical_trial_info){
  #input
  #drug_list- a vector of drugs 
  #clinical_trial_info- the disease specific table from  clinicaltrials.gov
  
  #output
  #return_list- a vector of true or false of drug list having a clinical trial presence 
  return_list <- c()
  for (i in 1:length(drug_list)){
    drug <- str_replace_all(drug_list[i], "[^[:alnum:]]", " ")
    drug <- str_squish(drug)
    test <- grep( drug, clinical_trial_info$Interventions, ignore.case= TRUE)
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


#this is a sig_filter function from SignatureSearchData. I adjsuted this function to work with the new LINCS data format. 
sig_filter_JLF<- function (meta, pert_type = "trt_cp", dose, time = "24 h") {
  #input
  #meta- this is metadata from LINCS 2020 data. 
  
  #output
  #meta- this is a filter metadata file for the subset of LINCS 2020
  meta %<>% dplyr::filter(pert_type == pert_type & pert_idose == 
                            dose & pert_itime == time)
  meta %<>% bind_cols(alt_id = paste(meta$cmap_name, meta$cell_iname, 
                                     sep = "_")) %>% bind_cols(pert_cell_factor = paste(meta$cmap_name, 
                                                                                        meta$cell_iname, meta$pert_type, sep = "__")) %>% distinct(alt_id, 
                                                                                                                                                   .keep_all = TRUE)
  return(meta)
}



#LINCS pathway & Drug targets

drug_analysis <- function(drug_name, target_list,  tpm_df, sample_vector, deseq_results, se, cell_lincs= "__GI1__trt_cp", result_path){
  #input
  #drug_name- character of the name of drug
  #target_list- the list of genes that are drug targets
  #tpm_df- a data frame with the transcript per million for the gene expresssion for tummor and control 
  #sample_vector- same order as the column names of the tpm_df, but the sample group (tumor or control)
  #deseq2_result- deseq2 result data frame with the symbol inculded for each gene
  #se- summarized experiment for LINCS level 5 data
  #cell_lincs- the cell line for the cancer in the following format "__GI1__trt_cp"
  #result_path- the path to the directory to store csv and plot images 
  
  #output
  #target different expression plots
  #pathway analysis for perturbation based on LINCS profiles
  
  #drug target expression
  tpm_sub <- tpm_df[rownames(tpm_df) %in% target_list,]
  tpm_sub$genes<- rownames(tpm_sub)
  tpm_sub_v2 <- tpm_sub %>%
    pivot_longer(!genes, names_to = "sample", values_to = "tpm")
  tpm_sub_v2$type <- rep(sample_vector, nrow(tpm_sub))
  name<- drug_name
  y_max<- max(tpm_sub_v2$tpm) +5
  
  
  deseq_sub <-  deseq_results[ deseq_results$Symbol %in% unique( tpm_sub$genes),]
  if (nrow(deseq_sub)>0 ){
  stat_table <- deseq_sub[,c(7,6)]
  stat_table$group1 <- rep("Primary Tumor" ,  nrow(deseq_sub))
  stat_table$group2 <- rep("Solid Tissue Normal" ,  nrow(deseq_sub))
  stat_table$padj <- ifelse(stat_table$padj < 0.05, formatC(stat_table$padj, format = "e", digits = 2), "")
  colnames(stat_table)<- c("genes", "padj" ,  "group1" ,"group2")
  
  
  file_name<- paste0(result_path, "/", drug_name, "_drug_target_expression_deseq2_padj.png")
  
  ggplot(tpm_sub_v2, aes(x=genes, y=tpm, color=type)) + geom_boxplot()+ stat_pvalue_manual(stat_table, x="genes" , label = "{padj}", y.position= y_max,  size = 3) + labs(title=name,x="Drug Targets", y = "Gene Expression (TPM)", color= "Sample Type") + scale_colour_manual(values =  c("#440154FF","#228C8DFF"), aesthetics = c("colour", "fill"))+theme_classic()
  ggsave(file_name, width = 20, height = 8)
  }
  #pathway enrichment
  lincs_name<- paste0(drug_name, cell_lincs)
  sub <- se[,grep(lincs_name, colnames(se))]
  up <- sub[sub >= 2  ]
  up<- names(up[order(-up)])
  down<- sub[ sub <= -2 ]
  down <- names(down[order(down)])
  #print(down)
  if (length(down) > 1 & length(up)> 1){
  
  set.seed(101)
  downset_pathway_results <- gost(query = down, 
                                  organism = "hsapiens", ordered_query = TRUE, 
                                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                                  user_threshold = 0.05, correction_method = "g_SCS", 
                                  domain_scope = "annotated", custom_bg = NULL, 
                                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
  downset_pathway_results<- downset_pathway_results$result
  
  set.seed(101)
  upset_pathway_results <- gost(query = up, 
                                organism = "hsapiens", ordered_query = TRUE, 
                                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                measure_underrepresentation = FALSE, evcodes = FALSE, 
                                user_threshold = 0.05, correction_method = "g_SCS", 
                                domain_scope = "annotated", custom_bg = NULL, 
                                numeric_ns = "", sources = NULL, as_short_link = FALSE)
  upset_pathway_results<- upset_pathway_results$result
  
  pathway_results<- rbind(upset_pathway_results, downset_pathway_results)
  pathway_results$set<- c(rep("Up", nrow(upset_pathway_results)), rep("Down", nrow(downset_pathway_results)))
  } 
  if (length(down) > 1 & length(up)  < 1){
    set.seed(101)
    downset_pathway_results <- gost(query = down, 
                                    organism = "hsapiens", ordered_query = TRUE, 
                                    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                    measure_underrepresentation = FALSE, evcodes = FALSE, 
                                    user_threshold = 0.05, correction_method = "g_SCS", 
                                    domain_scope = "annotated", custom_bg = NULL, 
                                    numeric_ns = "", sources = NULL, as_short_link = FALSE)
    downset_pathway_results<- downset_pathway_results$result
    pathway_results<-  downset_pathway_results
    pathway_results$set<-  rep("Down", nrow(downset_pathway_results))
  }
  if  (length(down) < 1 & length(up)  > 1){
    set.seed(101)
    upset_pathway_results <- gost(query = up, 
                                  organism = "hsapiens", ordered_query = TRUE, 
                                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                                  user_threshold = 0.05, correction_method = "g_SCS", 
                                  domain_scope = "annotated", custom_bg = NULL, 
                                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
    upset_pathway_results<- upset_pathway_results$result
    
    pathway_results<- upset_pathway_results
    pathway_results$set<- rep("Up", nrow(upset_pathway_results))
  }

  if ( length(down) > 1 | length(up)  > 1 ){
  file_name<- paste0(result_path, "/", drug_name, "_drug_pathways.csv")
  write.csv2(as.data.frame(pathway_results[,-14]), file_name )
  
  file_name<- paste0(result_path, "/", drug_name, "_drug_pathways.png")
  
  #making adjustments here
  #pathway_results$set <- factor(pathway_results$set, levels = c("up", "down"))
  
  if(nrow(pathway_results)<50){
    ggplot(pathway_results, aes(x=set, y=term_name, fill=p_value))+geom_tile()+scale_fill_viridis(direction=-1)+theme_classic()+ labs(title=name,x="Drug Gene Sets", y = "g:Profiler Gene Sets", fill= "p-value")
    
  }else{
    pathway_results<- pathway_results[pathway_results$source %in% c("KEGG","REAC","GO:MF" ),]
    ggplot(pathway_results, aes(x=set, y=term_name, fill=p_value)) + geom_tile()+ scale_fill_viridis(direction=-1)+ theme_classic() + labs(title=name,x="Drug Gene Sets", y = "g:Profiler Gene Sets", fill= "p-value")
  }
  ggsave(file_name, width = 12, height = 14)
  }
}

#clinical trials -testing
#To evaluate, if one of the methods was able to determine more candidates that were already in clinical trials for specific cancer, 
#a permutation test was conducted. For each cancer the drugs that were in the LINCS 2020 database, it was evaluated if it was FDA approved. 
#Then randomly the same number of drugs that were significant in the method being tested (i.e., DESeq2, limma, Transfer Learning) and 
#determine what fraction of the randomly selected drugs were in clinical trials for specific cancer. This was done 10,000 and a 
#one-tailed Wilcox test determined if the fraction of FDA-approved drugs that were higher than by randomly selected drugs.

clinical_trial_testing<- function(seed, cell_line,number_sig_drugs, fraction_clinical_drugs, CT_DATA){
  #input 
  #seed- the seed for the analysis for reproducible results
  #cell_line- name of cell line 
  #number_sig_drugs- the number of significant drugs from method
  #fraction_clinical_drugs- the fraction of drugs that are in clinical trial compared to all drug. Note this is decimal number
  #CT_DATA- the dataframe with the clinical trial data for the Cancer of interest.
  
  #Output
  #results for the one-tailed wilcox test 
  
  FDA_CT_table<- CT_DATA
  library(HDF5Array)
  lincs_samples<- HDF5Array("~/data/lincs_2020.h5", name="colnames")
  lincs_samples<- as.vector(lincs_samples)
  
  lincs_samples_GBM <- lincs_samples[grep(cell_line,lincs_samples )]
  library(stringr)
  gbm_lincs_drugs <- as.vector(str_split(lincs_samples_GBM,  "__", simplify = TRUE)[,1])
  
  gbm_lincs_drugs<- gbm_lincs_drugs[FDA_APPROVAL_CHECK(gbm_lincs_drugs)]
  set.seed(seed)
  fraction<-c()
  for (i in 1:1000){
    test <- gbm_lincs_drugs[ sample(x = 1:length(gbm_lincs_drugs),size = number_sig_drugs ,replace = FALSE)]
    value2 <- c()
    
    #tables<- list()
    for (k in 1:number_sig_drugs){
      res <- grep(test[k], FDA_CT_table$Interventions, ignore.case=TRUE )
      res <- as.logical(length(res) > 0 )
      value2[k] <- res
    }
    res_table <- cbind(test, value2)
    fraction[i]<- table(res_table[,2])["TRUE"]
    #tables[i]<- res_table
    
  }
  fraction <- ifelse(is.na(fraction), 0, fraction)
  fraction_adj <- fraction/number_sig_drugs

  plot_data<- as.data.frame(cbind(1:1000, fraction_adj))
  plot <- ggplot(plot_data, aes(x=fraction_adj)) + geom_histogram() + geom_vline(xintercept = fraction_clinical_drugs)
  
  print(plot)
  
  return(wilcox.test(fraction_adj,mu= fraction_clinical_drugs , alternative="less" ))

}

#Prism testing
# The PRISM study considered a cell line as sensitive to a treatment if the median-collapsed fold-change is less than 0.3. 
#To evaluate if one of the methods was able to identify a larger fraction of candidates that cancer cell lines (for specific cancer in question)
#were more sensitive than by random chance, a permutation test was conducted. For this, the same number of drugs were randomly selected. 
#The median of the fraction of cell lines sensitive (log2  fold change < 0.3) across the randomly selected drugs. A one-tail Wilcox 
#test was used to determine if the median fraction was higher than by random chance. 

PRISM_testing<- function(seed, cell_line, drugs, drug_senstive_precentage_all_drugs){
  #input 
  #seed- the seed for the analysis for reproducible results
  #cell_line- name of cell line 
  #drug_senstive_precentage_all_drugs- Prism results for the drug candidates signifcant in method
  
  #output
  #the results of the one-tail wilcox test 
  
  #check all drugs for fda approval
  library(HDF5Array)
  lincs_samples<- HDF5Array("~/data/lincs_2020.h5", name="colnames")
  lincs_samples<- as.vector(lincs_samples)
  
  lincs_samples_GBM <- lincs_samples[grep(cell_line,lincs_samples )]
  library(stringr)
  gbm_lincs_drugs <- as.vector(str_split(lincs_samples_GBM,  "__", simplify = TRUE)[,1])
  
  gbm_lincs_drugs<- gbm_lincs_drugs[FDA_APPROVAL_CHECK(gbm_lincs_drugs)]
  
  #get the median number of cell lines senestive to drugs for sigficant drug list 
  test_median <- median(drug_senstive_precentage_all_drugs$drug_fraction[ drug_senstive_precentage_all_drugs$drug_list %in% drugs ] )
  
  print(test_median)
  
  number_sig_drugs<- length(drugs)
  # do that again 10,000 times 
  set.seed(seed)
  fraction<-c()
  for (i in 1:10000){
    test <- gbm_lincs_drugs[ sample(x = 1:length(gbm_lincs_drugs),size = number_sig_drugs ,replace = FALSE)]
    fraction[i] <- median(drug_senstive_precentage_all_drugs$drug_fraction[ drug_senstive_precentage_all_drugs$drug_list %in% test ] )
    
  }
  #fraction <- ifelse(is.na(fraction), 0, fraction)
  #fraction_adj <- fraction
  #hist(fraction_adj)
  #print(fraction_adj)
  
  plot_data<- as.data.frame(cbind(1:10000, fraction))
  plot <- ggplot(plot_data, aes(x=fraction)) + geom_histogram() + geom_vline(xintercept = test_median)+ labs( x = "median fraction of sensitive cell lines") + annotate("text", x=(test_median + 0.1) , y=800, label= test_median)
  
  print(plot)
  
  return(wilcox.test(fraction, mu= test_median , alternative="less" ))
  #return(fraction_adj)
}

#create venn diagram for comparing genes, pathways, candidates, etc. 
venn_dia_methods <- function(limma_list, deseq2_list, tfl_list, file_name='~/output/liver_cancer/SR_liver_all_gene_venn_diagram.png' ){
  #input
  #limma_list, deseq2_list, tfl_list- list of the items for that method. order is important. 
  #file_name- path and name for the file with image of the venn_diagram
  
  #output- the venn diam image 
  
  myCol <- c("#440154FF" , "#31688EFF" ,"#35B779FF")
  venn.diagram(
    x = list(limma_list, deseq2_list, tfl_list),
    category.names = c("limma" , "DESeq2" , "Transfer Learning"),
    filename = file_name,
    output=TRUE,
    
    # Output features
    imagetype="png" ,
    height = 1000 , 
    width = 1000 , 
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
    #cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135),
    cat.dist = c(0.04, 0.04, 0.04),
    cat.fontfamily = "sans",
    rotation = 1
  )
}

#pathway analysis for the up and down gene signatures
method_up_down_gene_set_analysis<- function(up, down,result_path, method, name){
  #input
  #up- up gene vector
  #down- down gene vector 
  #result_path- the path to save the files 
  #mehtod- name of the method that found this signature
  #name- title for the plots
  
  #output
  #pathway results for both up and down genes 
  #saved as csv and a plot for quick veiw of results 
  
  library(gprofiler2)
  set.seed(101)
  downset_pathway_results <- gost(query = down, 
                                  organism = "hsapiens", ordered_query = TRUE, 
                                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                                  user_threshold = 0.05, correction_method = "g_SCS", 
                                  domain_scope = "annotated", custom_bg = NULL, 
                                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
  downset_pathway_results<- downset_pathway_results$result
  
  set.seed(101)
  upset_pathway_results <- gost(query = up, 
                                organism = "hsapiens", ordered_query = TRUE, 
                                multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                measure_underrepresentation = FALSE, evcodes = FALSE, 
                                user_threshold = 0.05, correction_method = "g_SCS", 
                                domain_scope = "annotated", custom_bg = NULL, 
                                numeric_ns = "", sources = NULL, as_short_link = FALSE)
  upset_pathway_results<- upset_pathway_results$result
  
  pathway_results<- rbind(upset_pathway_results, downset_pathway_results)
  pathway_results$set<- c(rep("Up", nrow(upset_pathway_results)), rep("Down", nrow(downset_pathway_results)))
  file_name<- paste0(result_path, "/", method, "_pathways.csv")
  write.csv2(as.data.frame(pathway_results[,-14]), file_name )
  
  file_name<- paste0(result_path, "/", method, "_pathways.png")
  
  #making adjustments here for case of more than 50 pathways enirched
  if(nrow(pathway_results)<50){
    ggplot(pathway_results, aes(x=set, y=term_name, fill=p_value))+geom_tile()+scale_fill_viridis(direction=-1)+theme_classic()+ labs(title=name,x="Drug Gene Sets", y = "g:Profiler Gene Sets", fill= "p-value")
    
  }else{
    pathway_results2<- pathway_results[pathway_results$source %in% c("KEGG","REAC","GO:MF" ),]
    ggplot(pathway_results2, aes(x=set, y=term_name, fill=p_value)) + geom_tile()+ scale_fill_viridis(direction=-1)+ theme_classic() + labs(title=name,x="Drug Gene Sets", y = "g:Profiler Gene Sets", fill= "p-value")
  }
  ggsave(file_name, width = 12, height = 14)
  return(as.data.frame(pathway_results[,-14]))
}

#functions from rrvgo
#rrvgo github (getGoTerm, loadOrgdb, getGoSize, reduceSimMatrix)
getGoTerm <- function(x) {
  sapply(x, function(x) tryCatch(GO.db::GOTERM[[x]]@Term, error=function(e) NA))
}
loadOrgdb <- function(orgdb) {
  if(!requireNamespace(orgdb, quietly = TRUE)) {
    stop("Bioconductor orgdb for ", orgdb, " not found. Consider installing it.",
         call. = FALSE)
  }
  eval(parse(text=paste0(orgdb, "::", orgdb)))
}
getGoSize <- function(terms, orgdb, keytype) {
  if(all(is(orgdb) != "OrgDb")) {
    orgdb <- loadOrgdb(orgdb)
  }
  
  # get all GO terms with genes associated
  go <- suppressMessages(
    AnnotationDbi::select(orgdb,
                          keytype=keytype,
                          columns=c("GO", "ONTOLOGY"),
                          keys=AnnotationDbi::keys(orgdb, keytype=keytype)))
  go <- go[!is.na(go$GO), ]
  go <- go[go$GO %in% terms, ]
  
  # count
  counts   <- table(go$GO)
  go <- go[go$GO %in% terms, ]
  empty    <- terms[!(terms %in% names(counts))]
  nocounts <- setNames(rep(0, length(empty)), empty)
  
  c(counts, nocounts)
}
reduceSimMatrix <- function (simMatrix, scores = NULL, threshold = 0.7, orgdb="org.Hs.eg.db", keytype = "ENTREZID") {
  if (!is.null(scores) && !all(rownames(simMatrix) %in% names(scores))) {
    stop("Scores vector does not contain all terms in the similarity matrix")
  }
  sizes <- getGoSize(rownames(simMatrix), orgdb, keytype)
  if (is.null(scores)) {
    message("No scores provided. Falling back to term's size")
    scores <- sizes
  }
  scores <- scores[match(rownames(simMatrix), names(scores))]
  orows <- match(rownames(simMatrix), names(scores))
  ocols <- match(colnames(simMatrix), names(scores))
  simMatrix <- simMatrix[orows, ocols]
  o <- rev(order(scores, sizes, na.last = FALSE))
  simMatrix <- simMatrix[o, o]
  cluster <- cutree(hclust(as.dist(1 - simMatrix)), h = threshold)
  clusterRep <- tapply(rownames(simMatrix), cluster, function(x) x[which.max(scores[x])])
  data.frame(go = rownames(simMatrix), cluster = cluster, parent = clusterRep[cluster], 
             parentSimScore = unlist(Map(seq_len(nrow(simMatrix)), 
                                         clusterRep[cluster], f = function(i, j) simMatrix[i, 
                                                                                           j])), score = scores[match(rownames(simMatrix), 
                                                                                                                      names(scores))], size = sizes[match(rownames(simMatrix), 
                                                                                                                                                          names(sizes))], term =getGoTerm(rownames(simMatrix)), 
             parentTerm = getGoTerm(clusterRep[cluster]))
}

#compare disease signatures between the different methods with go term semantic similarity based on go biological go terms and plot as a heatmap
go_term_heatmap<- function(limma_pathways, deseq2_pathways, TFL_pathways, list_type, threshold = 0.95, file_name){
  #input
  #limma_pathways, deseqq2_pathways, TFL_pathways- the gprofiler result data frame for each method with the last column (named set)
  #indicating if the pathway was from the up or down group
  #list_type- up or down regulated 
  #threshold- is the threshold of how far up or down the go ontology tree to go. default to 0.95
  #file name for the parent go term results
  
  #output 
  #heatmap
  
  #get terms
  limma_bp<- limma_pathways$term_id[ limma_pathways$source == "GO:BP" & limma_pathways$set == list_type]
  deseq2_bp <- deseq2_pathways$term_id[ deseq2_pathways$source == "GO:BP" & deseq2_pathways$set == list_type]
  TFL_bp <- TFL_pathways$term_id[ TFL_pathways$source == "GO:BP" & TFL_pathways$set == list_type]
  
  #run the go term semantic similarity 
  go1 <-  unique(c(TFL_bp, deseq2_bp, limma_bp ))
  go_gbm_up_sim <- mgoSim(go1, go1, semData=hsGO, measure="Wang", combine=NULL)
  
  #get the parent terms from rrvgo for the pathways
  res <- reduceSimMatrix(go_gbm_up_sim, threshold = threshold)
  res_v2<- res[match(colnames(go_gbm_up_sim), res$go),]
  
  
  
  #determine which pathway is enriched in the different methods for the heatmap annotation
  tfl_drugs <- grepl(paste(TFL_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  deseq2_drugs<- grepl(paste(deseq2_bp,collapse="|"), colnames(go_gbm_up_sim))
  limma_drugs <- grepl(paste(limma_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  
  #handing the case where there are no pathways enriched
  if(length(deseq2_bp) ==0){ deseq2_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  if(length(TFL_bp) ==0){ tfl_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  if(length(limma_bp) ==0){ limma_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  
  #pick the number of colors note limit is about 9 
  bp_color<- brewer.pal(n = length(unique(res_v2$parentTerm)), name = "Paired")
  names(bp_color)<- unique(res_v2$parentTerm)
  
  #create heatmap
  row_ha = HeatmapAnnotation(Transfer_Learning=tfl_drugs,DESeq2= deseq2_drugs, limma= limma_drugs,GO_BP_Group= res_v2$parentTerm , col = list(Transfer_Learning = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),DESeq2 = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),limma = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"), GO_BP_Group= bp_color ), annotation_name_gp= gpar(fontsize = 16,  fontface = "bold"), 
                             annotation_legend_param = list(limma = list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold")), 
                                                            Transfer_Learning = list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold")), 
                                                            DESeq2= list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold")),
                                                            GO_BP_Group = list(title_gp = gpar(fontsize = 15, fontface = "bold"), labels_gp = gpar(fontsize = 13, fontface = "bold"))))
  col_fun = colorRamp2(c(0,  1), c( "black", "yellow"))
  
  #SAVE THE PARENT MATCHES
  res_v3<- cbind(res_v2,tfl_drugs, deseq2_drugs, limma_drugs )
  write.csv(res_v3, file_name)
  heatmap <- Heatmap(go_gbm_up_sim, nam= "GO Term Similarity (Wang)", col = col_fun, show_column_names = FALSE,  show_row_names = FALSE, top_annotation = row_ha,  
          clustering_distance_rows= "euclidean",
          clustering_distance_columns=  "euclidean",
          clustering_method_rows = "ward.D2" ,
          clustering_method_columns="ward.D2", heatmap_legend_param = list( title_gp = gpar(fontsize = 15 , fontface = "bold"), labels_gp = gpar(fontsize = 12, fontface = "bold")))
  draw(heatmap, annotation_legend_side = "bottom")
  
  column_dend = hclust(dist(t(go_gbm_up_sim)), method = "ward.D2")
  
  heatmap_v2 <- Heatmap(matrix(nrow = 0, ncol = ncol(go_gbm_up_sim)),
          cluster_columns = column_dend, 
          show_column_names = FALSE,  
          show_row_names = FALSE, 
          top_annotation = row_ha) 
  draw(heatmap_v2, annotation_legend_side = "bottom")
  
}

#compare disease signatures between the different methods with go term semantic similarity based on go biological go terms and plot a tree plot
go_term_tree_plot<- function(limma_pathways, deseq2_pathways, TFL_pathways, list_type, num=80, k_num=5, colors= c( "#117733", "#661100",  "#0072B2", "#D55E00", "#AA4499"), shift = 4){
  #input
  #limma_pathways, deseqq2_pathways, TFL_pathways- the gprofiler result data frame for each method with the last column (named set)
  #indicating if the pathway was from the up or down group
  #list_type- up or down regulated 
  #num= the number of pathways to look at 
  #k_num- number of groups
  #colors- colors for the groups based of clustering
  
  #output 
  #treeplot
  
  
  limma_bp<- limma_pathways$term_id[ limma_pathways$source == "GO:BP" & limma_pathways$set == list_type]
  deseq2_bp <- deseq2_pathways$term_id[ deseq2_pathways$source == "GO:BP" & deseq2_pathways$set == list_type]
  TFL_bp <- TFL_pathways$term_id[ TFL_pathways$source == "GO:BP" & TFL_pathways$set == list_type]
  
  go1 <-  unique(c(TFL_bp, deseq2_bp, limma_bp ))
  
  
  go_gbm_up_sim <- mgoSim(go1, go1, semData=hsGO, measure="Wang", combine=NULL)
  
  limma_bp_up_term <- limma_pathways$term_name[ limma_pathways$source == "GO:BP" & limma_pathways$set == list_type]
  deseq2_bp_up_term <- deseq2_pathways$term_name[ deseq2_pathways$source == "GO:BP" & deseq2_pathways$set == list_type]
  TFL_bp_up_term <- TFL_pathways$term_name[ TFL_pathways$source == "GO:BP" & TFL_pathways$set == list_type]
  
  list <- unique(c(TFL_bp_up_term, deseq2_bp_up_term, limma_bp_up_term ))
  test_names <- paste0(strtrim(list, 25), " (", go1, ")")
  #print(test_names[1:5])
  #print(length(test_names))
  #dim(go_gbm_up_sim)
  colnames(go_gbm_up_sim)<- test_names 
  rownames(go_gbm_up_sim)<- test_names 
  
  med_up <- colMedians(go_gbm_up_sim)
  med_top<- colnames(go_gbm_up_sim)[order(-med_up)[1:num]]
  
  med_bottom<- test_names[! test_names %in% med_top]
  
  go_gbm_up_sim_sub <- go_gbm_up_sim[colnames(go_gbm_up_sim) %in% med_top, colnames(go_gbm_up_sim) %in% med_top]
  
  cluster_sub <- hclust(dist(go_gbm_up_sim_sub), method= "ward.D2")
  dend_sub  = as.dendrogram(cluster_sub )
  
  tfl_drugs <- grepl(paste(TFL_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  deseq2_drugs<- grepl(paste(deseq2_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  limma_drugs <- grepl(paste(limma_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  
  if(length(deseq2_bp) ==0){ deseq2_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  if(length(TFL_bp) ==0){ tfl_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  if(length(limma_bp) ==0){ limma_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  
  tfl_drugs2 <- tfl_drugs[ test_names %in% med_top]
  limma_drugs2 <- limma_drugs[ test_names %in% med_top]
  deseq2_drugs2<- deseq2_drugs[ test_names %in% med_top]
  
  test<- ifelse(tfl_drugs2== TRUE, "#440154FF", "#228C8DFF")
  print(test)
  test2<- ifelse(limma_drugs2== TRUE, "#440154FF", "#228C8DFF")
  test3<- ifelse(deseq2_drugs2== TRUE, "#440154FF", "#228C8DFF")
  #test<- ifelse(tfl_drugs== TRUE, "#440154FF", "#228C8DFF")
  #test2<- ifelse(limma_drugs== TRUE, "#440154FF", "#228C8DFF")
  #test3<- ifelse(deseq2_drugs== TRUE, "#440154FF", "#228C8DFF")
  
  par(oma=c(0.5,0.5,0.5,20), mar = c(4, 0.01, 0.01, 12))
  #dend2 <-prune( dend, med_bottom, reindex_dend=FALSE)
  dend_sub %>%
    set("labels_col",value =colors, k=k_num) %>%
    set("branches_k_color",value =colors,  k = k_num) %>%
    set("leaves_pch", 19)  %>% 
    set("nodes_cex", 0.6) %>% 
    set("branches_lwd", 3) %>% 
    set("labels_cex", 0.5) %>%
    plot( horiz=TRUE, axes=FALSE)
  #addjust eh y shift
  colored_bars(cbind(test, test2, test3), dend_sub ,rowLabels = c("Transfer", "limma", "DESeq2"), y_shift = shift, horiz = TRUE)
  
}

#plotting drug target and mechanism of action for all candidate results
drug_moa_target_plotting<- function(fda_approved_res, lincs_commpound_info=lincs_commpound_info, method, cancer, file_path_header ){
  #fda_approved_res- the signaturesearch results with only fda approved
  #lincs_commpound_info- the lincs compound info from the LINCS database
  #method- what method found these candidates
  #cancer- which cancer was this result for 
  #file_path_header- path and header for output files
  
  #output
  #several plots including alluvial and bar plots for the candidates
  file_starter<- paste0(file_path_header, "_", cancer, "_", method, "_", sep= "")
  
  test_fda <- fda_approved_res
  lincs_commpound_info[lincs_commpound_info$cmap_name %in% test_fda$pert,]
  
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
  
  summ_df<- summ_df[!duplicated(summ_df), ]
  
  TEST_LODES<- to_lodes_form(as.data.frame(summ_df),
                             axes = 1:3,
                             id = "INDEX")
  
  TEST_LODES$stratum = str_wrap(TEST_LODES$stratum, width = 20)
  TEST_LODES$x <- factor(TEST_LODES$x, levels= c("moa", "cmap_name", "target"))
  
  TEST_LODES$drug <- rep(NA, nrow(summ_df))
  drugs<- unique(summ_df$cmap_name)
  for (i in 1:length(drugs)){
    ind <- TEST_LODES$INDEX[TEST_LODES$stratum == drugs[i]]
    TEST_LODES$drug<- ifelse(TEST_LODES$INDEX %in% ind, drugs[i], TEST_LODES$drug)
  }
  
  #not need for this analysis but many for future.
  TEST_LODES$drug<- factor(TEST_LODES$drug, levels= test_fda$pert)
  
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
  
  
  counts_moa_v2<- counts_moa[counts_moa$Freq >0,]
  counts_moa_v2$Var1<- factor(counts_moa_v2$Var1, levels= counts_moa_v2$Var1[order(-counts_moa_v2$Freq)])
  if (nrow(counts_moa_v2) < 20){
    ggplot(counts_moa_v2, aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Mechanism of Action") + xlab("Number of Drugs") 
    
    file_name<- paste0(file_starter, "_", "barplot_moa.png")
    ggsave(file_name, width = 8, height=10, units= "in")
  }else{
    counts_moa_v2<- counts_moa_v2[order(-counts_moa_v2$Freq),]
    ggplot(counts_moa_v2[1:20,], aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Mechanism of Action") + xlab("Number of Drugs") 
    
    file_name<- paste0(file_starter, "_", "barplot_moa.png")
    ggsave(file_name, width = 8, height=10, units= "in")
  }
  
  
  counts$Var1<- factor(counts$Var1, levels= counts$Var1[order(-counts$Freq)])
  if( nrow(counts)< 20 ){
    #counts$Var1<- factor(counts$Var1, levels= counts$Var1[order(-counts$Freq),])
    ggplot(counts, aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Targets") + xlab("Number of Drugs") 
    file_name<- paste0(file_starter, "_", "barplot_target.png")
    ggsave(file_name, width = 8, height=10, units= "in")
  }else{
    
    counts<- counts[order(-counts$Freq),]
    
    ggplot(counts[1:20,], aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Targets") + xlab("Number of Drugs") 
    file_name<- paste0(file_starter, "_", "barplot_target.png")
    ggsave(file_name, width = 8, height=10, units= "in")
    
    
    
  }
  #add ploting function
  target_count <- counts
  
  
  #reduce to the top 20 drugs
  summ_df<- summ_df[summ_df$cmap_name %in% test_fda$pert[1:20],]
  
  summ_df_moa<- summ_df[,c(1,3)]
  summ_df_moa<- summ_df_moa[!duplicated(summ_df_moa), ]
  TEST_LODES_moa<- to_lodes_form(as.data.frame(summ_df_moa),
                                 axes = 1:2,
                                 id = "INDEX")
  
  TEST_LODES_moa$stratum = str_wrap(TEST_LODES_moa$stratum, width = 20)
  
  TEST_LODES_moa$drug <- rep(NA, nrow(summ_df_moa))
  drugs<- unique(summ_df_moa$cmap_name)
  for (i in 1:length(drugs)){
    ind <- TEST_LODES_moa$INDEX[TEST_LODES_moa$stratum == drugs[i]]
    TEST_LODES_moa$drug<- ifelse(TEST_LODES_moa$INDEX %in% ind, drugs[i], TEST_LODES_moa$drug)
  }
  
  #not need for this analysis but many for future.
  TEST_LODES_moa$drug<- factor(TEST_LODES_moa$drug, levels= test_fda$pert)
  
  #moa 
  moa_df <- TEST_LODES_moa[TEST_LODES_moa$x == "moa", 2:4]
  moa_df<- moa_df[!duplicated(moa_df), ]
  counts<- as.data.frame(table(moa_df$stratum))
  
  #counts<- as.data.frame(table(TEST_LODES_moa$stratum[TEST_LODES_moa$x == "moa"]))
  moa_order <- counts$Var1[ order(-counts$Freq)]
  counts_moa<- counts
  
  #factor the stratum column
  TEST_LODES_moa$stratum<- factor(TEST_LODES_moa$stratum, levels= unique(c(test_fda$pert, as.character(moa_order))))
  
  TEST_LODES_moa_com <- TEST_LODES_moa[ complete.cases(TEST_LODES_moa),]
  #TEST_LODES_com$index2 <- as.numeric(TEST_LODES_com$drug)
  #TEST_LODES_com_v2 <- TEST_LODES_com[!TEST_LODES_com$x == "cmap_name",]
  
  ggplot(TEST_LODES_moa_com,
         aes(x = x, stratum = stratum, alluvium = INDEX, label= stratum))+ geom_flow( aes(fill= drug), stat = "alluvium", color= "black") +
    geom_stratum(aes(fill= drug)) +geom_label(stat = "stratum", fill="white",   size= 5) + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.position="none", axis.text=element_text(size=14)) + scale_y_continuous(breaks=NULL) + scale_x_discrete(labels=c("cmap_name" = "Drug", "moa" = "Mechanism of Action")) +
    scale_fill_manual(
      values=viridis(n=length(drugs)),
      breaks=drugs,
      labels=drugs,
      na.value = NA
    )
  file_name<- paste0(file_starter, "_", "alluvial_moa.png")
  ggsave(file_name, width =10, height=20, units= "in")
  
  #add saving plot function
  
  summ_df_tar<- summ_df[,c(1,2)]
  #keep only the top targets
  summ_df_tar<- summ_df_tar[summ_df_tar[,2] %in% target_count$Var1[1:20], ]
  
  summ_df_tar<- summ_df_tar[!duplicated(summ_df_tar), ]
  TEST_LODES_tar<- to_lodes_form(as.data.frame(summ_df_tar),
                                 axes = 1:2,
                                 id = "INDEX")
  
  TEST_LODES_tar$stratum = str_wrap(TEST_LODES_tar$stratum, width = 20)
  
  TEST_LODES_tar$drug <- rep(NA, nrow(summ_df_tar))
  drugs<- unique(summ_df_tar$cmap_name)
  for (i in 1:length(drugs)){
    ind <- TEST_LODES_tar$INDEX[TEST_LODES_tar$stratum == drugs[i]]
    TEST_LODES_tar$drug<- ifelse(TEST_LODES_tar$INDEX %in% ind, drugs[i], TEST_LODES_tar$drug)
  }
  
  #not need for this analysis but many for future.
  TEST_LODES_tar$drug<- factor(TEST_LODES_tar$drug, levels= test_fda$pert)
  
  #targets
  target_df <- TEST_LODES[TEST_LODES$x == "target", 2:4]
  target_df<- target_df[!duplicated(target_df), ]
  counts<- as.data.frame(table(target_df$stratum))
  #counts
  target_order <- counts$Var1[ order(-counts$Freq)]
  
  #factor the stratum column
  TEST_LODES_tar$stratum<- factor(TEST_LODES_tar$stratum, levels= unique(c(test_fda$pert, as.character(target_order))))
  
  TEST_LODES_tar_com <- TEST_LODES_tar[ complete.cases(TEST_LODES_tar),]
  #TEST_LODES_com$index2 <- as.numeric(TEST_LODES_com$drug)
  #TEST_LODES_com_v2 <- TEST_LODES_com[!TEST_LODES_com$x == "cmap_name",]
  
  ggplot(TEST_LODES_tar_com,
         aes(x = x, stratum = stratum, alluvium = INDEX, label= stratum))+ geom_flow( aes(fill= drug), stat = "alluvium", color= "black") +
    geom_stratum(aes(fill= drug)) +geom_label(stat = "stratum", fill="white",  size= 8) + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.position="none", axis.text=element_text(size=14)) + scale_y_continuous(breaks=NULL) + scale_x_discrete(labels=c("cmap_name" = "Drug", "target" = "Targets")) +
    scale_fill_manual(
      values=viridis(n=length(drugs)),
      breaks=drugs,
      labels=drugs,
      na.value = NA
    )
  file_name<- paste0(file_starter, "_", "alluvial_target.png")
  ggsave(file_name, width = 10, height=20, units= "in")
}

#from the lsa r package (package not installed with docker image)
cosine <- function (x, y = NULL) {
  if (is.matrix(x) && is.null(y)) {
    co = array(0, c(ncol(x), ncol(x)))
    f = colnames(x)
    dimnames(co) = list(f, f)
    for (i in 2:ncol(x)) {
      for (j in 1:(i - 1)) {
        co[i, j] = cosine(x[, i], x[, j])
      }
    }
    co = co + t(co)
    diag(co) = 1
    return(as.matrix(co))
  }
  else if (is.vector(x) && is.vector(y)) {
    return(crossprod(x, y)/sqrt(crossprod(x) * crossprod(y)))
  }
  else if (is.vector(x) && is.matrix(y)) {
    co = vector(mode = "numeric", length = ncol(y))
    names(co) = colnames(y)
    for (i in 1:ncol(y)) {
      co[i] = cosine(x, y[, i])
    }
    return(co)
  }
  else {
    stop("argument mismatch. Either one matrix or two vectors needed as input.")
  }
}

candidate_cosine_data_heatmap<- function(lincs=se_assay, cell_line, deseq_drugs, limma_drugs, tfl_drugs ){
  #input
  #lincs- the summarized experiment for the lincs profiles (level 5)
  #cell_line- the name of the cell line
  #deseq_drug, limma_drugs, tfl_drugs - the list of drugs for the different methods
  
  #output
  # the cosine similarity matrix for the drug candidates, row and column should match
  
  se_assay_cell <- lincs[, grep( cell_line , colnames(lincs)) ]
  
  all_drugs<- c(deseq_drugs, limma_drugs, tfl_drugs)
  se_assay_cell_cand<- se_assay_cell[, grep(paste(all_drugs,collapse="|"), colnames(se_assay_cell))]
  
  cos <- cosine(as.matrix(se_assay_cell_cand))
  
  tfl_drugs <- grepl(paste(tfl_drugs,collapse="|"), colnames(cos)) 
  deseq2_drugs<- grepl(paste(deseq_drugs,collapse="|"), colnames(cos)) 
  limma_drugs <- grepl(paste(limma_drugs,collapse="|"), colnames(cos)) 
  
  row_ha = HeatmapAnnotation(Transfer_Learning=tfl_drugs,DESeq2= deseq2_drugs, limma= limma_drugs , col = list(Transfer_Learning = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),DESeq2 = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),limma = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF") ))
  id<- paste("__", cell_line, sep = "")
  rownames(cos) <- str_split(rownames(cos), id, simplify = TRUE)[,1]
  heatmap<- Heatmap(cos , nam= "Cosine Similarity of LINCS Profiles",  col = col_fun, show_column_names = FALSE, show_row_names = FALSE, top_annotation = row_ha, 
                    clustering_distance_rows= "euclidean",
                    clustering_distance_columns=  "euclidean",
                    clustering_method_rows = "ward.D2" ,
                    clustering_method_columns="ward.D2")
  print(heatmap)
  return(cos)
}


candidate_network_analysis<- function( cos_matrix, deseq_cand, limma_cand, tfl_cand, quantile= 9, header){
  #input
  #cos_matrix-  the cosine similarity matrix for the drug candidates, row and column should match
  #deseq_cand, limma_cand, tfl_cand - the list of drugs for the different methods
  #quantile- a number that indicates what quantile to cuttoff. oNly edges above this quantile will used in the plotting of the netwoork.
  # default is 9 so only the top 10% of edges are plotted. 
  #header- the file path and header for saving plots. 
  
  #output
  #html files with networks (one with lables and another without)
  #heatmap with the leiden communities and which candidates are in each. 
  
  #make the self match equal to zero
  cos_matrix[cos_matrix== 1]<- 0
  
  #add 1 to quantile because of 0 at index 1
  quantile <-quantile+1
  cutoff <- quantile(cos_matrix,probs = seq(0, 1, 0.1))[quantile]
  
  #set values below value to 0 
  cos_matrix[cos_matrix < cutoff ]<- 0 
  colnames(cos_matrix)<- rownames(cos_matrix)
  
  #remove candidate with only zero edge weights
  cos_matrix_v2<- cos_matrix[ !colSums(cos_matrix) == 0, !colSums(cos_matrix) == 0]
  
  #create network
  cos_cand_igraph <- graph_from_adjacency_matrix(cos_matrix_v2 ,mode= "undirected", weighted=TRUE)
  
  #cluster analysis 
  #leiden_clusters <- groups(cluster_leiden(cos_cand_igraph,n_iterations = 10,  objective_function ="modularity"))
  leiden_clusters <- cluster_leiden(cos_cand_igraph,n_iterations = 10,  objective_function ="modularity")
  
  #determine which candidate is in which cluster
  com <-rep(NA, nrow(cos_matrix_v2))
  for (i in 1: nrow(cos_matrix_v2)){
    com<- ifelse(colnames(cos_matrix_v2) %in% unlist(leiden_clusters[i]), i , com )
  }
  com
  print(com)
  
  tfl_drugs <- grepl(paste(tfl_cand,collapse="|"), colnames(cos_matrix_v2)) 
  deseq2_drugs<- grepl(paste(deseq_cand,collapse="|"), colnames(cos_matrix_v2))
  limma_drugs <- grepl(paste(limma_cand,collapse="|"), colnames(cos_matrix_v2))
  
  #shared by all methods or two methods
  all_drugs <- tfl_drugs & deseq2_drugs & limma_drugs
  
  deseq_tfl <- tfl_drugs & deseq2_drugs 
  deseq_tfl<- ifelse(all_drugs ==TRUE, FALSE, deseq_tfl)
  
  deseq_limma <- limma_drugs & deseq2_drugs 
  deseq_limma <- ifelse(all_drugs ==TRUE, FALSE, deseq_limma )
  
  limma_tfl<- limma_drugs & tfl_drugs 
  limma_tfl <- ifelse(all_drugs ==TRUE, FALSE, limma_tfl )
  #only indivudal methods
  tfl<- setdiff(tfl_cand,  c(limma_cand, deseq_cand) )
  tfl<- grepl(paste(tfl,collapse="|"), colnames(cos_matrix_v2))
  
  limma<- setdiff(limma_cand,  c(tfl_cand, deseq_cand) )
  limma<- grepl(paste(limma,collapse="|"), colnames(cos_matrix_v2))
  
  deseq<- setdiff(deseq_cand,  c(tfl_cand,limma_cand) )
  deseq<- grepl(paste(deseq,collapse="|"), colnames(cos_matrix_v2))
  
  #detemine colors for each drug node
  color <- rep(NA, length(tfl))
  color <- ifelse(tfl, "#440154FF", color)
  color <- ifelse(limma_tfl, "#443A83FF", color)
  color<- ifelse(limma, "#31688EFF", color)
  color<- ifelse(deseq_limma, "#21908CFF", color)
  color<- ifelse(deseq, "#35B779FF", color)
  color<- ifelse(deseq_tfl, "#8FD744FF", color)
  color<- ifelse(all_drugs, "#FDE725FF", color)
  #set colors for nodes and edges
  V(cos_cand_igraph)$color<- color
  E(cos_cand_igraph)$color <- "gray"
  
  #set groups for each node
  #this will match node color
  group <- rep(NA, length(tfl))
  group <- ifelse(tfl, "Transfer Learning", group)
  group <- ifelse(limma_tfl, "limma and Transfer Learning", group)
  group<- ifelse(limma, "limma", group)
  group<- ifelse(deseq_limma, "DESeq2 and limma", group)
  group<- ifelse(deseq, "DESeq2", group)
  group<- ifelse(deseq_tfl, "DESeq2 and Transfer Learning", group)
  group<- ifelse(all_drugs, "All Methods", group)
  
  #ploting with labels
  filename<- paste(header, "_network_with_labels.html")
  test.visn <- toVisNetworkData(cos_cand_igraph)
  test.visn$edges$value <- test.visn$edges$weight
  test.visn$nodes$group<- group
  test.visn$nodes$value<- rep(40, nrow(test.visn$nodes))
  lnodes <- data.frame(id = 1:7, label = c("Transfer Learning","limma and \nTransfer Learning", "limma", "DESeq2 and limma" ,  "DESeq2", "DESeq2 and \nTransfer Learning", "All Methods"), color = viridis(n=7), shape= "dot") 
  network <- visNetwork(test.visn$nodes, test.visn$edges) %>% 
    visIgraphLayout(layout = "layout.fruchterman.reingold")%>%
    visNodes(font = list(size = 25))
  visLegend(network, useGroups = FALSE, addNodes = lnodes, width= 0.3, zoom= TRUE) %>% visSave( file = filename)
  
  
  
  #ploting network with no labels 
  filename<- paste(header, "_network_without_labels.html")
  test.visn$nodes$label <- NULL
  network <- visNetwork(test.visn$nodes, test.visn$edges) %>% 
    visIgraphLayout(layout = "layout.fruchterman.reingold")%>%
    visNodes(font = list(size = 25))
  visLegend(network, useGroups = FALSE, addNodes = lnodes, width= 0.3, zoom= TRUE) %>% visSave( file = filename)
  
  #data wrangling the heatmap 
  names(com)<- colnames(cos_matrix_v2)
  
  method_list<- c("tfl_drugs", "deseq2_drugs", "limma_drugs")
  method<- c("Transfer\nLearning",  "DESeq2", "limma")
  current_method <- get(method_list[1])
  counts <- as.data.frame(table(com[current_method]))
  counts$method <- rep(method[1], nrow(counts))
  counts_final<- counts
  for(i in 2:3){
    current_method <- get(method_list[i])
    counts <- as.data.frame(table(com[current_method]))
    counts$method <- rep(method[i], nrow(counts))
    counts_final<- rbind(counts_final, counts)
  }
  counts_wider<- as.data.frame(pivot_wider(counts_final, names_from = method, values_from = Freq))
  rownames(counts_wider)<- counts_wider$Var1
  counts_wider<- counts_wider[,-1]
  counts_wider[is.na(counts_wider)]<- 0
  col_fun_v2 = colorRamp2(c( 0,max(counts_final$Freq)), c("white", "#31688EFF"))
  
  heatmap<- Heatmap(t(counts_wider), nam= "Number of Drugs from Method \nin Leiden Communities", col= col_fun_v2, #cell_fun = white.line,
                    clustering_distance_rows= "euclidean",
                    clustering_distance_columns=  "euclidean",
                    clustering_method_rows = "ward.D2" ,
                    clustering_method_columns="ward.D2",column_names_rot = 45, column_names_gp = grid::gpar(fontsize = 30, fontface= "bold"),
                    row_names_gp = grid::gpar(fontsize = 20, fontface= "bold"),
                    layer_fun = function(j, i, x, y, width, height, fill) {
                      v = pindex(t(counts_wider), i, j)
                      grid.text(sprintf("%.0f", v), x, y, gp = gpar(fontsize = 30, frontface= "bold"))
                      if(min(i)==1) {
                        grid.rect(gp = gpar(lwd = 2, fill = "transparent",col="black"))
                      }})
  draw(heatmap,legend_title_gp = gpar(fontsize = 30, fontface = "bold"))
}

#clinical trial summary bar plots
clinical_trial_summary_bar_plot <- function(summary_table, cancer, file ){
  #input
  #summary_table- table with infromation about how many candidates are in clinical trial vs not 
  #cancer - name of cancer
  #file- the file name to save the output bar plot 
  
  clinical_trial_gbm <- summary_table[summary_table$cancer == cancer,]
  t <- ggplot(clinical_trial_gbm ,aes(x=method,y=value,fill=type, label= value))+
    geom_bar(stat = "identity",color="white")+
    #facet_wrap(~cancer,ncol=1) +                                                              # Add values on top of bars
    geom_text(size = 10, position = position_stack(vjust = 0.5)) + scale_fill_manual(values= c("#FDE725FF","#21908CFF"), name="Legend") + ylab("Number of Drug Candidates")
  #dodger = position_dodge(width = 0.9)
  t + geom_text(aes(label=labels_v2),  position= position_stack(vjust = 1),  size= 20, face= "bold") + theme(text = element_text(size = 25,  face="bold"))
  ggsave(file, width= 12, height=7)
}

#Prism summary plot 
PRISM_methods_plotting<- function(deseq2_prism_path, tfl_prism_path, limma_prism_path, filename){
  #input
  #file paths to the prism data for each method. order is important
  #filename- the name of the file to save
  
  prism_deseq2_candidates <- read.csv( deseq2_prism_path)
  prism_tfl_candidates <- read.csv( tfl_prism_path)
  prism_limma_candidates <- read.csv(  limma_prism_path)
  
  #combine and remove duplicates and NAs
  prism_candidates<- rbind(prism_deseq2_candidates,prism_tfl_candidates ,  prism_limma_candidates)
  prism_candidates<- prism_candidates[!duplicated(prism_candidates),]
  prism_candidates <- prism_candidates[!is.na(prism_candidates$log_fold_change),]
  
  #create dataframe with drug medians for each drug and note methood 
  drug_list <- c(unique(prism_deseq2_candidates$Var2), unique(prism_tfl_candidates$Var2), unique(prism_limma_candidates$Var2))
  
  drug_median<- c()
  for (i in 1:length(drug_list)){
    cells <- prism_candidates$log_fold_change[prism_candidates$Var2 %in% drug_list[i]]
    drug_median[i] <- median(cells)
  }
  drug_prism_median_all_drugs <- data.frame( drug_list, drug_median)
  methods<- c(rep("DESeq2", length(unique(prism_deseq2_candidates$Var2))),  rep("Transfer Learning", length(unique(prism_tfl_candidates$Var2))),  rep("limma", length(unique(prism_limma_candidates$Var2))))
  drug_prism_median_all_drugs$methods<- methods
  
  #plot the results 
  t <- ggplot(drug_prism_median_all_drugs ,aes(x=methods, y=drug_median ,fill=methods))+
    #geom_point()+
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "grey", color = "black") + scale_fill_manual(values= c("#440154FF" ,"#21908CFF" ,"#FDE725FF"), name="Methods")  +
    ylab("Median log2fold change (PRISM)") + geom_hline(yintercept=0.3, linetype="dashed", color = "red") +theme(text = element_text(size = 25,  face="bold"))
  #t + geom_text(aes(label=labels_v2),vjust=- 15,  size= 10) +ylim(0, 90)
  ggsave(filename, width= 12, height=7)
}

#plotting and determine differences in centrality metrics in ppi networks 
network_cent_method_plotting<- function(limma_lists, deseq2_lists, tfl_lists, string_ppi_details= string_ppi_details, filepath){
  #limma_lists, deseq2_lists, tfl_lists- list of genes in disease signature
  #string_ppi_details- dataframe with network centality measurements for ppi string network
  #filepath- where to save the file
  
  #output
  #plots of the centrality metics and significance testing
  
  #unlist the genes
  limma_genes_v2<- unlist(limma_lists)
  deseq2_genes_v2<- unlist(deseq2_lists)
  topgenes<- unlist(tfl_lists)
  
  names(limma_genes_v2)<- NULL 
  names(deseq2_genes_v2) <- NULL 
  names(topgenes)<- NULL
  
  #convert entrez ids to symbols
  limma_genes_v2<- gconvert(as.character(limma_genes_v2), organism = "hsapiens", target = "HGNC", numeric_ns= "ENTREZGENE_ACC")
  deseq2_genes_v2<- gconvert( deseq2_genes_v2, organism = "hsapiens", target = "HGNC", numeric_ns= "ENTREZGENE_ACC")
  topgenes<- gconvert( topgenes, organism = "hsapiens", target = "HGNC", numeric_ns= "ENTREZGENE_ACC")
  
  #get centraility merics for all the genes and make dataframe 
  lima_ppi <- string_ppi_details[rownames(string_ppi_details) %in% limma_genes_v2$target,]
  deseq2_ppi <- string_ppi_details[rownames(string_ppi_details) %in% deseq2_genes_v2$target,]
  transfer_ppi <- string_ppi_details[rownames(string_ppi_details) %in% topgenes$target,]
  ppi_details_sig_3_methods <- rbind(lima_ppi ,deseq2_ppi , transfer_ppi)
  ppi_details_sig_3_methods$method <- c(rep("limma", nrow(lima_ppi)), rep("DESeq2", nrow(deseq2_ppi)), rep("Transfer Learning", nrow(transfer_ppi)))
  
  
  
  #run wilcox test on only unique genes from method and the different metrics 
  remove_genes<- unique(ppi_details_sig_3_methods$genes[duplicated(ppi_details_sig_3_methods$genes)])
  
  ppi_details_sig_3_methods_filtered <- ppi_details_sig_3_methods[! ppi_details_sig_3_methods$genes %in% remove_genes,]
  #print(ppi_details_sig_3_methods_filtered )
  
  
  my_comparisons <- list( c("Transfer Learning", "limma"), c("limma", "DESeq2"), c("DESeq2", "Transfer Learning") )
  #plot each metric 
  #ppi_details_sig_3_methods$degree_log <- log(ppi_details_sig_3_methods$degree)
  p1 <- ggplot(ppi_details_sig_3_methods_filtered, aes(x=method , y=degree_log, fill=method )) + 
    geom_violin()+ geom_boxplot(width = 0.1, fill = "grey", color = "black")+labs(title="PPI Degree",x="Genes from Disease Signature for Signature Reversion", y = "PPI log(Degree + 1)", color= "Method") + scale_colour_manual(values = c("#440154FF" ,"#21908CFF" ,"#FDE725FF"), aesthetics = c("colour", "fill"))+ theme_minimal()+ stat_compare_means(comparisons = my_comparisons, p.adjust.method = "bonf")+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = 10)
  file<- paste(filepath, "ppi_degree_methods.png", sep="")
  p1
  ggsave(file, bg = "white")
  
  #ppi_details_sig_3_methods$log_betweness <- log(ppi_details_sig_3_methods$betweenness)
  p2<- ggplot(ppi_details_sig_3_methods_filtered, aes(x=method, y=log_betweness , fill=method)) + 
    geom_violin()+ geom_boxplot(width = 0.1, fill = "grey", color = "black")+labs(title="PPI Betweeness",x="Genes from Disease Signature for Signature Reversion", y = "PPI log(Betweeness + 1)", color= "Method") + scale_colour_manual(values =  c("#440154FF" ,"#21908CFF" ,"#FDE725FF"), aesthetics = c("colour", "fill"))+ theme_minimal()+ 
    stat_compare_means(comparisons = my_comparisons, p.adjust.method = "bonf")+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = 20)
  file<- paste(filepath, "ppi_betweeness_methods.png", sep="")
  p2
  ggsave(file, bg = "white")
  
  #ppi_details_sig_3_methods$eign_cent_log <- log(ppi_details_sig_3_methods$eign_cent)
  p3<- ggplot(ppi_details_sig_3_methods_filtered, aes(x=method, y=eign_cent_log, fill=method)) + 
    geom_violin()+ geom_boxplot(width = 0.1, fill = "grey", color = "black")  +labs(title="PPI Eigenvector Centrality Scores",x="Genes from Disease Signature for Signature Reversion", y = "PPI log(Eigenvector Centrality Scores)", color= "Method") + scale_colour_manual(values =  c("#440154FF" ,"#21908CFF" ,"#FDE725FF"), aesthetics = c("colour", "fill"))+ theme_minimal()+ 
    stat_compare_means(comparisons = my_comparisons, method= "wilcox.test", p.adjust.method = "bonf")+ # Add pairwise comparisons p-value
    stat_compare_means(label.y = 10)
  
  file<- paste(filepath, "ppi_Eigenvector_Centrality_Scores_methods.png", sep="")
  p3
  ggsave(file, bg = "white")
  
  print(p1)
  print(p2)
  print(p3)
}

#from the customCMPdb github. not in bioconductor package
loadSDFwithName <- function(source="LINCS"){
  ah <- AnnotationHub()
  annot_path <- ah[["AH79563"]]
  conn <- DBI::dbConnect(SQLite(), annot_path)
  if(source == "LINCS"){
    lincs_annot <- dbReadTable(conn, "lincsAnnot")
    lincs_sdf_path <- ah[["AH79567"]]
    lincs_sdf <-ChemmineR::read.SDFset(lincs_sdf_path)
    cid(lincs_sdf) <- sdfid(lincs_sdf)
    ## make lincs cids as pert names
    brd_ids <- lincs_annot$lincs_id
    names(brd_ids) <- lincs_annot$pert_iname
    brd_uniq <- brd_ids[!duplicated(names(brd_ids))] # 19811
    brd_common <- brd_uniq[brd_uniq %in% cid(lincs_sdf)] # 19758
    res_sdf <- lincs_sdf[brd_common]
    cid(res_sdf) <- names(brd_common)
  }
  if(source == "CMAP2"){
    cmap_annot <- dbReadTable(conn, "cmapAnnot")
    cmap_sdf_path <- ah[["AH79566"]]
    cmap_sdf <-ChemmineR::read.SDFset(cmap_sdf_path)
    cid(cmap_sdf) <- sdfid(cmap_sdf)
    cmap_ids <- cmap_annot$cmap_id
    names(cmap_ids) <- cmap_annot$cmap_name
    res_sdf <- cmap_sdf[cmap_ids]
    cid(res_sdf) <- names(cmap_ids)
  }
  if(source == "DrugBank"){
    db_annot <- dbReadTable(conn, "drugBankAnnot")
    db_sdf_path <- ah[["AH79565"]]
    db_sdf <-ChemmineR::read.SDFset(db_sdf_path)
    cid(db_sdf) <- sdfid(db_sdf)
    db_ids <- db_annot$drugbank_id
    names(db_ids) <- db_annot$name
    db_common <- db_ids[db_ids %in% cid(db_sdf)] # 10569
    res_sdf <- db_sdf[db_common]
    cid(res_sdf) <- names(db_common)
  }
  if(source == "DrugAge"){
    da_annot <- dbReadTable(conn, "drugAgeAnnot")
    da_sdf_path <- ah[["AH79564"]]
    da_sdf <-ChemmineR::read.SDFset(da_sdf_path)
    cid(da_sdf) <- sdfid(da_sdf)
    da_ids <- da_annot$drugage_id
    names(da_ids) <- da_annot$compound_name
    da_uniq <- da_ids[!duplicated(names(da_ids))] # 420
    da_common <- da_uniq[da_uniq %in% cid(da_sdf)] # 223
    res_sdf <- da_sdf[da_common]
    cid(res_sdf) <- names(da_common)
  }
  dbDisconnect(conn)
  return(res_sdf)
}

drug_struc_heatmap<- function(limma_drugs, deseq2_drugs, TFL_drugs, file){
  
  drug_list <- unique(c(limma_drugs,deseq2_drugs,  TFL_drugs))
  IDs <- lincs_sdfset@ID[grepl(paste(drug_list,collapse="|"), lincs_sdfset@ID)]
  IDs<- intersect(IDs, drug_list)
  
  
  
  drug_sdfset <- lincs_sdfset[IDs]
  
  drug_res <- sapply(cid(drug_sdfset ), function(x) fmcsBatch(drug_sdfset[x], drug_sdfset , au=0, bu=0, timeout=10000000, matching.mode = "aromatic")[,"Tanimoto_Coefficient"])
  #return(drug_res)
  
  tfl_drugs3 <-  drug_sdfset@ID %in% TFL_drugs
  limma_drugs3 <- drug_sdfset@ID %in% limma_drugs
  deseq2_drugs3<- drug_sdfset@ID %in% deseq2_drugs
  
  row_ha = HeatmapAnnotation(Transfer_Learning=tfl_drugs3,DESeq2= deseq2_drugs3, limma= limma_drugs3 , col = list(Transfer_Learning = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),DESeq2 = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),limma = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF") ))
  col_fun = colorRamp2(c(0,  1), c( "black", "yellow"))
  #ha = rowAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 97:100), labels = TFL_bp_up[1:10]))
  
  png(file, width = 10, height =7,  units = "in", res = 120) 
  draw(Heatmap(drug_res, nam= "Tanimoto_Coefficient drug structure", col = col_fun, show_column_names = FALSE,  show_row_names = FALSE, top_annotation = row_ha,
               clustering_distance_rows= "euclidean",
               clustering_distance_columns=  "euclidean",
               clustering_method_rows = "ward.D2" ,
               clustering_method_columns="ward.D2"))
  dev.off()
  return(drug_res)
}


volcano_plots<- function(results, method, method_genes, transfer_learning_genes){
  #
  # results- results data frame from either deseq2 or limma 
  # method- "deseq2" or "limma" to indicate which method was used 
  # method_genes- disease signature genes used for signature reversion for method
  # transfer_learning_genes- disease signature genes used for signature reversion for transfer learning 
  
  #return 
  # a volcano plot with the disease signature genes for limma and transfer learning
  
  #adjust limma result names 
  if(method == "limma" ){
    #results<- results$limma
    
    colnames(results)[2]<- "log2FoldChange"
    colnames(results)[6]<- "padj"
    #limma symbol location
    colnames(results)[1]<- "Symbol"
  }
  
  #handle duplicated genes due to ids conversions
  aa <- results[order( -abs(results$log2FoldChange) ), ]
  results<- aa[ !duplicated(aa$Symbol), ]  
  
  # subset by LINCS genes 
  res_lincs_only <- results[!is.na(match_to_LINCS_genes(results$Symbol)), ]
  res_lincs_only$ids <- match_to_LINCS_genes( res_lincs_only$Symbol)
  
  #determine both genes
  both_genes<- intersect(unlist(transfer_learning_genes), unlist(method_genes))
  
  #limma or deseq2 only 
  method_genes_only <- unlist(method_genes)[!unlist(method_genes) %in% both_genes]
  
  #transfer learning only 
  tfl_genes_only <- unlist(transfer_learning_genes)[!unlist(transfer_learning_genes) %in% both_genes]
  
  res_lincs_only$neg_log_adj_p_value <- -1* log(res_lincs_only$padj)
  
  #if deseq2 the following colors based on the network colors and levels
  if(method == "DESeq2"){
    res_lincs_only$groups <- ifelse(res_lincs_only$ids %in% method_genes_only, "DESeq2", ifelse(res_lincs_only$ids %in% tfl_genes_only, "Transfer Learning", ifelse(res_lincs_only$ids %in% both_genes, "Both",  "NA")))
    if (length(both_genes) > 0){
      color_groups<- c("gray", "green4", "#FDE725FF", "#440154FF"  )
      group_levels<- c("NA", "Both", "Transfer Learning", "DESeq2")
    }else{
      color_groups<- c("gray", "#FDE725FF", "#440154FF"  )
      group_levels<- c("NA", "Transfer Learning", "DESeq2")
    }
  }else{
    #if limma the following colors
    res_lincs_only$groups <- ifelse(res_lincs_only$ids %in% method_genes_only, "limma", ifelse(res_lincs_only$ids %in% tfl_genes_only, "Transfer Learning", ifelse(res_lincs_only$ids %in% both_genes, "Both",  "NA")))
    if (length(both_genes) > 0){
      color_groups<- c("gray", "maroon", "#FDE725FF", "#21908CFF"  )
      group_levels<- c("NA", "Both", "Transfer Learning", "limma")
    }else{
      color_groups<- c("gray", "#FDE725FF", "#21908CFF"  )
      group_levels<- c("NA", "Transfer Learning", "limma")
    }
    
  }
  
  #rbind to change the order of the results
  res_lincs_only <- rbind(res_lincs_only[res_lincs_only$groups == "NA",], res_lincs_only[!res_lincs_only$groups == "NA",])
  res_lincs_only<- res_lincs_only[ complete.cases(res_lincs_only),]
  # Main plot
  res_lincs_only$groups <- factor(res_lincs_only$groups , levels= group_levels)
  
  p_value_cut<- -log(0.05)
  
  pmain <- ggplot(res_lincs_only, aes(x = log2FoldChange, y = neg_log_adj_p_value, color = groups))+
    geom_point( size =3) + scale_color_manual(values=color_groups)   +theme_bw() + xlab("log(Fold Change)") + ylab("-log(adj. p-value)") + guides(color=guide_legend(title="Method", override.aes = list(size = 10))) + geom_hline(yintercept=p_value_cut, linetype="dashed", color = "black")+ theme(legend.position="bottom", text = element_text(size = 35,  face="bold"))
  
  # Marginal densities along x axis
  xdens <- axis_canvas(pmain, axis = "x")+
    geom_density(data =res_lincs_only, aes(x = log2FoldChange, fill = groups),
                 alpha = 0.7, size = 0.3)+ scale_fill_manual(values=color_groups)           
  # Marginal densities along y axis
  # Need to set coord_flip = TRUE, if you plan to use coord_flip()
  ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
    geom_density(data = res_lincs_only, aes(x = neg_log_adj_p_value, fill = groups),
                 alpha = 0.7, size = 0.3)+
    coord_flip()+ scale_fill_manual(values=color_groups)     
  
  p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
  p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
  
  #change to histogram 
  #p3 <- ggplot(res_lincs_only, aes(x = log2FoldChange, fill = groups ))+ geom_density() +  scale_fill_manual(values=color_groups) + xlab("log Fold Change") + guides(color=guide_legend(title="Method"))
  #change color to fill, add points, and signficance line
 # p4 <- ggplot(res_lincs_only, aes(x = groups, y = neg_log_adj_p_value, fill = groups ))+ geom_histogram()+ geom_point(fill= "black")+ scale_fill_manual(values=color_groups) + geom_hline(yintercept=p_value_cut, linetype="dashed", color = "black")   + xlab("Methods") + theme(legend.position="none")
  
  #print(p3)
 # print(p4)
  print(ggdraw(p2))
  return(res_lincs_only )
  
  
}


compare_deseq2_limma <- function(deseq_results, limma_results, x_p1=3, y_p1=30, x_p2=400, y_p2=200){
  #comparing the log fold change and adjusted p-values between methods
  #deseq_results-  result output for DESeq2 
  #limma_results - result output for limma 
  #x_p1 - x location for plot 1 (fold change)
  #y_p1- y location for plot 1 (fold change)
  #x_p2 - x location for plot 2 (adjusted pvalue)
  #y_p2 - y location for plot 2 (adjusted pvalue)
  
  #return
  #plots with the spearman correlation between methods for both log fold change and adjsut pvalue 
  
  diff_metrics <- cbind(deseq_results$log2FoldChange, limma_results$logFC, -log(deseq_results$padj), -log(limma_results$adj.P.Val))
  diff_metrics<- diff_metrics[complete.cases(diff_metrics), ]
  
  #diff_metrics<- diff_metrics[diff_metrics$DESeq2_padj < 0.05 |diff_metrics$limma_padj < 0.05,]
  
  colnames(diff_metrics)<- c("DESeq2_logFC","limma_logFC" , "DESeq2_padj", "limma_padj")
  diff_metrics<- as.data.frame(diff_metrics)
  
  P1<-ggplot(diff_metrics, aes(x=limma_logFC, y=DESeq2_logFC) ) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+ stat_cor(method = "spearman", label.x = x_p1, label.y = y_p1) + xlab("limma log(Fold Change)") + ylab("DESeq2 log(Fold Change)") +theme(text = element_text(size = 20,  face="bold"))
  
  P2 <- ggplot(diff_metrics, aes(x=limma_padj, y=DESeq2_padj) ) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis") +
    theme_bw( )+ stat_cor(method = "spearman", label.x = x_p2, label.y = y_p2) + xlab("limma -log(adj. p-value)") + ylab("DESeq2 -log(adj. p-value)") +theme(text = element_text(size = 20,  face="bold"))
  # Add correlation coefficient
  print(P1)
  print(P2)
}



lincs_method_comparsion<- function(deseq_results_all,limma_results_all , tfl_results_all, cell_line= "GI1", x_p1=0, y_p1=2, x_p2=5, y_p2=10){
  #comparing the normalized connnectivity scores and false discovery rates between methods
  #deseq_results_all- the lincs signature search result output for DESeq2 (filter to the cell line of interest)
  #limma_results_all - the lincs signature search result output for limma (filter to the cell line of interest)
  #tfl_results_all- the lincs signature search result output for Transfer Learning (filter to the cell line of interest)
  #cell_line= "GI1" - cell line of interest
  #x_p1=0 - x location for plot 1 (NCS)
  #y_p1=2- y location for plot 1 (NCS)
  #x_p2=5 - x location for plot 2 (FDR)
  #y_p2=10 - y location for plot 2 (FDR)
  
  #return
  #plots with the spearman correlation between methods for both NCS and FDR
  
  #FOCUS ON ONLY THE CELL OF INTEREST
  deseq_results_all <- deseq_results_all[deseq_results_all$cell == "GI1",]
  limma_results_all <- limma_results_all[limma_results_all$cell == "GI1",]
  tfl_results_all  <- tfl_results_all[tfl_results_all$cell == "GI1",]
  
  #order the results for comparision
  deseq_results_all <- deseq_results_all[order(deseq_results_all$pert),]
  limma_results_all<- limma_results_all[order(limma_results_all$pert),]
  tfl_results_all<- tfl_results_all[order(tfl_results_all$pert),]
  
  #check and good
  #deseq_results_all$pert == limma_results_all$pert
  #deseq_results_all$pert == tfl_results_all$pert
  
  diff_metrics <- cbind(deseq_results_all$NCS, limma_results_all$NCS, tfl_results_all$NCS, -log(deseq_results_all$WTCS_FDR), -log(limma_results_all$WTCS_FDR), -log(tfl_results_all$WTCS_FDR))
  diff_metrics<- diff_metrics[complete.cases(diff_metrics), ]
  
  #diff_metrics<- diff_metrics[diff_metrics$DESeq2_padj < 0.05 |diff_metrics$limma_padj < 0.05,]
  
  colnames(diff_metrics)<- c("DESeq2_NCS","limma_NCS" , "TFL_NCS",  "DESeq2_FDR", "limma_FDR", "TFL_FDR")
  diff_metrics<- as.data.frame(diff_metrics)
  
  #deseq2 vs limma 
  P1<-ggplot(diff_metrics, aes(x=DESeq2_NCS, y=limma_NCS) ) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+ stat_cor(method = "spearman", label.x = x_p1, label.y = y_p1) + xlab("DESeq2 NCS") + ylab("limma NCS") +theme(text = element_text(size = 20,  face="bold"))
  
  P2 <- ggplot(diff_metrics, aes(x=DESeq2_FDR, y=limma_FDR) ) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+ stat_cor(method = "spearman", label.x = x_p2, label.y = y_p2) + xlab("DESeq2 -log(FDR)") + ylab("limma -log(FDR)") +theme(text = element_text(size = 20,  face="bold"))
  
  # Add correlation coefficient
  print(P1)
  print(P2)
  
  #deseq2 vs tfl 
  P1<-ggplot(diff_metrics, aes(x=DESeq2_NCS, y=TFL_NCS) ) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+ stat_cor(method = "spearman", label.x = x_p1, label.y = y_p1) + xlab("DESeq2 NCS") + ylab("Transfer Learning NCS") +theme(text = element_text(size = 20,  face="bold"))
  
  P2 <- ggplot(diff_metrics, aes(x=DESeq2_FDR, y=TFL_FDR) ) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+ stat_cor(method = "spearman", label.x = x_p2, label.y = y_p2) + xlab("DESeq2 -log(FDR)") + ylab("Transfer Learning -log(FDR)") +theme(text = element_text(size = 20,  face="bold"))
  
  # Add correlation coefficient
  print(P1)
  print(P2)
  
  #limma vs tfl 
  P1<-ggplot(diff_metrics, aes(x=limma_NCS, y=TFL_NCS) ) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+ stat_cor(method = "spearman", label.x = x_p1, label.y = y_p1) + xlab("limma NCS") + ylab("Transfer Learning NCS") +theme(text = element_text(size = 20,  face="bold"))
  
  P2 <- ggplot(diff_metrics, aes(x=limma_FDR, y=TFL_FDR) ) +
    geom_bin2d(bins = 60) +
    scale_fill_continuous(type = "viridis") +
    theme_bw()+ stat_cor(method = "spearman", label.x = x_p2, label.y = y_p2) + xlab("limma -log(FDR)") + ylab("Transfer Learning -log(FDR)") +theme(text = element_text(size = 20,  face="bold"))
  
  # Add correlation coefficient
  print(P1)
  print(P2)
  
} 

#plotting drug target and mechanism of action for all candidate results
drug_moa_target_plotting_ct_info <- function(fda_approved_res, lincs_commpound_info=lincs_commpound_info, method, cancer, file_path_header ){
  #fda_approved_res- the signaturesearch results with only fda approved
  #lincs_commpound_info- the lincs compound info from the LINCS database
  #method- what method found these candidates
  #cancer- which cancer was this result for 
  #file_path_header- path and header for output files
  
  #output
  #several plots including alluvial and bar plots for the candidates
  file_starter<- paste0(file_path_header, "_", cancer, "_", method, "_", sep= "")
  
  test_fda <- fda_approved_res
  lincs_commpound_info[lincs_commpound_info$cmap_name %in% test_fda$pert,]
  
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
    
    df<- expand.grid(test_fda$pert[i], targets, moa )
    
    summ_df<- rbind(summ_df, df)
    
  }
  
  
  summ_df<- summ_df[-1,]
  colnames(summ_df)<- c( "cmap_name", "target", "moa")
  
  summ_df<- summ_df[!duplicated(summ_df), ]
  #  return(summ_df)
  #}   
  TEST_LODES<- to_lodes_form(as.data.frame(summ_df),
                             axes = 1:3,
                             id = "INDEX")
  
  TEST_LODES$stratum = str_wrap(TEST_LODES$stratum, width = 20)
  TEST_LODES$x <- factor(TEST_LODES$x, levels= c("moa", "cmap_name", "target"))
  
  TEST_LODES$drug <- rep(NA, nrow(summ_df))
  drugs<- unique(summ_df$cmap_name)
  for (i in 1:length(drugs)){
    ind <- TEST_LODES$INDEX[TEST_LODES$stratum == drugs[i]]
    TEST_LODES$drug<- ifelse(TEST_LODES$INDEX %in% ind, drugs[i], TEST_LODES$drug)
  }
  
  #not need for this analysis but many for future.
  TEST_LODES$drug<- factor(TEST_LODES$drug, levels= test_fda$pert)
  
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
  
  
  counts_moa_v2<- counts_moa[counts_moa$Freq >0,]
  counts_moa_v2$Var1<- factor(counts_moa_v2$Var1, levels= counts_moa_v2$Var1[order(-counts_moa_v2$Freq)])
  if (nrow(counts_moa_v2) < 20){
    ggplot(counts_moa_v2, aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Mechanism of Action") + xlab("Number of Drugs")+ theme(text = element_text(size = 20,  face="bold"))
    
    file_name<- paste0(file_starter, "_", "barplot_moa_v2.png")
    ggsave(file_name, width = 8, height=10, units= "in")
  }else{
    counts_moa_v2<- counts_moa_v2[order(-counts_moa_v2$Freq),]
    ggplot(counts_moa_v2[1:20,], aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Mechanism of Action") + xlab("Number of Drugs")+ theme(text = element_text(size = 20,  face="bold"))
    
    file_name<- paste0(file_starter, "_", "barplot_moa_v2.png")
    ggsave(file_name, width = 8, height=10, units= "in")
  }
  
  
  counts$Var1<- factor(counts$Var1, levels= counts$Var1[order(-counts$Freq)])
  if( nrow(counts)< 20 ){
    #counts$Var1<- factor(counts$Var1, levels= counts$Var1[order(-counts$Freq),])
    ggplot(counts, aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Targets") + xlab("Number of Drugs")+ theme(text = element_text(size = 20,  face="bold"))
    file_name<- paste0(file_starter, "_", "barplot_target_v2.png")
    ggsave(file_name, width = 8, height=10, units= "in")
  }else{
    
    counts<- counts[order(-counts$Freq),]
    
    ggplot(counts[1:20,], aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Targets") + xlab("Number of Drugs")+ theme(text = element_text(size = 20,  face="bold"))
    file_name<- paste0(file_starter, "_", "barplot_target_v2.png")
    ggsave(file_name, width = 8, height=10, units= "in")
    
    
    
  }
  #add ploting function
  target_count <- counts
  
  
  #reduce to the top 20 drugs
  summ_df<- summ_df[summ_df$cmap_name %in% test_fda$pert[1:20],]
  
  summ_df_moa<- summ_df[,c(1,3)]
  summ_df_moa<- summ_df_moa[!duplicated(summ_df_moa), ]
  TEST_LODES_moa<- to_lodes_form(as.data.frame(summ_df_moa),
                                 axes = 1:2,
                                 id = "INDEX")
  
  TEST_LODES_moa$stratum = str_wrap(TEST_LODES_moa$stratum, width = 20)
  
  TEST_LODES_moa$drug <- rep(NA, nrow(summ_df_moa))
  drugs<- unique(summ_df_moa$cmap_name)
  for (i in 1:length(drugs)){
    ind <- TEST_LODES_moa$INDEX[TEST_LODES_moa$stratum == drugs[i]]
    TEST_LODES_moa$drug<- ifelse(TEST_LODES_moa$INDEX %in% ind, drugs[i], TEST_LODES_moa$drug)
  }
  
  #not need for this analysis but many for future.
  TEST_LODES_moa$drug<- factor(TEST_LODES_moa$drug, levels= test_fda$pert)
  ct_drugs<- test_fda$pert[as.logical(unlist(test_fda[,17]))]
  TEST_LODES_moa$ct<- TEST_LODES_moa$drug %in% ct_drugs
  TEST_LODES_moa$ct<- factor(TEST_LODES_moa$ct, levels= c(TRUE, FALSE))
  #moa
  moa_df <- TEST_LODES_moa[TEST_LODES_moa$x == "moa", 2:4]
  moa_df<- moa_df[!duplicated(moa_df), ]
  counts<- as.data.frame(table(moa_df$stratum))
  
  #counts<- as.data.frame(table(TEST_LODES_moa$stratum[TEST_LODES_moa$x == "moa"]))
  moa_order <- counts$Var1[ order(-counts$Freq)]
  counts_moa<- counts
  
  #factor the stratum column
  TEST_LODES_moa$stratum<- factor(TEST_LODES_moa$stratum, levels= unique(c(test_fda$pert, as.character(moa_order))))
  
  TEST_LODES_moa_com <- TEST_LODES_moa[ complete.cases(TEST_LODES_moa),]
  #TEST_LODES_com$index2 <- as.numeric(TEST_LODES_com$drug)
  #TEST_LODES_com_v2 <- TEST_LODES_com[!TEST_LODES_com$x == "cmap_name",]
  
  ct2<-  TEST_LODES_moa_com$ct 
  ggplot(TEST_LODES_moa_com,
         aes(x = x, stratum = stratum, alluvium = INDEX, label= stratum, fill= ct))+ geom_flow(aes(fill=ct), stat = "alluvium", color= "black") + geom_stratum() +geom_label(stat = "stratum", fill="white",   size= 6 ,  face="bold") + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.position="none", axis.text=element_text(size=25)) + scale_x_discrete(labels=c("cmap_name" = "Drug", "moa" = "Mechanism of Action")) +  scale_fill_manual(values=c("#0066CC", "#E0E0E0"))+ theme(text = element_text(size = 20,  face="bold"))
  
  file_name<- paste0(file_starter, "_", "alluvial_moa_v2.png")
  ggsave(file_name, width =10, height=20, units= "in")
  
  #add saving plot function
  
  summ_df_tar<- summ_df[,c(1,2)]
  #keep only the top targets
  summ_df_tar<- summ_df_tar[summ_df_tar[,2] %in% target_count$Var1[1:20], ]
  
  summ_df_tar<- summ_df_tar[!duplicated(summ_df_tar), ]
  TEST_LODES_tar<- to_lodes_form(as.data.frame(summ_df_tar),
                                 axes = 1:2,
                                 id = "INDEX")
  
  TEST_LODES_tar$stratum = str_wrap(TEST_LODES_tar$stratum, width = 20)
  
  TEST_LODES_tar$drug <- rep(NA, nrow(summ_df_tar))
  drugs<- unique(summ_df_tar$cmap_name)
  for (i in 1:length(drugs)){
    ind <- TEST_LODES_tar$INDEX[TEST_LODES_tar$stratum == drugs[i]]
    TEST_LODES_tar$drug<- ifelse(TEST_LODES_tar$INDEX %in% ind, drugs[i], TEST_LODES_tar$drug)
  }
  
  #not need for this analysis but many for future.
  TEST_LODES_tar$drug<- factor(TEST_LODES_tar$drug, levels= test_fda$pert)
  TEST_LODES_tar$ct<- TEST_LODES_tar$drug %in% ct_drugs
  TEST_LODES_tar$ct<- factor(TEST_LODES_tar$ct, levels= c(TRUE, FALSE))
  #targets
  target_df <- TEST_LODES[TEST_LODES$x == "target", 2:4]
  target_df<- target_df[!duplicated(target_df), ]
  counts<- as.data.frame(table(target_df$stratum))
  #counts
  target_order <- counts$Var1[ order(-counts$Freq)]
  
  #factor the stratum column
  TEST_LODES_tar$stratum<- factor(TEST_LODES_tar$stratum, levels= unique(c(test_fda$pert, as.character(target_order))))
  
  TEST_LODES_tar_com <- TEST_LODES_tar[ complete.cases(TEST_LODES_tar),]
  #TEST_LODES_com$index2 <- as.numeric(TEST_LODES_com$drug)
  #TEST_LODES_com_v2 <- TEST_LODES_com[!TEST_LODES_com$x == "cmap_name",]
  cts<- c(TRUE, FALSE)
  ggplot(TEST_LODES_tar_com,
         aes(x = x, stratum = stratum, alluvium = INDEX, label= stratum, fill= ct))+ geom_flow( aes(fill= ct), stat = "alluvium", color= "black") +
    geom_stratum() +geom_label(stat = "stratum", fill="white",  size= 10,  face="bold") + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.position="none", axis.text=element_text(size=25)) + scale_y_continuous(breaks=NULL) + scale_x_discrete(labels=c("cmap_name" = "Drug", "target" = "Targets")) +  scale_fill_viridis_d()+ theme(text = element_text(size = 20,  face="bold"))
  file_name<- paste0(file_starter, "_", "alluvial_target_v2.png")
  ggsave(file_name, width = 10, height=20, units= "in")
}



#plotting drug target and mechanism of action for all candidate results
drug_moa_target_plotting_all_info <- function(fda_approved_res, lincs_commpound_info=lincs_commpound_info, method, cancer, file_path_header ){
  #fda_approved_res- drug info with the methods and target info
  #lincs_commpound_info- the lincs compound info from the LINCS database
  #method- what method found these candidates
  #cancer- which cancer was this result for 
  #file_path_header- path and header for output files
  
  #output
  #several plots including alluvial and bar plots for the candidates
  file_starter<- paste0(file_path_header, "_", cancer, "_", method, "_", sep= "")
  
  test_fda <- fda_approved_res
  lincs_commpound_info[lincs_commpound_info$cmap_name %in% test_fda$pert,]
  
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
    
    df<- expand.grid(test_fda[i,1], targets, moa )
    
    summ_df<- rbind(summ_df, df)
    
  }
  
  
  summ_df<- summ_df[-1,]
  colnames(summ_df)<- c( "cmap_name", "target", "moa")
  
  summ_df<- summ_df[!duplicated(summ_df), ]
  #  return(summ_df)
  #}   
  TEST_LODES<- to_lodes_form(as.data.frame(summ_df),
                             axes = 1:3,
                             id = "INDEX")
  
  TEST_LODES$stratum = str_wrap(TEST_LODES$stratum, width = 20)
  TEST_LODES$x <- factor(TEST_LODES$x, levels= c("moa", "cmap_name", "target"))
  
  TEST_LODES$drug <- rep(NA, nrow(summ_df))
  
  drugs<- unique(summ_df$cmap_name)
  
  for (i in 1:length(drugs)){
    ind <- TEST_LODES$INDEX[TEST_LODES$stratum == drugs[i]]
    TEST_LODES$drug<- ifelse(TEST_LODES$INDEX %in% ind, drugs[i], TEST_LODES$drug)
  }
  TEST_LODES$method<- rep(NA, nrow(summ_df))
  for (i in 1:nrow(TEST_LODES)){
    TEST_LODES$method[i]<- test_fda$group[test_fda$pert %in% TEST_LODES$drug[i] ]
  }
  
  #not need for this analysis but many for future.
  TEST_LODES$drug<- factor(TEST_LODES$drug, levels= test_fda$pert)
  
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
  
  
  counts_moa_v2<- counts_moa[counts_moa$Freq >0,]
  counts_moa_v2$Var1<- factor(counts_moa_v2$Var1, levels= counts_moa_v2$Var1[order(-counts_moa_v2$Freq)])
  if (nrow(counts_moa_v2) < 20){
    ggplot(counts_moa_v2, aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Mechanism of Action") + xlab("Number of Drugs")+ theme(text = element_text(size = 20,  face="bold"))
    
    file_name<- paste0(file_starter, "_", "barplot_moa_all.png")
    ggsave(file_name, width = 8, height=10, units= "in")
  }else{
    counts_moa_v2<- counts_moa_v2[order(-counts_moa_v2$Freq),]
    ggplot(counts_moa_v2[1:20,], aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Mechanism of Action") + xlab("Number of Drugs")+ theme(text = element_text(size = 20,  face="bold"))
    
    file_name<- paste0(file_starter, "_", "barplot_moa_all.png")
    ggsave(file_name, width = 8, height=10, units= "in")
  }
  
  
  counts$Var1<- factor(counts$Var1, levels= counts$Var1[order(-counts$Freq)])
  if( nrow(counts)< 20 ){
    #counts$Var1<- factor(counts$Var1, levels= counts$Var1[order(-counts$Freq),])
    ggplot(counts, aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Targets") + xlab("Number of Drugs")+ theme(text = element_text(size = 20,  face="bold"))
    file_name<- paste0(file_starter, "_", "barplot_target_all.png")
    ggsave(file_name, width = 8, height=10, units= "in")
  }else{
    
    counts<- counts[order(-counts$Freq),]
    
    ggplot(counts[1:20,], aes(x= Freq, y= Var1)) + geom_bar(stat= "identity", fill="deepskyblue4", color= "black") + ylab("Targets") + xlab("Number of Drugs")+ theme(text = element_text(size = 20,  face="bold"))
    file_name<- paste0(file_starter, "_", "barplot_target_all.png")
    ggsave(file_name, width = 8, height=10, units= "in")
    
    
    
  }
  #add ploting function
  target_count <- counts
  
  
  #reduce to the top 20 drugs<= looking at all the drugs
  #summ_df<- summ_df[summ_df$cmap_name %in% test_fda$pert[1:20],]
  
  summ_df_moa<- summ_df[,c(1,3)]
  summ_df_moa<- summ_df_moa[!duplicated(summ_df_moa), ]
  TEST_LODES_moa<- to_lodes_form(as.data.frame(summ_df_moa),
                                 axes = 1:2,
                                 id = "INDEX")
  
  #TEST_LODES_moa$stratum = str_wrap(TEST_LODES_moa$stratum, width = 30)
  
  TEST_LODES_moa$drug <- rep(NA, nrow(summ_df_moa))
  drugs<- unique(summ_df_moa$cmap_name)
  for (i in 1:length(drugs)){
    ind <- TEST_LODES_moa$INDEX[TEST_LODES_moa$stratum == drugs[i]]
    TEST_LODES_moa$drug<- ifelse(TEST_LODES_moa$INDEX %in% ind, drugs[i], TEST_LODES_moa$drug)
  }
  
  #not need for this analysis but many for future.
  TEST_LODES_moa$drug<- factor(TEST_LODES_moa$drug, levels= test_fda$pert)
  ct_drugs<- test_fda$pert[as.logical(unlist(test_fda[,4]))]
  TEST_LODES_moa$ct<- TEST_LODES_moa$drug %in% ct_drugs
  TEST_LODES_moa$ct<- factor(TEST_LODES_moa$ct, levels= c(TRUE, FALSE))
  
  #add the methods info
  TEST_LODES_moa$method<- rep(NA, nrow(TEST_LODES_moa))
  for (i in 1:nrow(TEST_LODES_moa)){
    TEST_LODES_moa$method[i]<- test_fda$group[test_fda$pert %in% TEST_LODES_moa$drug[i] ]
  }
  
  
  
  #moa
  moa_df <- TEST_LODES_moa[TEST_LODES_moa$x == "moa", 2:4]
  moa_df<- moa_df[!duplicated(moa_df), ]
  counts<- as.data.frame(table(moa_df$stratum))
  
  #counts<- as.data.frame(table(TEST_LODES_moa$stratum[TEST_LODES_moa$x == "moa"]))
  moa_order <- counts$Var1[ order(-counts$Freq)]
  counts_moa<- counts
  
  #factor the stratum column
  TEST_LODES_moa$stratum<- factor(TEST_LODES_moa$stratum, levels= unique(c(test_fda$pert, as.character(moa_order))))
  
  TEST_LODES_moa_com <- TEST_LODES_moa[ complete.cases(TEST_LODES_moa),]
  #TEST_LODES_com$index2 <- as.numeric(TEST_LODES_com$drug)
  #TEST_LODES_com_v2 <- TEST_LODES_com[!TEST_LODES_com$x == "cmap_name",]
  
  ct2<-  TEST_LODES_moa_com$ct 
  print(TEST_LODES_moa_com)
  ggplot(TEST_LODES_moa_com,
         aes(x = x, stratum = stratum, alluvium = INDEX, label= stratum))+ geom_flow(aes(fill=method), stat = "alluvium", color= "gray") + geom_stratum(aes(fill=method), color="gray") + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.position="none", axis.text=element_text(size=40))+ theme(text = element_text(size = 20,  face="bold")) +scale_fill_manual(values=c("All Methods"="#FDE725FF","DESeq2 and Transfer Learning" ="#8FD744FF", "DESeq2"="#35B779FF", "DESeq2 and limma"= "#21908CFF", "limma"= "#31688EFF", "limma and Transfer Learning"= "#443A83FF", "Transfer Learning" = "#440154FF", "TRUE"= "maroon", "FALSE"= "white")) + theme(axis.text.y=element_blank(),axis.ticks.y=element_blank() )+  scale_x_discrete(expand = expansion(add = 2),labels=c("cmap_name" = "Drug", "moa" = "Mechanism of Action") )+
    ggrepel::geom_text_repel(
      aes(color= ct, label = ifelse(after_stat(x) == 1, as.character(after_stat(stratum)), NA)),
      stat = "stratum", size = 20, direction = "y", nudge_x = -.8, segment.size      = 4, segment.color = "gray", fontface= 2
    ) +  ggrepel::geom_text_repel(
      aes( label = ifelse(after_stat(x) == 2, as.character(after_stat(stratum)), NA)),
      stat = "stratum", size = 18, direction = "y", nudge_x = +1, segment.size      = 3, segment.color= "gray", fontface= 2
    )+ scale_color_manual(values=c("TRUE"="darkgray","FALSE" ="maroon")) 
  file_name<- paste0(file_starter, "_", "alluvial_moa_all.png")
  ggsave(file_name, width =49, height=25, units= "in")
  
  #add saving plot function
  
  summ_df_tar<- summ_df[,c(1,2)]
  #keep only the top targets
  summ_df_tar<- summ_df_tar[summ_df_tar[,2] %in% target_count$Var1[1:20], ]
  
  summ_df_tar<- summ_df_tar[!duplicated(summ_df_tar), ]
  TEST_LODES_tar<- to_lodes_form(as.data.frame(summ_df_tar),
                                 axes = 1:2,
                                 id = "INDEX")
  
  TEST_LODES_tar$stratum = str_wrap(TEST_LODES_tar$stratum, width = 20)
  
  TEST_LODES_tar$drug <- rep(NA, nrow(summ_df_tar))
  drugs<- unique(summ_df_tar$cmap_name)
  for (i in 1:length(drugs)){
    ind <- TEST_LODES_tar$INDEX[TEST_LODES_tar$stratum == drugs[i]]
    TEST_LODES_tar$drug<- ifelse(TEST_LODES_tar$INDEX %in% ind, drugs[i], TEST_LODES_tar$drug)
  }
  
  #not need for this analysis but many for future.
  TEST_LODES_tar$drug<- factor(TEST_LODES_tar$drug, levels= test_fda$pert)
  TEST_LODES_tar$ct<- TEST_LODES_tar$drug %in% ct_drugs
  TEST_LODES_tar$ct<- factor(TEST_LODES_tar$ct, levels= c(TRUE, FALSE))
  
  #add the methods info
  TEST_LODES_tar$method<- rep(NA, nrow(TEST_LODES_tar))
  for (i in 1:nrow(TEST_LODES_tar)){
    TEST_LODES_tar$method[i]<- test_fda$group[test_fda$pert %in% TEST_LODES_tar$drug[i] ]
  }
  
  
  #targets
  target_df <- TEST_LODES[TEST_LODES$x == "target", 2:4]
  target_df<- target_df[!duplicated(target_df), ]
  counts<- as.data.frame(table(target_df$stratum))
  #counts
  target_order <- counts$Var1[ order(-counts$Freq)]
  
  #factor the stratum column
  TEST_LODES_tar$stratum<- factor(TEST_LODES_tar$stratum, levels= unique(c(test_fda$pert, as.character(target_order))))
  
  TEST_LODES_tar_com <- TEST_LODES_tar[ complete.cases(TEST_LODES_tar),]
  #TEST_LODES_com$index2 <- as.numeric(TEST_LODES_com$drug)
  #TEST_LODES_com_v2 <- TEST_LODES_com[!TEST_LODES_com$x == "cmap_name",]
  cts<- c(TRUE, FALSE)
  #ggplot(TEST_LODES_tar_com,
  #    aes(x = x, stratum = stratum, alluvium = INDEX, label= stratum, fill= ct))+ geom_flow( aes(fill= method), stat = "alluvium", color= "black") +
  #geom_stratum(aes(fill= method)) +geom_label(stat = "stratum", fill=ct,  size= 10,  face="bold") + theme_bw() + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.ticks = element_blank(), legend.position="none", axis.text=element_text(size=25)) + scale_y_continuous(breaks=NULL) + scale_x_discrete(labels=c("cmap_name" = "Drug", "target" = "Targets")) +  scale_fill_viridis_d()+ theme(text = element_text(size = 20,  face="bold"))
  # file_name<- paste0(file_starter, "_", "alluvial_target_all.png")
  #ggsave(file_name, width = 20, height=10, units= "in")
}

drug_target_analysis <- function(drug_name, target_list,  tpm_df, sample_vector, deseq_results, result_path){
  #input
  #drug_name- character of the name of drug
  #target_list- the list of genes that are drug targets
  #tpm_df- a data frame with the transcript per million for the gene expression for tummor and control 
  #sample_vector- same order as the column names of the tpm_df, but the sample group (tumor or control)
  #deseq2_result- deseq2 result data frame with the symbol included for each gene
  #result_path- the path to the directory to store csv and plot images 
  
  #output
  #target different expression plots
  
  #drug target expression
  tpm_sub <- tpm_df[rownames(tpm_df) %in% target_list,]
  tpm_sub$genes<- rownames(tpm_sub)
  tpm_sub_v2 <- tpm_sub %>%
    pivot_longer(!genes, names_to = "sample", values_to = "tpm")
  tpm_sub_v2$type <- factor(rep(sample_vector, nrow(tpm_sub)))
  name<- drug_name
  y_max<- max(tpm_sub_v2$tpm) +5
  
  
  deseq_sub <-  deseq_results[ deseq_results$Symbol %in% unique( tpm_sub$genes),]
  if (nrow(deseq_sub)>0 ){
    stat_table <- deseq_sub[,c(6,5)]
    stat_table$group1 <- rep("Primary Tumor" ,  nrow(deseq_sub))
    stat_table$group2 <- rep("Solid Tissue Normal" ,  nrow(deseq_sub))
    stat_table$padj <- ifelse(stat_table$padj < 0.05, formatC(stat_table$padj, format = "e", digits = 2), "")
    colnames(stat_table)<- c("genes", "padj" ,  "group1" ,"group2")
    
    
    file_name<- paste0(result_path, "/", drug_name, "_drug_target_expression_deseq2_padj.png")
    
    ggplot(tpm_sub_v2, aes(x=genes, y=tpm, color= type)) + geom_boxplot()+ stat_pvalue_manual(stat_table, x="genes" , label = "{padj}", y.position= y_max,  size = 6) + labs(title=name,x="Drug Targets", y = "Gene Expression (TPM)", color= "Sample Type") + scale_colour_manual(values =  c("#440154FF","#228C8DFF"), aesthetics = c("colour", "fill"))+ theme(text = element_text(size = 30,  face="bold")) + theme(axis.text.x=element_text(angle=45, size=20, vjust = 0.5, hjust = 0.5))
    #return(tpm_sub_v2$type)
    ggsave(file_name, width = 20, height = 8)
  }
}

#these next four function are from the pzfx r package version 0.3.0 which is used to read and write "GraphPad Prism files
read_col <- function(col, strike_action="exclude", format="", col_name="") {
  if ("Title" %in% names(col)) {
    col_name <- paste(unlist(col[["Title"]]), collapse="")
  }
  subcol_lst <- list()
  for (i in seq_len(length(col))) {
    if (names(col)[i] == "Subcolumn") {
      this_subcol <- read_subcol(col[[i]], strike_action=strike_action)
      subcol_lst[[length(subcol_lst) + 1]] <- this_subcol
    }
  }
  
  if (length(subcol_lst) == 1) {
    col_names <- col_name
  } else if (format == "error") {
    col_names <- paste0(col_name, c("_X", "_ERROR"))
  } else if (format == "replicates") {
    col_names <- paste(col_name, seq_len(length(subcol_lst)), sep="_")
  } else if (format == "SDN") {
    col_names <- paste0(col_name, c("_MEAN", "_SD", "_N"))
  } else if (format == "SEN") {
    col_names <- paste0(col_name, c("_MEAN", "_SEM", "_N"))
  } else if (format == "CVN") {
    col_names <- paste0(col_name, c("_MEAN", "_CV", "_N"))
  } else if (format == "SD") {
    col_names <- paste0(col_name, c("_MEAN", "_SD"))
  } else if (format == "SE") {
    col_names <- paste0(col_name, c("_MEAN", "_SE"))
  } else if (format == "CV") {
    col_names <- paste0(col_name, c("_MEAN", "_CV"))
  } else if (format == "SD") {
    col_names <- paste0(col_name, c("_MEAN", "_SD"))
  } else if (format == "low-high") {
    col_names <- paste0(col_name, c("_MEAN", "_PLUSERROR", "_MINUSERROR"))
  } else if (format == "upper-lower-limits") {
    col_names <- paste0(col_name, c("_MEAN", "_UPPERLIMIT", "_LOWERLIMIT"))
  } else {
    stop("Sorry, don't know how to parse column format.")
  }
  
  names(subcol_lst) <- col_names
  max_len <- max(sapply(subcol_lst, length))
  long_subcol_lst <- lapply(subcol_lst, function(s) {
    length(s) <- max_len
    s
  })
  
  ret <- as.data.frame(long_subcol_lst, stringsAsFactors=FALSE)
  names(ret) <- col_names
  return(ret)
}

read_subcol <- function(subcol, strike_action="exclude") {
  strike_action <- tolower(strike_action)
  if (!strike_action %in% c("exclude", "keep", "star", "e", "k", "s")) {
    stop("strike_action must be one of c('exclude', 'keep', 'star', 'e', 'k', 's')")
  }
  vals <- rep(NA, length(subcol))
  for (i in seq_len(length(subcol))) {
    val <- unlist(subcol[[i]])
    if (is.null(val)) val <- NA
    if ("Excluded" %in% names(attributes(subcol[[i]]))) {
      if (attr(subcol[[i]], "Excluded") == "1") {
        if (strike_action %in% c("exclude", "e")) {
          val <- NA
        } else if (strike_action %in% c("keep", "k")) {
          val <- val
        } else if (strike_action %in% c("star", "s")) {
          val <- paste0(val, "*")
        }
      }
    }
    vals[i] <- val
  }
  if (!strike_action %in% c("star", "s")) {
    suppressWarnings(new_vals <- as.numeric(vals))
    if (all(is.na(new_vals) == is.na(vals))) vals <- new_vals
  }
  return(vals)
}
pzfx_tables<- function (path) {
  xml <- xml2::read_xml(path)
  table_nodes <- xml2::xml_find_all(xml, ".//*[name()='Table' or name()='HugeTable']")
  tables <- sapply(table_nodes, function(t) xml2::xml_text(xml2::xml_child(t, 
                                                                           ".//*[name()='Title']")))
  return(tables)
}

read_pzfx <- function (path, table = 1, strike_action = "exclude", date_x = "character"){
  date_x <- tolower(date_x)
  if (!date_x %in% c("numeric", "character", "both", "n", "c", 
                     "b")) {
    stop("date_x must be one of c('numeric', 'character', 'both', 'n', 'c', 'b')")
  }
  table_names <- pzfx_tables(path)
  if (is.numeric(table)) {
    if (table > length(table_names)) 
      stop("Table index out of range")
    this_idx <- table
  }
  else {
    table <- as.character(table)
    if (!table %in% table_names) 
      stop(sprintf("Can't find %s in prism file", table))
    this_idx <- which(table_names == table)
    if (length(this_idx) > 1) {
      warning(sprintf("Multiple tables named %s, returning the first one only", 
                      table))
      this_idx <- this_idx[1]
    }
  }
  xml <- xml2::read_xml(path)
  table_nodes <- xml2::xml_find_all(xml, ".//*[name()='Table' or name()='HugeTable']")
  this_table <- xml2::as_list(table_nodes[[this_idx]])
  if (!"Title" %in% names(this_table)) 
    stop("Can't work with this pzfx file, is it later than v6.0?")
  if (is.character(table) && table != this_table[["Title"]]) 
    stop("Can't work with this pzfx file, is it later than v6.0?")
  x_format <- ""
  if ("XFormat" %in% names(attributes(this_table))) {
    x_format <- attributes(this_table)$XFormat
  }
  y_format <- ""
  if ("YFormat" %in% names(attributes(this_table))) {
    y_format <- attributes(this_table)$YFormat
  }
  col_lst <- list()
  for (i in seq_len(length(this_table))) {
    if (names(this_table)[i] == "XColumn") {
      if (x_format == "date" & (date_x %in% c("numeric", 
                                              "n", "both", "b") | !"XAdvancedColumn" %in% names(this_table))) {
        this_col <- read_col(this_table[[i]], strike_action = strike_action, 
                             col_name = "X", format = "")
        if (date_x %in% c("both", "b")) {
          colnames(this_col) <- paste0(colnames(this_col), 
                                       "_1")
        }
      }
      else if (x_format == "date") {
        next
      }
      else {
        this_col <- read_col(this_table[[i]], strike_action = strike_action, 
                             col_name = "X", format = x_format)
      }
      if (nrow(this_col) > 0) {
        col_lst[[length(col_lst) + 1]] <- this_col
      }
    }
    else if (names(this_table)[i] == "XAdvancedColumn") {
      if (x_format == "date" & date_x %in% c("character", 
                                             "c", "both", "b")) {
        this_col <- read_col(this_table[[i]], strike_action = strike_action, 
                             col_name = "X", format = "")
        if (date_x %in% c("both", "b")) {
          colnames(this_col) <- paste0(colnames(this_col), 
                                       "_2")
        }
        if (nrow(this_col) > 0) {
          col_lst[[length(col_lst) + 1]] <- this_col
        }
      }
      else {
        next
      }
    }
    else if (names(this_table)[i] == "RowTitlesColumn") {
      this_col <- read_col(this_table[[i]], strike_action = strike_action, 
                           col_name = "ROWTITLE", format = "")
      if (nrow(this_col) > 0) {
        col_lst[[length(col_lst) + 1]] <- this_col
      }
    }
    else if (names(this_table)[i] == "YColumn") {
      this_col <- read_col(this_table[[i]], strike_action = strike_action, 
                           format = y_format)
      col_lst[[length(col_lst) + 1]] <- this_col
    }
  }
  if (length(col_lst) == 0) 
    return(data.frame())
  max_len <- max(sapply(col_lst, nrow))
  long_col_lst <- lapply(col_lst, function(c) {
    while (nrow(c) < max_len) {
      col_names <- colnames(c)
      c <- rbind(c, NA)
      colnames(c) <- col_names
    }
    c
  })
  ret <- Reduce("cbind", long_col_lst)
  return(ret)
}
