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

  #run deseq2 model
  dds <- DESeq(dds)
  
  #get results and save them and return it 
  dds_res <- results(dds, c(condition,B,A))
  res_df <- as.data.frame(dds_res)
  write.table(res_df, file=file_name, sep=",", row.names =TRUE)
  return(res_df)
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

#limma wrapper- smiliar to Transfer Learning wrapper, but adjusted to not plot the bxoplot diffLVs because it genes in this case it is just a ton plots
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
tau_scores<- function(results_df, tau_score_df){
  #download the tau values 
  #scores <- readRDS(tau_score_file)
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



FDA_APPROVAL_CHECK<- function(drug_list){
  fda_info <- read.table("~/data/fda_product_info_df.csv" )
  fda_approve <- fda_info[fda_info$market_status_v2 %in% c("Over-the-Counter", "Prescription" ),]
  return_list <- c()
  for (i in 1:length(drug_list)){
    test <- grep( drug_list[i], fda_approve$ActiveIngredient, ignore.case= TRUE)
    if (length(test)>0 ){
      return_list[i] <- TRUE
    }else{
      return_list[i] <- FALSE
    }
    #more than 3 character are needed for matching
    if (nchar(drug_list[i])< 3){
      return_list[i] <- FALSE
    }
  }
  names(return_list)<- drug_list
  return(return_list)
}

sig_filter_JLF<- function (meta, pert_type = "trt_cp", dose, time = "24 h") {
  meta %<>% dplyr::filter(pert_type == pert_type & pert_idose == 
                            dose & pert_itime == time)
  meta %<>% bind_cols(alt_id = paste(meta$cmap_name, meta$cell_iname, 
                                     sep = "_")) %>% bind_cols(pert_cell_factor = paste(meta$cmap_name, 
                                                                                        meta$cell_iname, meta$pert_type, sep = "__")) %>% distinct(alt_id, 
                                                                                                                                                   .keep_all = TRUE)
  return(meta)
}


clinical_trial_check<- function(drug_list, clinical_trial_info){
  return_list <- c()
  for (i in 1:length(drug_list)){
    test <- grep( drug_list[i], clinical_trial_info$Interventions, ignore.case= TRUE)
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
