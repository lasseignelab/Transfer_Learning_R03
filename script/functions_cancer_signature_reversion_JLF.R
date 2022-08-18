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
    test <- grep( drug_list[i], fda_approve$ActiveIngredient, ignore.case= TRUE)
    if (length(test)>0 ){
      return_list[i] <- TRUE
    }else{
      return_list[i] <- FALSE
    }
    #more than 4 character are needed for matching
    if (nchar(drug_list[i])< 4){
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
  downset_pathway_results <- gost(query = down, 
                                  organism = "hsapiens", ordered_query = TRUE, 
                                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                                  user_threshold = 0.05, correction_method = "g_SCS", 
                                  domain_scope = "annotated", custom_bg = NULL, 
                                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
  downset_pathway_results<- downset_pathway_results$result
  
  
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
clinical_trial_testing<- function(seed, cell_line,number_sig_drugs, fraction_clinical_drugs, CT_DATA){
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
PRISM_testing<- function(seed, cell_line, drugs, drug_senstive_precentage_all_drugs){
  
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


venn_dia_methods <- function(limma_list, deseq2_list, tfl_list, file_name='~/output/liver_cancer/SR_liver_all_gene_venn_diagram.png' ){
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


method_up_down_gene_set_analysis<- function(up, down,result_path, method, name){
  library(gprofiler2)
  downset_pathway_results <- gost(query = down, 
                                  organism = "hsapiens", ordered_query = TRUE, 
                                  multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                                  measure_underrepresentation = FALSE, evcodes = FALSE, 
                                  user_threshold = 0.05, correction_method = "g_SCS", 
                                  domain_scope = "annotated", custom_bg = NULL, 
                                  numeric_ns = "", sources = NULL, as_short_link = FALSE)
  downset_pathway_results<- downset_pathway_results$result
  
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
  
  #making adjustments here
  #pathway_results$set <- factor(pathway_results$set, levels = c("up", "down"))
  if(nrow(pathway_results)<50){
    ggplot(pathway_results, aes(x=set, y=term_name, fill=p_value))+geom_tile()+scale_fill_viridis(direction=-1)+theme_classic()+ labs(title=name,x="Drug Gene Sets", y = "g:Profiler Gene Sets", fill= "p-value")
    
  }else{
    pathway_results<- pathway_results[pathway_results$source %in% c("KEGG","REAC","GO:MF" ),]
    ggplot(pathway_results, aes(x=set, y=term_name, fill=p_value)) + geom_tile()+ scale_fill_viridis(direction=-1)+ theme_classic() + labs(title=name,x="Drug Gene Sets", y = "g:Profiler Gene Sets", fill= "p-value")
  }
  ggsave(file_name, width = 12, height = 14)
  return(as.data.frame(pathway_results[,-14]))
}

#functions from rrvgo
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

go_term_heatmap<- function(limma_pathways, deseq2_pathways, TFL_pathways, list_type, threshold = 0.95){
  
  #get terms
  limma_bp<- limma_pathways$term_id[ limma_pathways$source == "GO:BP" & limma_pathways$set == list_type]
  deseq2_bp <- deseq2_pathways$term_id[ deseq2_pathways$source == "GO:BP" & deseq2_pathways$set == list_type]
  TFL_bp <- TFL_pathways$term_id[ TFL_pathways$source == "GO:BP" & TFL_pathways$set == list_type]
  
  go1 <-  unique(c(TFL_bp, deseq2_bp, limma_bp ))
  go_gbm_up_sim <- mgoSim(go1, go1, semData=hsGO, measure="Wang", combine=NULL)
  
  
  res <- reduceSimMatrix(go_gbm_up_sim, threshold = threshold)
  res_v2<- res[match(colnames(go_gbm_up_sim), res$go),]
  
  tfl_drugs <- grepl(paste(TFL_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  deseq2_drugs<- grepl(paste(deseq2_bp,collapse="|"), colnames(go_gbm_up_sim))
  limma_drugs <- grepl(paste(limma_bp,collapse="|"), colnames(go_gbm_up_sim)) 
  
  if(length(deseq2_bp) ==0){ deseq2_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  if(length(TFL_bp) ==0){ tfl_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  if(length(limma_bp) ==0){ limma_drugs <- rep(FALSE,ncol(go_gbm_up_sim) )} 
  
  #bp_color<- rainbow(n= length(unique(res_v2$parentTerm)), alpha=0.8)
  bp_color<- brewer.pal(n = length(unique(res_v2$parentTerm)), name = "Paired")
  names(bp_color)<- unique(res_v2$parentTerm)
  
  row_ha = HeatmapAnnotation(Transfer_Learning=tfl_drugs,DESeq2= deseq2_drugs, limma= limma_drugs,GO_BP_Group= res_v2$parentTerm , col = list(Transfer_Learning = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),DESeq2 = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"),limma = c("TRUE" = "#440154FF", "FALSE" = "#228C8DFF"), GO_BP_Group= bp_color ))
  col_fun = colorRamp2(c(0,  1), c( "black", "yellow"))
  #ha = rowAnnotation(foo = anno_mark(at = c(1:4, 20, 60, 97:100), labels = TFL_bp_up[1:10]))
  
  Heatmap(go_gbm_up_sim, nam= "GO Term Similarity (Wang)", col = col_fun, show_column_names = FALSE,  show_row_names = FALSE, top_annotation = row_ha,  
          clustering_distance_rows= "euclidean",
          clustering_distance_columns=  "euclidean",
          clustering_method_rows = "ward.D2" ,
          clustering_method_columns="ward.D2")
  #file_name<- paste0(file_starter, "_", "Go_sem_heatmap.png")
  #ggsave(file_name, width = 20, height=20, units= "in")
}


go_term_tree_plot<- function(limma_pathways, deseq2_pathways, TFL_pathways, list_type, num=80, k_num=5, colors= c( "#117733", "#661100",  "#0072B2", "#D55E00", "#AA4499"), shift = 4){
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

drug_moa_target_plotting<- function(fda_approved_res, lincs_commpound_info=lincs_commpound_info, method, cancer, file_path_header ){
  
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
  #make the self match equal to zero
  
  cos_matrix[cos_matrix== 1]<- 0
  
  #add 1 to quantile because of 0 at index 1
  quantile <-quantile+1
  cutoff <- quantile(cos_matrix,probs = seq(0, 1, 0.1))[quantile]
  
  #set values below value to 0 
  cos_matrix[cos_matrix < cutoff ]<- 0 
  colnames(cos_matrix)<- rownames(cos_matrix)
  
  #remove candidate with zero edge weights
  cos_matrix_v2<- cos_matrix[ !colSums(cos_matrix) == 0, !colSums(cos_matrix) == 0]
  
  #create network
  cos_cand_igraph <- graph_from_adjacency_matrix(cos_matrix_v2 ,mode= "undirected", weighted=TRUE)
  
  #cluster analysis 
  leiden_clusters <- groups(cluster_leiden(cos_cand_igraph,n_iterations = 10,  objective_function ="modularity"))
  
  #determine which candidate is in which cluster
  com <-rep(NA, nrow(cos_matrix_v2))
  for (i in 1: nrow(cos_matrix_v2)){
    com<- ifelse(colnames(cos_matrix_v2) %in% unlist(leiden_clusters[i]), i , com )
  }
  com
  
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
  E(cos_cand_igraph)$color <- "black"
  
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
  method<- c("Transfer Learning",  "DESeq2", "limma")
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
                    clustering_method_columns="ward.D2",column_names_rot = 45, 
                    layer_fun = function(j, i, x, y, width, height, fill) {
                      v = pindex(t(counts_wider), i, j)
                      grid.text(sprintf("%.0f", v), x, y, gp = gpar(fontsize = 10))
                      if(min(i)==1) {
                        grid.rect(gp = gpar(lwd = 2, fill = "transparent",col="black"))
                      }})
  print(heatmap)
}






