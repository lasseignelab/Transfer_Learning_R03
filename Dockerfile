FROM bioconductor/bioconductor_docker:RELEASE_3_14

# location for mounting
#RUN mkdir /home/rstudio/data


RUN R -e 'BiocManager::install(c("DESeq2", "pasilla", "SummarizedExperiment", "TFEA.ChIP", "ensembldb", "EnsDb.Hsapiens.v75", "AnnotationDbi", "org.Hs.eg.db", "cogena", "BgeeDB", "projectR", "GEOquery", "recount3", "ComplexHeatmap", "factoextra", "limma"))'

RUN R -e 'install.packages(c("readr","tidyverse","ggforce", "plotly", "devtools", "ProliferativeIndex", "gprofiler2","gplots","RColorBrewer", "car", "VennDiagram", "DrInsight", "FactoMineR", "viridis", "ComplexUpset"))'

RUN R -e 'devtools::install_github("wgmao/PLIER")'
RUN R -e 'devtools::install_github("yduan004/signatureSearch", build_vignettes=FALSE)'
RUN R -e 'devtools::install_url("https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_1.3.0.tar.gz")'

RUN R -e 'BiocManager::install("signatureSearchData")'
