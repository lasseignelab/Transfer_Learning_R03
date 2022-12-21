
cd ../data/recount3/recount3_fix_download
#download the counts

wget http://duffel.rail.bio/recount3/human/data_sources/tcga/gene_sums/AD/PAAD/tcga.gene_sums.PAAD.G026.gz

wget http://duffel.rail.bio/recount3/human/data_sources/tcga/gene_sums/HC/LIHC/tcga.gene_sums.LIHC.G026.gz

wget http://duffel.rail.bio/recount3/human/data_sources/tcga/gene_sums/BM/GBM/tcga.gene_sums.GBM.G026.gz

wget http://duffel.rail.bio/recount3/human/data_sources/tcga/gene_sums/AD/LUAD/tcga.gene_sums.LUAD.G026.gz

wget http://duffel.rail.bio/recount3/human/data_sources/gtex/gene_sums/IN/BRAIN/gtex.gene_sums.BRAIN.G026.gz

wget http://duffel.rail.bio/recount3/human/data_sources/gtex/gene_sums/AS/PANCREAS/gtex.gene_sums.PANCREAS.G026.gz

wget http://duffel.rail.bio/recount3/human/data_sources/gtex/gene_sums/ER/LIVER/gtex.gene_sums.LIVER.G026.gz

wget http://duffel.rail.bio/recount3/human/data_sources/gtex/gene_sums/NG/LUNG/gtex.gene_sums.LUNG.G026.gz

#download the metadata
wget https://sciserver.org/public-data/recount3/data/human/data_sources/tcga/metadata/AD/PAAD/tcga.tcga.PAAD.MD.gz 

wget https://sciserver.org/public-data/recount3/data/human/data_sources/tcga/metadata/AD/LUAD/tcga.tcga.LUAD.MD.gz 

wget https://sciserver.org/public-data/recount3/data/human/data_sources/tcga/metadata/BM/GBM/tcga.tcga.GBM.MD.gz 

wget https://sciserver.org/public-data/recount3/data/human/data_sources/tcga/metadata/HC/LIHC/tcga.tcga.LIHC.MD.gz 

wget https://sciserver.org/public-data/recount3/data/human/data_sources/gtex/metadata/AS/PANCREAS/gtex.gtex.PANCREAS.MD.gz

wget https://sciserver.org/public-data/recount3/data/human/data_sources/gtex/metadata/ER/LIVER/gtex.gtex.LIVER.MD.gz

wget https://sciserver.org/public-data/recount3/data/human/data_sources/gtex/metadata/IN/BRAIN/gtex.gtex.BRAIN.MD.gz

wget https://sciserver.org/public-data/recount3/data/human/data_sources/gtex/metadata/NG/LUNG/gtex.gtex.LUNG.MD.gz

#unzip the files 
gunzip *.gz

#download the gif.gz file 
wget http://duffel.rail.bio/recount3/human/annotations/gene_sums/human.gene_sums.G026.gtf.gz