#!/bin/bash  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=jfisher7@uab.edu
#SBATCH --job-name=conda_packages
#SBATCH -n 1
#SBATCH --mem-per-cpu=20000 
#SBATCH --nodes=1
#SBATCH --time=0-02:00:00  
#SBATCH --share 
#SBATCH --partition=express
#SBATCH --error=%j.%N.err.txt
#SBATCH --output=%j.%N.conda.packages.out.txt

cd /data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/script/
module load intel/2017a
module load Anaconda3/2019.10
module load OpenSSL
source activate
conda activate SR_TAU_CELL

Rscript ./221221_conda_packages.R