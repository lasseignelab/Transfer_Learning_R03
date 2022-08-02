#!/bin/bash  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=jfisher7@uab.edu
#SBATCH --job-name=tau_test
#SBATCH -n 1
#SBATCH --mem-per-cpu=80000 
#SBATCH --nodes=1
#SBATCH --time=0-02:00:00  
#SBATCH --share 
#SBATCH --partition=express
#SBATCH --error=%j.%N.err.txt
#SBATCH --output=%j.%N.out.txt

cd /data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/script/SR_TAU_LINCS_2020
module load intel/2017a
module load R/4.1.0-foss-2021a

Rscript ./220125_TAU_calc_LINCS_step2_v2.R 1