#!/bin/bash  
#SBATCH --mail-type=ALL  
#SBATCH --mail-user=jfisher7@uab.edu
#SBATCH --job-name=tau_1PC3
#SBATCH -n 1
#SBATCH --mem-per-cpu=20000 
#SBATCH --nodes=1
#SBATCH --time=0-02:00:00  
#SBATCH --share 
#SBATCH --partition=express
#SBATCH --error=%j.%N.err.txt
#SBATCH --output=%j.%N.out.txt

cd /data/project/lasseigne_lab/JLF_scratch/Transfer_Learning_R03/script/SR_TAU_LINCS_2020
module load intel/2017a
module load Anaconda3/2019.10
module load OpenSSL
source activate
conda activate SR_TAU_CELL

Rscript ./220330_tau_calc_PC3.R $SLURM_ARRAY_TASK_ID