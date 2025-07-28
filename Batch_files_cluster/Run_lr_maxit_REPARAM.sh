#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="REPARAM_LR_MAXIt"
#SBATCH --time=02:00:00
#SBATCH --mem 500MB
#SBATCH --array=1-2420

module load R/4.1.0-foss-2021a 
Rscript ./LR_MAXIT_REPARAM.R 

