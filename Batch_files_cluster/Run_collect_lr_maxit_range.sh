#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="collect"
#SBATCH --partition=bdw
#SBATCH --time=00:01:00
#SBATCH --mem 1G
module load R
Rscript ./collect_lr_maxit_range.R 


