#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="NMF very long"
#SBATCH --time=96:00:00
#SBATCH --ntasks=101
#SBATCH --mem 10G
#SBATCH --array 1-100
module load R
Rscript ./Learning_curves_NMF_parallel_very_long.R
