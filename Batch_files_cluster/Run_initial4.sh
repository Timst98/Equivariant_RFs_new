#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="final"
#SBATCH --time=24:00:00
#SBATCH --mem 1G
#SBATCH --array=1-2000
module load R
Rscript ./LogS_scores.R 

