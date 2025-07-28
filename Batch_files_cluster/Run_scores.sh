#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="2Scores_all"
#SBATCH --time=03:00:00
#SBATCH --mem 1G
#SBATCH --array=1-4000
#SBATCH --partition=epyc2

module load R
Rscript ./Scores_all_2.R 

