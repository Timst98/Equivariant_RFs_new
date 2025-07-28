#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="Range_initial"
#SBATCH --time=00:40:00
#SBATCH --mem 500MB
#SBATCH --array=1-4000

module load R/4.1.0-foss-2021a 
Rscript ./rmses_initial_point2_140_range.R 

