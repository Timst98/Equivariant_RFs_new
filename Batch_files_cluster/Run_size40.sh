#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="size 40"
#SBATCH --time=96:00:00
#SBATCH --mem 1G
#SBATCH --array 1-2000
module load R
Rscript ./validation_ll_chemical_size40.R
