#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="collect"
#SBATCH --time=00:04:00
#SBATCH --mem 1G
module load R/4.1.0-foss-2021a
Rscript ./COLLECT_LOGS.R
