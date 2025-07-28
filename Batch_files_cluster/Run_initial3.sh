#!/bin/bash
#SBATCH --job-name="Initial_rmses"
#SBATCH --time=24:00:00
#SBATCH --mem 1G
#SBATCH --array=1-100
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err

module load R/4.1.0-foss-2021a
Rscript ./rmses_initial_adam2.R 
