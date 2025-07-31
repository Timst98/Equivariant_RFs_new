#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="NMF ind2"
#SBATCH --time=22:00:00
#SBATCH --cpus-per-task=50
#SBATCH --ntasks=1
#SBATCH --mem 50G
#SBATCH --array 1-1000

module load R/4.4.2-gfbf-2024a


Rscript ./Code/Learning_curves_NMF_inducing_points2.R
