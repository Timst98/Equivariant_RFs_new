#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="NMF ind2"
#SBATCH --time=96:00:00
#SBATCH --ntasks=128
#SBATCH --mem 30G
#SBATCH --array 1-1000
module load R
Rscript ./Code/Learning_curves_NMF_inducing_points2.R
