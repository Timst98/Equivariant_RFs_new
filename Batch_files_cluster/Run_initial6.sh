#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="Example 2"
#SBATCH --time=72:00:00
#SBATCH --mem 1G

module load R
Rscript ./first_example_plots_2.R 

