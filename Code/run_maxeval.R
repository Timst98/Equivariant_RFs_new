#!/bin/bash
#SBATCH --mail-user=tim.steinert@unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="myJob"
#SBATCH --time=96:00:00
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=all
#SBATCH --array=1-100
#SBATCH --cpus-per-task=4
#### Your shell commands below this line ####
module load R
R CMD BATCH --no-save --no-restore Rmses_time_maxeval.R