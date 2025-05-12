#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p node -n 1
#SBATCH -t 2:00:00
#SBATCH -J pig
#SBATCH --error /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/02.peak_analysis_pig.err
#SBATCH --output /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/02.peak_analysis_pig.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

module load bioinfo-tools
module load R_packages

R --no-save --quiet < RangeCalling_pig.R
