#!/bin/bash -l
#SBATCH -A uppmax2025-2-222
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 10:00:00
#SBATCH -J fitdistrplus
#SBATCH --error /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/dist_fitdistrplus.err
#SBATCH --output /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/dist_fitdistrplus.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load modules
module load bioinfo-tools
module load R_packages/4.2.1

R --no-save --quiet < data_distribution_fitdistrplus.R

