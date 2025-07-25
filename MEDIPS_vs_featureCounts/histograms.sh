#!/bin/bash -l
#SBATCH -A uppmax2025-2-302
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 10:00:00
#SBATCH -J histograms
#SBATCH --error /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/histo.err
#SBATCH --output /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/histo.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load modules
module load bioinfo-tools
module load R_packages

R --no-save --quiet < histograms.R
#R --no-save --quiet < histograms_null_distribution.R
