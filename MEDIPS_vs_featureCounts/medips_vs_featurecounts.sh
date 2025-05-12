#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core
#SBATCH -n 3
#SBATCH -t 10-00:00:00
#SBATCH -J medipsvsfeat
#SBATCH --error /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/medips_vs_featurecounts.err
#SBATCH --output /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/medips_vs_featurecounts.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load modules
module load bioinfo-tools
module load R_packages

R --no-save --quiet < plots_correlation_medips_featurecounts.R

