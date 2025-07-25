#!/bin/bash -l
#SBATCH -A uppmax2025-2-302
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J nullddswo
#SBATCH --error /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/null_distr_dds_wolf.err
#SBATCH --output /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/null_distr_dds_wolf.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load modules
module load bioinfo-tools
module load R_packages

R --no-save --quiet < wolf.dds.1000.R

