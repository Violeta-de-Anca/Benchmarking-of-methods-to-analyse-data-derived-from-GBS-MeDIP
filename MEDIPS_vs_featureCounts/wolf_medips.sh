#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core
#SBATCH -n 4
#SBATCH -t 10-00:00:00
#SBATCH -J medipswolf
#SBATCH --error /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/medips_wolf.err
#SBATCH --output /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/medips_wolf.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load modules
module load bioinfo-tools
module load R_packages/4.2.1

R --no-save --quiet < wolf_medips.R

