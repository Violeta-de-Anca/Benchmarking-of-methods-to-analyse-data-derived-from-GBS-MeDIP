#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10-00:00:00
#SBATCH -J MLbonfchicken
#SBATCH --error /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/null_distr_ML_bonfe_chicken.err
#SBATCH --output /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/null_distr_ML_bonfe_chicken.out
#SBATCH --mail-type=FAIL,BEGIN
#SBATCH --mail-user=violeta.deancaprado@ebc.uu.se

#load modules
module load bioinfo-tools
module load R_packages

R --no-save --quiet < chicken.1000.ML.bonfe.R

