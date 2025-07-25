#!/bin/bash -l
#SBATCH -A uppmax2025-2-302
#SBATCH -p core -n 1
#SBATCH -t 10:00:00
#SBATCH -J bedtools
#SBATCH --error /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/bedtools_all.err
#SBATCH --output /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/bedtools_all.out

module load bioinfo-tools
module load samtools/1.14
module load subread
module load sambamba BEDTools/2.31.1


working_dir=/proj/naiss2024-23-57/GBS_MeDIP_benchmark
aligned_dir=$working_dir/datasets
saf_dir=$working_dir/datasets/$1/merged
output_folder=$working_dir/benchmark_final_tables

#define variable with all the bams
bam_files=($aligned_dir/$1/*.MQ10.bam)

#for i in $aligned_dir/$1/*.MQ10.bam;do
#	samtools index -b $i
#done

tail -n +2 $saf_dir/final_sorted.saf | cut -f 2-4 > $saf_dir/bed_$1.bed

bedtools multicov -p -q 10 -bed $saf_dir/bed_$1.bed ${bam_files[@]/#/-bams } > $output_folder/${1}_bedtools_counts.bed
