#!/bin/bash -l
#SBATCH -A naiss2023-22-848
#SBATCH -p core -n 1
#SBATCH -t 10:00:00
#SBATCH -J pig.3
#SBATCH --error /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/pig_count_matrix.err
#SBATCH --output /proj/naiss2024-23-57/GBS_MeDIP_benchmark/log_files/pig_count_matrix.out

module load bioinfo-tools
module load samtools/1.14
module load subread
module load sambamba

working_dir=/proj/naiss2024-23-57/GBS_MeDIP_benchmark
aligned_dir=$working_dir/datasets/pig_sperm
count_dir=$working_dir/benchmark_final_tables
saf_dir=$working_dir/datasets/pig_sperm/merged

# First and additionally one can also filter by properly paired and quality map > 10
file_list=($aligned_dir/*.reheadered.bam)
for file in "${file_list[@]}"; do
	name=${file##*/}
	x=${name%.sorted.reheadered.bam}
	samtools view -f 2 -F 524 -q 10 -b -o $aligned_dir/$x.MQ10.bam $file
done

# Subset all the individuals to get the uniquely mapped reads and the multimap reads separatedly
#This only works if you have used bowtie2 as aligner, if you have used BWA-mem the flags for multimappers are XA:Z and SA:Z, so you will need to change this loop accordingly
#use this list if you haven't filter by MQ
#file_list=($aligned_dir/*.sorted.bam)

#use this list if you filter by MQ
file_list=($aligned_dir/*.MQ10.bam)

#Getting the counts from FeatureCounts
#For the multimapping both primary and secondary alignments will be considered in counting, -M allows multimapping to be counted
# -F specifies the format
# -p specifies that input data contain paired-end reads
# -B specifies that only fragments that have both ends successfully aligned will be considered for counting
# -d specifies minimum fragment length
# -D specifies maximum fragment length
# -o specifies output
# -T specifies number of threads used
#we can specify the samples using wildcards

featureCounts -F SAF -p -B --countReadPairs -o $count_dir/count.matrix.pig.bed -a $saf_dir/final_sorted.saf -d 0 -D 1500 -T 1 $aligned_dir/*.MQ10.bam


