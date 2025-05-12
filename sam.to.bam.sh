#!/bin/bash -l
#SBATCH -A naiss2023-22-108
#SBATCH -p core -n 4
#SBATCH --mem=80gb
#SBATCH -t 100:00:00
#SBATCH -J sam2bam
#SBATCH --error /proj/naiss2023-23-55/GBS_MeDIP_benchmark/sam2bam.err
#SBATCH --output /proj/naiss2023-23-55/GBS_MeDIP_benchmark/sam2bam.out

module load bioinfo-tools

module load Stacks/2.62
module load cutadapt/4.0
module load MultiQC
module load samtools/1.14
module load bowtie2/2.3.5.1
module load FastQC/0.11.9

#In the working directory folder should be a folder called barcodes and inside a tab separated file with barcode sequence and name
working_dir=/proj/naiss2023-23-55/GBS_MeDIP_benchmark

cd $working_dir/aligned

#for file in *.sam; do
#    samtools view -S -b $file > ${file/%sam/bam}
#done
#removing sam files
#rm *.sam

for file in *.bam; do
    samtools sort -m 768M $file > ${file/%bam/sorted.bam}
done
#removing only bam files and not sorted bam
#ls *.bam| grep -v .sorted.bam$| xargs rm

for i in *sorted.bam
do
  echo "Indexing: "$i
  samtools index $i $i.bai
done

for file in *.sorted.bam; do
  samtools depth $file | awk '{sum+=$3;cnt++} END {print sum/cnt" "sum}'
done

echo "finished SAMBAM"
