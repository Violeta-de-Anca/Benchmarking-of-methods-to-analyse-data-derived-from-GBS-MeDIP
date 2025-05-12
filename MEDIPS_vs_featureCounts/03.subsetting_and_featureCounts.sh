#!/bin/bash -l
#SBATCH -A naiss2023-22-108
#SBATCH -p core -n 1
#SBATCH -t 100:00:00
#SBATCH -J filterMQ10
#SBATCH --error /proj/naiss2024-23-57/trial.0.3.script/filter_MQ10.err
#SBATCH --output /proj/naiss2024-23-57/trial.0.3.script/filter_MQ10.out

module load bioinfo-tools
module load samtools/1.14
module load subread
module load sambamba

working_dir=/proj/naiss2024-23-57
aligned_dir=$working_dir/trial.0.3.script
count_dir=$working_dir/trial.0.3.script
saf_dir=$working_dir/trial.0.3.script

# First and additionally one can also filter by properly paired and quality map > 10
#file_list=($aligned_dir/*.sorted.bam)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.sorted.bam}
#	samtools view -f 2 -F 524 -q 10 -b -o $aligned_dir/$sample.MQ10.bam $file
#done

# Subset all the individuals to get the uniquely mapped reads and the multimap reads separatedly
#This only works if you have used bowtie2 as aligner, if you have used BWA-mem the flags for multimappers are XA:Z and SA:Z, so you will need to change this loop accordingly
#use this list if you haven't filter by MQ
#file_list=($aligned_dir/*.sorted.bam)

#use this list if you filter by MQ
#file_list=($aligned_dir/*.MQ10.bam)

#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.sorted.bam}
#	#get the uniquely mapped
#	samtools view -h $file | grep -v -e 'XS' | samtools view -b > $aligned_dir/$x.unique.bam
#	#get the multimapped
#	samtools view -h $file | awk '($1 ~ /^@/) || /\tXS:/' | samtools view -b > $aligned_dir/$x.multimap.bam
#done

#resort and index again
#file_list=($aligned_dir/*.unique.bam)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.unique.bam}
#	samtools sort -m 768M $file > $aligned_dir/$x.sorted.unique.bam
#done

#file_list=($aligned_dir/*.multimap.bam)
#for file in "${file_list[@]}"; do
#	name=${file##*/}
#	x=${name%.multimap.bam}
#	samtools sort -m 768M $file > $aligned_dir/$x.sorted.multimap.bam
#done

#file_list=($aligned_dir/*.bam)
#for file in "${file_list[@]}"; do
#	samtools index $file $file".bai"
#done

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

#multimapping run
#featureCounts -F SAF -p -B --countReadPairs -o $count_dir/count.matrix.multimapping -a $saf_dir/final_sorted.saf -M -d 0 -D 1500 -T 1 $aligned_dir/*.sorted.multimap.bam

#uniquely map run
#featureCounts -F SAF -p -B --countReadPairs -o $count_dir/count.matrix.unique -a $saf_dir/final_sorted.saf -d 0 -D 1500 -T 1 $aligned_dir/*.sorted.unique.bam

#Now we need to do a merge of both files, but first we need to add a column stating from
header="unique_mapped"

line_count=$(($(wc -l < "$count_dir/count.matrix.multimapping") - 2))
echo $header > $count_dir/false
for (( i=1; i<=line_count; i++ ))
do
    echo "FALSE" >> $count_dir/false
done
tail -n +2 $count_dir/count.matrix.multimapping > $count_dir/count.matrix.multimapping.noheader
paste $count_dir/false $count_dir/count.matrix.multimapping.noheader > $count_dir/add.col.count.matrix.multimapping


line_count=$(($(wc -l < "$count_dir/count.matrix.unique") - 2))
echo $header > $count_dir/true
for (( i=1; i<=line_count; i++ ))
do
    echo "TRUE" >> $count_dir/true
done
tail -n +2 $count_dir/count.matrix.unique > $count_dir/count.matrix.unique.noheader
paste $count_dir/true $count_dir/count.matrix.unique.noheader > $count_dir/add.col.count.matrix.unique

tail -n +2 $count_dir/add.col.count.matrix.unique >> $count_dir/add.col.count.matrix.multimapping

mv $count_dir/add.col.count.matrix.multimapping $count_dir/count_matrix_merged_unique_multimap.txt

rm $count_dir/true $count_dir/false $count_dir/add.col.count.matrix.unique $count_dir/count.matrix.multimapping.noheader $count_dir/count.matrix.unique.noheader


