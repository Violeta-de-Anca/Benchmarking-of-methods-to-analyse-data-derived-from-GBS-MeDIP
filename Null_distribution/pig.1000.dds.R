library(data.table)
library(edgeR)
library(doParallel)
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/")
pig_matrix=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/pig_sperm/design_matrix.txt")
load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/count.matrix.pig.rda")
list.p.values.loop.1000=list()

cl <- makeCluster(2)  # Create a cluster of workers
registerDoParallel(cl)  # Register the cluster with foreach
cluster.1.cpm.tmm.norm.rbc=pig_meth_countmatrix[,-c(1:6)]
colnames(cluster.1.cpm.tmm.norm.rbc)=sub("^/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/pig_sperm/","",colnames(cluster.1.cpm.tmm.norm.rbc))
colnames(cluster.1.cpm.tmm.norm.rbc)=sub(".MQ10.bam","",colnames(cluster.1.cpm.tmm.norm.rbc))
colnames(cluster.1.cpm.tmm.norm.rbc)=sub("HF_","",colnames(cluster.1.cpm.tmm.norm.rbc))
colnames(cluster.1.cpm.tmm.norm.rbc)=sub("LF_","",colnames(cluster.1.cpm.tmm.norm.rbc))
pign=pig_matrix$Individual
cluster.1.cpm.tmm.norm.rbc=cluster.1.cpm.tmm.norm.rbc[, ..pign]
win_1000=c(1:10000)
list.p.values.loop.1000 <- foreach(i = 1:1000, .combine = c, .packages = c("data.table","DESeq2","DEGreport")) %dopar% {
  print(i)
  
  # Sample rows of cluster.1.cpm.tmm.norm.rbc
  sampled_rows <- cluster.1.cpm.tmm.norm.rbc[sample(nrow(cluster.1.cpm.tmm.norm.rbc), 10000), ]
  shuffled=sample(colnames(sampled_rows))
  sampled_rows=sampled_rows[, ..shuffled]
  sampled_rows[]=lapply(sampled_rows, as.numeric)
  pig_matrix$Individual=colnames(sampled_rows)
  # Perform deseq2 pipeline on each row of cluster.1.cpm.tmm.norm.rbc
  dds=DESeqDataSetFromMatrix(countData = sampled_rows,colData = pig_matrix,design = ~ Group)
  dds=DESeq2::estimateSizeFactors(dds,type = "poscounts")
  dds=DESeq(dds)
  res=results(dds)
  deseq.pvalue=res$pvalue
  deseq.pvalue=data.frame(deseq.pvalue)
  colnames(deseq.pvalue)=c("p_value")
  list(deseq.pvalue)
}
list.p.values.loop.1000=unlist(list.p.values.loop.1000)
save(list.p.values.loop.1000,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.loop.pig.dds.rda")
stopCluster(cl)