library(data.table)
library("DESeq2")
library("DEGreport")
library(doParallel)
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/")
chicken_matrix=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/chichen_RBC/design_matrix.txt")
load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/count.matrix.chicken.rda")
list.p.values.loop.1000=list()

cl <- makeCluster(2)  # Create a cluster of workers
registerDoParallel(cl)  # Register the cluster with foreach
cluster.1.cpm.tmm.norm.rbc=chicken_meth_countmatrix[,-c(1:6)]
win_1000=c(1:10000)
list.p.values.loop.1000 <- foreach(i = 1:1000, .combine = function(...) rbindlist(list(...), use.names = F), .packages = c("data.table","DESeq2","DEGreport")) %dopar% {
  print(i)
  
  # Sample rows of cluster.1.cpm.tmm.norm.rbc
  sampled_rows <- cluster.1.cpm.tmm.norm.rbc[sample(nrow(cluster.1.cpm.tmm.norm.rbc), 10000), ]
  cluster.11.cpm.tmm.norm.rbc <- sapply(sampled_rows, sample, simplify = TRUE)
  
  # Perform deseq2 pipeline on each row of cluster.1.cpm.tmm.norm.rbc
  dds=DESeqDataSetFromMatrix(countData = cluster.11.cpm.tmm.norm.rbc,colData = df_design,design = ~ group)
  dds=DESeq(dds)
  res=results(dds)
  deseq.pvalue=res$pvalue
  deseq.pvalue=data.frame(deseq.pvalue)
  colnames(deseq.pvalue)=c("p_value")
  row.names(deseq.pvalue)=row.names(dds_mat)
  deseq.pvalue$corrected_p_value = p.adjust(((deseq.pvalue$p_value)), method = "bonferroni")
  list(deseq.pvalue)
  # Combine list.mann.whitney.RBC into a data table
  #list(top.quasi.RBC$table$FDR)
}
list.p.values.loop.1000=unlist(list.p.values.loop.1000)
save(list.p.values.loop.1000,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.loop.chicken_deseq2_bonferroni.rda")
stopCluster(cl)