library(data.table)
library(edgeR)
library(doParallel)
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/")
chicken_matrix=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/chichen_RBC/design_matrix.txt")
load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/count.matrix.chicken.rda")
list.p.values.loop.1000=list()

cl <- makeCluster(2)  # Create a cluster of workers
registerDoParallel(cl)  # Register the cluster with foreach
cluster.1.cpm.tmm.norm.rbc=chicken_meth_countmatrix[,-c(1:6)]
win_1000=c(1:10000)
list.p.values.loop.1000 <- foreach(i = 1:1000, .combine = function(...) rbindlist(list(...), use.names = F), .packages = c("data.table","edgeR")) %dopar% {
  print(i)
  
  # Sample rows of cluster.1.cpm.tmm.norm.rbc
  sampled_rows <- cluster.1.cpm.tmm.norm.rbc[sample(nrow(cluster.1.cpm.tmm.norm.rbc), 10000), ]
  cluster.11.cpm.tmm.norm.rbc <- sapply(sampled_rows, sample, simplify = TRUE)
  edgeR <- DGEList(counts=cluster.11.cpm.tmm.norm.rbc, genes=win_1000)
  rownames(edgeR$counts) <- rownames(edgeR$genes) <- win_1000
  edgeR <- calcNormFactors(edgeR)
  design=model.matrix(~chicken_matrix$Group)
  logCPM <- cpm(edgeR, log=TRUE)
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit, trend=TRUE)
  top.p.values.t.test.limma=topTable(fit, coef=design,number = 10000)
  top.p.values.t.test.limma$corrected_p_value = p.adjust(((top.p.values.t.test.limma$P.Value)), method = "bonferroni")
  # Combine list.mann.whitney.RBC into a data table
  list(top.p.values.t.test.limma$corrected_p_value)
}
list.p.values.loop.1000=unlist(list.p.values.loop.1000)
save(list.p.values.loop.1000,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.loop.chicken_ttest_bonfe.rda")
stopCluster(cl)
