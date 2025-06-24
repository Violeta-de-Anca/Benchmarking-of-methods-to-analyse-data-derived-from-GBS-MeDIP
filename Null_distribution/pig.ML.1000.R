library(data.table)
library(edgeR)
library(doParallel)
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/")
pig_matrix=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/chichen_RBC/design_matrix.txt")
load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/count.matrix.pig.rda")
list.p.values.loop.1000=list()

cl <- makeCluster(2)  # Create a cluster of workers
registerDoParallel(cl)  # Register the cluster with foreach
cluster.1.cpm.tmm.norm.rbc=pig_meth_countmatrix[,c(7:12,22:27)]
win_1000=c(1:10000)
list.p.values.loop.1000 <- foreach(i = 1:1000, .combine = c, .packages = c("data.table","edgeR")) %dopar% {
  print(i)
  
  sampled_rows <- cluster.1.cpm.tmm.norm.rbc[sample(nrow(cluster.1.cpm.tmm.norm.rbc), 10000), ]
  cluster.11.cpm.tmm.norm.rbc <- sapply(sampled_rows, sample, simplify = TRUE)
  # Perform Wilcoxon rank sum test on each row of cluster.1.cpm.tmm.norm.rbc
    edgeR.RBC <- DGEList(counts=(cluster.11.cpm.tmm.norm.rbc), genes=win_1000)
    rownames(edgeR.RBC$counts) <- rownames(edgeR.RBC$genes) <- win_1000
    edgeR.RBC$samples$lib.size <- colSums(edgeR.RBC$counts)
    edgeR.RBC <- calcNormFactors(edgeR.RBC)
    design=model.matrix(~pig_matrix$Group)
    edgeR.RBC <- estimateDisp(edgeR.RBC, design)
    model.fit.RBC=glmFit(edgeR.RBC, design)
    like.ratio.RBC = glmLRT(model.fit.RBC)
    top.RBC=like.ratio.RBC$table$PValue
    list(top.RBC)
    #top.RBC=topTags(like.ratio.RBC, n=10000, adjust.method = "fdr")
  
  # Combine list.mann.whitney.RBC into a data table
  #list(top.RBC$table$FDR)
}
list.p.values.loop.1000=unlist(list.p.values.loop.1000)
save(list.p.values.loop.1000,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.loop.pig_ML.rda")
#save(list.p.values.loop.1000,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.loop.pig_ML_FDR.rda")
stopCluster(cl)
