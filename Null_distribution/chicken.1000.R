library(data.table)
library(edgeR)
library(doParallel)
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/")
chicken_matrix=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/chichen_RBC/design_matrix.txt")
load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/chicken_normalized_counts.rda")
list.p.values.loop.1000=list()

cl <- makeCluster(1, timeout = 60*60)  # Create a cluster of workers
registerDoParallel(cl)  # Register the cluster with foreach
cluster.1.cpm.tmm.norm.rbc=chicken.MW[,-13]

for (i in 1:1000){
  print(i)
  sampled_rows <- cluster.1.cpm.tmm.norm.rbc[sample(nrow(cluster.1.cpm.tmm.norm.rbc), 10000), ]
  cluster.11.cpm.tmm.norm.rbc <- sapply(sampled_rows, sample, simplify = TRUE)
  list.mann.whitney.RBC <- apply(cluster.11.cpm.tmm.norm.rbc, 1, function(x) {wilcox.test(x~chicken_matrix$Group)$p.value})
  p.values.mann.whitney=unlist(list.mann.whitney.RBC)
  # p.values.mann.whitney.fdr = p.adjust(((p.values.mann.whitney)), method = "BH")
  # p.values.mann.whitney.bonfe =p.adjust(((p.values.mann.whitney)), method = "bonferroni")
  saveRDS(p.values.mann.whitney,file=paste0("chunckchicken_",i,".rds"))
  # saveRDS(p.values.mann.whitney.fdr,file=paste0("chunckchicken_fdr_",i,".rds"))
  # saveRDS(p.values.mann.whitney.bonfe,file=paste0("chunckchicken_bonfe_",i,".rds"))
}

merged_chicken=lapply(1:1000,function(i){readRDS(paste0("chunckchicken_",i,".rds"))})
save(merged_chicken,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.MW.loop.chicken.rda")
file.remove(list.files(pattern="chunckchicken_.*\\.rds"))

# merged_chicken_fdr=lapply(1:1000,function(i){readRDS(paste0("chunckchicken_fdr_",i,".rds"))})
# merged_chicken_bonfe=lapply(1:1000,function(i){readRDS(paste0("chunckchicken_bonfe_",i,".rds"))})
# 
# save(merged_chicken_fdr,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.MW.fdr.loop.chicken.rda")
# 
# save(merged_chicken_bonfe,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.MW.bonfe.loop.chicken.rda")
# 
# file.remove(list.files(pattern="chunckchicken_fdr_.*\\.rds"))
# file.remove(list.files(pattern="chunckchicken_bonfe_.*\\.rds"))

# 
# 
# list.p.values.loop.1000 <- foreach(i = 1:1000, .combine = function(...) rbindlist(list(...), use.names = F), .packages = "data.table") %dopar% {
#   print(i)
#   
#   # Sample rows of cluster.1.cpm.tmm.norm.rbc
#   cluster.1.cpm.tmm.norm.rbc <- sapply(cluster.1.cpm.tmm.norm.rbc, sample, simplify = TRUE)
#   
#   # Perform Wilcoxon rank sum test on each row of cluster.1.cpm.tmm.norm.rbc
#   list.mann.whitney.RBC <- apply(cluster.1.cpm.tmm.norm.rbc, 1, function(x) {
#     counts.window.RBC <- data.table(x)
#     setnames(counts.window.RBC, "V1")
#     out <- wilcox.test(counts.window.RBC$V1~chicken_matrix$Group) #change here the design matrix
#     data.table::data.table(p_value = out$p.value)
#   })
#   
#   # Combine list.mann.whitney.RBC into a data table
#   df.RBC <- rbindlist(list.mann.whitney.RBC)
#   
#   # Remove unnecessary columns and rows
#   df.RBC <- df.RBC[!is.na(p_value)]
#   
#   list(df.RBC)
# }
# list.p.values.loop.1000=unlist(list.p.values.loop.1000)
# save(list.p.values.loop.1000,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.loop.chicken.rda")
# stopCluster(cl)
