library(data.table)
library(edgeR)
library(doParallel)
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/")
wolf_matrix=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/wolf/design_matrix.txt")
load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/wolf_normalized_counts.rda")
list.p.values.loop.1000=list()

# cl <- makeCluster(1, timeout = 60*60)  # Create a cluster of workers
# registerDoParallel(cl)  # Register the cluster with foreach
cluster.1.cpm.tmm.norm.rbc=wolf.MW[,-7]

for(i in 1:1000){
  print(i)
  sampled_rows <- cluster.1.cpm.tmm.norm.rbc[sample(nrow(cluster.1.cpm.tmm.norm.rbc), 10000), ]
  cluster.11.cpm.tmm.norm.rbc <- sapply(sampled_rows, sample, simplify = TRUE)
  list.mann.whitney.RBC <- apply(cluster.11.cpm.tmm.norm.rbc, 1, function(x) {wilcox.test(x~wolf_matrix$Group)$p.value})
  p.values.mann.whitney=unlist(list.mann.whitney.RBC)
  saveRDS(p.values.mann.whitney,file=paste0("chunckwolf_",i,".rds"))
  # p.values.mann.whitney.fdr = p.adjust(((p.values.mann.whitney)), method = "BH")
  # p.values.mann.whitney.bonfe =p.adjust(((p.values.mann.whitney)), method = "bonferroni")
  # saveRDS(p.values.mann.whitney.fdr,file=paste0("chunckwolf_fdr_",i,".rds"))
  # saveRDS(p.values.mann.whitney.bonfe,file=paste0("chunckwolf_bonfe_",i,".rds"))

}

merged_wolf=lapply(1:1000,function(i){readRDS(paste0("chunckwolf_",i,".rds"))})
save(merged_wolf,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.MW.loop.wolf.rda")
file.remove(list.files(pattern="chunckwolf_.*\\.rds"))




# merged_wolf_fdr=lapply(1:1000,function(i){readRDS(paste0("chunckwolf_fdr_",i,".rds"))})
# merged_wolf_bonfe=lapply(1:1000,function(i){readRDS(paste0("chunckwolf_bonfe_",i,".rds"))})
# 
# save(merged_wolf_fdr,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.MW.fdr.loop.wolf.rda")
# 
# save(merged_wolf_bonfe,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.MW.bonfe.loop.wolf.rda")
# 
# file.remove(list.files(pattern="chunckwolf_fdr_.*\\.rds"))
# file.remove(list.files(pattern="chunckwolf_bonfe_.*\\.rds"))


# list.p.values.loop.1000 <- foreach(i = 1:500, .combine = function(...) rbindlist(list(...), use.names = F), .packages = "data.table") %dopar% {
#   print(i)
#   
#   # Sample rows of cluster.1.cpm.tmm.norm.rbc
#   cluster.1.cpm.tmm.norm.rbc <- sapply(cluster.1.cpm.tmm.norm.rbc, sample, simplify = TRUE)
#   
#   # Perform Wilcoxon rank sum test on each row of cluster.1.cpm.tmm.norm.rbc
#   list.mann.whitney.RBC <- apply(cluster.1.cpm.tmm.norm.rbc, 1, function(x) {wilcox.test(x~wolf_matrix$Group)$p.value})
#    # counts.window.RBC <- data.table(x)
#    # setnames(counts.window.RBC, "V1")
#    # out= apply((edgeR.MW), 1, function(x){wilcox.test(x~df_design$group)$p.value})
#    # out <- wilcox.test(counts.window.RBC$V1~wolf_matrix$Group) #change here the design matrix
#    # data.table::data.table(p_value = out$p.value)
#  # })
#   
#   # Combine list.mann.whitney.RBC into a data table
#   df.RBC <- rbindlist(list.mann.whitney.RBC)
#   
#   # Remove unnecessary columns and rows
#   df.RBC <- df.RBC[!is.na(p_value)]
#   
#   list(df.RBC)
# }

#list.p.values.loop.1000=unlist(list.p.values.loop.1000)
#save(list.p.values.loop.1000,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.loop.wolf.rda")

