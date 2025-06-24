library(data.table)
library(edgeR)
library(doParallel)
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/")
pig_matrix=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/datasets/pig_sperm/design_matrix.txt")
load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/pig_normalized_counts.rda")
list.p.values.loop.1000=list()
 
# # cl <- makeCluster(1, timeout = 60*60)  # Create a cluster of workers
# # registerDoParallel(cl)  # Register the cluster with foreach
cluster.1.cpm.tmm.norm.rbc=pig.MW[,-27]
 
 for (i in 1:1000){
 	print(i)
 	sampled_rows <- cluster.1.cpm.tmm.norm.rbc[sample(nrow(cluster.1.cpm.tmm.norm.rbc), 10000), ]
 	cluster.11.cpm.tmm.norm.rbc <- sapply(sampled_rows, sample, simplify = TRUE)
 	list.mann.whitney.RBC <- apply(cluster.11.cpm.tmm.norm.rbc, 1, function(x) {wilcox.test(x~pig_matrix$Group)$p.value})
 	p.values.mann.whitney=unlist(list.mann.whitney.RBC)
 	saveRDS(p.values.mann.whitney,file=paste0("chunckpig_",i,".rds"))
 	# p.values.mann.whitney.fdr = p.adjust(((p.values.mann.whitney)), method = "BH")
 	# p.values.mann.whitney.bonfe =p.adjust(((p.values.mann.whitney)), method = "bonferroni")
 	# saveRDS(p.values.mann.whitney.fdr,file=paste0("chunckpig_fdr_",i,".rds"))
 	# saveRDS(p.values.mann.whitney.bonfe,file=paste0("chunckpig_bonfe_",i,".rds"))
 }

merged_pig=lapply(1:1000,function(i){readRDS(paste0("chunckpig_",i,".rds"))})
save(merged_pig,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.MW.loop.pig.rda")
file.remove(list.files(pattern="chunckpig_.*\\.rds"))


# merged_pig_fdr=lapply(1:1000,function(i){readRDS(paste0("chunckpig_fdr_",i,".rds"))})
# merged_pig_bonfe=lapply(1:1000,function(i){readRDS(paste0("chunckpig_bonfe_",i,".rds"))})
# 
# save(merged_pig_fdr,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.MW.fdr.loop.pig.rda")
# 
# save(merged_pig_bonfe,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.MW.bonfe.loop.pig.rda")
# 
# file.remove(list.files(pattern="chunckpig_fdr_.*\\.rds"))
# file.remove(list.files(pattern="chunckpig_bonfe_.*\\.rds"))

# list.p.values.loop.1000 <- foreach(i = 1:500, .combine = function(...) rbindlist(list(...), use.names = F), .packages = "data.table") %dopar% {
#   print(i)
#   
#   # Sample rows of cluster.1.cpm.tmm.norm.rbc
#   cluster.1.cpm.tmm.norm.rbc <- sapply(cluster.1.cpm.tmm.norm.rbc, sample, simplify = TRUE)
#   
#   # Perform Wilcoxon rank sum test on each row of cluster.1.cpm.tmm.norm.rbc
#   list.mann.whitney.RBC <- apply(cluster.1.cpm.tmm.norm.rbc, 1, function(x) {wilcox.test(x~~pig_matrix$Group)$p.value})
# 
# #    counts.window.RBC <- data.table(x)
# #    setnames(counts.window.RBC, "V1")
# #    out <- wilcox.test(counts.window.RBC$V1~pig_matrix$Group) #change here the design matrix
# #    data.table::data.table(p_value = out$p.value)
# #  })
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
#save(list.p.values.loop.1000,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/list.1000.p.values.loop.pig.rda")

