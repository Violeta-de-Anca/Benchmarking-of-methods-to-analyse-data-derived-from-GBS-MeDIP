setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged/")
library(fitdistrplus)
library(data.table)
library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grDevices)
library(gridGraphics)
library(grid)
#to try and see which distribution works best, also you need to say that is discrete
#descdist()

#NON normalized data - chicken ####
#load data, but it has to be in a vector, let's see per individual####
# load("chicken_RBC/meth_countmatrix_chicken.rda")
# meth_countmatrix=meth_countmatrix[rowSums(meth_countmatrix[,-c(1:3)])>=1,]
# for(i in seq_len(ncol(meth_countmatrix[,-c(1:3)]))){
#   dir_out="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/"
#   file_names=file.path(dir_out, paste0("chicken_individual_distribution",i,".tiff"))
#   tiff(file_names,width = 1500, height = 1600, units = "px",res = 300)
#   f=as.numeric(meth_countmatrix[,-c(1:3)][,i])
#   descdist(f,discrete = T)
#   dev.off()
# }
# 
# #now let's do it for each window, it gives error if the window is all 0 ####
# meth_countmatrix=meth_countmatrix[rowSums(meth_countmatrix[,-c(1:3)])>=1,]
# for(i in 1:100){
#   dir_out="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/"
#   file_names=file.path(dir_out, paste0("chicken_window_distribution",i,".tiff"))
#   tiff(file_names,width = 1500, height = 1600, units = "px",res = 300)
#   f=as.numeric(meth_countmatrix[,-c(1:3)][i,])
#   descdist(f,discrete = T)
#   dev.off()
# }
# 
# #turns out, by individual and by window the distribution changes?
# # so i am going to use fitdist() as it will give a numerial value ####
# AIC_chicken_windows=list()
# models=c("pois","nbinom","norm")
# for(i in models){
#   AIC_chicken_windows[[i]]=list()
#   for (b in 1:nrow(meth_countmatrix)){
#     print(b)
#     f=as.numeric(meth_countmatrix[,-c(1:3)][b,])
#     AIC_chicken_windows[[i]][[b]]=fitdist(f,i)
#     
#   }
# }
# save(AIC_chicken_windows,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/AIC_non_normalized_chicken.rda")
# load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/AIC_non_normalized_chicken.rda")
# models=c("pois","nbinom","norm")
# #extract the AIC
# results_AIC_chicken=matrix(NA,nrow = 209452, ncol = length(models))
# colnames(results_AIC_chicken)=models
# for (i in seq_along(models)){
#   model_name=models[i]
#   results_AIC_chicken[,i]=sapply(AIC_chicken_windows[[model_name]], function(x){
#     x$aic
#   }
#     )
# }
# save(results_AIC_chicken,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/AIC_only_no_norm_chicken.rda")
# #extract the BIC
# results_BIC_chicken=matrix(NA,nrow = 209452, ncol = length(models))
# colnames(results_BIC_chicken)=models
# for (i in seq_along(models)){
#   model_name=models[i]
#   results_BIC_chicken[,i]=sapply(AIC_chicken_windows[[model_name]], function(x){
#     x$bic
#   }
#   )
# }
# save(results_BIC_chicken,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/BIC_non_norm_chicken.rda")

# # do the AIC for all non-normalized data ####
# load("chicken_RBC/meth_countmatrix_chicken.rda")
# chicken_meth_countmatrix=meth_countmatrix
# chicken_meth_countmatrix=chicken_meth_countmatrix[rowSums(chicken_meth_countmatrix[,-c(1:3)])>0,]
# load("pig_sperm/meth_countmatrix_pig.rda")
# pig_meth_countmatrix=meth_countmatrix
# pig_meth_countmatrix=pig_meth_countmatrix[rowSums(pig_meth_countmatrix[,-c(1:3)])>0,]
# load("wolf/meth_countmatrix_wolf.rda")
# wolf_meth_countmatrix=meth_countmatrix
# wolf_meth_countmatrix=wolf_meth_countmatrix[rowSums(wolf_meth_countmatrix[,-c(1:3)])>0,]
# 
# data_norm=c("wolf","pig","chicken")
# # AIC_windows_norm=list()
# models=c("norm", "pois", "nbinom")
# dis_windows_nonnorm=list()
# for(i in models){
#    for (a in data_norm){
#      list_name=paste0(a,"_",i)
#      dis_windows_nonnorm[[list_name]]=list() #in here i need to put the name with a function for pig and wolf
#      for (b in 1:nrow(get(paste0(a,"_meth_countmatrix")))){
#        print(b)
#        print(i)
#        print(a)
#        f=as.numeric(get(paste0(a,"_meth_countmatrix"))[,-c(1:3)][b,])
#        res=tryCatch({
#          fitdist(f,i,discrete = T)
#        },error=function(e){NA})
#        dis_windows_nonnorm[[list_name]][[b]]=res
#      }
#    }
#  }
# save(dis_windows_nonnorm,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/AIC_non_normalized_chicken_pig_wolf.rda")
load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/AIC_non_normalized_chicken_pig_wolf.rda")
data_norm=c("wolf","pig","chicken")
models=c("norm", "pois", "nbinom")
#extract the aic of chicken, pig and wolf
for (a in data_norm) {
  for (i in seq_along(models)){
    print(a)
    print(i)
    model_name=models[i]
    if (length(dis_windows_nonnorm[[paste0(a,"_",model_name)]])==0) {
      results_aic=NA
    }else{
      results_aic=vapply(dis_windows_nonnorm[[paste0(a,"_",model_name)]], function(x){
        if (is.list(x) && !is.null(x$aic) && !is.na(x$aic)){
          x$aic}else{
            NA
          }
      },FUN.VALUE = numeric(1))
    }
    assign(paste0(a,"_",model_name,"_aic"),results_aic)
  }

}

# do the visualization, right now everything is on vectors ####
#they have to be in a matrix
chicken_aIc=matrix(data = c(chicken_nbinom_aic,chicken_norm_aic,chicken_pois_aic),nrow = length(chicken_nbinom_aic),ncol = 3)
colnames(chicken_aIc)=c("negative_binomial","normal","poisson")
wolf_aic=matrix(data = c(wolf_nbinom_aic,wolf_norm_aic,wolf_pois_aic),nrow = length(wolf_nbinom_aic),ncol = 3)
colnames(wolf_aic)=c("negative_binomial","normal","poisson")
pig_aic=matrix(data = c(pig_nbinom_aic,pig_norm_aic,pig_pois_aic),nrow = length(pig_nbinom_aic),ncol = 3)
colnames(pig_aic)=c("negative_binomial","normal","poisson")

#check by row which is the model with the lowest AIC
lowest.chicken=apply(chicken_aIc,1,function(row) names(row)[which.min(row)])
table(lowest.chicken)
length(lowest.chicken)

lowest.wolf=apply(wolf_aic,1,function(row) names(row)[which.min(row)])
table(lowest.wolf)
length(lowest.wolf)

lowest.pig=apply(pig_aic,1,function(row) names(row)[which.min(row)])
table(lowest.pig)
length(lowest.pig)

# AIC visualization ####
# tiff("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/all_AIC_non_normalized_heatmaps.tiff",
#      width = 500,height = 2000,res = 200,units = "px")
all_vec=c(chicken_nbinom_aic,chicken_norm_aic,chicken_pois_aic)
max(all_vec,na.rm = T)
median(all_vec,na.rm = T)
color_AIC=colorRamp2(breaks = c(min(all_vec,na.rm = T),median(all_vec,na.rm = T)),
                     colors = c("#990000","#fff3f3"))
chicken_heat=Heatmap(chicken_aIc,col = color_AIC,
        row_title = "Normalized count of methylated windows",
        column_title = "AIC of different distributions",
        na_col = "black",
        border=T,
        column_gap = unit(0.05,"mm"),
        column_split = c("Negative Binomial","Normal","Poisson"),
        column_title_side = "bottom",
        cluster_rows=F,
        cluster_columns = F,
        show_row_names=F,
        show_column_dend = F,
        width = unit(10,"mm")*5,
        column_names_rot = 35,
        name = "AIC")
#wolf heatmap
all_vec=c(wolf_nbinom_aic,wolf_norm_aic,wolf_pois_aic)
max(all_vec,na.rm = T)
median(all_vec,na.rm = T)
color_AIC=colorRamp2(breaks = c(min(all_vec,na.rm = T),median(all_vec,na.rm = T)),
                     colors = c("#990000","#fff3f3"))
wolf_heat=Heatmap(wolf_aic,col = color_AIC,
        row_title = "Normalized count of methylated windows",
        column_title = "AIC of different distributions",
        na_col = "black",border_gp = gpar(col="black"),
        column_gap = unit(0.05,"mm"),
        column_split = c("Negative Binomial","Normal","Poisson"),
        column_title_side = "bottom",
        cluster_rows=F,
        show_row_names=F,
        show_column_dend = F,
        width = unit(10,"mm")*5,
        column_names_rot = 35,
        name = "AIC")

#pig heatmap
all_vec=c(pig_nbinom_aic,pig_norm_aic,pig_pois_aic)
max(all_vec,na.rm = T)
median(all_vec,na.rm = T)
color_AIC=colorRamp2(breaks = c(min(all_vec,na.rm = T),median(all_vec,na.rm = T)),
                     colors = c("#990000","#fff3f3"))
pig_heat=Heatmap(pig_aic,col = color_AIC,
        row_title = "Normalized count of methylated windows",
        column_title = "AIC of different distributions",
        na_col = "black",border_gp = gpar(col="black"),
        column_gap = unit(0.05,"mm"),
        column_split = c("Negative Binomial","Normal","Poisson"),
        column_title_side = "bottom",
        cluster_rows=F,
        show_row_names=F,
        show_column_dend = F,
        width = unit(10,"mm")*5,
        column_names_rot = 35,
        name = "AIC")

#heat_list=pig_heat+wolf_heat+chicken_heat
tiff("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_non_normalized_AIC_heatmaps.tiff",
     width = 800,height = 2000,res = 200,units = "px")
pig_heat
dev.off()

tiff("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_non_normalized_AIC_heatmaps.tiff",
     width = 800,height = 2000,res = 200,units = "px")
wolf_heat
dev.off()

tiff("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_non_normalized_AIC_heatmaps.tiff",
     width = 800,height = 2000,res = 200,units = "px")
chicken_heat
dev.off()

# NORMALIZED DATA - Chicken ####
# Now do the same but with the normalized counts as the distribution can change
# but we will need to test every major distribution as normaloized it is continous
# load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/chicken_normalized_counts.rda")
# AIC_chicken_windows_norm=list()
 models=c("norm", "pois", "nbinom", "unif","logis")
# for(i in models){
#   AIC_chicken_windows_norm[[i]]=list()
#   for (b in 1:nrow(chicken.MW)){
#     print(b)
#     print(i)
#     f=as.numeric(chicken.MW[,-c(13)][b,])
#     res=tryCatch({
#       fitdist(f,i,discrete = F)
#     },error=function(e){NULL})
#     AIC_chicken_windows_norm[[i]][[b]]=res
# 
#   }
# }
# save(AIC_chicken_windows_norm,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/AIC_normalized_chicken.rda")
load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/AIC_normalized_chicken.rda")
#extract the AIC
# for (i in seq_along(models)){
#   print(i)
#   model_name=models[i]
#   if (length(AIC_chicken_windows_norm[[model_name]])==0) {
#     results_AIC=NA
#   }else{
#     results_AIC=vapply(AIC_chicken_windows_norm[[model_name]], function(x){
#       if (!is.null(x$aic)){
#         x$aic}else{
#           NA
#         }
#     },FUN.VALUE = numeric(1))
#     
#   }
# assign(paste0("chicken_norm_",model_name,"_AIC"),results_AIC)}
# 
# save(chicken_norm_logis_AIC,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_logis_AIC.rda")
# save(chicken_norm_nbinom_AIC,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_nbinom_AIC.rda")
# save(chicken_norm_norm_AIC,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_norm_AIC.rda")
# save(chicken_norm_pois_AIC,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_pois_AIC.rda")
# save(chicken_norm_unif_AIC,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_unif_AIC.rda")
# 
load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_logis_AIC.rda")
load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_nbinom_AIC.rda")
load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_norm_AIC.rda")
load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_pois_AIC.rda")
load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_unif_AIC.rda")


#extract the BIC
# for (i in seq_along(models)){
#   print(i)
#   model_name=models[i]
#   if (length(AIC_chicken_windows_norm[[model_name]])==0) {
#     results_bic=NA
#   }else{
#     results_bic=vapply(AIC_chicken_windows_norm[[model_name]], function(x){
#       if (!is.null(x$bic)){
#         x$bic}else{
#           NA
#         }
#     },FUN.VALUE = numeric(1))
#     
#   }
# assign(paste0("chicken_norm_",model_name,"_bic"),results_bic)}
# 
# save(chicken_norm_logis_bic,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_logis_bic.rda")
# save(chicken_norm_nbinom_bic,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_nbinom_bic.rda")
# save(chicken_norm_norm_bic,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_norm_bic.rda")
# save(chicken_norm_pois_bic,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_pois_bic.rda")
# save(chicken_norm_unif_bic,file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_unif_bic.rda")

load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_logis_bic.rda")
load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_nbinom_bic.rda")
load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_norm_bic.rda")
load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_pois_bic.rda")
load(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_norm_unif_bic.rda")


# NORMALIZED DATA - wolf and pigs ####
# Now do the same but with the normalized counts as the distribution can change
# but we will need to test every major distribution as normaloized it is continous
# load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/pig_normalized_counts.rda")
# load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/wolf_normalized_counts.rda")
data_norm=c("wolf","pig")
# AIC_windows_norm=list()
models=c("norm", "pois", "exp", "gamma", "nbinom", "geom", "unif","logis")
# for(i in models){
#   for (a in data_norm){
#     list_name=paste0(a,"_",i)
#     AIC_windows_norm[[list_name]]=list() #in here i need to put the name with a function for pig and wolf
#     for (b in 1:nrow(get(paste0(a,".MW")))){
#       print(b)
#       print(i)
#       print(a)
#       f=as.numeric(get(paste0(a,".MW"))[,-ncol(get(paste0(a,".MW")))][b,])
#       res=tryCatch({
#         fitdist(f,i,discrete = F)
#       },error=function(e){NA})
#       AIC_windows_norm[[list_name]][[b]]=res
#     }
#   }
# }
# save(AIC_windows_norm,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/AIC_normalized_pig_wolf.rda")
# load("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/AIC_normalized_pig_wolf.rda")

#extract the AIC of both pig and wolf
# for (a in data_norm) {
#   for (i in seq_along(models)){
#     print(a)
#     print(i)
#     model_name=models[i]
#     if (length(AIC_windows_norm[[paste0(a,"_",model_name)]])==0) {
#       results_AIC=NA
#     }else{
#       results_AIC=vapply(AIC_windows_norm[[paste0(a,"_",model_name)]], function(x){
#         if (!is.null(x$aic)){
#           x$aic}else{
#             NA
#           }
#       },FUN.VALUE = numeric(1))
#     }
#     assign(paste0(a,"_",model_name,"_AIC"),results_AIC)
#   }
#   
# }
# 
# all_vec=list(wolf_logis_AIC=wolf_logis_AIC,wolf_nbinom_AIC=wolf_nbinom_AIC,wolf_norm_AIC=wolf_norm_AIC,wolf_pois_AIC=wolf_pois_AIC,wolf_unif_AIC=wolf_unif_AIC)
# save(wolf_logis_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_logis_AIC.rda")
# save(wolf_nbinom_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_nbinom_AIC.rda")
# save(wolf_norm_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_norm_AIC.rda")
# save(wolf_pois_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_pois_AIC.rda")
# save(wolf_unif_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_unif_AIC.rda")
# save(pig_logis_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_logis_AIC.rda")
# save(pig_nbinom_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_nbinom_AIC.rda")
# save(pig_norm_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_norm_AIC.rda")
# save(pig_pois_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_pois_AIC.rda")
# save(pig_unif_AIC,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_unif_AIC.rda")

load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_logis_AIC.rda")
load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_nbinom_AIC.rda")
load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_norm_AIC.rda")
load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_pois_AIC.rda")
load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_unif_AIC.rda")
load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_logis_AIC.rda")
load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_nbinom_AIC.rda")
load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_norm_AIC.rda")
load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_pois_AIC.rda")
load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_unif_AIC.rda")

#extract the bic of both pig and wolf
# for (a in data_norm) {
#   for (i in seq_along(models)){
#     print(a)
#     print(i)
#     model_name=models[i]
#     if (length(AIC_windows_norm[[paste0(a,"_",model_name)]])==0) {
#       results_bic=NA
#     }else{
#       results_bic=vapply(AIC_windows_norm[[paste0(a,"_",model_name)]], function(x){
#         if (!is.null(x$bic)){
#           x$bic}else{
#             NA
#           }
#       },FUN.VALUE = numeric(1))
#     }
#     assign(paste0(a,"_",model_name,"_bic"),results_bic)
#   }
#   
# }
# 
# all_vec=list(wolf_logis_bic=wolf_logis_bic,wolf_nbinom_bic=wolf_nbinom_bic,wolf_norm_bic=wolf_norm_bic,wolf_pois_bic=wolf_pois_bic,wolf_unif_bic=wolf_unif_bic)
# save(wolf_logis_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_logis_bic.rda")
# save(wolf_nbinom_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_nbinom_bic.rda")
# save(wolf_norm_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_norm_bic.rda")
# save(wolf_pois_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_pois_bic.rda")
# save(wolf_unif_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_unif_bic.rda")
# save(pig_logis_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_logis_bic.rda")
# save(pig_nbinom_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_nbinom_bic.rda")
# save(pig_norm_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_norm_bic.rda")
# save(pig_pois_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_pois_bic.rda")
# save(pig_unif_bic,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_unif_bic.rda")
# 
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_logis_bic.rda")
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_nbinom_bic.rda")
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_norm_bic.rda")
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_pois_bic.rda")
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_unif_bic.rda")
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_logis_bic.rda")
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_nbinom_bic.rda")
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_norm_bic.rda")
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_pois_bic.rda")
# load(file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_unif_bic.rda")

# 
# #now let's do it for each window, it gives error if the window is all 0 ####
# load("wolf/meth_countmatrix_wolf.rda")
# meth_countmatrix=meth_countmatrix[rowSums(meth_countmatrix[,-c(1:3)])>=1,]
# for(i in 1:100){
#   print(i)
#   dir_out="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/"
#   file_names=file.path(dir_out, paste0("wolf_window_distribution",i,".tiff"))
#   tiff(file_names,width = 1500, height = 1600, units = "px",res = 300)
#   f=as.numeric(meth_countmatrix[,-c(1:3)][i,])
#   descdist(f,discrete = T)
#   dev.off()
# }
# 
# load("pig_sperm/meth_countmatrix_pig.rda")
# meth_countmatrix=meth_countmatrix[rowSums(meth_countmatrix[,-c(1:3)])>=1,]
# for(i in 1:100){
#   print(i)
#   dir_out="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/"
#   file_names=file.path(dir_out, paste0("pig_window_distribution",i,".tiff"))
#   tiff(file_names,width = 1500, height = 1600, units = "px",res = 300)
#   f=as.numeric(meth_countmatrix[,-c(1:3)][i,])
#   descdist(f,discrete = T)
#   dev.off()
# }
# 
# 
# 
# do the visualization, right now everything is on vectors ####
#they have to be in a matrix
chicken_aIc_norm=matrix(data = c(chicken_norm_logis_AIC[1:80360],chicken_norm_nbinom_AIC,chicken_norm_norm_AIC[1:80360],chicken_norm_pois_AIC,chicken_norm_unif_AIC[1:80360]),nrow = 80360,ncol = 5)
colnames(chicken_aIc_norm)=c("logis","negative_binomial","normal","poisson","uniform")
wolf_norm_aic_norm=matrix(data = c(wolf_logis_AIC[1:22188],wolf_nbinom_AIC,wolf_norm_AIC[1:22188],
                                   wolf_pois_AIC,wolf_unif_AIC[1:22188]),nrow = 22188,ncol = 5)
colnames(wolf_norm_aic_norm)=c("logis","negative_binomial","normal","poisson","uniform")
pig_aic_norm=matrix(data = c(pig_logis_AIC,pig_nbinom_AIC,pig_norm_AIC,pig_pois_AIC,pig_unif_AIC),nrow = 109439,ncol = 5)
colnames(pig_aic_norm)=c("logis","negative_binomial","normal","poisson","uniform")


#check by row which is the model with the lowest AIC
lowest.chicken=apply(chicken_aIc_norm,1,function(row) names(row)[which.min(row)])
table(lowest.chicken)
length(lowest.chicken)

lowest.wolf=apply(wolf_norm_aic_norm,1,function(row) names(row)[which.min(row)])
table(lowest.wolf)
length(lowest.wolf)

lowest.pig=apply(pig_aic_norm,1,function(row) names(row)[which.min(row)])
table(lowest.pig)
length(lowest.pig)

# AIC visualization ####
par(mfrow=c(1,3))
tiff("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/all_AIC_heatmaps.tiff",
     width = 500,height = 2000,res = 200,units = "px")
all_vec=c(chicken_norm_logis_AIC,chicken_norm_nbinom_AIC,chicken_norm_norm_AIC,chicken_norm_pois_AIC,chicken_norm_unif_AIC)
max(all_vec,na.rm = T)
median(all_vec,na.rm = T)
color_AIC=colorRamp2(breaks = c(min(all_vec,na.rm = T),median(all_vec,na.rm = T)),
                     colors = c("#990000","#fff3f3"))
chicken_heat=Heatmap(chicken_aIc_norm,col = color_AIC,
        row_title = "Normalized count of methylated windows",
        column_title = "AIC of different distributions",
        na_col = "black",
        border=T,
        column_gap = unit(0.05,"mm"),
        column_split = c("Poisson","Normal","Uniform","Logarithmic","Negative Binomial"),
        column_title_side = "bottom",
        cluster_rows=F,
        cluster_columns = F,
        show_row_names=F,
        show_column_dend = F,
        width = unit(10,"mm")*5,
        column_names_rot = 35,
        name = "AIC")
#wolf heatmap
all_vec=c(wolf_logis_AIC,wolf_nbinom_AIC,wolf_norm_AIC,wolf_pois_AIC,wolf_unif_AIC)
max(all_vec,na.rm = T)
median(all_vec,na.rm = T)
color_AIC=colorRamp2(breaks = c(min(all_vec,na.rm = T),median(all_vec,na.rm = T)),
                     colors = c("#990000","#fff3f3"))
wolf_heat=Heatmap(wolf_norm_aic_norm,col = color_AIC,
        row_title = "Normalized count of methylated windows",
        column_title = "AIC of different distributions",
        na_col = "black",border_gp = gpar(col="black"),
        column_gap = unit(0.05,"mm"),
        column_split = c("Poisson","Normal","Uniform","Logarithmic","Negative Binomial"),
        column_title_side = "bottom",
        cluster_rows=F,
        cluster_columns = F,
        show_row_names=F,
        show_column_dend = F,
        width = unit(10,"mm")*5,
        column_names_rot = 35,
        name = "AIC")

#pig heatmap
all_vec=c(pig_logis_AIC,pig_nbinom_AIC,pig_norm_AIC,pig_pois_AIC,pig_unif_AIC)
max(all_vec,na.rm = T)
median(all_vec,na.rm = T)
color_AIC=colorRamp2(breaks = c(min(all_vec,na.rm = T),median(all_vec,na.rm = T)),
                     colors = c("#990000","#fff3f3"))
pig_heat=Heatmap(pig_aic_norm,col = color_AIC,
        row_title = "Normalized count of methylated windows",
        column_title = "AIC of different distributions",
        na_col = "black",border_gp = gpar(col="black"),
        column_gap = unit(0.05,"mm"),
        column_split = c("Poisson","Normal","Uniform","Logarithmic","Negative Binomial"),
        column_title_side = "bottom",
        cluster_rows=F,
        cluster_columns = F,
        show_row_names=F,
        show_column_dend = F,
        width = unit(10,"mm")*5,
        column_names_rot = 35,
        name = "AIC")


tiff("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/pig_AIC_heatmaps.tiff",
     width = 800,height = 2000,res = 200,units = "px")
pig_heat
dev.off()

tiff("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/wolf_AIC_heatmaps.tiff",
     width = 800,height = 2000,res = 200,units = "px")
wolf_heat
dev.off()

tiff("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/distribution_plots/chicken_AIC_heatmaps.tiff",
     width = 800,height = 2000,res = 200,units = "px")
chicken_heat
dev.off()







