setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged")
library(data.table)
library(ggplot2)
library(cleaner)
library(gridExtra)
library(cowplot)
library(tidyr)
library(dplyr)

#read medips
wolf_medips=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged/wolf/meth_countmatrix_wolf.txt")
wolf_medips$location=paste0(wolf_medips$chr,":",wolf_medips$start,"-",wolf_medips$stop)
colnames(wolf_medips)=sub(".*/","",colnames(wolf_medips))
colnames(wolf_medips)=sub("_MeDIP\\.MQ10\\.bam\\.counts$","",colnames(wolf_medips))
colnames(wolf_medips)=sub("X","",colnames(wolf_medips))
#wolf_medips$test="MEDIPS"
wolf_medips=pivot_longer(wolf_medips,cols = 4:9, names_to = "individual",values_to = "medips_counts")

pig_medips=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged/pig_sperm/meth_countmatrix_pig.txt")
pig_medips$location=paste0(pig_medips$chr,":",pig_medips$start,"-",pig_medips$stop)
colnames(pig_medips)=sub(".*/","",colnames(pig_medips))
colnames(pig_medips)=sub("\\.MQ10\\.bam\\.counts$","",colnames(pig_medips))
colnames(pig_medips)=sub("X","",colnames(pig_medips))
#pig_medips$test="MEDIPS"
pig_medips=pivot_longer(pig_medips,cols = 4:29, names_to = "individual",values_to = "medips_counts")

chicken_medips=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/merged/chicken_RBC/meth_countmatrix_chicken1.txt")
chicken_medips$location=paste0(chicken_medips$chr,":",chicken_medips$start,"-",chicken_medips$stop)
colnames(chicken_medips)=sub(".*/","",colnames(chicken_medips))
colnames(chicken_medips)=sub("\\.MQ10\\.bam\\.counts$","",colnames(chicken_medips))
colnames(chicken_medips)=sub("X","",colnames(chicken_medips))
#chicken_medips$test="MEDIPS"
chicken_medips=pivot_longer(chicken_medips,cols = 4:15, names_to = "individual",values_to = "medips_counts")

#read featurecounts
wolf_featurecounts=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/count.matrix.wolf.bed")
wolf_featurecounts$location=paste0(wolf_featurecounts$Chr,":",wolf_featurecounts$Start,"-",wolf_featurecounts$End)
colnames(wolf_featurecounts)=sub(".*/","",colnames(wolf_featurecounts))
colnames(wolf_featurecounts)=sub("_MeDIP\\.MQ10\\.bam$","",colnames(wolf_featurecounts))
#wolf_featurecounts$test="featureCounts"
wolf_names=colnames(wolf_featurecounts)
wolf_names=wolf_names[7:12]
wolf_featurecounts=pivot_longer(wolf_featurecounts,cols = 7:12, names_to = "individual",values_to = "featurecounts_counts")

pig_featurecounts=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/count.matrix.pig.bed")
pig_featurecounts$location=paste0(pig_featurecounts$Chr,":",pig_featurecounts$Start,"-",pig_featurecounts$End)
colnames(pig_featurecounts)=sub(".*/","",colnames(pig_featurecounts))
colnames(pig_featurecounts)=sub("\\.MQ10\\.bam$","",colnames(pig_featurecounts))
#pig_featurecounts$test="featureCounts"
pig_names=colnames(pig_featurecounts)
pig_names=pig_names[7:32]
pig_featurecounts=pivot_longer(pig_featurecounts,cols = 7:32, names_to = "individual",values_to = "featurecounts_counts")

chicken_featurecounts=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/count.matrix.chicken")
chicken_featurecounts$location=paste0(chicken_featurecounts$Chr,":",chicken_featurecounts$Start,"-",chicken_featurecounts$End)
colnames(chicken_featurecounts)=sub(".*/","",colnames(chicken_featurecounts))
colnames(chicken_featurecounts)=sub("\\.MQ10\\.bam$","",colnames(chicken_featurecounts))
#chicken_featurecounts$test="featureCounts"
chicken_names=colnames(chicken_featurecounts)
chicken_names=chicken_names[7:18]
chicken_featurecounts=pivot_longer(chicken_featurecounts,cols = 7:18, names_to = "individual",values_to = "featurecounts_counts")

#bedtools files #####
wolf_bedtools=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/wolf_bedtools_counts.bed")
wolf_bedtools$location=paste0(wolf_bedtools$V1,":",wolf_bedtools$V2,"-",wolf_bedtools$V3)
colnames(wolf_bedtools)[4:9]=wolf_names
wolf_bedtools=pivot_longer(wolf_bedtools,cols = 4:9, names_to = "individual",values_to = "bedtools_counts")

pig_bedtools=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/pig_sperm_bedtools_counts.bed")
pig_bedtools$location=paste0(pig_bedtools$V1,":",pig_bedtools$V2,"-",pig_bedtools$V3)
colnames(pig_bedtools)[4:29]=pig_names
pig_bedtools=pivot_longer(pig_bedtools,cols = 4:29, names_to = "individual",values_to = "bedtools_counts")

chicken_bedtools=fread("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/chichen_RBC_bedtools_counts.bed")
chicken_bedtools$location=paste0(chicken_bedtools$V1,":",chicken_bedtools$V2,"-",chicken_bedtools$V3)
colnames(chicken_bedtools)[4:15]=chicken_names
chicken_bedtools=pivot_longer(chicken_bedtools,cols = 4:15, names_to = "individual",values_to = "bedtools_counts")

#merge the two datasets, medips and featurecounts
wolf_total=inner_join(wolf_featurecounts,wolf_medips,by=c("location","individual"))
wolf_total$specie="Wolf"
wolf_total=wolf_total[,c(7,8,9,13,14)]
#wolf_total=pivot_longer(wolf_total,cols = c(7:12,18:23), names_to = "individuals",values_to = "counts")

pig_total=inner_join(pig_featurecounts,pig_medips,by=c("location","individual"))
pig_total$specie="Pig"
pig_total=pig_total[,c(7,8,9,13,14)]
#pig_total=pivot_longer(pig_total,cols = c(7:32,37:62), names_to = "individuals",values_to = "counts")

chicken_total=inner_join(chicken_featurecounts,chicken_medips,by=c("location","individual"))
chicken_total$specie="Chicken"
chicken_total=chicken_total[,c(7,8,9,13,14)]
#chicken_total=pivot_longer(chicken_total,cols = ,names_to = "individuals",values_to = "counts")
total=bind_rows(chicken_total,pig_total,wolf_total)

#merge the two datasets, medips and bedtools
wolf_bed_medips=inner_join(wolf_bedtools,wolf_medips,by=c("location","individual"))
wolf_bed_medips$specie="Wolf"
wolf_bed_medips=wolf_bed_medips[,c(7,8,9,10,6,11)]
#wolf_bed_medips=pivot_longer(wolf_bed_medips,cols = c(7:12,18:23), names_to = "individuals",values_to = "counts")

pig_bed_medips=inner_join(pig_bedtools,pig_medips,by=c("location","individual"))
pig_bed_medips$specie="Pig"
pig_bed_medips=pig_bed_medips[,c(7,8,9,10,6,11)]
#pig_bed_medips=pivot_longer(pig_bed_medips,cols = c(7:32,37:62), names_to = "individuals",values_to = "counts")

chicken_bed_medips=inner_join(chicken_bedtools,chicken_medips,by=c("location","individual"))
chicken_bed_medips$specie="Chicken"
chicken_bed_medips=chicken_bed_medips[,c(7,8,9,10,6,11)]
#chicken_bed_medips=pivot_longer(chicken_bed_medips,cols = ,names_to = "individuals",values_to = "counts")
bed_medips=bind_rows(chicken_bed_medips,pig_bed_medips,wolf_bed_medips)

#merge the two datasets, bedtools and featurecounts #####
wolf_bed_feature=inner_join(wolf_featurecounts,wolf_bedtools,by=c("location","individual"))
wolf_bed_feature$specie="Wolf"
wolf_bed_feature=wolf_bed_feature[,c(10:13,9,14)]
#wolf_bed_feature=pivot_longer(wolf_bed_feature,cols = c(7:12,18:23), names_to = "individuals",values_to = "counts")

pig_bed_feature=inner_join(pig_featurecounts,pig_bedtools,by=c("location","individual"))
pig_bed_feature$specie="Pig"
pig_bed_feature=pig_bed_feature[,c(10:13,9,14)]
#pig_bed_feature=pivot_longer(pig_bed_feature,cols = c(7:32,37:62), names_to = "individuals",values_to = "counts")

chicken_bed_feature=inner_join(chicken_featurecounts,chicken_bedtools,by=c("location","individual"))
chicken_bed_feature$specie="Chicken"
chicken_bed_feature=chicken_bed_feature[,c(10:13,9,14)]
#chicken_bed_feature=pivot_longer(chicken_bed_feature,cols = ,names_to = "individuals",values_to = "counts")
bed_feature=bind_rows(chicken_bed_feature,pig_bed_feature,wolf_bed_feature)

# MEDIPS vs featureCounts ####
#do chicken ####
# tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/plotchicken_medips_vs_featurecounts.tiff",width = 1500,height = 1500,units = "px",res = 100)
# ggplot(chicken_total, aes(x=featurecounts_counts,y=medips_counts))+
#   geom_point(alpha=1,size = 7,color="#E69F00")+labs(x="featureCounts",y="MEDIPS")+theme_minimal()+
#   #scale_fill_manual(values = c("Chicken"="#E69F00"))+
#   theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50),
#         legend.text = element_text(size = 50),
#         legend.title = element_text(size = 50))
# dev.off()
# #do wolf #####
# tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/plotwolf_medips_vs_featurecounts.tiff",width = 1500,height = 1500,units = "px",res = 100)
# ggplot(wolf_total, aes(x=featurecounts_counts,y=medips_counts))+
#   geom_point(size = 7,color="#009E73")+labs(x="featureCounts",y="MEDIPS")+theme_minimal()+
#   #scale_color_manual(values = c("#009E73"))+
#   theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50))
# dev.off()
# #do pig ####
# tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/plotpig_medips_vs_featurecounts.tiff",width = 1500,height = 1500,units = "px",res = 100)
# ggplot(pig_total, aes(x=featurecounts_counts,y=medips_counts,color=specie))+
#   geom_point(alpha=1,size = 7,color="#56B4E9")+labs(x="featureCounts",y="MEDIPS")+theme_minimal()+
#   #scale_fill_manual(values = c("Pig"="#56B4E9"))+
#   theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50),
#         legend.text = element_text(size = 50),
#         legend.title = element_text(size = 50))
# dev.off()

#bedtools vs MEDIPS #####
#do chicken ####
tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/plotchicken_bedtools_vs_medips.tiff",width = 1500,height = 1500,units = "px",res = 100)
ggplot(chicken_bed_medips, aes(x=bedtools_counts,y=medips_counts))+
  geom_point(alpha=1,size = 7,color="#E69F00")+labs(x="BEDTools",y="MEDIPS")+theme_minimal()+
  #scale_fill_manual(values = c("Chicken"="#E69F00"))+
  theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.title = element_text(size = 50))
dev.off()
#do wolf #####
tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/plotwolf_bedtools_vs_medips.tiff",width = 1500,height = 1500,units = "px",res = 100)
ggplot(wolf_bed_medips, aes(x=bedtools_counts,y=medips_counts))+
  geom_point(size = 7,color="#009E73")+labs(x="BEDTools",y="MEDIPS")+theme_minimal()+
  #scale_color_manual(values = c("#009E73"))+
  theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50))
dev.off()
#do pig ####
tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/plotpig_bedtools_vs_medips.tiff",width = 1500,height = 1500,units = "px",res = 100)
ggplot(pig_bed_medips, aes(x=bedtools_counts,y=medips_counts,color=specie))+
  geom_point(alpha=1,size = 7,color="#56B4E9")+labs(x="BEDTools",y="MEDIPS")+theme_minimal()+
  #scale_fill_manual(values = c("Pig"="#56B4E9"))+
  theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.title = element_text(size = 50))
dev.off()

#bedtools vs featureCounts #####
#do chicken ####
tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/plotchicken_bedtools_vs_featurecounts.tiff",width = 1500,height = 1500,units = "px",res = 100)
ggplot(chicken_bed_feature, aes(x=featurecounts_counts,y=bedtools_counts))+
  geom_point(alpha=1,size = 7,color="#E69F00")+labs(x="featureCounts",y="BEDTools")+theme_minimal()+
  #scale_fill_manual(values = c("Chicken"="#E69F00"))+
  theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.title = element_text(size = 50))
dev.off()
#do wolf #####
tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/plotwolf_bedtools_vs_featurecounts.tiff",width = 1500,height = 1500,units = "px",res = 100)
ggplot(wolf_bed_feature, aes(x=featurecounts_counts,y=bedtools_counts))+
  geom_point(size = 7,color="#009E73")+labs(x="featureCounts",y="BEDTools")+theme_minimal()+
  #scale_color_manual(values = c("#009E73"))+
  theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50))
dev.off()
#do pig ####
tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/plotpig_bedtools_vs_featurecounts.tiff",width = 1500,height = 1500,units = "px",res = 100)
ggplot(pig_bed_feature, aes(x=featurecounts_counts,y=bedtools_counts,color=specie))+
  geom_point(alpha=1,size = 7,color="#56B4E9")+labs(x="featureCounts",y="BEDTools")+theme_minimal()+
  #scale_fill_manual(values = c("Pig"="#56B4E9"))+
  theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50),
        legend.text = element_text(size = 50),
        legend.title = element_text(size = 50))
dev.off()

# all together #####
# graph_total=ggplot(total, aes(x=featurecounts_counts,y=medips_counts,color=specie))+
#   geom_point(alpha=0.5,size = 20)+labs(x="featureCounts",y="MEDIPS",color="Specie")+theme_minimal()+
#   scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
#   theme(text = element_text(size = 50),axis.title = element_text(size = 50),axis.text = element_text(size = 50),
#         legend.text = element_text(size = 50),
#         legend.title = element_text(size = 50))
# 
# tiff(filename = "plot_medips_vs_featurecounts.tiff",width = 6000,height = 6000,units = "px",res = 100)
# print(graph_total)
# dev.off()

























