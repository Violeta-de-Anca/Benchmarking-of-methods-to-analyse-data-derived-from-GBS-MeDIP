setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/")
library(ggplot2)
library(cleaner)
library(gridExtra)
library(cowplot)
library(dplyr)
load("combination.rda")
# wolf.combin=combination%>%filter(grepl("Wolf",dataset))
# save(wolf.combin,file="wolf.combination.rda")
load("wolf.combination.rda")
# chicken.combin=combination%>%filter(grepl("Chicken",dataset))
# save(chicken.combin,file="chicken.combination.rda")
load("chicken.combination.rda")
# pig.combin=combination%>%filter(grepl("Pig",dataset))
# save(pig.combin,file="pig.combination.rda")
load("pig.combination.rda")

pp_plot_w=ggplot(wolf.combin,aes(x=sorted_p_values,y=expected_p_values,color=method))+
  geom_point()+geom_abline(slope = 1,intercept = 0,linetype="dashed")+labs(x="P-values",y="Theoretical quantiles")+
  theme_minimal()+ scale_color_manual(values = c("Mann-Whitney"="#1f78b4","Maximum-likelihood"="#dca60d","Quasi-likelihood"="#31a354","Moderated T-test"="#e6550d","DESeq2"="#756bb1"))+
  theme(axis.title = element_text(size = 40),axis.text = element_text(size = 40),  legend.title = element_text(size = 40),   legend.text = element_text(size = 40),plot.title = element_text(size = 40),strip.text = element_text(size = 40))
 tiff(filename = "p-p_plot_wolf.tiff",width = 2500,height = 1500,units = "px",res = 100)
 print(pp_plot_w)
 dev.off()

 pp_plot_p=ggplot(pig.combin,aes(x=sorted_p_values,y=expected_p_values,color=method))+
   geom_point()+geom_abline(slope = 1,intercept = 0,linetype="dashed")+labs(x="P-values",y="Theoretical quantiles")+
   theme_minimal()+ scale_color_manual(values = c("Mann-Whitney"="#1f78b4","Maximum-likelihood"="#dca60d","Quasi-likelihood"="#31a354","Moderated T-test"="#e6550d","DESeq2"="#756bb1"),guide = guide_legend(override.aes = list(size = 20)),labels = c("Mann-Whitney" ="Mann-Whitney","Maximum-likelihood"="EdgeR - ML","Quasi-likelihood"="EdgeR - QL", "Moderated T-test"="Limma - Moderated T-test","DESeq2"="DESeq2"))+
   theme(axis.title = element_text(size = 40),axis.text = element_text(size = 40),  legend.title = element_text(size = 40),   legend.text = element_text(size = 40),plot.title = element_text(size = 40),strip.text = element_text(size = 40))
 tiff(filename = "p-p_plot_pig.tiff",width = 2500,height = 1500,units = "px",res = 100)
 print(pp_plot_p)
 dev.off()

pp_plot_c=ggplot(chicken.combin,aes(x=sorted_p_values,y=expected_p_values,color=method))+
   geom_point()+geom_abline(slope = 1,intercept = 0,linetype="dashed")+labs(x="P-values",y="Theoretical quantiles")+
   theme_minimal()+ scale_color_manual(values = c("Mann-Whitney"="#1f78b4","Maximum-likelihood"="#dca60d","Quasi-likelihood"="#31a354","Moderated T-test"="#e6550d","DESeq2"="#756bb1"))+
  theme(axis.title = element_text(size = 40),axis.text = element_text(size = 40),
 legend.title = element_text(size = 40),   legend.text = element_text(size = 40),
 plot.title = element_text(size = 40),strip.text = element_text(size = 40))
 tiff(filename = "p-p_plot_chicken.tiff",width = 2500,height = 1500,units = "px",res = 100)
 print(pp_plot_c)
 dev.off()

#pp_plot=ggplot(combination,aes(x=sorted_p_values,y=expected_p_values,color=method,shape=dataset))+
#  geom_point()+geom_abline(slope = 1,intercept = 0,linetype="dashed")+labs(x="P-values",y="Theoretical quantiles")+
#  theme_minimal()+ scale_color_manual(values = c("Mann-Whitney"="blue","Maximum-likelihood"="red","Quasi-likelihood"="green","Moderated T-test"="orange"))+
#  scale_shape_manual(values = c(16,17,18))
#tiff(filename = "p-p_plot.tiff",width = 6000,height = 6000,units = "px",res = 100)
#print(pp_plot)
#dev.off()

 #do histograms as another metric #####
tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/histogram_null_distr_dog_wolf_withDESeq2.tiff",
      width = 4500,height = 5500,units = "px",res = 200)
ggplot(wolf.combin,aes(x=sorted_p_values,fill=method))+ facet_grid(method ~.)+
   geom_histogram(bins = 100,color="black")+
   labs(title = "Distribution of p-values\nderived from a random distribution of\nthe dog and wolf dataset",
      x = "P-values",
      y = "Frequency",
      fill= "Statistical method")+
  geom_hline(yintercept = 0.01*nrow(wolf.combin),linetype="dashed",color="black")+
     scale_fill_manual(values = c("Mann-Whitney"="#1f78b4","Maximum-likelihood"="#dca60d",
                                   "Quasi-likelihood"="#31a354","Moderated T-test"="#e6550d","DESeq2"="#756bb1"))+
   theme_minimal()+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white",color = NA),
                         plot.background = element_rect(fill = "white",colour = NA),
                         axis.title = element_text(size = 40),axis.text = element_text(size = 40),
                         legend.title = element_text(size = 40),   legend.text = element_text(size = 40),
                         plot.title = element_text(size = 40),strip.text = element_text(size = 40))
dev.off()
tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/histogram_null_pig_withDESeq2.tiff",
     width = 4500,height = 5500,units = "px",res = 200)
   ggplot(pig.combin,aes(x=sorted_p_values,fill=method))+ facet_grid(method ~.)+
     geom_histogram(position = "identity",bins = 100,color="black")+
     geom_hline(yintercept = 0.01*nrow(pig.combin),linetype="dashed",color="black")+
     labs(title = "Distribution of p-values\nderived from a random distribution of\nthe pig dataset",
          x = "P-values",
          y = "Frequency",
          fill= "Statistical method")+
     scale_fill_manual(values = c("Mann-Whitney"="#1f78b4","Maximum-likelihood"="#dca60d",
                                   "Quasi-likelihood"="#31a354","Moderated T-test"="#e6550d","DESeq2"="#756bb1"))+
   theme_minimal()+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white",color = NA),
                         plot.background = element_rect(fill = "white",colour = NA),
                         axis.title = element_text(size = 40),axis.text = element_text(size = 40),
                         legend.title = element_text(size = 40),   legend.text = element_text(size = 40),
                         plot.title = element_text(size = 40),strip.text = element_text(size = 40))
dev.off()
tiff(filename = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/histogram_null_distr_chicken_withDESeq2.tiff",
        width = 4500,height = 5500,units = "px",res = 200)
   ggplot(chicken.combin,aes(x=sorted_p_values,fill=method))+ 
     geom_hline(yintercept = 0.01*nrow(chicken.combin),linetype="dashed",color="black")+
     geom_histogram(position = "identity",bins = 100,color="black")+facet_grid(method ~.)+
     labs(title = "Distribution of p-values\nderived from a random distribution of\nthe chicken dataset",
          x = "P-values",
          y = "Frequency",
          fill= "Statistical method")+
     scale_fill_manual(values = c("Mann-Whitney"="#1f78b4","Maximum-likelihood"="#dca60d",
                                   "Quasi-likelihood"="#31a354","Moderated T-test"="#e6550d","DESeq2"="#756bb1"))+
   theme_minimal()+theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white",color = NA),
                         plot.background = element_rect(fill = "white",colour = NA),
                         axis.title = element_text(size = 40),axis.text = element_text(size = 40),
                         legend.title = element_text(size = 40),   legend.text = element_text(size = 40),
                         plot.title = element_text(size = 40),strip.text = element_text(size = 40))
   dev.off()
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 