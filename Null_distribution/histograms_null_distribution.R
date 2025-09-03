setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/")
library(ggplot2)
library(cleaner)
library(gridExtra)
library(cowplot)

#Lets do the distribution of null p-values against expectation
#mann-whitney
load(file = "list.1000.p.values.MW.loop.wolf.rda")
mw.wolf=unlist(merged_wolf)

# mw.wolf=sort(mw.wolf)
#max(table(cut(mw.wolf,breaks=30)))
load(file = "list.1000.p.values.MW.loop.pig.rda")
mw.pig=unlist(merged_pig)

load(file = "list.1000.p.values.MW.loop.chicken.rda")
mw.chicken=unlist(merged_chicken)

#ML
load(file = "list.1000.p.values.loop.wolf_ML.rda")
ML.wolf=unlist(list.p.values.loop.1000)

#max(table(cut(ML.wolf,breaks=30)))
load(file = "list.1000.p.values.loop.pig_ML.rda")
ML.pig=unlist(list.p.values.loop.1000)

#max(table(cut(ML.pig,breaks=30)))
load(file = "list.1000.p.values.loop.chicken_ML.rda")
ML.chicken=unlist(list.p.values.loop.1000)

#QL
load(file = "list.1000.p.values.loop.wolf_QL.rda")
QL.wolf=unlist(list.p.values.loop.1000)

#max(table(cut(QL.wolf,breaks=30)))
load(file = "list.1000.p.values.loop.pig_QL.rda")
QL.pig=unlist(list.p.values.loop.1000)

load(file = "list.1000.p.values.loop.chicken_QL.rda")
QL.chicken=unlist(list.p.values.loop.1000)

#ttest
load(file = "list.1000.p.values.loop.wolf_ttest.rda")
ttest.wolf=unlist(list.p.values.loop.1000)

#max(table(cut(ttest.wolf,breaks=30)))
load(file = "list.1000.p.values.loop.pig_ttest.rda")
ttest.pig=unlist(list.p.values.loop.1000)

load(file = "list.1000.p.values.loop.chicken_ttest.rda")
ttest.chicken=unlist(list.p.values.loop.1000)

#dds
load(file = "list.1000.p.values.loop.wolf_dds.rda")
dds.wolf=unlist(list.p.values.loop.1000)

#max(table(cut(dds.wolf,breaks=30)))
load(file = "list.1000.p.values.loop.pig.dds.rda")
dds.pig=unlist(list.p.values.loop.1000)

load(file = "list.1000.p.values.loop.chicken_deseq2.rda")
dds.chicken=unlist(list.p.values.loop.1000)


list_p_values=list(mw.wolf,mw.pig,mw.chicken,ML.wolf,ML.pig,ML.chicken,QL.wolf,QL.pig,QL.chicken,ttest.wolf,ttest.pig,ttest.chicken,dds.pig,dds.wolf,dds.chicken)
save(list_p_values,file = "list_p_values_null.rda")
sorted_p_values=list()
expected_p_values=list()
for (i in 1:length(list_p_values)){
  sorted_p_values[[i]]=sort(list_p_values[[i]])
  expected_p_values[[i]]=seq(1/length(sorted_p_values[[i]]),1,length.out=length(sorted_p_values[[i]]))
}
method=c(rep("Mann-Whitney",3),rep("Maximum-likelihood",3),rep("Quasi-likelihood",3),rep("Moderated T-test",3),rep("DESeq2",3))
datasets=c(rep(c("Wolf","Pig","Chicken"),times=4),c("Wolf","Pig","Chicken"))
combination=data.frame()
for (i in 1:length(sorted_p_values)){
  print(i)
  temp_data=data.frame(
    sorted_p_values=sorted_p_values[[i]],
    expected_p_values=expected_p_values[[i]],
    method=method[i],
    dataset=datasets[i]
  )
  combination=rbind(combination,temp_data)
}
save(combination,file="combination.rda")
#load("combination.rda")
#pp_plot=ggplot(combination,aes(x=sorted_p_values,y=expected_p_values,color=method,shape=dataset))+
  geom_point()+geom_abline(slope = 1,intercept = 0,linetype="dashed")+labs(x="P-values",y="Theoretical quantiles")+
  theme_minimal()+ scale_color_manual(values = c("Mann-Whitney"="#1f78b4","Maximum-likelihood"="#dca60d","Quasi-likelihood"="#31a354","Moderated T-test"="#e6550d","DESeq2"="#756bb1"))+
  scale_shape_manual(values = c(16,17,18))

#print(pp_plot)

#get the false positive rate under a null distribution of p-values
percentage.w.mw=mean(mw.wolf<0.05,na.rm = T)*100
percentage.p.mw=mean(mw.pig<0.05,na.rm = T)*100
percentage.c.mw=mean(mw.chicken<0.05,na.rm = T)*100

percentage.w.ml=mean(ML.wolf<0.05)*100
percentage.p.ml=mean(ML.pig<0.05)*100
percentage.c.ml=mean(ML.chicken<0.05)*100

percentage.w.ql=mean(QL.wolf<0.05)*100
percentage.p.ql=mean(QL.pig<0.05)*100
percentage.c.ql=mean(QL.chicken<0.05)*100

percentage.w.ttest=mean(ttest.wolf<0.05)*100
percentage.p.ttest=mean(ttest.pig<0.05)*100
percentage.c.ttest=mean(ttest.chicken<0.05)*100

percentage.w.dds=mean(dds.wolf<0.05,na.rm = T)*100
percentage.p.dds=mean(dds.pig<0.05,na.rm = T)*100
percentage.c.dds=mean(dds.chicken<0.05,na.rm = T)*100
print(percentage.c.dds)
#############################################################################
#distribution of null p-values
#mann-whitney
load(file = "list.1000.p.values.MW.loop.wolf.rda")
merged_wolf=unlist(merged_wolf)
mw.wolf=merged_wolf
#mw.wolf=sort(mw.wolf)
#max(table(cut(mw.wolf,breaks=30)))
load(file = "list.1000.p.values.MW.loop.pig.rda")
merged_pig=unlist(merged_pig)
mw.pig=merged_pig
load(file = "list.1000.p.values.MW.loop.chicken.rda")
merged_chicken=unlist(merged_chicken)
mw.chicken=merged_chicken
vector.mw=data.frame(value=c(mw.chicken,mw.pig,mw.wolf),
                         Group=rep(c("Chicken","Pig","Wolf"),each=length(mw.chicken)))
max_y_value.mw=max(table(cut(vector.mw$value,breaks=30)))
mw=ggplot(vector.mw,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",title="Distribution of null p-values under Mann-Whitney test")+
  #ylim(0,3115500)+
  scale_y_continuous(limits = c(0,max_y_value.mw),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                        axis.title = element_text(size = 50),
                        axis.text = element_text(size = 50),
                        legend.text = element_text(size = 50),
                        legend.title = element_text(size = 50),
                        plot.title = element_text(size = 50),
                        plot.title.position = "plot")
save(mw,file = "mw_all.rda")

#distribution of null p-values
#ML
load(file = "list.1000.p.values.loop.wolf_ML.rda")
merged_wolf=(list.p.values.loop.1000)
ML.wolf=merged_wolf
#max(table(cut(ML.wolf,breaks=30)))
load(file = "list.1000.p.values.loop.pig_ML.rda")
merged_pig=(list.p.values.loop.1000)
ML.pig=merged_pig
max(table(cut(ML.pig,breaks=30)))
load(file = "list.1000.p.values.loop.chicken_ML.rda")
merged_chicken=(list.p.values.loop.1000)
ML.chicken=merged_chicken
vector.ML=data.frame(value=c(ML.chicken,ML.pig,ML.wolf),
                     Group=rep(c("Chicken","Pig","Wolf"),each=length(ML.chicken)))
max_y_value.ML=max(table(cut(vector.ML$value,breaks=30)))
ML=ggplot(vector.ML,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",title = "Distribution of null p-values under Maximum Likelihood")+
  ylim(0,2917203)+
  #scale_y_continuous(limits = c(0,max_y_value.ML),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                        axis.title = element_text(size = 50),
                        axis.text = element_text(size = 50),
                        legend.text = element_text(size = 50),
                        legend.title = element_text(size = 50),
                        plot.title = element_text(size = 50),
                        plot.title.position = "plot")
save(ML,file = "ML_all.rda")

#QL
load(file = "list.1000.p.values.loop.wolf_QL.rda")
merged_wolf=(list.p.values.loop.1000)
QL.wolf=merged_wolf
#max(table(cut(QL.wolf,breaks=30)))
load(file = "list.1000.p.values.loop.pig_QL.rda")
merged_pig=(list.p.values.loop.1000)
QL.pig=merged_pig
load(file = "list.1000.p.values.loop.chicken_QL.rda")
merged_chicken=(list.p.values.loop.1000)
QL.chicken=merged_chicken
max(table(cut(QL.chicken,breaks=30)))ttest.chicken
vector.QL=data.frame(value=c(QL.chicken,QL.pig,QL.wolf),
                     Group=rep(c("Chicken","Pig","Wolf"),each=length(QL.chicken)))
max_y_value.QL=max(table(cut(vector.QL$value,breaks=30)))
QL=ggplot(vector.QL,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",title = "Distribution of null p-values under Quasi Likelihood")+
  ylim(0,5392811)+
  #scale_y_continuous(limits = c(0,max_y_value.QL),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                        axis.title = element_text(size = 50),
                        axis.text = element_text(size = 50),
                        legend.text = element_text(size = 50),
                        legend.title = element_text(size = 50),
                        plot.title = element_text(size = 50),
                        plot.title.position = "plot")
save(QL,file = "QL_all.rda")

#ttest
load(file = "list.1000.p.values.loop.wolf_ttest.rda")
merged_wolf=(list.p.values.loop.1000)
ttest.wolf=merged_wolf
#max(table(cut(ttest.wolf,breaks=30)))
load(file = "list.1000.p.values.loop.pig_ttest.rda")
merged_pig=(list.p.values.loop.1000)
ttest.pig=merged_pig
load(file = "list.1000.p.values.loop.chicken_ttest.rda")
merged_chicken=(list.p.values.loop.1000)
ttest.chicken=merged_chicken
max(table(cut(ttest.chicken,breaks=30)))
vector.ttest=data.frame(value=c(ttest.chicken,ttest.pig,ttest.wolf),
                     Group=rep(c("Chicken","Pig","Wolf"),each=length(ttest.chicken)))
max_y_value.ttest=max(table(cut(vector.ttest$value,breaks=30)))
ttest=ggplot(vector.ttest,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",title = "Distribution of null p-values under Moderated T-test")+
  ylim(0,1425792)+
  #scale_y_continuous(limits = c(0,max_y_value.ttest),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                        axis.title = element_text(size = 50),
                        axis.text = element_text(size = 50),
                        legend.text = element_text(size = 50),
                        legend.title = element_text(size = 50),
                        plot.title = element_text(size = 50),
                        plot.title.position = "plot")
save(ttest,file = "ttest_all.rda")

legend.all=get_legend(ttest+theme(legend.position = "bottom"))
ttest=ttest+theme(legend.position = "none")
QL=QL+theme(legend.position = "none")
ML=ML+theme(legend.position = "none")
mw=mw+theme(legend.position = "none")
combined_all=plot_grid(mw,ML,QL,ttest,ncol = 2,align = "v")
final_all=plot_grid(combined_all,legend.all,ncol = 1,rel_heights = c(1,0.1))
save(final_all,file = "plot_null_all.rda")

tiff(filename = "plot_all_nocorrection_null.tiff",width = 6000,height = 6000,units = "px",res = 100)
print(final_all)
dev.off()

################################################################################
#load all mann whitney FDR
load(file = "list.1000.p.values.MW.fdr.loop.wolf.rda")
merged_wolf_fdr=unlist(merged_wolf_fdr)
mw.fdr.wolf=merged_wolf_fdr
max(table(cut(mw.fdr.wolf,breaks=30)))
load(file = "list.1000.p.values.MW.fdr.loop.pig.rda")
merged_pig_fdr=unlist(merged_pig_fdr)
mw.fdr.pig=merged_pig_fdr
load(file = "list.1000.p.values.MW.fdr.loop.chicken.rda")
merged_chicken_fdr=unlist(merged_chicken_fdr)
mw.fdr.chicken=merged_chicken_fdr
vector.mw.fdr=data.frame(value=c(mw.fdr.chicken,mw.fdr.pig,mw.fdr.wolf),
                            Group=rep(c("Chicken","Pig","Wolf"),each=length(mw.fdr.chicken)))
max_y_value.mw.fdr=max(table(cut(vector.mw.fdr$value,breaks=30)))
mw.fdr=ggplot(vector.mw.fdr,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",main="Mann-Whitney with FDR correction")+
  ylim(0,3531880)+
  #scale_y_continuous(limits = c(0,max_y_value.mw.fdr),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                        axis.title = element_text(size = 50),
                        axis.text = element_text(size = 50),
                        legend.text = element_text(size = 50),
                        legend.title = element_text(size = 50),
                        plot.title = element_text(size = 50),
                        plot.title.position = "plot")
save(mw.fdr,file = "mw.fdr_all.rda")

#load all mann whitney bonfe
load(file = "list.1000.p.values.MW.bonfe.loop.wolf.rda")
merged_wolf_bonfe=unlist(merged_wolf_bonfe)
mw.bonfe.wolf=merged_wolf_bonfe
#max(table(cut(mw.bonfe.wolf,breaks=30)))
load(file = "list.1000.p.values.MW.bonfe.loop.pig.rda")
merged_pig_bonfe=unlist(merged_pig_bonfe)
mw.bonfe.pig=merged_pig_bonfe
max(table(cut(mw.bonfe.pig,breaks=30)))
load(file = "list.1000.p.values.MW.bonfe.loop.chicken.rda")
merged_chicken_bonfe=unlist(merged_chicken_bonfe)
mw.bonfe.chicken=merged_chicken_bonfe
vector.mw.bonfe=data.frame(value=c(mw.bonfe.chicken,mw.bonfe.pig,mw.bonfe.wolf),
                         Group=rep(c("Chicken","Pig","Wolf"),each=length(mw.bonfe.chicken)))
max_y_value.mw.bonfe=max(table(cut(vector.mw.bonfe$value,breaks=30)))
mw.bonfe=ggplot(vector.mw.bonfe,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",main="Mann-Whitney with bonfe correction")+
  ylim(0,9185243)+xlim(-0.1,1.1)+
  #scale_y_continuous(limits = c(0,max_y_value.mw.bonfe),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                        axis.title = element_text(size = 50),
                        axis.text = element_text(size = 50),
                        legend.text = element_text(size = 50),
                        legend.title = element_text(size = 50),
                        plot.title = element_text(size = 50),
                        plot.title.position = "plot")
save(mw.bonfe,file = "mw.bonfe_all.rda")

# load t test bonferroni
load(file="list.1000.p.values.loop.chicken_ttest_bonfe.rda")
ttest.bonfe.chicken=list.p.values.loop.1000
load(file="list.1000.p.values.loop.wolf_ttest_bonfe.rda")
ttest.bonfe.wolf=list.p.values.loop.1000
load(file="list.1000.p.values.loop.pig_ttest_bonfe.rda")
ttest.bonfe.pig=list.p.values.loop.1000
max(table(cut(ttest.bonfe.pig,breaks=30)))

#ttest bonfe
vector.ttest.bonfe=data.frame(value=c(ttest.bonfe.chicken,ttest.bonfe.pig,ttest.bonfe.wolf),
                  Group=rep(c("Chicken","Pig","Wolf"),each=length(ttest.bonfe.chicken)))
max_y_value.ttest.bonfe=max(table(cut(vector.ttest.bonfe$value,breaks=100)))
max_y_value.ttest.bonfe=max(table((vector.ttest.bonfe$value)))
ttest.bonfe=ggplot(vector.ttest.bonfe,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",main="Moderated t-test with Bonferroni correction")+
  ylim(0,9999872)+
  #scale_y_continuous(limits = c(0,max_y_value),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                axis.title = element_text(size = 50),
                axis.text = element_text(size = 50),
                legend.text = element_text(size = 50),
                legend.title = element_text(size = 50),
                plot.title = element_text(size = 50),
                plot.title.position = "plot")
save(ttest.bonfe,file = "ttest.bonfe_all.rda")

# load t test FDR
load(file="list.1000.p.values.loop.chicken_ttest_FDR.rda")
ttest.fdr.chicken=list.p.values.loop.1000
load(file="list.1000.p.values.loop.wolf_ttest_FDR.rda")
ttest.fdr.wolf=list.p.values.loop.1000
load(file="list.1000.p.values.loop.pig_ttest_FDR.rda")
ttest.fdr.pig=list.p.values.loop.1000
max(table(cut(ttest.fdr.pig,breaks=30)))

#ttest fdr
vector.ttest.fdr=data.frame(value=c(ttest.fdr.chicken,ttest.fdr.pig,ttest.fdr.wolf),
                  Group=rep(c("Chicken","Pig","Wolf"),each=length(ttest.fdr.chicken)))
max_y_value.ttest.fdr=max(table(cut(vector.ttest.fdr$value,breaks=30)))
ttest.fdr=ggplot(vector.ttest.fdr,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",main="Moderated t-test with FDR correction")+
  ylim(0,9999982)+
  #  scale_y_continuous(limits = c(0,max_y_value),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                axis.title = element_text(size = 50),
                axis.text = element_text(size = 50),
                legend.text = element_text(size = 50),
                legend.title = element_text(size = 50),
                plot.title = element_text(size = 50),
                plot.title.position = "plot")
save(ttest.fdr,file = "ttest.fdr_all.rda")

# load all quasi likelihood FDR
load(file = "list.1000.p.values.loop.chicken_QL_FDR.rda")
ql.fdr.chicken=list.p.values.loop.1000
max(table(cut(ql.fdr.chicken,breaks=30)))
load(file = "list.1000.p.values.loop.wolf_QL_FDR.rda")
ql.fdr.wolf=list.p.values.loop.1000
load(file = "list.1000.p.values.loop.pig_QL_FDR.rda")
ql.fdr.pig=list.p.values.loop.1000

#ql fdr
vector.ql.fdr=data.frame(value=c(ql.fdr.chicken,ql.fdr.pig,ql.fdr.wolf),
                            Group=rep(c("Chicken","Pig","Wolf"),each=length(ql.fdr.chicken)))
max_y_value.ql.fdr=max(table(cut(vector.ql.fdr$value,breaks=30)))
ql.fdr=ggplot(vector.ql.fdr,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",main="Quasi-likelihood with FDR correction")+
  ylim(0,5621381)+
  #scale_y_continuous(limits = c(0,max_y_value),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                axis.title = element_text(size = 50),
                axis.text = element_text(size = 50),
                legend.text = element_text(size = 50),
                legend.title = element_text(size = 50),
                plot.title = element_text(size = 50),
                plot.title.position = "plot")
save(ql.fdr,file = "ql.fdr_all.rda")

# load all quasi likelihood bonferroni
load(file = "list.1000.p.values.loop.chicken_QL_bonfe.rda")
ql.bonfe.chicken=list.p.values.loop.1000
load(file = "list.1000.p.values.loop.wolf_QL_bonfe.rda")
ql.bonfe.wolf=list.p.values.loop.1000
max(table(cut(ql.bonfe.wolf,breaks=30)))
load(file = "list.1000.p.values.loop.pig_QL_bonfe.rda")
ql.bonfe.pig=list.p.values.loop.1000

#ql bonfe
vector.ql.bonfe=data.frame(value=c(ql.bonfe.chicken,ql.bonfe.pig,ql.bonfe.wolf),
                              Group=rep(c("Chicken","Pig","Wolf"),each=length(ql.bonfe.chicken)))
max_y_value.ql.bonfe=max(table(cut(vector.ql.bonfe$value,breaks=30)))
ql.bonfe=ggplot(vector.ql.bonfe,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",main="Quasi-likelihood with bonferroni correction")+
  ylim(0,10000001)+
  #scale_y_continuous(limits = c(0,max_y_value),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                axis.title = element_text(size = 50),
                axis.text = element_text(size = 50),
                legend.text = element_text(size = 50),
                legend.title = element_text(size = 50),
                plot.title = element_text(size = 50),
                plot.title.position = "plot")
save(ql.bonfe,file = "ql.bonfe_all.rda")

# load all maximum likelihood bonferroni
load(file = "list.1000.p.values.loop.chicken_ML_bonferroni.rda")
ml.bonfe.chicken=list.p.values.loop.1000
load(file = "list.1000.p.values.loop.wolf_ML_bonferroni.rda")
ml.bonfe.wolf=list.p.values.loop.1000
load(file = "list.1000.p.values.loop.pig_ML_bonferroni.rda")
ml.bonfe.pig=list.p.values.loop.1000

#ml bonfe
vector.ml.bonfe=data.frame(value=c(ml.bonfe.chicken,ml.bonfe.pig,ml.bonfe.wolf),
                              Group=rep(c("Chicken","Pig","Wolf"),each=length(ml.bonfe.chicken)))
max_y_value.ml.bonfe=max(table(cut(vector.ml.bonfe$value,breaks=30)))
y.final=max(vector.ml.bonfe$value)+0.1
ml.bonfe=ggplot(vector.ml.bonfe,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",main="Maximum-likelihood with Bonferroni correction")+
  ylim(0,10000001)+xlim(-0.1,y.final)+
  #scale_y_continuous(limits = c(0,max_y_value),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 50),
                axis.title = element_text(size = 50),
                axis.text = element_text(size = 50),
                legend.text = element_text(size = 50),
                legend.title = element_text(size = 50),
                plot.title = element_text(size = 50),
                plot.title.position = "plot")
save(ml.bonfe,file = "ml.bonfe_all.rda")

#load all maximum likelihood FDR
load(file = "list.1000.p.values.loop.chicken_ML_FDR.rda")
ml.fdr.chicken=list.p.values.loop.1000
load(file = "list.1000.p.values.loop.wolf_ML_FDR.rda")
ml.fdr.wolf=list.p.values.loop.1000
load(file = "list.1000.p.values.loop.pig_ML_FDR.rda")
ml.fdr.pig=list.p.values.loop.1000
max(table(cut(ml.fdr.pig,breaks=30)))

#ml fdr
vector.ml.fdr=data.frame(value=c(ml.fdr.chicken,ml.fdr.pig,ml.fdr.wolf),
                            Group=rep(c("Chicken","Pig","Wolf"),each=length(ml.fdr.chicken)))
max_y_value.ml.fdr=max(table(cut(vector.ml.fdr$value,breaks=30)))
ml.fdr=ggplot(vector.ml.fdr,aes(x=value,fill=Group))+
  geom_histogram(position = "identity",alpha=0.5,bins=30)+
  scale_fill_manual(values = c("Chicken"="#E69F00","Pig"="#56B4E9","Wolf"="#009E73"))+
  labs(x="P-values",y="Frequency",main="Maximum-likelihood with FDR correction")+
  ylim(0,3959416)+
  #scale_y_continuous(limits = c(0,max_y_value),expand=expansion(mult=c(0,0.05)))+
  theme_minimal()+theme(text = element_text(size = 40),
                axis.title = element_text(size = 50),
                axis.text = element_text(size = 50),
                legend.text = element_text(size = 50),
                legend.title = element_text(size = 50),
                plot.title = element_text(size = 50),
                plot.title.position = "plot")
save(ml.fdr,file = "ml.fdr_all.rda")

#merge everything
load("ml.fdr_all.rda")
load("ml.bonfe_all.rda")
load("ql.bonfe_all.rda")
load("ql.fdr_all.rda")
load("ttest.fdr_all.rda")
load("ttest.bonfe_all.rda")
load("mw.bonfe_all.rda")
load("mw.fdr_all.rda")

# large_font=theme(text = element_text(size = 20),
#                  axis.title = element_text(size = 20),
#                  axis.text = element_text(size = 20),
#                  legend.text = element_text(size = 20),
#                  legend.title = element_text(size = 20),
#                  plot.title = element_text(size = 20))
legend.all=get_legend(ttest.bonfe+theme(legend.position = "bottom"))
ttest.bonfe=ttest.bonfe+theme(legend.position = "none")
ttest.fdr=ttest.fdr+theme(legend.position = "none")
ql.bonfe=ql.bonfe+theme(legend.position = "none")
ql.fdr=ql.fdr+theme(legend.position = "none")
ml.bonfe=ml.bonfe+theme(legend.position = "none")
ml.fdr=ml.fdr+theme(legend.position = "none")
mw.bonfe=mw.bonfe+theme(legend.position = "none")
mw.fdr=mw.fdr+theme(legend.position = "none")
combined_all=plot_grid(mw.fdr,mw.bonfe,ml.fdr,ml.bonfe,ql.fdr,ql.bonfe,ttest.fdr,ttest.bonfe,ncol = 2,align = "v")
final_all=plot_grid(combined_all,legend.all,ncol = 1,rel_heights = c(1,0.1))
save(final_all,file = "final_plot_null_all.rda")
# tiff(filename = "plot_all_butMW_null.tiff",width = 1000,height = 10000,units = "px",res = 500)
# print(final_all)
# dev.off()

tiff(filename = "plot_all_butMW_null.tiff",width = 6000,height = 6000,units = "px",res = 100)
print(final_all)
dev.off()








