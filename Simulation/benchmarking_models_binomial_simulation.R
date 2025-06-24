library(stringr)
library(edgeR)
library(dplyr)
library(PLNmodels)
library(corrplot)
library(factoextra)
library(FactoMineR)
library(ade4)
library(gridExtra)
setwd("/proj/naiss2024-23-57/GBS_MeDIP_benchmark/bin")
nind<-50
#first 3 are high DMRs, next 3 are low DMRs, next 3 are no DMRs
wind_nb<-c(300,300,300,300,300,200,100,50,50)
coverage<-c(0.1,0.5,0.7,0.9,1,5,10,20,50)
sd_cov<-c(10,10,10,10,1,1,1,2,5)
groups<-c('case','control')
p_groups<-list(c(0.8,0.7,0.67,0.65,0.63,0.6,0.57,0.5,0.4),c(0.1,0.2,0.3,0.33,0.35,0.37,0.4, 0.44,0.5))

Mat<-matrix(nrow=1,ncol=nind)
Design<-list()
for(a in 1:length(groups)){
  df_design<-data.frame(coverage=NA,sd_cov=NA,prob=NA,test=NA,wind_id=NA)
  prob<-p_groups[[a]]
  for(i in prob){
    for(j in 1:length(wind_nb)){
      cov<-ceiling(rnorm(wind_nb[j],coverage[j],sd_cov[j]))
      cov[cov<0]<-0
      ind<-list()
      for(k in 1:nind){ind[[k]]<-rbinom(wind_nb[j],cov,i)}
      ind<-do.call(cbind,ind)
      Mat<-rbind(Mat,ind)
      x<-data.frame(coverage=cov,sd_cov=sd_cov[j],prob=i,test=groups[[a]],wind_id=NA)
      df_design<-rbind(df_design,x)
    }
  }
  df_design$wind_id<-c('wind_0',paste0('wind',seq(1,sum(wind_nb)*length(prob),1)))
  Design[[a]]<-df_design
}
Mat<-Mat[-1,]
Design<-do.call(rbind,Design)
Design<-Design[!is.na(Design$coverage),]

Mat<-t(Mat)
colnames(Mat)<-Design$test
Mat_group<-list()
for(i in 1:length(groups)){Mat_group[[i]]<-Mat[,colnames(Mat)==groups[i]]}
df_Mat<-do.call(rbind,Mat_group)
colnames(df_Mat)<-paste0('wind',seq(1,ncol(df_Mat),by=1))
rownames(df_Mat)<-unlist(lapply(groups,function(x){paste0(x,seq(1:nind))}))
df_Mat<-as.data.frame(df_Mat)
df_design<-data.frame(group=unlist(lapply(groups,function(x){paste0(x,seq(1:nind))})))
df_design$group=str_replace(df_design$group,"[0-9]","")
df_design$group=str_replace(df_design$group,"[0-9]","")
rownames(df_design)<-rownames(df_Mat)
#in Design you have all the conditions (high DMR, low DMR, no DMR 9 times),
desing_win=data.frame(unique(Design$wind_id),c(rep("high DMR",8100),rep("low DMR",7200),rep("nonDMR",1800)),c(rep(1,11400),rep(0,5700)))
colnames(desing_win)=c("wind_id","classification","binary")
DMR=(Design[Design$prob>=0.6|Design$prob<=0.37,])
case=DMR[DMR$test=="case",]
control=DMR[DMR$test=="control",]
DMR=left_join(case,control,by='wind_id')
nonDMR=Design[Design$prob<=0.57&Design$prob>=0.4,]
case=nonDMR[nonDMR$test=="case",]
control=nonDMR[nonDMR$test=="control",]
nonDMR=left_join(case,control,by='wind_id')

##### EdgeR methods, normal dispersion ######
edgeR.RBC <- DGEList(counts=t(df_Mat), genes=colnames(df_Mat))
rownames(edgeR.RBC$counts) <- rownames(edgeR.RBC$genes) <- colnames(df_Mat)
edgeR.RBC$samples$lib.size <- colSums(edgeR.RBC$counts)
edgeR.RBC <- calcNormFactors(edgeR.RBC)
design=model.matrix(~df_design$group)
edgeR.RBC <- estimateDisp(edgeR.RBC, design)
#Maximum-likelihood approach for calculating the coefficients
model.fit.RBC=glmFit(edgeR.RBC, design)
# Now do likelihood ratio tests to see which windows are differentialy methylated
like.ratio.RBC = glmLRT(model.fit.RBC)
#FDR
top.RBC=topTags(like.ratio.RBC, n=ncol(df_Mat), adjust.method = "fdr")
# getting the windows that are marked as significant
ML.sig=desing_win[desing_win$wind_id%in%row.names(top.RBC$table)[top.RBC$table$FDR<0.01],]
# getting the windows that are marked as non-significant
ML.nonsig=desing_win[desing_win$wind_id%in%row.names(top.RBC$table)[top.RBC$table$FDR>0.01],]
# get the true positives
TP_ML=nrow(semi_join(ML.sig,DMR,by='wind_id'))
# get the false positives
FP_ML=nrow(ML.sig)-TP_ML
# get the true negatives
TN_ML=nrow(semi_join(ML.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_ML=nrow(ML.nonsig)-TN_ML

#instead of FDR do bonferroni
top.RBC.bonfe=topTags(like.ratio.RBC, n=ncol(df_Mat), adjust.method = "bonferroni")
# getting the windows that are marked as significant
ML.bonfe.sig=desing_win[desing_win$wind_id%in%row.names(top.RBC.bonfe$table)[top.RBC.bonfe$table$FWER<0.01],]
# getting the windows that are marked as non-significant
ML.bonfe.nonsig=desing_win[desing_win$wind_id%in%row.names(top.RBC.bonfe$table)[top.RBC.bonfe$table$FWER>0.01],]
# get the true positives
TP_ML.bonfe=nrow(semi_join(ML.bonfe.sig,DMR,by='wind_id'))
# get the false positives
FP_ML.bonfe=nrow(ML.bonfe.sig)-TP_ML.bonfe
# get the true negatives
TN_ML.bonfe=nrow(semi_join(ML.bonfe.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_ML.bonfe=nrow(ML.bonfe.nonsig)-TN_ML.bonfe

### Fit the quasi-likelihood Negative Binomial model ####
model.quasi.fit.RBC=glmQLFit(edgeR.RBC, design)
#Now do quasi-likelihood ratio tests to see which ROIS are DE
quasi.like.ratio.RBC = glmQLFTest(model.quasi.fit.RBC, coef = 2)
#Do the False Discovery Rate
top.quasi.RBC=topTags(quasi.like.ratio.RBC, n=ncol(df_Mat), adjust.method = "fdr")
# getting the windows that are marked as significant
QL.sig=desing_win[desing_win$wind_id%in%row.names(top.quasi.RBC$table)[top.quasi.RBC$table$FDR<0.01],]
# getting the windows that are marked as non-significant
QL.nonsig=desing_win[desing_win$wind_id%in%row.names(top.quasi.RBC$table)[top.quasi.RBC$table$FDR>0.01],]
# get the true positives
TP_QL=nrow(semi_join(QL.sig,DMR,by='wind_id'))
# get the false positives
FP_QL=nrow(QL.sig)-TP_QL
# get the true negatives
TN_QL=nrow(semi_join(QL.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_QL=nrow(QL.nonsig)-TN_QL

#do bonferroni instead of BH
#Do the False Discovery Rate
top.quasi.bonfe.RBC=topTags(quasi.like.ratio.RBC, n=ncol(df_Mat), adjust.method = "bonferroni")
# getting the windows that are marked as significant
QL.bonfe.sig=desing_win[desing_win$wind_id%in%row.names(top.quasi.bonfe.RBC$table)[top.quasi.bonfe.RBC$table$FWER<0.01],]
# getting the windows that are marked as non-significant
QL.bonfe.nonsig=desing_win[desing_win$wind_id%in%row.names(top.quasi.bonfe.RBC$table)[top.quasi.bonfe.RBC$table$FWER>0.01],]
# get the true positives
TP_QL.bonfe=nrow(semi_join(QL.bonfe.sig,DMR,by='wind_id'))
# get the false positives
FP_QL.bonfe=nrow(QL.bonfe.sig)-TP_QL.bonfe
# get the true negatives
TN_QL.bonfe=nrow(semi_join(QL.bonfe.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_QL.bonfe=nrow(QL.bonfe.nonsig)-TN_QL.bonfe

##### EdgeR methods, robust dispersion ######
edgeR.robust.RBC <- DGEList(counts=t(df_Mat), genes=colnames(df_Mat))
rownames(edgeR.robust.RBC$counts) <- rownames(edgeR.robust.RBC$genes) <- colnames(df_Mat)
edgeR.robust.RBC$samples$lib.size <- colSums(edgeR.robust.RBC$counts)
edgeR.robust.RBC <- calcNormFactors(edgeR.robust.RBC)
design=model.matrix(~df_design$group)
edgeR.robust.RBC <- estimateDisp(edgeR.robust.RBC, design, robust=TRUE)
#Maximum-likelihood approach for calculating the coefficients
model.fit.robust.RBC=glmFit(edgeR.robust.RBC, design)
# Now do likelihood ratio tests to see which windows are differentialy methylated
like.ratio.robust.RBC = glmLRT(model.fit.robust.RBC)
#FDR
top.robust.RBC=topTags(like.ratio.robust.RBC, n=ncol(df_Mat), adjust.method = "fdr")
# getting the windows that are marked as significant
ML.rob.sig=desing_win[desing_win$wind_id%in%row.names(top.robust.RBC$table)[top.robust.RBC$table$FDR<0.01],]
# getting the windows that are marked as non-significant
ML.rob.nonsig=desing_win[desing_win$wind_id%in%row.names(top.robust.RBC$table)[top.robust.RBC$table$FDR>0.01],]
# get the true positives
TP_ML.rob=nrow(semi_join(ML.rob.sig,DMR,by='wind_id'))
# get the false positives
FP_ML.rob=nrow(ML.rob.sig)-TP_ML.rob
# get the true negatives
TN_ML.rob=nrow(semi_join(ML.rob.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_ML.rob=nrow(ML.rob.nonsig)-TN_ML.rob

#instead do bonferroni
top.robust.RBC.bonfe=topTags(like.ratio.robust.RBC, n=ncol(df_Mat), adjust.method = "bonferroni")
# getting the windows that are marked as significant
ML.rob.bonfesig=desing_win[desing_win$wind_id%in%row.names(top.robust.RBC.bonfe$table)[top.robust.RBC.bonfe$table$FWER<0.01],]
# getting the windows that are marked as non-significant
ML.rob.bonfenonsig=desing_win[desing_win$wind_id%in%row.names(top.robust.RBC.bonfe$table)[top.robust.RBC.bonfe$table$FWER>0.01],]
# get the true positives
TP_ML.rob.bonfe=nrow(semi_join(ML.rob.bonfesig,DMR,by='wind_id'))
# get the false positives
FP_ML.rob.bonfe=nrow(ML.rob.bonfesig)-TP_ML.rob.bonfe
# get the true negatives
TN_ML.rob.bonfe=nrow(semi_join(ML.rob.bonfenonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_ML.rob.bonfe=nrow(ML.rob.bonfenonsig)-TN_ML.rob.bonfe

### Fit the quasi-likelihood Negative Binomial model robust ####
model.quasi.fit.robust.RBC=glmQLFit(edgeR.robust.RBC, design)
#Now do quasi-likelihood ratio tests to see which ROIS are DE
quasi.like.ratio.robust.RBC = glmQLFTest(model.quasi.fit.robust.RBC, coef = 2)
#Do the False Discovery Rate
top.quasi.robust.RBC=topTags(quasi.like.ratio.robust.RBC, n=ncol(df_Mat), adjust.method = "fdr")
# getting the windows that are marked as significant
QL.rob.sig=desing_win[desing_win$wind_id%in%row.names(top.quasi.robust.RBC$table)[top.quasi.robust.RBC$table$FDR<0.01],]
# getting the windows that are marked as non-significant
QL.rob.nonsig=desing_win[desing_win$wind_id%in%row.names(top.quasi.robust.RBC$table)[top.quasi.robust.RBC$table$FDR>0.01],]
# get the true positives
TP_QL.rob=nrow(semi_join(QL.rob.sig,DMR,by='wind_id'))
# get the false positives
FP_QL.rob=nrow(QL.rob.sig)-TP_QL.rob
# get the true negatives
TN_QL.rob=nrow(semi_join(QL.rob.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_QL.rob=nrow(QL.rob.nonsig)-TN_QL.rob

#do bonferroni
top.quasirob.bonfeust.RBC=topTags(quasi.like.ratio.robust.RBC, n=ncol(df_Mat), adjust.method = "bonferroni")
# getting the windows that are marked as significant
QLrob.bonfe.sig=desing_win[desing_win$wind_id%in%row.names(top.quasirob.bonfeust.RBC$table)[top.quasirob.bonfeust.RBC$table$FWER<0.01],]
# getting the windows that are marked as non-significant
QLrob.bonfe.nonsig=desing_win[desing_win$wind_id%in%row.names(top.quasirob.bonfeust.RBC$table)[top.quasirob.bonfeust.RBC$table$FWER>0.01],]
# get the true positives
TP_QLrob.bonfe=nrow(semi_join(QLrob.bonfe.sig,DMR,by='wind_id'))
# get the false positives
FP_QLrob.bonfe=nrow(QLrob.bonfe.sig)-TP_QLrob.bonfe
# get the true negatives
TN_QLrob.bonfe=nrow(semi_join(QLrob.bonfe.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_QLrob.bonfe=nrow(QLrob.bonfe.nonsig)-TN_QLrob.bonfe

##### moderated t test #####
edgeR <- DGEList(counts=t(df_Mat), genes=colnames(df_Mat))
rownames(edgeR$counts) <- rownames(edgeR$genes) <- colnames(df_Mat)
edgeR <- calcNormFactors(edgeR)
design=model.matrix(~df_design$group)
logCPM <- cpm(edgeR, log=TRUE)
fit <- lmFit(logCPM, design)
fit <- eBayes(fit, trend=TRUE)
top.p.values.t.test.limma=topTable(fit, coef=ncol(design),number = ncol(df_Mat))
top.p.values.t.test.limma$corrected_p_value = p.adjust(((top.p.values.t.test.limma$P.Value)), method = "BH")
# getting the windows that are marked as significant
ttest.sig=desing_win[desing_win$wind_id%in%row.names(top.p.values.t.test.limma)[top.p.values.t.test.limma$corrected_p_value<0.01],]
# getting the windows that are marked as non-significant
ttest.nonsig=desing_win[desing_win$wind_id%in%row.names(top.p.values.t.test.limma)[top.p.values.t.test.limma$corrected_p_value>0.01],]
# get the true positives
TP_ttest=nrow(semi_join(ttest.sig,DMR,by='wind_id'))
# get the false positives
FP_ttest=nrow(ttest.sig)-TP_ttest
# get the true negatives
TN_ttest=nrow(semi_join(ttest.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_ttest=nrow(ttest.nonsig)-TN_ttest

#bonferroni
top.p.values.ttest.bonfe.limma=top.p.values.t.test.limma
top.p.values.ttest.bonfe.limma$corrected_p_value = p.adjust(((top.p.values.ttest.bonfe.limma$P.Value)), method = "bonferroni")
# getting the windows that are marked as significant
ttest.bonfe.sig=desing_win[desing_win$wind_id%in%row.names(top.p.values.ttest.bonfe.limma)[top.p.values.ttest.bonfe.limma$corrected_p_value<0.01],]
# getting the windows that are marked as non-significant
ttest.bonfe.nonsig=desing_win[desing_win$wind_id%in%row.names(top.p.values.ttest.bonfe.limma)[top.p.values.ttest.bonfe.limma$corrected_p_value>0.01],]
# get the true positives
TP_ttest.bonfe=nrow(semi_join(ttest.bonfe.sig,DMR,by='wind_id'))
# get the false positives
FP_ttest.bonfe=nrow(ttest.bonfe.sig)-TP_ttest.bonfe
# get the true negatives
TN_ttest.bonfe=nrow(semi_join(ttest.bonfe.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_ttest.bonfe=nrow(ttest.bonfe.nonsig)-TN_ttest.bonfe

##### mann whitney #####
edgeR.MW <- DGEList(counts=t(df_Mat), genes=colnames(df_Mat))
rownames(edgeR.MW$counts) <- rownames(edgeR.MW$genes) <- colnames(df_Mat)
edgeR.MW <- calcNormFactors(edgeR.MW)
efective= as.data.frame(edgeR.MW$samples$lib.size*edgeR.MW$samples$norm.factors)

edgeR.MW=data.frame(mapply(`*`,(df_Mat),(efective)))
p.values.mann.whitney=apply((edgeR.MW), 2, function(x){wilcox.test(x~df_design$group)$p.value})
p.values.mann.whitney=as.data.frame(p.values.mann.whitney)
colnames(p.values.mann.whitney)=c("p_value")
p.values.mann.whitney$corrected_p_value = p.adjust(((p.values.mann.whitney$p_value)), method = "BH")
# getting the windows that are marked as significant
MW.sig=desing_win[desing_win$wind_id%in%row.names(p.values.mann.whitney)[p.values.mann.whitney$corrected_p_value<0.01],]
# getting the windows that are marked as non-significant
MW.nonsig=desing_win[desing_win$wind_id%in%row.names(p.values.mann.whitney)[p.values.mann.whitney$corrected_p_value>=0.01],]
# get the true positives
TP_MW=nrow(semi_join(MW.sig,DMR,by='wind_id'))
# get the false positives
FP_MW=nrow(MW.sig)-TP_MW
# get the true negatives
TN_MW=nrow(semi_join(MW.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_MW=nrow(MW.nonsig)-TN_MW

#bonferroni
p.values.mann.bonfe.whitney=p.values.mann.whitney
p.values.mann.bonfe.whitney$corrected_p_value = p.adjust(((p.values.mann.bonfe.whitney$p_value)), method = "bonferroni")
# getting the windows that are marked as significant
MW.bonfe.sig=desing_win[desing_win$wind_id%in%row.names(p.values.mann.bonfe.whitney)[p.values.mann.bonfe.whitney$corrected_p_value<0.01],]
# getting the windows that are marked as non-significant
MW.bonfe.nonsig=desing_win[desing_win$wind_id%in%row.names(p.values.mann.bonfe.whitney)[p.values.mann.bonfe.whitney$corrected_p_value>=0.01],]
# get the true positives
TP_MW.bonfe=nrow(semi_join(MW.bonfe.sig,DMR,by='wind_id'))
# get the false positives
FP_MW.bonfe=nrow(MW.bonfe.sig)-TP_MW.bonfe
# get the true negatives
TN_MW.bonfe=nrow(semi_join(MW.bonfe.nonsig,nonDMR,by='wind_id'))
# get the false negatives
FN_MW.bonfe=nrow(MW.bonfe.nonsig)-TN_MW.bonfe

# ##### PCA+BCA approach #####
# # for the sake of consistency i am going to normalize the counts with the TMM normalization factor + CPM
# edgeR.BCA <- DGEList(counts=t(df_Mat), genes=colnames(df_Mat))
# rownames(edgeR.BCA$counts) <- rownames(edgeR.BCA$genes) <- colnames(df_Mat)
# edgeR.BCA <- calcNormFactors(edgeR.BCA)
# efective= as.data.frame(edgeR.BCA$samples$lib.size*edgeR.BCA$samples$norm.factors)
# edgeR.BCA=data.frame(mapply(`*`,(df_Mat),efective))
# row.names(edgeR.BCA)=row.names(df_Mat)
# #PCA to have the eigthteen vectors and to see how your data looks like
# pca_TMM_CPM<-dudi.pca(edgeR.BCA,scannf=F,nf=5) #keep the first 5 axis
# fviz_pca_biplot(pca_TMM_CPM,col.var='contrib')
# # now the eigthteen vectors is passed to the BCA, with this BCA you get the windows that contributes the most to the difference between the two groups
# mod_dud_TMM_CPM<-bca(pca_TMM_CPM,as.factor(df_design$group),scannf=F,nf=1)
# # now that we have the contribution of each window to each of the groups, we rank them by order of contribution
# contribution_TMM_CPM<-mod_dud_TMM_CPM$co
# rownames(contribution_TMM_CPM)[order(contribution_TMM_CPM[,1])]
# barplot(abs(contribution_TMM_CPM[,1]))
# #permutation test to have the threshold tailored to each experiment
# threshold_ini<-0
# threshold<-10
# perm<-list()
# while(abs(threshold_ini-threshold)>2){
#   threshold_ini<-threshold
#   for(i in 1:1000){
#     mod<-bca(pca_TMM_CPM,sample(as.factor(df_design$group)),scannf=F,nf=1)
#     perm[[i]]<-mod$co
#   }
#   perm<-do.call(rbind,perm)
#   
# }
# #now that you got the permutation where you have the null population, set the threshold at the upper quantile for 0.01 (so 1%) 
# threshold=abs(quantile(perm[,1],prob=0.01/length(unique(desing_win$wind_id))))
# # getting the windows that are marked as significant
# contribution_TMM_CPM=data.frame(wind_id=row.names(contribution_TMM_CPM),contribution=contribution_TMM_CPM$Comp1)
# sign.PCA=as.data.frame((contribution_TMM_CPM[abs(contribution_TMM_CPM$contribution)>threshold,]))
# # getting the windows that are marked as non-significant
# PCA.nonsig=as.data.frame((contribution_TMM_CPM[abs(contribution_TMM_CPM$contribution)<threshold,]))
# # get the true positives
# TP_PCA=nrow(semi_join(sign.PCA,DMR,by='wind_id'))
# # get the false positives
# FP_PCA=nrow(sign.PCA)-TP_PCA
# # get the true negatives
# TN_PCA=nrow(semi_join(PCA.nonsig,nonDMR,by='wind_id'))
# # get the false negatives
# FN_PCA=nrow(PCA.nonsig)-TN_PCA


##### Getting accuracy ######
acc.ML=(TP_ML+TN_ML)/(TP_ML+TN_ML+FN_ML+FP_ML)
acc.ML.bonfe=(TP_ML.bonfe+TN_ML.bonfe)/(TP_ML.bonfe+TN_ML.bonfe+FN_ML.bonfe+FP_ML.bonfe)
acc.QL=(TP_QL+TN_QL)/(TP_QL+TN_QL+FN_QL+FP_QL)
acc.QL.bonfe=(TP_QL.bonfe+TN_QL.bonfe)/(TP_QL.bonfe+TN_QL.bonfe+FN_QL.bonfe+FP_QL.bonfe)
acc.ML.rob=(TP_ML.rob+TN_ML.rob)/(TP_ML.rob+TN_ML.rob+FN_ML.rob+FP_ML.rob)
acc.ML.rob.bonfe=(TP_ML.rob.bonfe+TN_ML.rob.bonfe)/(TP_ML.rob.bonfe+TN_ML.rob.bonfe+FN_ML.rob.bonfe+FP_ML.rob.bonfe)
acc.QL.rob=(TP_QL.rob+TN_QL.rob)/(TP_QL.rob+TN_QL.rob+FN_QL.rob+FP_QL.rob)
acc.ttest=(TP_ttest+TN_ttest)/(TP_ttest+TN_ttest+FN_ttest+FP_ttest)
acc.ttest.bonfe=(TP_ttest.bonfe+TN_ttest.bonfe)/(TP_ttest.bonfe+TN_ttest.bonfe+FN_ttest.bonfe+FP_ttest.bonfe)
acc.MW=(TP_MW+TN_MW)/(TP_MW+TN_MW+FN_MW+FP_MW)
acc.MW.bonfe=(TP_MW.bonfe+TN_MW.bonfe)/(TP_MW.bonfe+TN_MW.bonfe+FN_MW.bonfe+FP_MW.bonfe)
#acc.PCA=(TP_PCA+TN_PCA)/(TP_PCA+FP_PCA+TN_PCA+FN_PCA)
acc.total=data.frame("maximum-likelihood"=acc.ML,"quasi-likelihood"=acc.QL,"maximum-likelihood+observational.weigths"=acc.ML.rob,"quasi-likelihood+observational.weigths"=acc.QL.rob,"moderated t-test"=acc.ttest,"Mann-whitney"=acc.MW,"maximum-likelihood bonferroni"=acc.ML.bonfe,"quasi-likelihood bonferroni"=acc.QL.bonfe,"maximum-likelihood+observational.weigths bonferroni"=acc.ML.rob.bonfe,"moderated t-test bonferroni"=acc.ttest.bonfe,"Mann-whitney bonferroni"=acc.MW.bonfe)
write.table(acc.total,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/accuracy_simulation.txt",quote = F,sep = "\t",row.names = F)

##### Getting precision #####
pre.ML=TP_ML/(TP_ML+FP_ML)
pre.ML.bonfe=TP_ML.bonfe/(TP_ML.bonfe+FP_ML.bonfe)
pre.QL=TP_QL/(TP_QL+FP_QL)
pre.QL.bonfe=TP_QL.bonfe/(TP_QL.bonfe+FP_QL.bonfe)
pre.ML.rob=TP_ML.rob/(TP_ML.rob+FP_ML.rob)
pre.QL.rob=TP_QL.rob/(TP_QL.rob+FP_QL.rob)
pre.ttest=TP_ttest/(TP_ttest+FP_ttest)
pre.ttest.bonfe=TP_ttest.bonfe/(TP_ttest.bonfe+FP_ttest.bonfe)
pre.MW=TP_MW/(FP_MW+TP_MW)
pre.MW.bonfe=TP_MW.bonfe/(FP_MW.bonfe+TP_MW.bonfe)
# pre.PCA=TP_PCA/(FP_PCA+TP_PCA)
pre.total=data.frame("maximum-likelihood"=pre.ML,"quasi-likelihood"=pre.QL,"maximum-likelihood+observational.weigths"=pre.ML.rob,"quasi-likelihood+observational.weigths"=pre.QL.rob,"moderated t-test"=pre.ttest,"Mann-whitney"=pre.MW,"maximum-likelihood bonferroni"=pre.ML.bonfe,"quasi-likelihood bonferroni"=pre.QL.bonfe,"moderated t-test bonferroni"=pre.ttest.bonfe,"Mann-whitney bonferroni"=pre.MW.bonfe)
write.table(pre.total,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/precision_simulation.txt",quote = F,sep = "\t",row.names = F)

##### Getting recall #####
rec.ML=TP_ML/(TP_ML+FN_ML)
rec.ML.bonfe=TP_ML.bonfe/(TP_ML.bonfe+FN_ML.bonfe)
rec.QL=TP_QL/(TP_QL+FN_QL)
rec.QL.bonfe=TP_QL.bonfe/(TP_QL.bonfe+FN_QL.bonfe)
rec.ML.rob=TP_ML.rob/(TP_ML.rob+FN_ML.rob)
rec.QL.rob=TP_QL.rob/(TP_QL.rob+FN_QL.rob)
rec.ttest=TP_ttest/(TP_ttest+FN_ttest)
rec.ttest.bonfe=TP_ttest.bonfe/(TP_ttest.bonfe+FN_ttest.bonfe)
rec.MW=TP_MW/(TP_MW+FN_MW)
rec.MW.bonfe=TP_MW.bonfe/(TP_MW.bonfe+FN_MW.bonfe)
# rec.PCA=TP_PCA/(TP_PCA+FN_PCA)
rec.total=data.frame("maximum-likelihood"=rec.ML,"quasi-likelihood"=rec.QL,"maximum-likelihood+observational.weigths"=rec.ML.rob,"quasi-likelihood+observational.weigths"=rec.QL.rob,"moderated t-test"=rec.ttest,"Mann-whitney"=rec.MW,"maximum-likelihood bonferroni"=rec.ML.bonfe,"quasi-likelihood bonferroni"=rec.QL.bonfe,"moderated t-test bonferroni"=rec.ttest.bonfe,"Mann-whitney bonferroni"=rec.MW.bonfe)

write.table(rec.total,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/recall_simulation.txt",quote = F,sep = "\t",row.names = F)

##### ROC plots #####
library(pROC)
# the direction of the arrow symbolizes where is the significance, < means the lower the p-value the more significance
names(top.RBC$table)[names(top.RBC$table) == 'genes'] <- 'wind_id'
final.ml=left_join(desing_win,top.RBC$table,by='wind_id')
final.ml=final.ml[order(final.ml$FDR),]
roc.ml=roc(final.ml$binary,final.ml$FDR,direction=">")
auc.ml=auc(final.ml$binary,final.ml$FDR,direction=">")

names(top.RBC.bonfe$table)[names(top.RBC.bonfe$table) == 'genes'] <- 'wind_id'
final.ml.bonfe=left_join(desing_win,top.RBC.bonfe$table,by='wind_id')
final.ml.bonfe=final.ml.bonfe[order(final.ml.bonfe$FWER),]
roc.ml.bonfe=roc(final.ml.bonfe$binary,final.ml.bonfe$FWER,direction=">")
auc.ml.bonfe=auc(final.ml.bonfe$binary,final.ml.bonfe$FWER,direction=">")

names(top.quasi.RBC$table)[names(top.quasi.RBC$table) == 'genes'] <- 'wind_id'
final.ql=left_join(desing_win,top.quasi.RBC$table,by='wind_id')
roc.ql=roc(final.ql$binary,final.ql$FDR,direction=">")
auc.ql=auc(final.ql$binary,final.ql$FDR,direction=">")

names(top.quasi.bonfe.RBC$table)[names(top.quasi.bonfe.RBC$table) == 'genes'] <- 'wind_id'
final.ql.bonfe=left_join(desing_win,top.quasi.bonfe.RBC$table,by='wind_id')
roc.ql.bonfe=roc(final.ql.bonfe$binary,final.ql.bonfe$FWER,direction=">")
auc.ql.bonfe=auc(final.ql.bonfe$binary,final.ql.bonfe$FWER,direction=">")

names(top.robust.RBC$table)[names(top.robust.RBC$table) == 'genes'] <- 'wind_id'
final.ml.rob=left_join(desing_win,top.robust.RBC$table,by='wind_id')
roc.ml.rob=roc(final.ml.rob$binary,final.ml.rob$FDR,direction=">")
auc.ml.rob=auc(final.ml.rob$binary,final.ml.rob$FDR,direction=">")

names(top.quasi.robust.RBC$table)[names(top.quasi.robust.RBC$table) == 'genes'] <- 'wind_id'
final.ql.rob=left_join(desing_win,top.quasi.robust.RBC$table,by='wind_id')
roc.ql.rob=roc(final.ql.rob$binary,final.ql.rob$FDR,direction=">")
auc.ql.rob=auc(final.ql.rob$binary,final.ql.rob$FDR,direction=">")

p.values.mann.whitney$wind_id=row.names(p.values.mann.whitney)
final.mw=left_join(desing_win,p.values.mann.whitney,by='wind_id')
roc.mw=roc(final.mw$binary,final.mw$corrected_p_value,direction=">")
auc.mw=auc(final.mw$binary,final.mw$corrected_p_value,direction=">")

p.values.mann.bonfe.whitney$wind_id=row.names(p.values.mann.bonfe.whitney)
final.mw.bonfe=left_join(desing_win,p.values.mann.bonfe.whitney,by='wind_id')
roc.mw.bonfe=roc(final.mw.bonfe$binary,final.mw.bonfe$corrected_p_value,direction=">")
auc.mw.bonfe=auc(final.mw.bonfe$binary,final.mw.bonfe$corrected_p_value,direction=">")

top.p.values.t.test.limma$wind_id=row.names(top.p.values.t.test.limma)
final.ttest=left_join(desing_win,top.p.values.t.test.limma,by='wind_id')
roc.ttest=roc(final.ttest$binary,final.ttest$corrected_p_value,direction=">")
auc.ttest=auc(final.ttest$binary,final.ttest$corrected_p_value,direction=">")

top.p.values.ttest.bonfe.limma$wind_id=row.names(top.p.values.ttest.bonfe.limma)
final.ttest.bonfe=left_join(desing_win,top.p.values.ttest.bonfe.limma,by='wind_id')
roc.ttest.bonfe=roc(final.ttest.bonfe$binary,final.ttest.bonfe$corrected_p_value,direction=">")
auc.ttest.bonfe=auc(final.ttest.bonfe$binary,final.ttest.bonfe$corrected_p_value,direction=">")

# final.PCA=left_join(desing_win,contribution_TMM_CPM,by='wind_id')
# roc.PCA=roc(final.PCA$binary,final.PCA$contribution)
# auc.PCA=auc(final.PCA$binary,final.PCA$contribution,direction="<")
# auc.total=data.frame("Maximum likelihood with FDR correction"=auc.ml,"Maximum likelihood with Bonferroni correction"=auc.ml.bonfe,"Maximum likelihood with observational weights with FDR correction"=auc.ml.rob,"Quasi-likelihood with FDR correction"=auc.ql,"Quasi-likelihood with Bonferroni correction"=auc.ql.bonfe,"Quasi-likelihood with observational weights with FDR correction"=auc.ql.rob,"Moderated t-test with FDR correction"=auc.ttest,"Moderated t-test with Bonferroni correction"=auc.ttest.bonfe,"Mann Withney with FDR correction"=auc.mw,"Mann Withney with Bonferroni correction"=auc.mw.bonfe,"PCA + BCA with FDR correction"=auc.PCA)
# write.table(auc.total,file = "/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/auc.final.table.txt",sep = "\t",row.names = F,quote = F)

# tiff(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/roc_plot_simulation.tiff", width=1200, height=1200, units="px", res=150)
# par(xpd = T,mar=c(5,5,4,100)+0.1)
# plot(roc.ml,main="ROC curves",xlim=c(1,0),col="#800000",lwd=6,cex.lab=2,cex.axis=2)
# plot(roc.ql,add=TRUE,col="#f58231",lwd=6)
# plot(roc.ttest,add=TRUE,col="#808000",lwd=6)
# plot(roc.mw,add=TRUE,col="#469990",lwd=6)
# # plot(roc.PCA,add=TRUE,col="#000075",lwd=6)
# plot(roc.ml.bonfe,add=T,col="#a9a9a9",lwd=6)
# plot(roc.mw.bonfe,add=T,col="#fabed4",lwd=6)
# plot(roc.ttest.bonfe,add=T,col="#f032e6",lwd=6)
# legend(x="right", inset = c(-0.4,-0.4),xpd = T,
#        cex=1,bty = 'n',text.width=0.1,ncol=1,xjust=0,lwd=4,
#        y.intersp =0,x.intersp = 0.1,legend = c("Maximum-likelihood and with\nobservational weights\nwith Benjamini-Honchberg correction","Quasi-likelihood and with\nobservational weights\nwith Benjamini-Honchberg\nand with Bonferroni correction","Moderated t test with\nBenjamini-Honchberg correction","Mann-whitney with\nBenjamini-Honchberg correction","PCA+BCA","Maximum-likelihood and\nwith observational weights\nwith Bonferroni correction","Mann-whitney with\nBonferroni correction","Moderated t test with\nBonferroni correction"), 
#        col = c("#800000","#f58231", "#808000", "#469990","#000075","#a9a9a9","#fabed4","#f032e6"), lty = c(1,1,1,1,1,1,1,1))

#legend(x=0.6,y=1,cex=1,bty = 'n',text.width=0.1,ncol=2,xjust=0,lwd=6,y.intersp =0,x.intersp = 0.1,legend = c("Maximum-likelihood and with\nobservational weights\nwith Benjamini-Honchberg correction","Quasi-likelihood and with\nobservational weights\nwith Benjamini-Honchberg\nand with Bonferroni correction","Moderated t test with\nBenjamini-Honchberg correction","Mann-whitney with\nBenjamini-Honchberg correction","PCA+BCA","Maximum-likelihood and\nwith observational weights\nwith Bonferroni correction","Mann-whitney with\nBonferroni correction","Moderated t test with\nBonferroni correction"), col = c("#800000","#f58231", "#808000", "#469990","#000075","#a9a9a9","#fabed4","#f032e6"), lty = c(1,1,1,1,1,1,1,1))
# plot(1, type = "n", xlab = "",
#      ylab = "", xlim = c(0, 5), 
#      ylim = c(0, 5))
# dev.off()
# without PCA+BCA #####
tiff(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/roc_plot_simulation_noPCABCA.tiff", width=1200, height=1200, units="px", res=150)
par(xpd = T,mar=c(5,5,4,100)+0.1)
plot(roc.ml,main="ROC curves",xlim=c(1,0),col="#800000",lwd=6,cex.lab=2,cex.axis=2)
plot(roc.ql,add=TRUE,col="#f58231",lwd=6)
plot(roc.ttest,add=TRUE,col="#808000",lwd=6)
plot(roc.mw,add=TRUE,col="#469990",lwd=6)
plot(roc.ml.bonfe,add=T,col="#a9a9a9",lwd=6)
plot(roc.mw.bonfe,add=T,col="#fabed4",lwd=6)
plot(roc.ttest.bonfe,add=T,col="#f032e6",lwd=6)
legend(x="right", inset = c(-0.4,-0.4),xpd = T,
       cex=1,bty = 'n',text.width=0.1,ncol=1,xjust=0,lwd=4,
       y.intersp =0,x.intersp = 0.1,legend = c("EdgeR - ML with FDR","EdgeR - QL with FDR","Limma with FDR","Mann-whitney with FDR","EdgeR - ML with Bonferroni","Mann-whitney with Bonferroni","Limma with Bonferroni"), 
       col = c("#800000","#f58231", "#808000", "#469990","#a9a9a9","#fabed4","#f032e6"), lty = c(1,1,1,1,1,1,1))
dev.off()
tiff(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/legend_roc_plot_simulationnoPCABCA.tiff", width=1200, height=1200, units="px", res=150)
par(xpd = T,mar=c(5,5,4,10)+0.1)
plot(1, type = "n", xlab = "",
     ylab = "", xlim = c(0, 0),
     ylim = c(0, 0))
legend(x="center", inset = c(0.5,0.5),xpd = T,
       cex=1,bty = 'n',text.width=0.5,ncol=1,xjust=0,lwd=4,
       y.intersp =3,x.intersp = 1,legend = c("EdgeR - ML with FDR","EdgeR - QL with FDR","Limma with FDR","Mann-whitney with FDR","EdgeR - ML with Bonferroni","Mann-whitney with Bonferroni","Limma with Bonferroni"), 
       col = c("#800000","#f58231", "#808000", "#469990","#a9a9a9","#fabed4","#f032e6"), lty = c(1,1,1,1,1,1,1))
dev.off()

# let's plot it in ggplot so we can see the all curves#####
ro_data=function(roc){
  data.frame(
    Specificity=rev(roc$specificities),
    Sensitivity=rev(roc$sensitivities),
    Model="name"
  )
}
roc_list=list(roc_ml=roc.ml,roc_ml_bonfe=roc.ml.bonfe,roc_ql=roc.ql,roc_ql_bonfe=roc.ql.bonfe,roc_mw=roc.mw,roc_mw_bonfe=roc.mw.bonfe,roc_ttest=roc.ttest,roc_ttest_bonfe=roc.ttest.bonfe)
roc_dfs=lapply(roc_list,ro_data)
roc_dfs=lapply(names(roc_list),function(name){
  roc_data
}
  )
roc_dfs=do.call(rbind,roc_dfs)


#########################################
tiff(file="/proj/naiss2024-23-57/GBS_MeDIP_benchmark/benchmark_final_tables/legend_roc_plot_simulation.tiff", width=1200, height=1200, units="px", res=150)
par(xpd = T,mar=c(5,5,4,10)+0.1)
plot(1, type = "n", xlab = "",
     ylab = "", xlim = c(0, 5),
     ylim = c(0, 5))
legend(x="center", inset = c(0,0),xpd = T,
       cex=1,bty = 'n',text.width=0.5,ncol=1,xjust=0.5,lwd=4,
       y.intersp =3,x.intersp = 1,legend = c("Maximum-likelihood\nwith FDR correction","Quasi-likelihood\nwith FDR correction","Moderated t-test with\nFDR correction","Mann-whitney with\nFDR correction","PCA+BCA","Maximum-likelihood\nwith Bonferroni correction","Mann-whitney with\nBonferroni correction","Moderated t-test with\nBonferroni correction"), 
       col = c("#800000","#f58231", "#808000", "#469990","#000075","#a9a9a9","#fabed4","#f032e6"), lty = c(1,1,1,1,1,1,1,1))
dev.off()

#venn diagram
#install.packages("VennDiagram")
library(VennDiagram)
overlapping=full_join(sign.PCA, MW.sig,by="wind_id")
nrow(overlapping)
sum(complete.cases(overlapping))
nrow(overlapping)
length(complete.cases(overlapping)==FALSE)
grid.newpage()
draw.pairwise.venn(area1=1907,area2=1757,cex=4,cat.cex=4,cat.dist=c(-0.07,-0.07),cross.area=1757,lwd=c(0.4,0.4),category=c("Mann-Whitney with\nBenjamini-Honchberg\ncorrection","PCA+BCA"),col=c("#42d4f4","#fabed4"),fill=c("#42d4f4","#fabed4"))
venn.pca.mw
sum(complete.cases(full_join(MW.sig, sign.PCA,by="wind_id")))
sum(complete.cases(full_join(MW.sig, MW.bonfe.sig,by="wind_id")))
sum(complete.cases(full_join(MW.bonfe.sig, sign.PCA,by="wind_id")))
o=full_join(MW.bonfe.sig, sign.PCA,by="wind_id")
o=full_join(o, sign.PCA,by="wind_id")
sum(complete.cases(o))
tiff(file="venn_diagram_simulation.tiff", width=1200, height=1200, units="px", res=1000)
grid.newpage()                                        
draw.triple.venn(1907,1757, 1698,1757,1698,1698, 1698,euler.d=T,scaled=T,cex=4,sep.dist=10,ext.text=T,ext.pos=c(2),ext.dist=3,cat.dist=c(-0.07,-0.07,-0.07),lwd=c(0.4,0.4,0.4),cat.cex=4,category=c("Mann-Whitney with\nBenjamini-Honchberg\ncorrection","PCA+BCA","Mann-Whitney with\nBonferroni\ncorrection"),col=c("#42d4f4","#fabed4","#aaffc3"),fill=c("#42d4f4","#fabed4","#aaffc3"),cat.col=c("#42d4f4","#fabed4","#808000"))
dev.off()