library(limma)
library(lumi)
library(gplots)
library(stringr)
library(qvalue)

RUV2 = function(Y, ctl, k, Z=matrix(rep(1, ncol(Y))),first=NULL,second=NULL,old=NULL)
{
    # Project onto the orthogonal complement of Z
    # DC - this just means regress out the intercept
    RZY = Y - Y%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)
    # Perform SVD
    W = svd(RZY[ctl,])$v
    # Keep the first k factors
    W = W[,1:k]
    # Fit using Limma and return
    mod1 = cbind(Z,W)
    gammahat = (Y %*% mod1 %*% solve(t(mod1) %*% mod1))[,2:dim(mod1)[2]]
    Yfit = Y - gammahat %*% t(W)
    Ysub1 = Yfit[,c(first,second)]
    X1 = rep(c(0,1),times=c(length(first),length(second)))
    fit1 = lmFit(Ysub1,cbind(X1,Z[1:dim(Ysub1)[2]]))
    fit1 = eBayes(fit1)
    fit1 = fit1$p.value[,1]
    fit1.q = qvalue(fit1)$qvalues
    fit1 = length(which(fit1.q < 0.05))
    
    Ysub2 = Yfit[,c(first,old)]
    X2 = rep(c(0,1),times=c(length(first),length(old)))
    fit2 = lmFit(Ysub2,cbind(X2,Z[1:dim(Ysub2)[2]]))
    fit2 = eBayes(fit2)
    fit2 = fit2$p.value[,1]
    fit2.q = qvalue(fit2)$qvalues
    fit2 = length(which(fit2.q < 0.05))
    return(list(fit1,fit2))
}


#Set working directory...
#setwd("E:/Research Efforts/TF Knockdowns/Knockdowns/RNAi/Expression Results/Kd_Arrays")

#Load annotation files...
classy = c('numeric',rep('factor',times=9),'numeric','factor',
           rep('numeric',times=2),rep('factor',times=5))
covariates = read.csv("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/CombiningInfo.csv",colClasses=classy)
probereport = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt",header=T)
detection = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/allthree_detection.txt")
allprobes = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_Stranded.txt",header=T)
targets = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/TargetSummary.txt")
tfs = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/TFcensus.txt",header=T,sep="\t",fill=NA)
#covariates = read.csv("~/home/Kd_Arrays/Analysis/RawExprs/CombiningInfo.csv",colClasses=classy)
#probereport = read.table("~/home/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt",header=T)
#detection = read.table("~/home/Kd_Arrays/Analysis/RawExprs/allthree_detection.txt")
#allprobes = read.table("~/home/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_Stranded.txt",header=T)
#targets = read.table("~/home/Kd_Arrays/Analysis/Annotations/TargetSummary.txt")
#tfs = read.table("~/home/Kd_Arrays/CombinedBinding/Annotations/TFcensus.txt",header=T,sep="\t",fill=NA)
tfs = tfs[,c(1,2,6)]
tfs = tfs[tfs[,1] != 'x',]
tfs = tfs[tfs[,1] != 'c',]
ilmns = probereport[match(targets[,3],probereport[,7]),4]
ilmns = cbind(as.character(targets[,3]),as.character(ilmns))

#Read in arrays, remove unexpressed probes...
expr_trans = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/allthree_trans_expr.txt")
#expr_trans = read.table("~/home/Kd_Arrays/Analysis/RawExprs/allthree_trans_expr.txt")
arrays = colnames(expr_trans)
arrays = str_sub(arrays,2,-1)
probes = rownames(expr_trans)
goodprobes = intersect(probereport$probeID,probes)
good.ind = which(probes %in% goodprobes)
seventwos = grep("72",covariates$X48or72)
covariates_seventwos = covariates[seventwos,]



expr_trans_good = expr_trans[good.ind,]
detection_good = detection[good.ind,]
expr_trans_seventwos = expr_trans_good[,seventwos]
detection_seventwos = detection_good[,seventwos]
byebye = c("X86_NS3_4_72","X156_NS2_1B_72")
byebye.ind = match(byebye,colnames(expr_trans_seventwos))
expr_trans_seventwo = expr_trans_seventwos[,-byebye.ind]
detection_seventwo = detection_seventwos[,-byebye.ind]
covariates_seventwo = covariates_seventwos[-byebye.ind,]
firstbatch = grep("10/17/2011",covariates_seventwo$NFDate)
secondbatch = grep("3/26/2012",covariates_seventwo$NFDate)
oldbatch = grep("3/29/2011",covariates_seventwo$NFDate)
good_probes = rownames(expr_trans_good)
detect.ind = which(rowSums(detection_seventwo<0.01) > 1)
expr_trans_best = expr_trans_seventwo[detect.ind,]
expr_quant_best = lumiN(as.matrix(expr_trans_best),method="quantile")
detection_best = detection_seventwo[detect.ind,]
expressed.ind = which(rowSums(detection_best<0.01) == dim(detection_best)[2])
ns = grep("NS",colnames(expr_quant_best))
ns.first = intersect(ns,firstbatch)
ns.second = intersect(ns,secondbatch)
ns.old = intersect(ns,oldbatch)

var1 = apply(expr_quant_best[expressed.ind,firstbatch],1,var)
var2 = apply(expr_quant_best[expressed.ind,secondbatch],1,var)
var3 = apply(expr_quant_best[expressed.ind,oldbatch],1,var)
means = apply(expr_quant_best[expressed.ind,],1,mean)
vars.ind = intersect(order(var1)[1:2000],order(var2)[1:2000])
var.ind = intersect(order(var3)[1:2000],vars.ind)
ctl = expressed.ind[var.ind]

dereport = matrix(0,50,2)
for(i in 1:50){
    testruv = RUV2(Y=as.matrix(expr_quant_best),ctl=ctl,k=i,first=ns.first,second=ns.second,old=ns.old)
    dereport[i,1] = testruv[[1]][1]
    dereport[i,2] = testruv[[2]][1]
}
plot(dereport[,1])
points(dereport[,2],pch=19)

k=8
Y = expr_quant_best
Z=matrix(rep(1, ncol(Y)))
RZY = Y - Y%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)
W = svd(RZY[ctl,])$v
W = W[,1:k]
mod1 = cbind(Z,W)
gammahat = (Y %*% mod1 %*% solve(t(mod1) %*% mod1))[,2:dim(mod1)[2]]
Yfit = Y - gammahat %*% t(W)

test = paste(colnames(Yfit)[133:201],"_Old",sep="")

splitter = function(expr,firstbatch=NULL,secondbatch=NULL,oldbatch=NULL,first.ind=NULL,second.ind=NULL,old.ind=NULL){
    expr_firstbatch = expr[,firstbatch]
    expr_secondbatch = expr[,secondbatch]
    expr_oldbatch = expr[,oldbatch]
    first = expr_firstbatch[first.ind,]
    second = expr_secondbatch[second.ind,]
    old = expr_oldbatch[old.ind,]
    
    return(list(first = first,second = second,old = old))
}

detection_firstbatch = detection_best[,firstbatch]
detection_secondbatch = detection_best[,secondbatch]
detection_oldbatch = detection_best[,oldbatch]
detect_firstbatch.ind = which(rowSums(detection_firstbatch<0.01) > 1)
detect_secondbatch.ind = which(rowSums(detection_secondbatch<0.01) > 1)
detect_oldbatch.ind = which(rowSums(detection_oldbatch<0.01) > 1)

Ysplit = splitter(expr=Yfit,firstbatch=firstbatch,secondbatch=secondbatch,oldbatch=oldbatch,first.ind=detect_firstbatch.ind,second.ind=detect_secondbatch.ind,old.ind=detect_oldbatch.ind)
Yquant = splitter(expr=expr_quant_best,firstbatch=firstbatch,secondbatch=secondbatch,oldbatch=oldbatch,first.ind=detect_firstbatch.ind,second.ind=detect_secondbatch.ind,old.ind=detect_oldbatch.ind)

nfdate.col = c("gray","red","blue")[unlist(covariates_seventwo$NFDate)]
heatmap.2(cor(expr_quant_best,method="spearman"),distfun=function(x) as.dist(1-abs(x)),trace="none",ColSideColors=nfdate.col)
heatmap.2(cor(Yfit,method="spearman"),distfun=function(x) as.dist(1-abs(x)),trace="none",ColSideColors=nfdate.col)

#write.table(Ysplit$first,"expr_firstbatch_ruv2.txt")
#rite.table(Ysplit$second,"expr_secondbatch_ruv2.txt")
#write.table(Ysplit$old,"expr_oldbatch_ruv2.txt")
#write.table(Yquant$old,"expr_oldbatch.txt")

RUV2.test = function(Y, ctl, k, Z=matrix(rep(1, ncol(Y))),irf=NULL,sp1=NULL,ns=NULL)
{
  # Project onto the orthogonal complement of Z
  # DC - this just means regress out the intercept
  RZY = Y - Y%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)
  # Perform SVD
  W = svd(RZY[ctl,])$v
  # Keep the first k factors
  W = W[,1:k]
  # Fit using Limma and return
  mod1 = cbind(Z,W)
  gammahat = (Y %*% mod1 %*% solve(t(mod1) %*% mod1))[,2:dim(mod1)[2]]
  Yfit = Y - gammahat %*% t(W)
  X = rep(c(0,1),times=c(3,6))
  Ysub.irf.old = Yfit[,c(irf[1:3],ns[1:6])]
  fit1 = lmFit(Ysub.irf.old,cbind(X,Z[1:dim(Ysub.irf.old)[2]]))
  fit1 = eBayes(fit1)
  fit1 = fit1$p.value[,1]
  fit1.q = qvalue(fit1)$qvalues
  fit1s = which(fit1.q < 0.05)
  
  Ysub.irf.new = Yfit[,c(irf[4:6],ns[7:12])]
  fit2 = lmFit(Ysub.irf.new,cbind(X,Z[1:dim(Ysub.irf.new)[2]]))
  fit2 = eBayes(fit2)
  fit2 = fit2$p.value[,1]
  fit2.q = qvalue(fit2)$qvalues
  fit2s = which(fit2.q < 0.05)
  
  irf.match = length(intersect(fit1s,fit2s))
  irf.cor = cor(-log10(fit1),-log10(fit2),use="pairwise.complete.obs")
  irf.cor.alt = cor(as.vector(Ysub.irf.old),as.vector(Ysub.irf.new),use="pairwise.complete.obs")
  
  Ysub.sp1.old = Yfit[,c(sp1[1:3],ns[1:6])]
  fit3 = lmFit(Ysub.sp1.old,cbind(X,Z[1:dim(Ysub.sp1.old)[2]]))
  fit3 = eBayes(fit3)
  fit3 = fit3$p.value[,1]
  fit3.q = qvalue(fit3)$qvalues
  fit3s = which(fit3.q < 0.05)
  
  Ysub.sp1.new = Yfit[,c(sp1[4:6],ns[13:18])]
  fit4 = lmFit(Ysub.sp1.new,cbind(X,Z[1:dim(Ysub.sp1.new)[2]]))
  fit4 = eBayes(fit4)
  fit4 = fit4$p.value[,1]
  fit4.q = qvalue(fit4)$qvalues
  fit4s = which(fit4.q < 0.05)
  
  sp1.match = length(intersect(fit3s,fit4s))
  sp1.cor = cor(-log10(fit3),-log10(fit4),use="pairwise.complete.obs")
  sp1.cor.alt = cor(as.vector(Ysub.sp1.old),as.vector(Ysub.sp1.new),use="pairwise.complete.obs")
  return(list(irf.match,sp1.match,irf.cor,sp1.cor,irf.cor.alt,sp1.cor.alt))
}

oldIRF = c("X68_IRF5_C","X110_IRF5_B","X97_IRF5_A")
newIRF = c("X31_IRF5_1","X62_IRF5_3","X142_IRF5_2")
oldSP1 = c("X46_SP1_C","X107_SP1_E","X19_SP1_D")
newSP1 = c("X26_SP1_2","X82_SP1_1","X110_SP1_3")
oldNS = c("X76_NS4_J","X88_NS2_E","X50_NS3_H","X36_NS_P3_1","X104_NS_P6_2","X1_NS1_D")
newIRFNS = c("X7_NS_P2_3_72","X44_NS_P1_2_72","X71_NS1_2_72","X107_NS2_3_72","X111_NS3_3_72","X127_NS4_2_72")
newSP1NS = c("X22_NS1_1A_72","X25_NS_P1_1B_72","X52_NS2_4_72","X75_NS_P2_4_72","X108_NS3_1A_72","X136_NS4_1B_72")

irf.ind = c(match(oldIRF,colnames(expr_quant_best)),match(newIRF,colnames(expr_quant_best)))
sp1.ind = c(match(oldSP1,colnames(expr_quant_best)),match(newSP1,colnames(expr_quant_best)))
nss.ind = c(match(oldNS,colnames(expr_quant_best)),match(newIRFNS,colnames(expr_quant_best)),match(newSP1NS,colnames(expr_quant_best)))

dereport = matrix(0,50,6)
for(i in 1:50){
  testruv = RUV2.test(Y=as.matrix(expr_quant_best),ctl=ctl,k=i,irf=irf.ind,sp1=sp1.ind,ns=nss.ind)
  #    dereport[i] = testruv
  dereport[i,1] = testruv[[1]][1]
  dereport[i,2] = testruv[[2]][1]
  dereport[i,3] = testruv[[3]][1]
  dereport[i,4] = testruv[[4]][1]
  dereport[i,5] = testruv[[5]][1]
  dereport[i,6] = testruv[[6]][1]
  print(i)
}

par(mfrow=c(2,2))
plot(dereport[,1],main = "DE O/L (Old vs New) at FDR 0.05",xlab="k",ylab="IRF5")
abline(v=10)
plot(dereport[,2],main = "DE O/L (Old vs New) at FDR 0.05",xlab="k",ylab="SP1")
abline(v=10)
plot(dereport[,3],main = "Cor(-Log10(P-values)) (Old vs New)",xlab="k",ylab="IRF5")
abline(v=10)
plot(dereport[,4],main = "Cor(-Log10(P-values)) (Old vs New)",xlab="k",ylab="SP1")
abline(v=10)
plot(dereport[,5],main = "Correlation (Old vs New)",xlab="k",ylab="IRF5")
abline(v=10)
plot(dereport[,6],main = "Correlation (Old vs New)",xlab="k",ylab="SP1")
abline(v=10)