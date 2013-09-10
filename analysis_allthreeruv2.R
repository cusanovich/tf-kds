library('limma')
library('lumi')
library('gplots')
library('stringr')
library('qvalue')
library('lattice')
library('grid')
library('gridBase')

#Load arrays, annotation files...
expr_trans = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/allthree_trans_expr.txt")
outdir = "/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/QCPlots/"
classygeo = c('numeric',rep('factor',times=8),rep('numeric',times=2),
              rep('factor',times=2))
classy = c('numeric',rep('factor',times=9),'numeric','factor',
           rep('numeric',times=2),rep('factor',times=5))
covariates = read.csv("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/CombiningInfo.csv",colClasses=classy)
covariatesgeo = read.csv("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/CombiningInfoGEO.csv",colClasses=classygeo)
probereport = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt",header=T)
detection = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/allthree_detection.txt")
targets = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/TargetSummary.txt")
tfs = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/TFcensus.txt",header=T,sep="\t",fill=NA)
#expr_trans = read.table("~/home/Kd_Arrays/Analysis/RawExprs/allthree_trans_expr.txt")
#covariates = read.csv("~/home/Kd_Arrays/Analysis/RawExprs/CombiningInfo.csv",colClasses=classy)
#probereport = read.table("~/home/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt",header=T)
#detection = read.table("~/home/Kd_Arrays/Analysis/RawExprs/allthree_detection.txt")
#targets = read.table("~/home/Kd_Arrays/Analysis/Annotations/TargetSummary.txt")
#tfs = read.table("~/home/Kd_Arrays/CombinedBinding/Annotations/TFcensus.txt",header=T,sep="\t",fill=NA)

#Define functions...
RUV2 = function(Y,ctl,k){
  Z = matrix(rep(1, ncol(Y)))
  # Regress out the intercept
  RZY = Y - Y%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)
  # Perform SVD
  W = svd(RZY[ctl,])$v
  # Keep the first k factors
  W = W[,1:k]
  # Regress out k factors
  mod1 = cbind(Z,W)
  gammahat = (Y %*% mod1 %*% solve(t(mod1) %*% mod1))[,2:dim(mod1)[2]]
  Yfit = Y - gammahat %*% t(W)
}

splitter = function(expr,firstbatch=NULL,secondbatch=NULL,oldbatch=NULL,
                    first.ind=NULL,second.ind=NULL,old.ind=NULL){
  expr_firstbatch = expr[,firstbatch]
  expr_secondbatch = expr[,secondbatch]
  expr_oldbatch = expr[,oldbatch]
  first = expr_firstbatch[first.ind,]
  second = expr_secondbatch[second.ind,]
  old = expr_oldbatch[old.ind,]
  return(list(first = first,second = second,old = old))
}

ruvfit = function(Yfit,group1=NULL,group2=NULL){
    Ysub = Yfit[,c(group1,group2)]
    X = rep(c(0,1),times=c(length(group1),length(group2)))
    Z = matrix(rep(1, ncol(Yfit)))
    fit = lmFit(Ysub,cbind(X,Z[1:dim(Ysub)[2]]))
    fit.ebayes = eBayes(fit)
    fit.p = fit.ebayes$p.value[,1]
    fit.q = qvalue(fit.p)$qvalues
    fit.count = length(which(fit.q < 0.05))
    return(list(p=fit.p,q=fit.q,count=fit.count))
}

ruvtest = function(Y,ctl,k,group1=NULL,group2=NULL,group3=NULL,nsgroup1=NULL,
                   nsgroup2=NULL,nsgroup3=NULL){
  Yfit = RUV2(Y=Y,ctl=ctl,k=k)
  if(k == 0){
    Yfit = Y
  }
  fit1 = ruvfit(Yfit,group1=group1,group2=group2)
  fit2 = ruvfit(Yfit,group1=group1,group2=group3)
  fit3 = ruvfit(Yfit,group1=group2,group2=group3)
  fit4 = ruvfit(Yfit,group1=nsgroup1,group2=nsgroup2)
  fit5 = ruvfit(Yfit,group1=nsgroup1,group2=nsgroup3)
  fit6 = ruvfit(Yfit,group1=nsgroup2,group2=nsgroup3)
  return(list(fit1,fit2,fit3,fit4,fit5,fit6))
}

ruvtest2 = function(Y, ctl, k,irf=NULL,sp1=NULL,ns=NULL){
  Yfit = RUV2(Y=Y,ctl=ctl,k=k)
  if(k == 0){
    Yfit = Y
  }
  fit1 = ruvfit(Yfit,group1=irf[1:3],group2=ns[1:6])
  fit2 = ruvfit(Yfit,group1=irf[4:6],group2=ns[7:12])
  fit3 = ruvfit(Yfit,group1=sp1[1:3],group2=ns[1:6])
  fit4 = ruvfit(Yfit,group1=sp1[4:6],group2=ns[13:18])
  irf.match = sum(fit1$q < 0.05 & fit2$q < 0.05)
  irf.cor = cor(-log10(fit1$p),-log10(fit2$p),use="pairwise.complete.obs")
  irf.cor.alt = cor(as.vector(Yfit[,c(irf[1:3],ns[1:6])]),
                    as.vector(Yfit[,c(irf[4:6],ns[7:12])]),use="pairwise.complete.obs")
  sp1.match = sum(fit3$q < 0.05 & fit4$q < 0.05)
  sp1.cor = cor(-log10(fit3$p),-log10(fit4$p),use="pairwise.complete.obs")
  sp1.cor.alt = cor(as.vector(Yfit[,c(sp1[1:3],ns[1:6])]),
                    as.vector(Yfit[,c(sp1[4:6],ns[13:18])]),use="pairwise.complete.obs")
  return(list(irf.match,sp1.match,irf.cor,sp1.cor,irf.cor.alt,sp1.cor.alt))
}

#RLE plot maker
rlemaker = function(expr){
  exprmeds = apply(expr,1,median)
  exprrle = expr - exprmeds
  return(exprrle)
}


tfs = tfs[,c(1,2,6)]
tfs = tfs[tfs[,1] != 'x',]
tfs = tfs[tfs[,1] != 'c',]
ilmns = probereport[match(targets[,3],probereport[,7]),4]
ilmns = cbind(as.character(targets[,3]),as.character(ilmns))

#Remove unexpressed probes...
arrays = colnames(expr_trans)
arrays = str_sub(arrays,2,-1)
probes = rownames(expr_trans)
goodprobes = intersect(probereport$probeID,probes)
good.ind = which(probes %in% goodprobes)
seventwos = grep("72",covariates$X48or72)
covariates_seventwos = covariates[seventwos,]
covariatesgeo_seventwos = covariatesgeo[seventwos,]

expr_trans_good = expr_trans[good.ind,]
detection_good = detection[good.ind,]
expr_trans_seventwos = expr_trans_good[,seventwos]
detection_seventwos = detection_good[,seventwos]
byebye = c("X86_NS3_4_72","X156_NS2_1B_72")
byebye.ind = match(byebye,colnames(expr_trans_seventwos))
expr_trans_seventwo = expr_trans_seventwos[,-byebye.ind]
detection_seventwo = detection_seventwos[,-byebye.ind]
covariates_seventwo = covariates_seventwos[-byebye.ind,]
covariatesgeo_seventwo = covariatesgeo_seventwos[-byebye.ind,]
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

detection_firstbatch = detection_best[,firstbatch]
detection_secondbatch = detection_best[,secondbatch]
detection_oldbatch = detection_best[,oldbatch]
detect_firstbatch.ind = which(rowSums(detection_firstbatch<0.01) > 1)
detect_secondbatch.ind = which(rowSums(detection_secondbatch<0.01) > 1)
detect_oldbatch.ind = which(rowSums(detection_oldbatch<0.01) > 1)

Yfit = RUV2(Y=expr_quant_best,ctl=ctl,k=8)
Ysplit = splitter(expr=Yfit,firstbatch=firstbatch,secondbatch=secondbatch,
                  oldbatch=oldbatch,first.ind=detect_firstbatch.ind,
                  second.ind=detect_secondbatch.ind,old.ind=detect_oldbatch.ind)
Yquant = splitter(expr=expr_quant_best,firstbatch=firstbatch,
                  secondbatch=secondbatch,oldbatch=oldbatch,
                  first.ind=detect_firstbatch.ind,
                  second.ind=detect_secondbatch.ind,old.ind=detect_oldbatch.ind)
fit_rle = rlemaker(Yfit)
quant_rle = rlemaker(expr_quant_best)

#write.table(Ysplit$first,"expr_firstbatch_ruv2.txt")
#rite.table(Ysplit$second,"expr_secondbatch_ruv2.txt")
#write.table(Ysplit$old,"expr_oldbatch_ruv2.txt")
#write.table(Yquant$old,"expr_oldbatch.txt")

nfdate.col = c("indianred","dodgerblue2","goldenrod")[unlist(covariates_seventwo$NFDate)]
#covariates_pc = covariates_seventwo[,-15]
covariates_pc = covariatesgeo_seventwo[,-12]
pc.q = prcomp(t(expr_quant_best))
covtable = matrix(NA,10,5)
for(i in 3:12){
  for(j in 1:5){
#    print(i)
#    print(j)
    currlm = lm(pc.q$x[,j] ~ covariates_pc[,i])
    covtable[(i-2),j] = pf(summary(currlm)$f[1],summary(currlm)$f[2],summary(currlm)$f[3],lower.tail=F)
  }
}
colnames(covtable) = paste("PC",1:5,sep="")
rownames(covtable) = colnames(covariates_pc)[3:12]
covtable = rbind(summary(pc.q)$importance[2,1:5],covtable)
rownames(covtable)[1] = "PCImport"
write.table(covtable,paste0(outdir,"Quantile_PCA_Covariate_Correlations.txt"),sep="\t")

pc.ruv = prcomp(t(Yfit))
covtable = matrix(NA,10,5)
for(i in 3:12){
  for(j in 1:5){
#    print(i)
#    print(j)
    currlm = lm(pc.ruv$x[,j] ~ covariates_pc[,i])
    covtable[(i-2),j] = pf(summary(currlm)$f[1],summary(currlm)$f[2],summary(currlm)$f[3],lower.tail=F)
  }
}
colnames(covtable) = paste("PC",1:5,sep="")
rownames(covtable) = colnames(covariates_pc)[3:12]
covtable = rbind(summary(pc.ruv)$importance[2,1:5],covtable)
rownames(covtable)[1] = "PCImport"
write.table(covtable,paste0(outdir,"RUV2_PCA_Covariate_Correlations.txt"),sep="\t")

dereport = matrix(0,26,6)
for(i in 0:25){
    testerruv = ruvtest(Y=as.matrix(expr_quant_best),ctl=ctl,k=i,
                        group1=secondbatch,group2=firstbatch,group3=oldbatch,
                        nsgroup1=ns.second,nsgroup2=ns.first,nsgroup3=ns.old)
    dereport[i+1,1] = testerruv[[1]]$count
    dereport[i+1,2] = testerruv[[2]]$count
    dereport[i+1,3] = testerruv[[3]]$count
    dereport[i+1,4] = testerruv[[4]]$count
    dereport[i+1,5] = testerruv[[5]]$count
    dereport[i+1,6] = testerruv[[6]]$count
    print(i)
}

oldIRF = c("X68_IRF5_C","X110_IRF5_B","X97_IRF5_A")
newIRF = c("X31_IRF5_1","X62_IRF5_3","X142_IRF5_2")
oldSP1 = c("X46_SP1_C","X107_SP1_E","X19_SP1_D")
newSP1 = c("X26_SP1_2","X82_SP1_1","X110_SP1_3")
oldNS = c("X76_NS4_J","X88_NS2_E","X50_NS3_H","X36_NS_P3_1","X104_NS_P6_2",
          "X1_NS1_D")
newIRFNS = c("X7_NS_P2_3_72","X44_NS_P1_2_72","X71_NS1_2_72","X107_NS2_3_72",
             "X111_NS3_3_72","X127_NS4_2_72")
newSP1NS = c("X22_NS1_1A_72","X25_NS_P1_1B_72","X52_NS2_4_72","X75_NS_P2_4_72",
             "X108_NS3_1A_72","X136_NS4_1B_72")

irf.ind = c(match(oldIRF,colnames(expr_quant_best)),
            match(newIRF,colnames(expr_quant_best)))
sp1.ind = c(match(oldSP1,colnames(expr_quant_best)),
            match(newSP1,colnames(expr_quant_best)))
nss.ind = c(match(oldNS,colnames(expr_quant_best)),
            match(newIRFNS,colnames(expr_quant_best)),
            match(newSP1NS,colnames(expr_quant_best)))

dereport2 = matrix(0,26,6)
for(i in 0:25){
  testerruv2 = ruvtest2(Y=as.matrix(expr_quant_best),ctl=ctl,k=i,irf=irf.ind,
                        sp1=sp1.ind,ns=nss.ind)
  dereport2[i+1,1] = testerruv2[[1]][1]
  dereport2[i+1,2] = testerruv2[[2]][1]
  dereport2[i+1,3] = testerruv2[[3]][1]
  dereport2[i+1,4] = testerruv2[[4]][1]
  dereport2[i+1,5] = testerruv2[[5]][1]
  dereport2[i+1,6] = testerruv2[[6]][1]
  print(i)
}
mypalette = colorRampPalette(c("cornflowerblue", "white", "darkblue"), space = "Lab")


pdf(paste0(outdir,"FigS1alt.pdf"))
heatmap.2(cor(expr_quant_best,method="spearman"),
          distfun=function(x) as.dist(1-abs(x)),trace="none",
          ColSideColors=nfdate.col,col=mypalette(20),main="Quantile Normalized",labRow="",labCol="")
dev.off()
pdf(paste0(outdir,"FigS3alt.pdf"))
heatmap.2(cor(Yfit,method="spearman"),
          distfun=function(x) as.dist(1-abs(x)),trace="none",
          ColSideColors=nfdate.col,col=mypalette(20),main="Quantile+RUV-2",labRow="",labCol="")
dev.off()

pdf(paste0(outdir,"FigSPCs.pdf"))
interestingpcs = c(3,4,5,9,12,10,11)
renew = c(1:23)[match(as.character(1:23),levels(covariates_pc[,3]))]
covariates_pc[,3] = reorder(covariates_pc[,3],renew[covariates_pc[,3]])
renew = c(1:12)[match(as.character(1:12),levels(covariates_pc[,4]))]
covariates_pc[,4] = reorder(covariates_pc[,4],renew[covariates_pc[,4]])
renew = c(1:3)[match(levels(covariates_pc[,12]),c("3/29/11","10/17/11","3/26/12"))]
covariates_pc[,12] = reorder(covariates_pc[,12],renew[covariates_pc[,12]])
xer=0
for(i in 1:5){
  for(j in 1:7){
    sig = ifelse(covtable[(interestingpcs[j]-1),i] < 0.05,"Significant","Not Significant")
    if(xer%%4 == 0){
      plot.new()
      gl <- grid.layout(nrow=2, ncol=2)
      vp <- viewport(layout.pos.col=1, layout.pos.row=1)
      pushViewport(viewport(layout=gl))
      grid.text(sig, vp = viewport(layout.pos.row = 1, layout.pos.col = 1),gp=gpar(col="white"))
    }
    if(xer%%4 == 1){
      vp <- viewport(layout.pos.col=2, layout.pos.row=1)
    }
    if(xer%%4 == 2){
      vp <- viewport(layout.pos.col=1, layout.pos.row=2)
    }
    if(xer%%4 == 3){
      vp <- viewport(layout.pos.col=2, layout.pos.row=2)
    }
    pushViewport(vp)
        if(j < 6){
#      stripchart(pc.ruv$x[,i] ~ covariates_pc[,interestingpcs[j]],vertical=T,
#                 method="jitter",col=nfdate.col,cex=1.5,pch=20,
#                 ylab=paste0("PC",i," - ",round(summary(pc.ruv)$importance[2,i]*100,3),"% of Variance"))
      print(stripplot(pc.ruv$x[,i] ~ covariates_pc[,interestingpcs[j]],
                      jitter.data = T,col=nfdate.col,cex=1.5,pch=20,main=list(sig),
                      xlab=list(label=names(covariates_pc)[interestingpcs[j]]),
                      scales = list(x = list(cex = 0.75)),
                      ylab=list(label=paste0("PC",i," - ",round(summary(pc.ruv)$importance[2,i]*100,3),"% of Variance"))),
            newpage=F)
      if(xer%%4 == 0){
        
      }
    }else{
      par(new=TRUE, fig=gridFIG())
      plot(covariates_pc[,interestingpcs[j]],pc.ruv$x[,i],pch=20,cex=1.5,
           col=nfdate.col,
           ylab=paste0("PC",i," - ",round(summary(pc.ruv)$importance[2,i]*100,3),"% of Variance"),
           xlab=colnames(covariates_pc)[interestingpcs[j]],las=1,main=sig)
      abline(lm(pc.ruv$x[,i] ~ covariates_pc[,interestingpcs[j]]),col="mediumorchid3",lty=2,lwd=3)
    }
    popViewport()
    xer = xer + 1
  }
}
dev.off()

pdf(paste0(outdir,"FigS4.pdf"))
par(oma=c(2, 2, 2, 2) + 0.1)
par(mar=c(2, 2, 2, 2) + 0.1)
par(mfrow=c(2,1))
boxplot(quant_rle,col=nfdate.col,outline=F,ylim=c(-0.4,0.4),
        main="A. Quantile Normalized RLEs",xaxt="n",las=1)
abline(h=.2,lwd=2,lty="dashed",col="mediumorchid3")
abline(h=-.2,lwd=2,lty="dashed",col="mediumorchid3")
abline(h=0,lwd=2)
boxplot(quant_rle,col=nfdate.col,outline=F,add=T,xaxt="n",yaxt="n")
legend("bottom",legend=c("Batch 1","Batch 2", "Batch 3"),
       fill= c("indianred","dodgerblue2","goldenrod"),horiz=T,bg="white")
boxplot(fit_rle,col=nfdate.col,outline=F,ylim=c(-0.4,0.4),
        main="B. Quantile Normalized + RUV-2 RLEs",xaxt="n",las=1)
abline(h=.2,lwd=2,lty="dashed",col="mediumorchid3")
abline(h=-.2,lwd=2,lty="dashed",col="mediumorchid3")
abline(h=0,lwd=2)
boxplot(fit_rle,col=nfdate.col,outline=F,add=T,xaxt="n",yaxt="n")
legend("bottom",legend=c("Batch 1","Batch 2", "Batch 3"),
       fill= c("indianred","dodgerblue2","goldenrod"),horiz=T,bg="white")
dev.off()

pdf(paste0(outdir,"FigS2.pdf"))
par(mfrow=c(3,2))
plot(0:25,dereport[,1],col="indianred",type="b",lwd=2,xlim=c(0,10),
     ylim=c(0,max(dereport[,1:3])),main="A. Genes DE Between Batches vs. RUV-2 k",
     las=1,xlab="k",ylab="No. Genes DE")
lines(0:25,dereport[,2],col="dodgerblue2",lwd=2,type="b")
lines(0:25,dereport[,3],col="goldenrod",lwd=2,type="b")
legend("topright",legend=c("1v2","1v3","2v3"),
       fill=c("indianred","dodgerblue2","goldenrod"))
plot(0:25,dereport[,4],col="indianred",type="b",lwd=2,xlim=c(0,10),
     ylim=c(0,max(dereport[,1:3])),
     main="A. Genes DE Between NS Batches vs. RUV-2 k",las=1,xlab="k",
     ylab="No. Genes DE")
lines(0:25,dereport[,5],col="dodgerblue2",lwd=2,type="b")
lines(0:25,dereport[,6],col="goldenrod",lwd=2,type="b")
legend("topright",legend=c("1v2","1v3","2v3"),
       fill=c("indianred","dodgerblue2","goldenrod"))
plot(0:25,dereport2[,1],type="b",col="dodgerblue2",
     main = "C. IRF5 Common DE Genes at FDR 0.05",xlab="k",ylab="No. Genes DE",
     las=1)
abline(v=8,col="indianred",lty="dashed")
plot(0:25,dereport2[,2],type="b",col="dodgerblue2",
     main = "D. SP1 Common DE Genes at FDR 0.05",xlab="k",ylab="No. Genes DE",
     las=1)
abline(v=8,col="indianred",lty="dashed")
plot(0:25,dereport2[,3],type="b",col="dodgerblue2",
     main = "E. IRF5 Correlation of -Log10(P-values)",xlab="k",
     ylab="Pearson Correlation",las=1)
abline(v=8,col="indianred",lty="dashed")
plot(0:25,dereport2[,4],type="b",col="dodgerblue2",
     main = "F. SP1 Correlation of -Log10(P-values)",xlab="k",
     ylab="Pearson Correlation",las=1)
abline(v=8,col="indianred",lty="dashed")
#plot(dereport2[,5],type="b",col="dodgerblue2",main = "Correlation (Old vs New)",xlab="k",ylab="IRF5")
#abline(v=8,col="indianred",lty="dashed")
#plot(dereport2[,6],type="b",col="dodgerblue2",main = "Correlation (Old vs New)",xlab="k",ylab="SP1")
#abline(v=8,col="indianred",lty="dashed")
dev.off()

pdf(paste0(outdir,"FigS5.pdf"),height=5,width=8.5)
par(mfrow=c(1,2))
plot(pc.q$x,pch=20,cex=2,col=nfdate.col,main="A. Quantile Normalization PCA",
     xlab=paste0("PC1 - ",round(summary(pc.q)$importance[2,1]*100,3),"% of Variance"),
     ylab=paste0("PC2 - ",round(summary(pc.q)$importance[2,2]*100,3),"% of Variance"),
     las=1)
plot(pc.ruv$x,pch=20,cex=2,col=nfdate.col,main="B. Quantile+RUV-2 Normalization PCA",
     xlab=paste0("PC1 - ",round(summary(pc.ruv)$importance[2,1]*100,3),"% of Variance"),
     ylab=paste0("PC2 - ",round(summary(pc.ruv)$importance[2,2]*100,3),"% of Variance"),
     las=1)
dev.off()