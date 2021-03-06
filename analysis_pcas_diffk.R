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

RUV2 = function(Y,ctl,k){
  Z = matrix(rep(1, ncol(Y)))
  RZY = Y - Y%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)
  W = svd(RZY[ctl,])$v
  W = W[,1:k]
  mod1 = cbind(Z,W)
  gammahat = (Y %*% mod1 %*% solve(t(mod1) %*% mod1))[,2:dim(mod1)[2]]
  Yfit = Y - gammahat %*% t(W)
}

tfs = tfs[,c(1,2,6)]
tfs = tfs[tfs[,1] != 'x',]
tfs = tfs[tfs[,1] != 'c',]
ilmns = probereport[match(targets[,3],probereport[,7]),4]
ilmns = cbind(as.character(targets[,3]),as.character(ilmns))

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

var1 = apply(expr_quant_best[expressed.ind,firstbatch],1,var)
var2 = apply(expr_quant_best[expressed.ind,secondbatch],1,var)
var3 = apply(expr_quant_best[expressed.ind,oldbatch],1,var)
means = apply(expr_quant_best[expressed.ind,],1,mean)
vars.ind = intersect(order(var1)[1:2000],order(var2)[1:2000])
var.ind = intersect(order(var3)[1:2000],vars.ind)
ctl = expressed.ind[var.ind]

nfdate.col = c("indianred","dodgerblue2","goldenrod")[unlist(covariates_seventwo$NFDate)]
covariates_pc = covariatesgeo_seventwo[,-12]
pc.q = prcomp(t(expr_quant_best))

interestingpcs = c(3,4,5,9,12,10,11)
renew = c(1:23)[match(as.character(1:23),levels(covariates_pc[,3]))]
covariates_pc[,3] = reorder(covariates_pc[,3],renew[covariates_pc[,3]])
renew = c(1:12)[match(as.character(1:12),levels(covariates_pc[,4]))]
covariates_pc[,4] = reorder(covariates_pc[,4],renew[covariates_pc[,4]])
renew = c(1:3)[match(levels(covariates_pc[,12]),c("3/29/11","10/17/11","3/26/12"))]
covariates_pc[,12] = reorder(covariates_pc[,12],renew[covariates_pc[,12]])
for(z in 1:11){
  print(z)
  Yfit = RUV2(Y=expr_quant_best,ctl=ctl,k=z)
  pc.ruv = prcomp(t(Yfit))
  covtable = matrix(NA,10,5)
  for(i in 3:12){
    for(j in 1:5){
      currlm = lm(pc.ruv$x[,j] ~ covariates_pc[,i])
      covtable[(i-2),j] = pf(summary(currlm)$f[1],summary(currlm)$f[2],summary(currlm)$f[3],lower.tail=F)
    }
  }
  colnames(covtable) = paste("PC",1:5,sep="")
  rownames(covtable) = colnames(covariates_pc)[3:12]
  covtable = rbind(summary(pc.ruv)$importance[2,1:5],covtable)
  rownames(covtable)[1] = "PCImport"
  write.table(covtable,paste0(outdir,"k",z,"RUV2_PCA_Covariate_Correlations.txt"),sep="\t")
  pdf(paste0(outdir,"k",z,"PCs.pdf"))
  xer=0
  for(i in 1:5){
    for(j in 1:7){
      sig = ifelse(covtable[(interestingpcs[j]-1),i] < 0.05,"Significant","Not Significant")
      if(sig == "Significant"){
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
          print(stripplot(pc.ruv$x[,i] ~ covariates_pc[,interestingpcs[j]],
                          jitter.data = T,col=nfdate.col,cex=1.5,pch=20,main=list(sig),
                          xlab=list(label=names(covariates_pc)[interestingpcs[j]]),
                          scales = list(x = list(cex = 0.75)),
                          ylab=list(label=paste0("PC",i," - ",round(summary(pc.ruv)$importance[2,i]*100,3),"% of Variance"))),
                newpage=F)
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
  }
  dev.off()
}