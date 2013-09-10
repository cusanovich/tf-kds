library('stringr')
library('lumi')
expr_trans = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/allthree_trans_expr.txt")
classygeo = c('numeric',rep('factor',times=8),rep('numeric',times=2),
              rep('factor',times=2))
classy = c('numeric',rep('factor',times=9),'numeric','factor',
           rep('numeric',times=2),rep('factor',times=5))
covariates = read.csv("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/CombiningInfo.csv",colClasses=classy)
covariatesgeo = read.csv("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/CombiningInfoGEO.csv",colClasses=classygeo)
probereport = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt",
                         header=T)
detection = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/allthree_detection.txt")

RUV2 = function(Y,ctl,k){
  Z = matrix(rep(1, ncol(Y)))
  RZY = Y - Y%*%Z%*%solve(t(Z)%*%Z)%*%t(Z)
  W = svd(RZY[ctl,])$v
  W = W[,1:k]
  mod1 = cbind(Z,W)
  gammahat = (Y %*% mod1 %*% solve(t(mod1) %*% mod1))[,2:dim(mod1)[2]]
  Yfit = Y - gammahat %*% t(W)
}


arrays = colnames(expr_trans)
arrays = str_sub(arrays,2,-1)
probes = rownames(expr_trans)
goodprobes = intersect(probereport$probeID,probes)
good.ind = which(probes %in% goodprobes)
seventwos = grep("72",covariatesgeo$X48or72)
covariates_seventwos = covariates[seventwos,]
covariatesgeo_seventwos = covariatesgeo[seventwos,]

expr_trans_good = expr_trans[good.ind,]
detection_good = detection[good.ind,]
expr_trans_seventwos = expr_trans_good[,seventwos]
detection_seventwos = detection_good[,seventwos]
expr_trans_geos = expr_trans[,seventwos]
detection_geos = detection[,seventwos]
byebye = c("X86_NS3_4_72","X156_NS2_1B_72")
byebye.ind = match(byebye,colnames(expr_trans_seventwos))
expr_trans_seventwo = expr_trans_seventwos[,-byebye.ind]
detection_seventwo = detection_seventwos[,-byebye.ind]
covariates_seventwo = droplevels(covariates_seventwos[-byebye.ind,])
expr_trans_geo = expr_trans_geos[,-byebye.ind]
expr_trans_geo = round(2^expr_trans_geo,1)
detection_geo = detection_geos[,-byebye.ind]
detection_geo = round(detection_geo,5)
covariatesgeo_seventwo = droplevels(covariatesgeo_seventwos[-byebye.ind,])
colnames(expr_trans_geo) = as.character(covariatesgeo_seventwo[,2])
geoexprtable = matrix(NA,47170,403)
geoexprtable = cbind(expr_trans_geo[,1],detection_geo[,1])
collingnames = c("ID_REF",colnames(expr_trans_geo)[1],"Detection Pval")
for(i in 2:dim(expr_trans_geo)[2]){
  geoexprtable = cbind(geoexprtable,expr_trans_geo[,i])
  collingnames = c(collingnames,colnames(expr_trans_geo)[i])
  geoexprtable = cbind(geoexprtable,detection_geo[,i])
  collingnames = c(collingnames,"Detection Pval")
}
geoexprtable = cbind(rownames(expr_trans_geo),geoexprtable)
colnames(geoexprtable) = collingnames
write.table(geoexprtable,"/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/cusanovichGEOmatrixnon-normalized.txt",sep="\t",quote=F,row.names=F)

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

var1 = apply(expr_quant_best[expressed.ind,firstbatch],1,var)
var2 = apply(expr_quant_best[expressed.ind,secondbatch],1,var)
var3 = apply(expr_quant_best[expressed.ind,oldbatch],1,var)
means = apply(expr_quant_best[expressed.ind,],1,mean)
vars.ind = intersect(order(var1)[1:2000],order(var2)[1:2000])
var.ind = intersect(order(var3)[1:2000],vars.ind)
ctl = expressed.ind[var.ind]

Yfit = RUV2(Y=expr_quant_best,ctl=ctl,k=8)

exprgenes = rownames(detection_best)[which(rowSums(detection_best[,ns]<0.01) == length(ns))]
matchers = levels(covariatesgeo_seventwo[,8])
for(i in 1:length(matchers)){
  currexpr = detection_best[,grep(matchers[i],colnames(detection_best))]
  currgenes = rownames(currexpr)[which(rowSums(currexpr<0.01) == dim(currexpr)[2])]
  exprgenes = union(exprgenes,currgenes)
}

geonormtable = as.matrix(sort(exprgenes))
rownames(geonormtable) = sort(exprgenes)
normcolling = c("ID_REF")
normind = match(sort(exprgenes),rownames(Yfit))
for(i in 1:dim(Yfit)[2]){
  geonormtable = cbind(geonormtable,Yfit[normind,i])
  normcolling = c(normcolling,as.character(covariatesgeo_seventwo[i,2]))
  geonormtable = cbind(geonormtable,detection_best[normind,i])
  normcolling = c(normcolling,"Detection Pval")
}
colnames(geonormtable) = normcolling
write.table(geonormtable,"/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/cusanovichGEOmatrixnormalized.txt",sep="\t",quote=F,row.names=F)


metatable = as.matrix(covariatesgeo_seventwo[,2])
metatable = cbind(metatable,paste0("GM19238_Kd_",covariatesgeo_seventwo[,2]))
metatable = cbind(metatable,paste0("LCL GM19238, Kd ",covariatesgeo_seventwo[,8],", 72hrs"))
metatable = cbind(metatable,"Homo sapiens")
metatable = cbind(metatable,as.character(covariatesgeo_seventwo[,3]))
metatable = cbind(metatable,as.character(covariatesgeo_seventwo[,5]))
metatable = cbind(metatable,as.character(covariatesgeo_seventwo[,8]))
metatable = cbind(metatable,as.character(covariatesgeo_seventwo[,10]))
metatable = cbind(metatable,as.character(covariatesgeo_seventwo[,11]))
metatable = cbind(metatable,as.character(covariatesgeo_seventwo[,13]))
metatable = cbind(metatable,"RNA")
metatable = cbind(metatable,"Biotin")
metatable = cbind(metatable,"GPL10558")
colnames(metatable) = c("Sample name","title","source name","organism","characteristics: ChipNumber","characteristics: TransfectionPlate","characteristics: TargetGene","characteristics: KdLevel","characteristics: RIN","characteristics: TransfectionDate","molecule","label","platform")
write.table(metatable,"/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/RawExprs/Metatable.txt",sep="\t",quote=F,row.names=F)
