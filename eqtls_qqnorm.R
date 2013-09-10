classy = c(rep('character',times=2),rep('numeric',times=70))
noqq = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/gene_expression_levels_noqq.txt",header=T,colClasses=classy)
tss = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Centipede/Annotation/hg18_HT12ensemblTSScombinedsorted.bed")
#noqq = read.table("~/home/Kd_Arrays/Trans_eQTLs/Expression/gene_expression_levels_noqq.txt",header=T,colClasses=classy)
#tss = read.table("~/home/Kd_Arrays/Centipede/Annotation/hg18_HT12ensemblTSScombinedsorted.bed")

noqq.ind = match(tss$V4,noqq$gene)
tssupdate = tss[!is.na(noqq.ind),]
noqq.ind = noqq.ind[!is.na(noqq.ind)]
noqq.tss = noqq[noqq.ind,]

genelen = noqq.tss[,3]
noqq.rpkm = (noqq.tss[,4:72]/as.numeric(genelen))*1000000
expressed = rowSums(noqq.rpkm) > 0
noqq.genes = noqq.tss[expressed,1]
noqq.expr = noqq.rpkm[expressed,]
noqq.qq = apply(noqq.expr,2,function(x){qqnorm(x,plot.it=F)$x})
#pc.q = prcomp(t(noqq.rpkm))
#plot(pc.q$x,pch=20,cex=2,main="Quantile Normalization PCA",
#     xlab=paste0("PC1 - ",round(summary(pc.q)$importance[2,1]*100,3),"% of Variance"),
#     ylab=paste0("PC2 - ",round(summary(pc.q)$importance[2,2]*100,3),"% of Variance"),
#     las=1)
noqq.out = cbind(noqq.genes,noqq.qq)
colnames(noqq.out)[1] = "YRI"
colnames(noqq.out) = as.character(colnames(noqq.out))
write.table(t(noqq.out),"/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/gene_expression_qqnorm.txt",row.names=T,col.names=F,quote=F,sep="\t")
write.table(cov(noqq.qq),"/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Expression/gene_expression_cov.txt",row.names=F,col.names=F,quote=F,sep="\t")