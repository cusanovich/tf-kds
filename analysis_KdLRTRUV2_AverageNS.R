library('plyr')
library('qvalue')
library('topGO')
library('gplots')
Args <- commandArgs(TRUE)
#Args = c("EP300","expr_secondbatch_a3ruv2k8")
#Args = c("HCST","expr_firstbatch_ruv2")
print(Args[1])
print(Args[2])
resultsbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/"

print('Loading expression file...')
expr = as.matrix(read.table(paste("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Exprs/",Args[2],".txt",sep="")))
#expr = as.matrix(read.table(paste("~/home/Kd_Arrays/Analysis/Exprs/",Args[2],".txt",sep="")))


print('Loading NS arrays...')
batches = c("oldbatch","firstbatch","secondbatch")
currbatch = strsplit(Args[2],"_")[[1]][2]
batches = batches[-match(currbatch,batches)]
nsa = expr[,grep("NS",colnames(expr))]
nsb = as.matrix(read.table(paste("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Exprs/",gsub(currbatch,batches[1],Args[2]),".txt",sep="")))
nsc = as.matrix(read.table(paste("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Exprs/",gsub(currbatch,batches[2],Args[2]),".txt",sep="")))
#nsb = as.matrix(read.table(paste("~/home/Kd_Arrays/Analysis/Exprs/",gsub(currbatch,batches[1],Args[2]),".txt",sep="")))
#nsc = as.matrix(read.table(paste("~/home/Kd_Arrays/Analysis/Exprs/",gsub(currbatch,batches[2],Args[2]),".txt",sep="")))
nsb = nsb[,grep("NS",colnames(nsb))]
nsc = nsc[,grep("NS",colnames(nsc))]
nsclass = list("NS1","NS2","NS3","NS4",c("NS_P1","NS_P3"),c("NS_P2","NS_P6"))
commongenes = intersect(rownames(nsa),rownames(nsb))
commongenes = intersect(commongenes,rownames(nsc))
for(i in 1:6){
  currns = nsclass[[i]]
  if(length(currns)<2){
    nsm = cbind(nsa[match(commongenes,rownames(nsa)),grep(currns,colnames(nsa))],nsb[match(commongenes,rownames(nsb)),grep(currns,colnames(nsb))])
    nsm = cbind(nsm,nsc[match(commongenes,rownames(nsc)),grep(currns,colnames(nsc))])
    nsave = apply(nsm,1,mean)
  }else{
    nsind = c()
    nsind[1] = ifelse(length(grep(currns[1],colnames(nsa))) > 0,grep(currns[1],colnames(nsa)),grep(currns[2],colnames(nsa)))
    nsind[2] = ifelse(length(grep(currns[1],colnames(nsb))) > 0,grep(currns[1],colnames(nsb)),grep(currns[2],colnames(nsb)))
    nsind[3] = ifelse(length(grep(currns[1],colnames(nsc))) > 0,grep(currns[1],colnames(nsc)),grep(currns[2],colnames(nsc)))
    nsm = cbind(nsa[match(commongenes,rownames(nsa)),nsind[1]],nsb[match(commongenes,rownames(nsb)),nsind[2]])
    nsm = cbind(nsm,nsc[match(commongenes,rownames(nsc)),nsind[3]])
    nsave = apply(nsm,1,mean)
  }
  if(i==1){
    ns_expr = as.matrix(nsave)
  }else{
    ns_expr = cbind(ns_expr,nsave)
  }
}
colnames(ns_expr) = c("NS1","NS2","NS3","NS4","NS_P1","NS_P2")

print('Loading detection file...')
detection = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Exprs/Detection_Scores_All3.txt")
#detection = read.table("~/home/Kd_Arrays/Analysis/Exprs/Detection_Scores_All3.txt")


print('Cleaning up a bit...')
nsout.ind = grep("NS",colnames(expr))
nsprobe.ind = match(rownames(ns_expr),rownames(expr))
expr = expr[,-nsout.ind]
expr = cbind(expr[nsprobe.ind,],ns_expr)
tfs.ind = grep(Args[1],colnames(expr))
nss.ind = grep("NS",colnames(expr))
svastable = expr[,c(nss.ind,tfs.ind)]
probings = as.factor(rownames(expr))
detect.ind = match(probings,rownames(detection))
detcol.ind = c(match(colnames(expr)[tfs.ind],colnames(detection)),grep("NS",colnames(detection)))
detecters = detection[detect.ind,detcol.ind]
detgood = which(rowSums(detecters<0.01) > 1)
det.ns = which(rowSums(detecters[,4:dim(detecters)[2]] < 0.01) == 18)
det.tf = which(rowSums(detecters[,1:3] < 0.01) == 3)
detgood = union(det.ns,det.tf)
svatable = svastable[detgood,]
probing = droplevels(probings[detgood])

probereport = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt",header=T)
#probereport = read.table("~/home/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt",header=T)
probereport.ind = match(probing,probereport[,4])
probereportupdate = probereport[probereport.ind,]
probereportupdate = probereportupdate[!is.na(probereport.ind),]
genes = probereportupdate[,7:8]
probes.ind = match(probereportupdate[,4],probing)
probing = droplevels(probing[probes.ind])

expr_batch = svatable[probes.ind,]

print('Building dataframe...')
ns.ind = grep("NS",colnames(expr_batch))
nstable = data.frame()
for(i in 1:length(ns.ind)){
  temp = cbind(probing,data.frame(expr_batch[,ns.ind[i]]))
  names(temp) = c("V1","V2")
  temp$V3 = as.factor("NS")
  temp$V4 = colnames(expr_batch)[ns.ind[i]]
  nstable = rbind(nstable,temp)
}
dim(nstable)


tf.ind = grep(Args[1],colnames(expr_batch))
tftable = data.frame()
for(j in 1:length(tf.ind)){
  temp = cbind(probing,data.frame(expr_batch[,tf.ind[j]]))
  names(temp) = c("V1","V2")
  temp$V3 = as.factor(Args[1])
  temp$V4 = colnames(expr_batch)[tf.ind[j]]
  tftable = rbind(tftable,temp)
}
tftable = rbind(nstable,tftable)
rm(temp)
dim(tftable)

print('Running models...')
lnmodels = dlply(tftable, "V1", function(lmfit) lm(lmfit$V2 ~ lmfit$V3))
lnresults = ldply(lnmodels, logLik)
rm(lnmodels)
nullmodels = dlply(tftable, "V1", function(lmfit) lm(lmfit$V2 ~ 1))
nullresults = ldply(nullmodels, logLik)
rm(nullmodels)
probeorder = levels(probing)
probe.ind = match(probeorder,probing)

print('Calculating P-values...')
lrt = 2*(lnresults[,1]-nullresults[,1])
plrt = pchisq(lrt,df=1,lower.tail=FALSE)
qlrt = qvalue(plrt)$qvalues
fivepercers = which(qlrt <= 0.05)
fiveperc = length(fivepercers)
fivecol = rep("black",times=length(qlrt))
fivecol[fivepercers] = "goldenrod"
fivecex = rep(0.5,times=length(qlrt))
fivecex[fivepercers] = 1

# Calculate stats for MA plot and Volcano plot
meanns = apply(expr_batch[,ns.ind],1,mean)
meantf = apply(expr_batch[,tf.ind],1,mean)
a = apply(expr_batch[,c(ns.ind,tf.ind)],1,mean)
m = meantf - meanns
meanns = meanns[probe.ind]
meantf = meantf[probe.ind]
a = a[probe.ind]
m = m[probe.ind]
tfprobes = which(genes[probe.ind,2] == Args[1])

# Calculate stats for QQplot
obsp <- -log10(plrt)
N <- length(obsp)
null <- -log10(ppoints(N))
MAX <- max(c(obsp,null))
c95 <- rep(0,N)
c05 <- rep(0,N)
for(i in 1:N){
  c95[i] <- qbeta(0.95,i,N-i+1)
  c05[i] <- qbeta(0.05,i,N-i+1)
}
xx <- c(sort(null),sort(null,decreasing=T))
yy <- c(sort(-log10(c95)),sort(-log10(c05),decreasing=T))


print('Making PDFs...')
# Print Figures to PDF
pdf(paste(resultsbin,Args[1],"_",strsplit(Args[2],split="_")[[1]][2],"_Results.pdf",sep=""))
par(mfrow=c(2,2))

# Histogram of P-values
hist(plrt,main=paste(Args[1],"\nNo. Probes = ",N),xlab="P-values")

# QQPlot
plot(1,1, type="n", xlim=c(0,max(null)+1),ylim=c(0,max(obsp)+1),xlab="Expected", ylab="Observed",main=Args[1])
polygon(xx, yy, col="gray",border="gray")
lines(null,null,col="dodgerblue2",lwd=2)
points(sort(null),sort(obsp),pch=20,cex=0.5)
if(length(tfprobes) > 0){
  for(i in 1:length(tfprobes)){
    points(sort(null)[which(sort(obsp) == obsp[tfprobes[i]])],obsp[tfprobes[i]],
           pch=20,cex=1.2,col="indianred")
  }
}

# MA plot
plot(a,m,pch=20,ylim=c(min(-1,min(m)),max(1,max(m))),cex=fivecex,main=Args[1],
     xlab="A",ylab="M",col=fivecol)
if(length(tfprobes) > 0){
  points(a[tfprobes],m[tfprobes],pch=20,col="indianred",cex=1.5)
}
abline(h=0,lty="dashed",lwd=2,col="dodgerblue2")

# Volcano plot
plot(m,-log10(plrt),xlim=c(min(-1,min(m)),max(1,max(m))),cex=0.5,pch=20,
     main=paste(Args[1],"\n5% FDR = ",fiveperc," probes",sep=""),
     xlab="Log2(Fold-Change)",ylab="-Log10(P-value)")
if(length(tfprobes) > 0){
  points(m[tfprobes],-log10(plrt[tfprobes]),pch=20,col="indianred",cex=1.5)
}
if(fiveperc > 0){
#  abline(h=-log10(plrt[which(qlrt == max(qlrt[qlrt <= 0.01]))[1]]),lty="dashed",lwd=2,col="blue")
  abline(h=-log10(plrt[which(qlrt == max(qlrt[qlrt <= 0.05]))[1]]),lty="dashed",
         lwd=2,col="dodgerblue2")
}

par(mfrow=c(1,1))
heatmap.2(cor(svatable,method="spearman"),main=paste(Args[1]," QN",sep=""),
          trace="none",distfun=function(x) as.dist(1-abs(x)))
dev.off()


print('Making pngs...')
# Print Figures to png
png(paste(resultsbin,Args[1],"_",strsplit(Args[2],split="_")[[1]][2],"_Results.png",
          sep=""))
par(mfrow=c(2,2))

# Histogram of P-values
hist(plrt,main=paste(Args[1],"\nNo. Probes = ",N),xlab="P-values")

# QQPlot
plot(1,1, type="n", xlim=c(0,max(null)+1),ylim=c(0,max(obsp)+1),xlab="Expected",
     ylab="Observed",main=Args[1])
polygon(xx, yy, col="gray",border="gray")
lines(null,null,col="dodgerblue2",lwd=2)
points(sort(null),sort(obsp),pch=20,cex=0.5)
if(length(tfprobes) > 0){
  for(i in 1:length(tfprobes)){
    points(sort(null)[which(sort(obsp) == obsp[tfprobes[i]])],obsp[tfprobes[i]],
           pch=20,cex=1.2,col="indianred")
  }
}

# MA plot
plot(a,m,pch=20,ylim=c(min(-1,min(m)),max(1,max(m))),cex=fivecex,main=Args[1],xlab="A",
     ylab="M",col=fivecol)
if(length(tfprobes) > 0){
  points(a[tfprobes],m[tfprobes],pch=20,col="indianred",cex=1.5)
}
abline(h=0,lty="dashed",lwd=2,col="dodgerblue2")

# Volcano plot
plot(m,-log10(plrt),xlim=c(min(-1,min(m)),max(1,max(m))),cex=0.5,pch=20,
     main=paste(Args[1],"\n5% FDR = ",fiveperc," probes",sep=""),
     xlab="Log2(Fold-Change)",ylab="-Log10(P-value)")
if(length(tfprobes) > 0){
  points(m[tfprobes],-log10(plrt[tfprobes]),pch=20,col="indianred",cex=1.5)
}
if(fiveperc > 0){
#  abline(h=-log10(plrt[which(qlrt == max(qlrt[qlrt <= 0.01]))[1]]),lty="dashed",lwd=2,col="blue")
  abline(h=-log10(plrt[which(qlrt == max(qlrt[qlrt <= 0.05]))[1]]),lty="dashed",
         lwd=2,col="dodgerblue2")
}
dev.off()



print('Writing results...')
likes = cbind(probeorder,data.frame(lnresults[,1]))
likes = cbind(likes, nullresults[,1])
names(likes) = c("Probe","LogLik","nullLogLik")
write.table(likes,paste(resultsbin,Args[1],"_",strsplit(Args[2],split="_")[[1]][2],
                        "_logLiks.txt",sep=""))


p.ind = order(plrt)
results = cbind(data.frame(plrt[p.ind]),qlrt[p.ind])
results = cbind(results,m[p.ind])
results = cbind(results,probeorder[p.ind])
results = cbind(results,genes[probe.ind[p.ind],])
names(results) = c("Pvalue","Qvalue","Log2FC","ProbeID","ENSGID","Symbol")
write.table(results,paste(resultsbin,Args[1],"_",strsplit(Args[2],split="_")[[1]][2],
                          "_Pvalues.txt",sep=""))

counters = cbind(c("TotalProbes","0.05FDR"),c(N,fiveperc))
write.table(counters,paste(resultsbin,Args[1],"_",strsplit(Args[2],split="_")[[1]][2],
                           "_counts.txt",sep=""),col.names=F,row.names=F,quote=F,
            sep="\t")

print('Running GO Analysis...')
genes4GO = results$ENSGID
siggenes = unique(genes4GO[1:fiveperc])
universe = rep(0,times=length(unique(genes4GO)))
universe[1:length(siggenes)] = 1
universe = as.factor(universe)
names(universe) = unique(genes4GO)

print('MF Categories first...')
GOdataMF <- new("topGOdata", ontology = "MF", allGenes=universe, geneSel=siggenes,nodeSize = 5,
              annot = annFUN.org,mapping = "org.Hs.eg",ID = "ensembl")
resultFisherMF <- runTest(GOdataMF, algorithm = "classic", statistic = "fisher")
allResMF <- GenTable(GOdataMF, classicFisher = resultFisherMF,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(score(resultFisherMF)))

print('Then BP categories...')
GOdataBP <- new("topGOdata", ontology = "BP", allGenes=universe, geneSel=siggenes,nodeSize = 5,
              annot = annFUN.org,mapping = "org.Hs.eg",ID = "ensembl")
resultFisherBP <- runTest(GOdataBP, algorithm = "classic", statistic = "fisher")
allResBP <- GenTable(GOdataBP, classicFisher = resultFisherBP,orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = length(score(resultFisherBP)))

allRes = rbind(allResMF,allResBP)
allRes.e = allRes[which(allRes$Significant > 0),]
pvals = allRes.e$classicFisher
allRes.e$AdjustedP = p.adjust(pvals,method="BH")
allRes.clean = allRes.e[which(allRes.e$Significant > 1),]
allRes.o = allRes.clean[order(allRes.clean$AdjustedP,allRes.clean$classicFisher,allRes.clean$Significant,allRes.clean$Term),]
if(dim(allRes.o)[1] == 0){allRes.o = "Wha-wha!  You lose!"}
write.table(allRes.o,paste(resultsbin,Args[1],"_",strsplit(Args[2],split="_")[[1]][2],"_GOResults.txt",sep=""))
print('Done.')