library('plyr')
library('gplots')
library('beanplot')
library('plyr')
library('qvalue')
library('topGO')
library('gplots')
#Args <- commandArgs(TRUE)
#Args = c("EP300","expr_secondbatch_a3ruv2k8")
Args = c("IRF4","expr_oldbatch_ruv2")
print(Args[1])
print(Args[2])
resultsbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/"

print('Loading expression file...')
expr = as.matrix(read.table(paste("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Exprs/",Args[2],".txt",sep="")))
#expr = as.matrix(read.table(paste("~/home/Kd_Arrays/Analysis/Exprs/",Args[2],".txt",sep="")))
print('Loading detection file...')
detection = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Exprs/Detection_Scores_All3.txt")
#detection = read.table("~/home/Kd_Arrays/Analysis/Exprs/Detection_Scores_All3.txt")
probereport = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt",header=T)
#probereport = read.table("~/home/Kd_Arrays/Analysis/Annotations/HT-12v4R2_Probes_inhg19EnsemblGenes_NoGM19238SNPs_NoChrY_Stranded_OneProbePerGene_alt.txt",header=T)


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

probereport.ind = match(probing,probereport[,4])
probereportupdate = probereport[probereport.ind,]
probereportupdate = probereportupdate[!is.na(probereport.ind),]
genes = probereportupdate[,7:8]
probes.ind = match(probereportupdate[,4],probing)
probing = droplevels(probing[probes.ind])

expr_batches = svatable[probes.ind,]
rownames(expr_batches) = probereport[match(rownames(expr_batches),probereport[,4]),7]

de_threshold = 0.05
windowsize = 10000
windowname = paste0(windowsize/1000,'kb')
resultsbin = "~/home/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/"
bindingmatrix = paste0("~/home/Kd_Arrays/CombinedBinding/Binding/allunionbindingresults",windowname,".txt")

effectspdf = paste0(resultsbin,windowname,"_union_EffectsPlots.pdf")
resultsmatrix = as.matrix(read.table(bindingmatrix,sep="\t"))

all.pvals <- list.files(path = resultsbin,pattern="Pvalues.txt")
pvals <- llply(paste(resultsbin,all.pvals,sep=""), read.table)
namers = c()
for(i in 1:length(all.pvals)){
  namers[i] = strsplit(all.pvals[i],"batch")[[1]][1]
}
names(pvals) = namers

names.clean = c()
for(i in 1:length(namers)){
  names.clean[i] = strsplit(namers[i],"_")[[1]][1]
}

master = list()
for(j in 1:3){
  master[[j]] = pvals[[1]][,c(4:6,j)]
  names(master)[j] = colnames(pvals[[1]])[j]
  names(master[[j]])[4] = namers[1]
  for(i in 2:length(namers)){
    newbies = setdiff(as.character(pvals[[i]][,4]),master[[j]][,1])
    master[[j]] = merge(master[[j]],pvals[[i]][,c(4,j)],by="ProbeID",all=T)
    master[[j]][,1] = as.character(master[[j]][,1])
    master[[j]][,2] = as.character(master[[j]][,2])
    master[[j]][,3] = as.character(master[[j]][,3])
    names(master[[j]])[(3+i)] = namers[i]
    if(length(newbies) == 0){
      next
    }
    updates = match(newbies,as.character(pvals[[i]][,4]))
    places = match(newbies,master[[j]][,1])
    master[[j]][places,2] = as.character(pvals[[i]][updates,5])
    master[[j]][places,3] = as.character(pvals[[i]][updates,6])
  }
  
  for(k in 4:dim(master[[j]])[2]){
    curcol = master[[j]][,k]
    master[[j]][which(is.na(curcol)),k] = 5
  }
}
master[[2]] = master[[2]][match(master[[1]][,2],master[[2]][,2]),]
master[[3]] = master[[3]][match(master[[1]][,2],master[[3]][,2]),]

irf4.ind = grep("IRF4",colnames(master[[2]]))
irf4bound = rownames(resultsmatrix)[resultsmatrix[,match("IRF4",colnames(resultsmatrix))]>0]
irf4de = master[[2]][,2][master[[2]][,irf4.ind]<0.05]
irf4deandbound = intersect(irf4bound,irf4de)
irf4up = master[[3]][,2][master[[3]][,irf4.ind]>0]
irf4repressed = intersect(irf4up,irf4deandbound)
irf4activated = setdiff(irf4deandbound,irf4up)

probe.ind = match(irf4bound,rownames(expr_batches))
probe.ind = probe.ind[!is.na(probe.ind)]
tfprobeensg = probereport[match(Args[1],probereport[,8]),7]
tfprobe.ind = match(tfprobeensg,rownames(expr_batches))
probe.ind = setdiff(probe.ind,tfprobe.ind)
expr_batch = expr_batches[probe.ind,]
repressed.ind = match(irf4repressed,rownames(expr_batch))
activated.ind = match(irf4activated,rownames(expr_batch))
coling = rep("black",times=length(probe.ind))
coling[repressed.ind] = "indianred"
coling[activated.ind] = "dodgerblue2"
#also need genes that are de
# Calculate stats for MA plot and Volcano plot
tf.ind = grep(Args[1],colnames(expr_batch))
ns.ind = grep("NS",colnames(expr_batch))
meanns = apply(expr_batch[,ns.ind],1,mean)
meantf = apply(expr_batch[,tf.ind],1,mean)
a = apply(expr_batch[,c(ns.ind,tf.ind)],1,mean)
m = meantf - meanns

pdf("IRF4_inferring_repression.pdf")
plot(a,m,pch=20,ylim=c(-2,2),cex=1.5,main=Args[1],las=1,xlab="A",ylab="M",
     col=coling)
abline(h=0,lty="dashed",lwd=2,col="goldenrod")
dev.off()