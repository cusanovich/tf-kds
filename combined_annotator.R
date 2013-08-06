library('plyr')
library('stringr')
windowsize="1kb"
resultsbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/"
outdir = "/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Results/"
grepbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/"
de_threshold = 0.05
factors = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt")
resultsmatrix = as.matrix(read.table(paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/allbindingresults",windowsize,".txt"),sep="\t"))
resultsmatrix = resultsmatrix[which(rowSums(resultsmatrix) > 0),]
resultsmatrix = resultsmatrix[,which(colSums(resultsmatrix) > 0)]
resultsbinary = resultsmatrix>0
resultsbinary = resultsbinary + 0

factors = factors[match(colnames(resultsmatrix),factors$V1),]

all.pvals <- list.files(path = resultsbin,pattern="Pvalues.txt")
pvals <- llply(paste(resultsbin,all.pvals,sep=""), read.table)
namers = c()
for(i in 1:length(all.pvals)){
namers[i] = strsplit(all.pvals[i],"batch")[[1]][1]
}
names(pvals) = namers

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
}
master[[2]] = master[[2]][match(master[[1]][,2],master[[2]][,2]),]
master[[3]] = master[[3]][match(master[[1]][,2],master[[3]][,2]),]

bound = list()
bound.bound = list()
bound.de = list()
bound.t = c()
down = list()
down.bound = list()
down.de = list()
binaries = list()
binary.bound = list()
binary.de = list()
down.t = c()
for(i in 4:dim(master[[2]])[2]){
  currgene = names(master[[2]])[i]
  matchgene = strsplit(currgene,"_")[[1]][1]
  commongenes = intersect(master[[2]][which(master[[2]][,i]< 5),2],rownames(resultsmatrix))
  commonind = match(commongenes,master[[2]][,2])
  currmaster = master[[2]][commonind,c(2,i)]
  currresultsmatrix = resultsmatrix[match(commongenes,rownames(resultsmatrix)),]
  currresultsbinary = resultsbinary[match(commongenes,rownames(resultsbinary)),]
  winners = unique(currmaster[which(currmaster[,2]<= de_threshold),1])
  winnerfactors = factors$V3[factors$V3 %in% winners]
  winnerfactors = setdiff(winnerfactors,factors[match(matchgene,factors[,1]),3])
  winnerfactornames = as.character(factors[match(winnerfactors,factors$V3),1])
  matrixcol = grep(matchgene,colnames(currresultsmatrix))
  downstreamcol = matrixcol
  for(k in 1:length(winnerfactornames)){
    downstreamcol = c(downstreamcol,match(winnerfactornames[k],
                                          colnames(currresultsmatrix)))
  }
  if(length(matrixcol) == 0 & length(downstreamcol) == 1 & is.na(downstreamcol[1])){
    print(paste0("Whoa! ",matchgene," doesn't bind anything!!!"))
    next
  }
  if(length(winnerfactornames) == 0){
    downstreamcol = downstreamcol[1]
  }
  if(length(matrixcol) > 1){
    print(paste0("Whoa! ",matchgene,
                 " matches too many columns in the binding matrix!!!"))
    next
  }
  boundgenes = rownames(currresultsmatrix)[which(currresultsmatrix[,matrixcol] > 0)]
  downstreambound = rownames(currresultsmatrix)[which(rowSums(as.matrix(currresultsmatrix[,downstreamcol])) > 0)]
  winnerbound = intersect(winners,boundgenes)
  justbound = setdiff(boundgenes,winners)
  downwinnerbound = intersect(winners,downstreambound)
  justdownbound = setdiff(downstreambound,winners)
  bound[[length(bound)+1]] = currresultsmatrix[match(winnerbound,rownames(currresultsmatrix)),matrixcol]
  bound.de[[length(bound.de)+1]] = currresultsmatrix[match(winnerbound,rownames(currresultsmatrix)),matrixcol]
  bound[[length(bound)+1]] = currresultsmatrix[match(justbound,rownames(currresultsmatrix)),matrixcol]
  bound.bound[[length(bound.bound)+1]] = currresultsmatrix[match(justbound,rownames(currresultsmatrix)),matrixcol]
  write.table(paste(winnerbound,collapse="|"),
              paste0(grepbin,windowsize,"/",currgene,".DEandBound.txt"),row.names=F,
              col.names=F,quote=F)
  write.table(paste(justbound,collapse="|"),
              paste0(grepbin,windowsize,"/",currgene,".Bound.txt"),row.names=F,col.names=F,
              quote=F)
  write.table(paste(downwinnerbound,collapse="|"),
              paste0(grepbin,windowsize,"/",currgene,".DEandBoundPlus.txt"),row.names=F,
              col.names=F,quote=F)
  write.table(paste(justdownbound,collapse="|"),
              paste0(grepbin,windowsize,"/",currgene,".BoundPlus.txt"),row.names=F,
              col.names=F,quote=F)
  write.table(paste(colnames(currresultsmatrix)[matrixcol],collapse="|"),
              paste0(grepbin,windowsize,"/",currgene,".BindingTF.txt"),row.names=F,
              col.names=F,quote=F)
  write.table(paste(colnames(currresultsmatrix)[downstreamcol],collapse="|"),
              paste0(grepbin,windowsize,"/",currgene,".DownstreamTFs.txt"),row.names=F,
              col.names=F,quote=F)
  if(length(boundgenes) < 1){
    bound.t = c(bound.t,NA)
    next
  }
  if(var(c(bound[[length(bound)-1]],bound[[length(bound)]])) != 0){
    bound.t = c(bound.t,t.test(bound[[length(bound)-1]],bound[[length(bound)]])$statistic)
  }else{
    bound.t = c(bound.t,0)
  }
  if(length(downstreamcol) == 1){
    down[[length(down)+1]] = currresultsmatrix[match(downwinnerbound,rownames(currresultsmatrix)),downstreamcol]
    binaries[[length(binaries)+1]] = currresultsbinary[match(downwinnerbound,rownames(currresultsbinary)),downstreamcol]
    down.de[[length(down.de)+1]] = currresultsmatrix[match(downwinnerbound,rownames(currresultsmatrix)),downstreamcol]
    binary.de[[length(binary.de)+1]] = currresultsbinary[match(downwinnerbound,rownames(currresultsbinary)),downstreamcol]
    down[[length(down)+1]] = currresultsmatrix[match(justdownbound,rownames(currresultsmatrix)),downstreamcol]
    binaries[[length(binaries)+1]] = currresultsbinary[match(justdownbound,rownames(currresultsbinary)),downstreamcol]
    down.bound[[length(down.bound)+1]] = currresultsmatrix[match(justdownbound,rownames(currresultsmatrix)),downstreamcol]
    binary.bound[[length(binary.bound)+1]] = currresultsbinary[match(justdownbound,rownames(currresultsbinary)),downstreamcol]
    if(var(c(down[[length(down)-1]],down[[length(down)]])) != 0){
      down.t = c(down.t,t.test(down[[length(down)-1]],down[[length(down)]])$statistic)
    }else{
      down.t = c(down.t,0)
    }
    next
  }
  down[[length(down)+1]] = rowSums(currresultsmatrix[match(downwinnerbound,rownames(currresultsmatrix)),downstreamcol])
  binaries[[length(binaries)+1]] = rowSums(currresultsbinary[match(downwinnerbound,rownames(currresultsbinary)),downstreamcol])
  down.de[[length(down.de)+1]] = rowSums(currresultsmatrix[match(downwinnerbound,rownames(currresultsmatrix)),downstreamcol])
  binary.de[[length(binary.de)+1]] = rowSums(currresultsbinary[match(downwinnerbound,rownames(currresultsbinary)),downstreamcol])  
  down[[length(down)+1]] = rowSums(currresultsmatrix[match(justdownbound,rownames(currresultsmatrix)),downstreamcol])
  binaries[[length(binaries)+1]] = rowSums(currresultsbinary[match(justdownbound,rownames(currresultsbinary)),downstreamcol])
  down.bound[[length(down.bound)+1]] = rowSums(currresultsmatrix[match(justdownbound,rownames(currresultsmatrix)),downstreamcol])
  binary.bound[[length(binary.bound)+1]] = rowSums(currresultsbinary[match(justdownbound,rownames(currresultsbinary)),downstreamcol])  
  if(var(c(down[[length(down)-1]],down[[length(down)]])) != 0){
    down.t = c(down.t,t.test(down[[length(down)-1]],down[[length(down)]])$statistic)
  }else{
    down.t = c(down.t,0)
  }
}

u.bound = signif(wilcox.test(unlist(bound.de),unlist(bound.bound),na.rm=T)$p.value,4)
u.down = signif(wilcox.test(unlist(down.de),unlist(down.bound),na.rm=T)$p.value,4)
u.binary = signif(wilcox.test(unlist(binary.de),unlist(binary.bound),na.rm=T)$p.value,4)

pdf(paste0(outdir,windowsize,"NumTFsBindingDEGenes.pdf"))
par(mfrow=c(2,2))
boxplot(unlist(bound.de),unlist(bound.bound),na.rm=T,log="y",outline=F,
        col=c("indianred","dodgerblue2"),las=1,ylab="No. of Binding Events",
        main=paste0("Kd Binding Only\nP-value = ",format(u.bound)),
        names=c("DE","Bound Only"),notch=T,boxlwd=3,medlwd=4,cex.lab=1.5)
print(length(unlist(bound.bound)))
print(mean(unlist(bound.bound)))
boxplot(unlist(down.de),unlist(down.bound),na.rm=T,log="y",outline=F,
        col=c("indianred","dodgerblue2"),las=1,ylab="No. of Binding Events",
        main=paste0("Kd & Downstream Binding\nP-value = ",format(u.down)),
        names=c("DE","Bound Only"),notch=T,boxlwd=3,medlwd=4,cex.lab=1.5)
boxplot(unlist(binary.de),unlist(binary.bound),na.rm=T,log="y",outline=F,
        col=c("indianred","dodgerblue2"),las=1,ylab="No. of Factors Binding",
        main=paste0("Factor Variety\nP-value = ",format(u.binary)),
        names=c("DE","Bound Only"),notch=T,boxlwd=3,medlwd=4,cex.lab=1.5)
dev.off()