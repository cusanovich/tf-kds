library('plyr')
library('stringr')
resultsbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_Results/"
outdir = "/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Results/"
de_threshold = 0.05
factors = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt")
resultsmatrix = as.matrix(read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/allbindingresults10kb.txt",sep="\t"))
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
  currmaster = master[[2]][,c(2,i)]
  winners = unique(currmaster[which(currmaster[,2]<= de_threshold),1])
  winnerfactors = factors$V3[factors$V3 %in% winners]
  winnerfactors = setdiff(winnerfactors,factors[match(matchgene,factors[,1]),3])
  winnerfactornames = as.character(factors[match(winnerfactors,factors$V3),1])
  matrixcol = grep(matchgene,colnames(resultsmatrix))
  downstreamcol = matrixcol
  for(k in 1:length(winnerfactornames)){
    downstreamcol = c(downstreamcol,match(winnerfactornames[k],
                                          colnames(resultsmatrix)))
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
  boundgenes = rownames(resultsmatrix)[which(resultsmatrix[,matrixcol] > 0)]
  downstreambound = rownames(resultsmatrix)[which(rowSums(as.matrix(resultsmatrix[,downstreamcol])) > 0)]
  winnerbound = intersect(winners,boundgenes)
  justbound = setdiff(boundgenes,winners)
  downwinnerbound = intersect(winners,downstreambound)
  justdownbound = setdiff(downstreambound,winners)
  bound[[length(bound)+1]] = resultsmatrix[match(winnerbound,rownames(resultsmatrix)),matrixcol]
  bound.de[[length(bound.de)+1]] = resultsmatrix[match(winnerbound,rownames(resultsmatrix)),matrixcol]
  bound[[length(bound)+1]] = resultsmatrix[match(justbound,rownames(resultsmatrix)),matrixcol]
  bound.bound[[length(bound.bound)+1]] = resultsmatrix[match(justbound,rownames(resultsmatrix)),matrixcol]
  write.table(paste(winnerbound,collapse="|"),paste("/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/",currgene,".DEandBound.txt",sep=""),
              row.names=F,col.names=F,quote=F)
  write.table(paste(justbound,collapse="|"),paste("/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/",currgene,".Bound.txt",sep=""),
              row.names=F,col.names=F,quote=F)
  write.table(paste(downwinnerbound,collapse="|"),paste("/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/",currgene,".DEandBoundPlus.txt",sep=""),
              row.names=F,col.names=F,quote=F)
  write.table(paste(justdownbound,collapse="|"),paste("/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/",currgene,".BoundPlus.txt",sep=""),
              row.names=F,col.names=F,quote=F)
  write.table(paste(colnames(resultsmatrix)[matrixcol],collapse="|"),
              paste("/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/",currgene,".BindingTF.txt",sep=""),
              row.names=F,col.names=F,quote=F)
  write.table(paste(colnames(resultsmatrix)[downstreamcol],collapse="|"),
              paste("/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/",currgene,".DownstreamTFs.txt",sep=""),
              row.names=F,col.names=F,quote=F)
  if(length(boundgenes) < 1){
    bound.t = c(bound.t,NA)
    next
  }
  bound.t = c(bound.t,t.test(bound[[length(bound)-1]],bound[[length(bound)]])$statistic)
  if(length(downstreamcol) == 1){
    down[[length(down)+1]] = resultsmatrix[match(downwinnerbound,rownames(resultsmatrix)),downstreamcol]
    binaries[[length(binaries)+1]] = resultsbinary[match(downwinnerbound,rownames(resultsbinary)),downstreamcol]
    down.de[[length(down.de)+1]] = resultsmatrix[match(downwinnerbound,rownames(resultsmatrix)),downstreamcol]
    binary.de[[length(binary.de)+1]] = resultsbinary[match(downwinnerbound,rownames(resultsbinary)),downstreamcol]
    down[[length(down)+1]] = resultsmatrix[match(justdownbound,rownames(resultsmatrix)),downstreamcol]
    binaries[[length(binaries)+1]] = resultsbinary[match(justdownbound,rownames(resultsbinary)),downstreamcol]
    down.bound[[length(down.bound)+1]] = resultsmatrix[match(justdownbound,rownames(resultsmatrix)),downstreamcol]
    binary.bound[[length(binary.bound)+1]] = resultsbinary[match(justdownbound,rownames(resultsbinary)),downstreamcol]
    down.t = c(down.t,t.test(down[[length(down)-1]],down[[length(down)]])$statistic)
#    write.table(paste(winnerbound,collapse="|"),paste("../../GenomeAnnotations/Grepers/",currgene,".DEandBound.txt",sep=""),
#                row.names=F,col.names=F,quote=F)
#    write.table(paste(justbound,collapse="|"),paste("../../GenomeAnnotations/Grepers/",currgene,".Bound.txt",sep=""),
#                row.names=F,col.names=F,quote=F)
#    write.table(paste(downwinnerbound,collapse="|"),paste("../../GenomeAnnotations/Grepers/",currgene,".DEandBoundPlus.txt",sep=""),
#                row.names=F,col.names=F,quote=F)
#    write.table(paste(justdownbound,collapse="|"),paste("../../GenomeAnnotations/Grepers/",currgene,".BoundPlus.txt",sep=""),
#                row.names=F,col.names=F,quote=F)
    next
  }
  down[[length(down)+1]] = rowSums(resultsmatrix[match(downwinnerbound,rownames(resultsmatrix)),downstreamcol])
  binaries[[length(binaries)+1]] = rowSums(resultsbinary[match(downwinnerbound,rownames(resultsbinary)),downstreamcol])
  down.de[[length(down.de)+1]] = rowSums(resultsmatrix[match(downwinnerbound,rownames(resultsmatrix)),downstreamcol])
  binary.de[[length(binary.de)+1]] = rowSums(resultsbinary[match(downwinnerbound,rownames(resultsbinary)),downstreamcol])  
  down[[length(down)+1]] = rowSums(resultsmatrix[match(justdownbound,rownames(resultsmatrix)),downstreamcol])
  binaries[[length(binaries)+1]] = rowSums(resultsbinary[match(justdownbound,rownames(resultsbinary)),downstreamcol])
  down.bound[[length(down.bound)+1]] = rowSums(resultsmatrix[match(justdownbound,rownames(resultsmatrix)),downstreamcol])
  binary.bound[[length(binary.bound)+1]] = rowSums(resultsbinary[match(justdownbound,rownames(resultsbinary)),downstreamcol])  
  down.t = c(down.t,t.test(down[[length(down)-1]],down[[length(down)]])$statistic)
}

u.bound = wilcox.test(unlist(bound.de),unlist(bound.bound),na.rm=T)$p.value
u.down = wilcox.test(unlist(down.de),unlist(down.bound),na.rm=T)$p.value
u.binary = wilcox.test(unlist(binary.de),unlist(binary.bound),na.rm=T)$p.value

pdf("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Results/NumTFsBindingDEGenes.pdf")
par(mfrow=c(1,1))
boxplot(unlist(bound.de),unlist(bound.bound),na.rm=T,log="y",
        col=c("tomato","dodgerblue2"),las=1,ylab="No. of Binding Events",
        main=paste0("Kd Binding Only\nP-value = ",format(u.bound)),
        names=c("DE","Bound Only"),pch=20,outcol=c("tomato","dodgerblue2"))
boxplot(unlist(down.de),unlist(down.bound),na.rm=T,log="y",
        col=c("tomato","dodgerblue2"),las=1,ylab="No. of Binding Events",
        main=paste0("Kd & Downstream Binding\nP-value = ",format(u.down)),
        names=c("DE","Bound Only"),pch=20,outcol=c("tomato","dodgerblue2"))
boxplot(unlist(binary.de),unlist(binary.bound),na.rm=T,log="y",
        col=c("tomato","dodgerblue2"),las=1,ylab="No. of Factors Binding",
        main=paste0("Factor Variety\nP-value = ",format(u.binary)),
        names=c("DE","Bound Only"),pch=20,outcol=c("tomato","dodgerblue2"))
dev.off()
