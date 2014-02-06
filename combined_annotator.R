library('plyr')
library('stringr')
library('beanplot')
library('Hmisc')
#library('vioplot')
source('./config.R')

outdir = "/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Results/"
grepbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/Union/"
pdfout = paste0(outdir,windowname,"_union_NumTFsBindingDEGenes.pdf")
factors = as.matrix(read.table(factored))
resultsmatrix = as.matrix(read.table(bindingmatrix,sep="\t"))
resultsmatrix = resultsmatrix[which(rowSums(resultsmatrix) > 0),]
resultsmatrix = resultsmatrix[,which(colSums(resultsmatrix) > 0)]
resultsbinary = resultsmatrix>0
resultsbinary = resultsbinary + 0

factors = factors[match(colnames(resultsmatrix),factors[,1]),]

all.pvals <- list.files(path = resultsbin,pattern="Pvalues.txt")
pvals <- llply(paste(resultsbin,all.pvals,sep=""), read.table)
namers = c()
matchers = c()
for(i in 1:length(all.pvals)){
namers[i] = strsplit(all.pvals[i],"batch")[[1]][1]
matchers[i] = strsplit(all.pvals[i],"_")[[1]][1]
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
down = list()
down.bound = list()
down.de = list()
binaries = list()
binary.bound = list()
binary.de = list()
funcers = rep(0,length(rownames(resultsmatrix)))
winning = rep(0,length(rownames(resultsmatrix)))
nonfuncers = rep(0,length(rownames(resultsmatrix)))
mastercol = c()
expressed = c()
for(i in 4:dim(master[[2]])[2]){
  currgene = names(master[[2]])[i]
  matchgene = strsplit(currgene,"_")[[1]][1]
  commongenes = intersect(master[[2]][which(master[[2]][,i]< 5),2],rownames(resultsmatrix))
  commonind = match(commongenes,master[[2]][,2])
  currmaster = master[[2]][commonind,c(2,i)]
  currresultsmatrix = resultsmatrix[match(commongenes,rownames(resultsmatrix)),]
  currresultsbinary = resultsbinary[match(commongenes,rownames(resultsbinary)),]
  winners = unique(currmaster[which(currmaster[,2]<= de_threshold),1])
  winnerfactors = factors[,3][factors[,3] %in% winners]
  winnerfactors = setdiff(winnerfactors,factors[match(matchgene,factors[,1]),3])
  winnerfactornames = as.character(factors[match(winnerfactors,factors[,3]),1])
  expressed = union(expressed,rownames(currresultsmatrix))
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
  mastercol = union(mastercol,downstreamcol)
  boundgenes = rownames(currresultsmatrix)[which(currresultsmatrix[,matrixcol] > 0)]
  downstreambound = rownames(currresultsmatrix)[which(rowSums(as.matrix(currresultsmatrix[,downstreamcol])) > 0)]
  winnerbound = intersect(winners,boundgenes)
  justbound = setdiff(boundgenes,winners)
  downwinnerbound = intersect(winners,downstreambound)
  downwinnerrows = match(downwinnerbound,rownames(resultsmatrix))
  funcers[downwinnerrows] = funcers[downwinnerrows] + 1
  winning[match(winners,rownames(resultsmatrix))] = winning[match(winners,rownames(resultsmatrix))] + 1
  justdownbound = setdiff(downstreambound,winners)
  downboundrows = match(justdownbound,rownames(resultsmatrix))
  nonfuncers[downboundrows] = nonfuncers[downboundrows] + 1
  bound[[length(bound)+1]] = currresultsmatrix[match(winnerbound,rownames(currresultsmatrix)),matrixcol]
  bound.de[[length(bound.de)+1]] = currresultsmatrix[match(winnerbound,rownames(currresultsmatrix)),matrixcol]
  bound[[length(bound)+1]] = currresultsmatrix[match(justbound,rownames(currresultsmatrix)),matrixcol]
  bound.bound[[length(bound.bound)+1]] = currresultsmatrix[match(justbound,rownames(currresultsmatrix)),matrixcol]
  if(length(winnerbound) > 0){
    write.table(paste(winnerbound,collapse="|"),
                paste0(grepbin,windowname,"/",currgene,".DEandBound.txt"),
                row.names=F,col.names=F,quote=F)
  }
  if(length(justbound) > 0){
    write.table(paste(justbound,collapse="|"),
                paste0(grepbin,windowname,"/",currgene,".Bound.txt"),
                row.names=F,col.names=F,quote=F)
  }
  if(length(downwinnerbound) > 0){
    write.table(paste(downwinnerbound,collapse="|"),
                paste0(grepbin,windowname,"/",currgene,".DEandBoundPlus.txt"),
                row.names=F,col.names=F,quote=F)
  }
  if(length(justdownbound) > 0){
    write.table(paste(justdownbound,collapse="|"),
                paste0(grepbin,windowname,"/",currgene,".BoundPlus.txt"),
                row.names=F,col.names=F,quote=F)
  }
  if(length(matrixcol) > 0){
    write.table(paste(colnames(currresultsmatrix)[matrixcol],collapse="|"),
                paste0(grepbin,windowname,"/",currgene,".BindingTF.txt"),
                row.names=F,col.names=F,quote=F)
  }
  if(length(downstreamcol) > 0){
    write.table(paste(colnames(currresultsmatrix)[downstreamcol],collapse="|"),
                paste0(grepbin,windowname,"/",currgene,".DownstreamTFs.txt"),
                row.names=F,col.names=F,quote=F)
  }
  if(length(downstreamcol) == 1){
    down[[length(down)+1]] = currresultsmatrix[match(downwinnerbound,rownames(currresultsmatrix)),downstreamcol]
    binaries[[length(binaries)+1]] = currresultsbinary[match(downwinnerbound,rownames(currresultsbinary)),downstreamcol]
    down.de[[length(down.de)+1]] = currresultsmatrix[match(downwinnerbound,rownames(currresultsmatrix)),downstreamcol]
    names(down.de)[length(down.de)] = matchgene
    binary.de[[length(binary.de)+1]] = currresultsbinary[match(downwinnerbound,rownames(currresultsbinary)),downstreamcol]
    down[[length(down)+1]] = currresultsmatrix[match(justdownbound,rownames(currresultsmatrix)),downstreamcol]
    binaries[[length(binaries)+1]] = currresultsbinary[match(justdownbound,rownames(currresultsbinary)),downstreamcol]
    down.bound[[length(down.bound)+1]] = currresultsmatrix[match(justdownbound,rownames(currresultsmatrix)),downstreamcol]
    names(down.bound)[length(down.bound)] = matchgene
    binary.bound[[length(binary.bound)+1]] = currresultsbinary[match(justdownbound,rownames(currresultsbinary)),downstreamcol]
    next
  }
  down[[length(down)+1]] = rowSums(currresultsmatrix[match(downwinnerbound,rownames(currresultsmatrix)),downstreamcol])
  binaries[[length(binaries)+1]] = rowSums(currresultsbinary[match(downwinnerbound,rownames(currresultsbinary)),downstreamcol])
  down.de[[length(down.de)+1]] = rowSums(currresultsmatrix[match(downwinnerbound,rownames(currresultsmatrix)),downstreamcol])
  names(down.de)[length(down.de)] = matchgene
  binary.de[[length(binary.de)+1]] = rowSums(currresultsbinary[match(downwinnerbound,rownames(currresultsbinary)),downstreamcol])  
  down[[length(down)+1]] = rowSums(currresultsmatrix[match(justdownbound,rownames(currresultsmatrix)),downstreamcol])
  binaries[[length(binaries)+1]] = rowSums(currresultsbinary[match(justdownbound,rownames(currresultsbinary)),downstreamcol])
  down.bound[[length(down.bound)+1]] = rowSums(currresultsmatrix[match(justdownbound,rownames(currresultsmatrix)),downstreamcol])
  names(down.bound)[length(down.bound)] = matchgene
  binary.bound[[length(binary.bound)+1]] = rowSums(currresultsbinary[match(justdownbound,rownames(currresultsbinary)),downstreamcol])
}

u.bound = signif(wilcox.test(unlist(bound.de),unlist(bound.bound),na.rm=T)$p.value,4)
u.down = signif(wilcox.test(unlist(down.de),unlist(down.bound),na.rm=T)$p.value,4)
u.binary = signif(wilcox.test(unlist(binary.de),unlist(binary.bound),na.rm=T)$p.value,4)

pdf(pdfout)
par(mfrow=c(2,2))
boxplot(unlist(bound.de),unlist(bound.bound),na.rm=T,log="y",outline=F,
        col=c("indianred","dodgerblue2"),las=1,ylab="No. of Binding Events",
        main=paste0("Kd Binding Only\nP-value = ",format(u.bound)),
        names=c("DE","Bound Only"),notch=T,boxlwd=3,medlwd=4,cex.lab=1.5)
#plot(seq(0,2,length.out=6),
#     seq(0,max(c(log10(unlist(bound.de)),log10(unlist(bound.bound)))),length.out=6),
#     cex.lab=1.5,type="n",xlab="",ylab="Log10(No. of Binding Events)",las=1)
#vioplot(log10(unlist(bound.bound)), horizontal=FALSE, at=1.5,col="dodgerblue2",add=T)
#vioplot(log10(unlist(bound.de)), horizontal=FALSE, at=0.5,col="indianred",add=T)
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
de.counts = unlist(llply(down.de,length))
bound.counts = unlist(llply(down.bound,length))
countering = de.counts/(de.counts+bound.counts)
de.meds = unlist(llply(down.de,median))
bound.meds = unlist(llply(down.bound,median))
de.tops = unlist(llply(down.de,quantile,0.95))
bound.tops = unlist(llply(down.bound,quantile,0.95))
de.bottoms = unlist(llply(down.de,quantile,0.05))
bound.bottoms = unlist(llply(down.bound,quantile,0.05))
frac.order = order(countering)
expressind = match(expressed,rownames(resultsmatrix))
binding = rowSums(resultsmatrix[expressind,mastercol])
par(mfrow=c(1,1))
boxplot(down.bound[frac.order],na.rm=T,log="y",outline=F,las=2,col=adjustcolor("dodgerblue2",0.5),border=adjustcolor("dodgerblue2",0.8),
        cex.lab=0.5)
boxplot(down.de[frac.order],na.rm=T,log="y",outline=F,las=2,col=adjustcolor("indianred",0.5),border=adjustcolor("indianred",0.8),
        cex.lab=0.5,add=T)
boxplot(down.bound[frac.order],na.rm=T,log="y",outline=F,las=2,col=adjustcolor("dodgerblue2",0.6),border="dodgerblue2",
        cex.lab=1.5,varwidth=T,names=namers[frac.order])
boxplot(down.de[frac.order],na.rm=T,log="y",outline=F,las=2,col=adjustcolor("indianred",0.6),border="indianred",
        add=T,varwidth=T,names=namers[frac.order])
plot(de.tops[frac.order],type="b",col="indianred",pch=20)
points(de.meds[frac.order],type="b",col="indianred",pch=20)
points(bound.tops[frac.order],type="b",col="dodgerblue2",pch=20)
points(bound.meds[frac.order],type="b",col="dodgerblue2",pch=20)
plot(de.meds[frac.order],type="b",col="indianred",pch=20)
points(bound.meds[frac.order],type="b",col="dodgerblue2",pch=20)
plot(de.tops,type="n")
xs = c(1:56,56:1)
ys = c(de.tops[frac.order],de.bottoms[frac.order])
ysb = c(bound.tops[frac.order],bound.bottoms[frac.order])
polygon(xs,ys,col=adjustcolor("indianred",0.5),border=NA)
polygon(xs,ysb,col=adjustcolor("dodgerblue2",0.5),border=NA)
nonzeroes = which(binding > 0)
bindings = cut2(binding[nonzeroes],g=20)
fracfunc1 = by(funcers[expressind][nonzeroes],bindings,function(x){length(which(x > 0))/length(x)})
fracfunc2 = by(funcers[expressind][nonzeroes],bindings,function(x){length(which(x > 1))/length(x)})
fracfunc3 = by(funcers[expressind][nonzeroes],bindings,function(x){length(which(x > 4))/length(x)})
fracfunc4 = by(funcers[expressind][nonzeroes],bindings,function(x){length(which(x > 9))/length(x)})
fracfunc5 = by(funcers[expressind][nonzeroes],bindings,function(x){length(which(x > 19))/length(x)})
smoothScatter(log10(binding[nonzeroes]),funcers[expressind][nonzeroes])
smoothScatter(log10(binding[nonzeroes]),nonfuncers[expressind][nonzeroes])
hist(winning[expressind],breaks=50)
hist(binding,breaks=50)
plot(fracfunc1,type="b",col=1,ylim=c(min(fracfunc1,fracfunc2,fracfunc3,fracfunc4,fracfunc5),1))
points(fracfunc2,type="b",col=2)
points(fracfunc3,type="b",col=3)
points(fracfunc4,type="b",col=4)
points(fracfunc5,type="b",col=5)
print(countering[frac.order[1]])
print(countering[frac.order[length(frac.order)]])
print(length(countering))
print(min(funcers[expressind][nonzeroes]))
print(length(which(funcers[expressind][nonzeroes] == 0)))
print(levels(bindings))
print(range(binding))
print(median(binding))
dev.off()