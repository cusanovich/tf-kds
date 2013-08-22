library('plyr')
library('gplots')
library('beanplot')
source('./config.R')

effectspdf = paste0(resultsbin,windowname,"_EffectsPlots.pdf")
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

effects = list()
directions = list()
ups = c()
bounddirections = list()
boundups = c()
for(i in 4:dim(master[[3]])[2]){
  matchgene = strsplit(names(master[[2]])[i],"_")[[1]][1]
  commongenes = intersect(master[[2]][which(master[[2]][,i]< 5),2],rownames(resultsmatrix))
  commonind = match(commongenes,master[[2]][,2])
  currmaster = list()
  currmaster[[1]] = master[[1]][commonind,]
  currmaster[[2]] = master[[2]][commonind,]
  currmaster[[3]] = master[[3]][commonind,]
  currresultsmatrix = resultsmatrix[match(commongenes,rownames(resultsmatrix)),]
  matrixcol = match(matchgene,colnames(currresultsmatrix))
  directions[[i-3]] = currmaster[[3]][currmaster[[2]][,i] < 0.05,i]
  ups[i-3] = length(which(currmaster[[3]][currmaster[[2]][,i] < 0.05,i] > 0))/length(currmaster[[3]][currmaster[[2]][,i] < 0.05,i])
  effects[[i-3]] = abs(currmaster[[3]][currmaster[[2]][,i] < 0.05,i])  
  if(is.na(matrixcol)){next}
  boundgenes = rownames(currresultsmatrix)[currresultsmatrix[,matrixcol] > 0]
  currbound = match(boundgenes,currmaster[[1]][,2])
  currbound = currbound[!is.na(currbound)]
  bounddirections[[length(bounddirections)+1]] = currmaster[[3]][currbound[currmaster[[2]][currbound,i] < 0.05],i]
  names(bounddirections)[length(bounddirections)] = matchgene
  boundups = c(boundups,length(which(currmaster[[3]][currbound[currmaster[[2]][currbound,i] < 0.05],i] > 0))/length(currmaster[[3]][currbound[currmaster[[2]][currbound,i] < 0.05],i]))
}

densedirs = llply(directions,density)
denseffects = llply(effects,density)

pdf(effectspdf,height=11,width=8.5)
par(mgp=c(3,1,0))
par(mar=c(5, 5, 4, 2) + 0.1)
par(mfrow=c(2,1))
par(oma=c(0,4,0,4) + 0.1)
boxplot(directions,outline=F,las=2,names=names.clean,axes=F)
axis(1, at = 1:59, labels = names.clean,las=2,cex.axis=0.4)
axis(2,las=2)
box()
abline(h=0,lwd=3)
boxplot(directions,outline=F,add=T,col="dodgerblue2",axes=F)
beanplot(directions,ylab="Log2(Fold-Change)",names=names.clean,las=2,cex.lab=1.5,axes=F,beanlines="median",what=c(1,1,1,0),col="dodgerblue2")
axis(1, at = 1:length(directions), labels = names.clean,las=2,cex.axis=0.4)
axis(2,las=2)
box()
beanplot(bounddirections,ylab="Log2(Fold-Change) of Direct Targets",names=names(bounddirections),las=2,cex.lab=1.5,axes=F,beanlines="median",what=c(1,1,1,0),col="dodgerblue2")
axis(1, at = 1:length(bounddirections), labels = names(bounddirections),las=2,cex.axis=0.4)
axis(2,las=2)
box()
plot(densedirs[[1]],xlab="Log2(Fold-Change)",ylim=c(0,4.5),las=1,lwd=2,cex.lab=2,col="dodgerblue2")
for(i in 2:length(densedirs)){
  lines(densedirs[[i]],col="dodgerblue2",lwd=2)
}
boxplot(effects,outline=F,col="indianred",las=2,cex.lab=2,axes=F,ylab="Log2(Fold-Change)",cex.lab=1.5)
axis(1, at = 1:59, labels = names.clean,las=2,cex.axis=0.4)
axis(2,las=2)
box()
beanplot(effects,ylab="Log2(Fold-Change)",names=names.clean,col="indianred",las=2,log="",cex.lab=1.5,beanlines="median",what=c(1,1,1,0),overallline="median",axes=F)
axis(1, at = 1:59, labels = names.clean,las=2,cex.axis=0.4)
axis(2,las=2)
box()
plot(denseffects[[1]],xlab="Absolute Log2(Fold-Change)",ylim=c(0,13),las=1,lwd=2,cex.lab=2,col="indianred")
for(i in 2:length(denseffects)){
  lines(denseffects[[i]],col="indianred",lwd=2)
}
boxplot(ups,notch=T,outpch=20,boxlwd=3,medlwd=4,las=1,cex=2,cex.lab=1.5,col="indianred",outcol="indianred",ylab="Fraction of Genes Repressed")
abline(h=0.5,lty=3)
boxplot(ups,notch=T,outpch=20,boxlwd=3,medlwd=4,add=T,axes=F,col="indianred",outcol="indianred")
plot(seq(-.5,.5,length.out=6),seq(0,4.5,length.out=6),cex.lab=1.5,type="n",xlab="Log2(Fold-Change)",ylab="Density",las=1)
abline(v=0,lty="dashed")
polygon(density(unlist(directions)),col="dodgerblue2")
legend("topright","Repressed Genes",bty="n")
legend("topleft","Enhanced Genes",bty="n")

boxplot(boundups,notch=T,outpch=20,boxlwd=3,medlwd=4,las=1,cex=2,cex.lab=1.5,col="indianred",outcol="indianred",ylab="Fraction of Direct Targets Repressed")
abline(h=0.5,lty=3)
boxplot(boundups,notch=T,outpch=20,boxlwd=3,medlwd=4,add=T,axes=F,col="indianred",outcol="indianred")
plot(seq(-.5,.5,length.out=6),seq(0,4.5,length.out=6),cex.lab=1.5,type="n",xlab="Log2(Fold-Change) of Direct Targets",ylab="Density",las=1)
abline(v=0,lty="dashed")
polygon(density(unlist(bounddirections)),col="dodgerblue2")
legend("topright","Repressed Direct Targets",bty="n")
legend("topleft","Enhanced Direct Targets",bty="n")
dev.off()
print(names(bounddirections)[which(boundups == min(boundups))])
print(min(boundups))
print(names(bounddirections)[which(boundups == max(boundups))])
print(max(boundups))
print(boundups)
print(names(bounddirections))