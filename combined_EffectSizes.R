library('plyr')
library('gplots')
library('beanplot')
source('./config.R')

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

effects = list()
directions = list()
ups = c()
bounddirections = list()
boundups = c()
frac = c()
denom = c()
for(i in 4:dim(master[[3]])[2]){
  matchgene = strsplit(names(master[[2]])[i],"_")[[1]][1]
  commongenes = intersect(master[[2]][which(master[[2]][,i]< 5),2],rownames(resultsmatrix))
  commongenes = setdiff(commongenes,master[[2]][match(matchgene,master[[2]][,3]),2])
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
  frac = c(frac,paste(length(which(currmaster[[3]][currbound[currmaster[[2]][currbound,i] < 0.05],i] > 0)),length(currmaster[[3]][currbound[currmaster[[2]][currbound,i] < 0.05],i]),sep="/"))
  denom = c(denom,length(currmaster[[3]][currbound[currmaster[[2]][currbound,i] < 0.05],i]))
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
indish = order(colSums(master[[2]][,4:62] < 0.05))
inders = c(59,16,44,2,32,46,3,28,19,58,11,57,24,35,50,12,14,34,41,48,23,25,36,56,
           10,31,49,7,26,45,40,18,9,55,4,37,29,38,53,42,5,33,54,22,8,13,47,27,1,
           43,6,30,51,21,17,15,39,20,52)
qers = c(quantile(unlist(effects),0.25),
         quantile(unlist(effects),0.75))
boxplot(effects[indish],outline=F,col="mediumorchid3",las=2,cex.lab=2,axes=F,
        ylab="Log2(Fold-Change)",cex.lab=1.5)
polygon(c(-10,100,100,-10),c(qers[1],qers[1],qers[2],qers[2]),col="gray",
        border="gray")
boxplot(effects[indish],outline=F,col="mediumorchid3",las=2,cex.lab=2,axes=F,add=T)
axis(1, at = 1:59, labels = names.clean[indish],las=2,cex.axis=0.4)
axis(2,las=2)
box()
boxplot(effects[inders],outline=F,col="dodgerblue2",las=2,cex.lab=2,axes=F,
        ylab="Log2(Fold-Change)",cex.lab=1.5)
polygon(c(-10,100,100,-10),c(qers[1],qers[1],qers[2],qers[2]),col="gray",
        border="gray")
boxplot(effects[inders],outline=F,col="dodgerblue2",las=2,cex.lab=2,axes=F,add=T)
axis(1, at = 1:59, labels = names.clean[inders],las=2,cex.axis=0.4)
axis(2,las=2)
box()
beanplot(effects[indish],ylab="Log2(Fold-Change)",names=names.clean,col="indianred",las=2,log="",cex.lab=1.5,beanlines="median",what=c(1,1,1,0),overallline="median",axes=F)
axis(1, at = 1:59, labels = names.clean[indish],las=2,cex.axis=0.4)
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
stripchart(ups,vertical=T,cex=1.5,pch=20,col="indianred",method="jitter",las=1,ylab="Fraction of Genes Repressed")
abline(h=0.5,lty=3)
lines(c(0.9,1.1),c(median(ups),median(ups)),lwd=4)
kernels = c("gaussian", "epanechnikov", "rectangular",
           "triangular", "biweight",
           "cosine", "optcosine")
adjusters = c(1)
plot(seq(-.5,.5,length.out=6),seq(0,4.5,length.out=6),cex.lab=1.5,type="n",xlab="Log2(Fold-Change)",ylab="Density",las=1)
abline(v=0,lty="dashed")
polygon(density(unlist(directions)),col="dodgerblue2")
legend("topright","Repressed Genes",bty="n")
legend("topleft","Enhanced Genes",bty="n")
boxplot(boundups,notch=T,outpch=20,boxlwd=3,medlwd=4,las=1,cex=2,cex.lab=1.5,col="indianred",outcol="indianred",ylab="Fraction of Direct Targets Repressed")
abline(h=0.5,lty=3)
boxplot(boundups,notch=T,outpch=20,boxlwd=3,medlwd=4,add=T,axes=F,col="indianred",outcol="indianred")
for(i in 1:length(adjusters)){
  plot(seq(-2,2,length.out=6),seq(0,4.5,length.out=6),cex.lab=1.5,type="n",xlab="Log2(Fold-Change) of Direct Targets",ylab="Density",las=1)
  abline(v=0,lty="dashed")
  polygon(density(unlist(bounddirections),adjust=adjusters[i]),col="dodgerblue2")
  legend("topright","Repressed Direct Targets",bty="n")
  legend("topleft","Enhanced Direct Targets",bty="n")
}
hist(unlist(bounddirections),xlim=c(-2,2),col="dodgerblue2",breaks=75,las=1,xlab="Log2(Fold-Change) of Direct Targets")
box()
abline(v=0,lty="dashed")
legend("topright","Repressed Direct Targets",bty="n")
legend("topleft","Enhanced Direct Targets",bty="n")
par(oma=c(0,8,0,8) + 0.1)
stripchart(boundups,vertical=T,cex=2,pch=20,col="indianred",method="jitter",las=1,ylab="Fraction of Genes Repressed")
abline(h=0.5,lty=2)
lines(c(0.9,1.1),c(median(boundups),median(boundups)),lwd=3)
dev.off()
directioner = cbind(names(bounddirections),frac,boundups)
write.table(directioner,paste0(resultsbin,windowname,"directions.txt"),row.names=F,
            col.names=F, quote=F, sep="\t")
print(names(bounddirections)[which(boundups == min(boundups))])
print(min(boundups))
print(names(bounddirections)[which(boundups == max(boundups))])
print(max(boundups))
print(boundups)
print(names(bounddirections))
print(frac)
print(median(denom))