library('plyr')
library('gplots')
library('beanplot')
library('stringr')
resultsbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/"
kdlevels = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/TargetKdLevels.txt")
targetgenes = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Annotations/TargetSummary.txt")

all.pvals <- list.files(path = resultsbin,pattern="Pvalues.txt")
pvals <- llply(paste(resultsbin,all.pvals,sep=""), read.table)
namers = c()
for(i in 1:length(all.pvals)){
namers[i] = strsplit(all.pvals[i],"batch")[[1]][1]
}
names(pvals) = namers

all.counts <- list.files(path = resultsbin,pattern="counts.txt")
counts <- llply(paste(resultsbin,all.counts,sep=""), read.table)
namers = c()
for(i in 1:length(all.counts)){
namers[i] = strsplit(all.counts[i],"batch")[[1]][1]
}
names(counts) = namers
names.clean = c()
for(i in 1:length(namers)){
names.clean[i] = strsplit(namers[i],"_")[[1]][1]
}

all.pdfs <- list.files(path = resultsbin,pattern="_Results.pdf")
batchers = c()
for(i in 1:length(all.pdfs)){
batchers[i] = strsplit(all.pdfs[i],"_")[[1]][2]
}
batchers= as.factor(batchers)

avekds = apply(kdlevels[,2:4],1,mean)
varkds = apply(kdlevels[,2:4],1,var)
#Build a master list of matrices that contains all P-values, Q-values, and Log2(FC) for each knockdown.  Genes DE are coded with a 5 instead of 'NA'.
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

master.overlap = rowSums(master[[2]][,4:dim(master[[2]])[2]] < 0.05)
lessthanfive = master[[2]][which(master.overlap < 5),1]

knockdowncors = matrix(NA,1,5)
for(i in 1:dim(targetgenes)[1]){
  inder = grep(targetgenes[i,1],kdlevels[,1])
  indy = match(targetgenes[i,3],master[[3]][,2])
  if(length(inder) > 0 & !is.na(indy)){
    newline = c(as.character(targetgenes[i,1]),as.character(kdlevels[inder,1]),
                avekds[inder],master[[3]][indy,match(kdlevels[inder,1],colnames(master[[3]]))],
                as.numeric(counts[[match(kdlevels[inder,1],namers)]][2,2]))
    knockdowncors = rbind(knockdowncors,newline)
  }
}
knockdowncors = knockdowncors[-1,]
knockdowncors = knockdowncors[knockdowncors[,4] < 5,]

summar = matrix(NA,length(counts),3)
for(i in 1:length(counts)){
	summar[i,1] = namers[i]
	summar[i,2] = as.numeric(counts[[i]][1,2])
	summar[i,3] = as.numeric(counts[[i]][2,2])
}
names(summar) = c("Factor","TotalProbes","SigProbes")
write.table(summar,paste(resultsbin,"SummaryCounts.txt",sep=""),row.names=F,quote=F,sep="\t")

batchcol = c("indianred","goldenrod","dodgerblue2")[unlist(batchers)]
pdf(paste(resultsbin,"SummaryPlots.pdf",sep=""),height=11,width=8.5)
barplot(as.numeric(summar[,2]),names.arg=names.clean,las=2,cex.names=0.8,col=batchcol,ylim=c(0,(max(as.numeric(summar[,2]))+3000)),
    main="Total No. of Probes",ylab="No. Probes",legend.text=c("Batch 1","Batch 2","Batch 3"),
    args.legend = list(fill=c("indianred","dodgerblue2","goldenrod")))
#legend(x="topright",legend=c("Batch 1", "Batch 2", "Batch 3"),fill=c("blue","red","black"))
degenes = as.numeric(summar[,3])
inder = order(degenes)
bp = barplot(degenes[inder],horiz=T,names.arg=names.clean[inder],xlim=c(0,(max(degenes)+3000)),
    las=2,cex.names=0.8,col=batchcol[inder],main="No. DE Genes",xlab="No. Genes",
    legend.text=c("Batch 1","Batch 2","Batch 3"),args.legend = list(fill=c("goldenrod","dodgerblue2","indianred")))
text(x=degenes[inder],y=bp,labels=degenes[inder],cex=0.8,pos=4)
bp2 = barplot(degenes[inder],horiz=T,names.arg=names.clean[inder],xlim=c(0,(max(degenes)+1100)),
              las=1,cex.names=0.8,col="indianred",main="No. DE Genes",xlab="No. Genes",)
text(x=degenes[inder],y=bp2,labels=degenes[inder],cex=0.8,pos=4)
bp3 = barplot(degenes[inder],ylim=c(0,(max(degenes)+1100)),las=1,cex.names=0.8,
              col="indianred",main="No. DE Genes",xlab="No. Genes",)
par(mfrow=c(2,2))
# par(oma=c(4,4,4,4) + 0.1)
# par(mgp=c(3,1,0))
# par(mar=c(5, 5, 4, 2) + 0.1)
a = avekds
b = log10(as.numeric(summar[,3]))
plot(a,b,las=1,xlab="%Knockdown (qPCR)",type="n",cex.lab=1.5,cex.axis=1.5,
     ylab="Log10(No. Genes DE)",main=cor.test(a,b)$p.value)
abline(lm(b ~ a),lty=2,lwd=2)
points(a,b,pch=20,col="dodgerblue2",cex=2)
x = round(cor(a,b)^2,3)
legend("topleft",legend=bquote(R^2 == .(x)),bty="n",cex=1.5)
a = varkds
plot(a,b,las=1,xlab="%Knockdown (qPCR) Variance",cex.lab=1.5,cex.axis=1.5,type="n",
     ylab="Log10(No. Genes DE)",main=cor.test(a,b)$p.value)
abline(lm(b ~ a),lty=2,lwd=2)
points(a,b,pch=20,col="mediumseagreen",cex=2)
x = round(cor(a,b)^2,3)
legend("topleft",legend=bquote(R^2 == .(x)),bty="n",cex=1.5)
a = (1-(2^as.numeric(knockdowncors[,4])))
b = log10(as.numeric(knockdowncors[,5]))
plot(a,b,las=1,xlab="%Knockdown (Microarray)",cex.lab=1.5,cex.axis=1.5,type="n",
     ylab="Log10(No. Genes DE)",main=cor.test(a,b)$p.value)
abline(lm(b ~ a),lty=2,lwd=2)
points(a,b,pch=20,col="goldenrod",cex=2)
x = round(cor(a,b)^2,3)
legend("topleft",legend=bquote(R^2 == .(x)),bty="n",cex=1.5)
a = as.numeric(knockdowncors[,3])
b = (1-(2^as.numeric(knockdowncors[,4])))
plot(a,b,las=1,xlab="%Knockdown (qPCR)",cex.lab=1.5,cex.axis=1.5,type="n",
     ylab="%Knockdown (Microarray)",main=cor.test(a,b)$p.value)
abline(lm(b ~ a),lty=2,lwd=2)
points(a,b,pch=20,col="mediumorchid3",cex=2)
x = round(cor(a,b)^2,3)
legend("topleft",legend=bquote(R^2 == .(x)),bty="n",cex=1.5)
dev.off()

results.m = matrix(NA,1,4)
for(i in 1:(length(names(pvals))-1)){
	for(j in (i+1):length(names(pvals))){
		genes1 = pvals[[i]]
		genes2 = pvals[[j]]
		siggenes1 = genes1[which(genes1$Qvalue < 0.05),4]
		siggenes2 = genes2[which(genes2$Qvalue < 0.05),4]
		q = length(intersect(siggenes1,siggenes2)) - 1
		m = length(siggenes1)
		n = length(intersect(genes1$ProbeID,genes2$ProbeID)) - m
		k = length(siggenes2)
		p = phyper(q,m,n,k, lower.tail=FALSE)
		newline = c(names(pvals)[i],names(pvals)[j],(q+1),p)
		results.m = rbind(results.m,newline)
	}
}
names(results.m) = c("Factor1","Factor2","Overlap","Pvalue")
results.m = results.m[-1,]
write.table(results.m,paste(resultsbin,"SummaryOverlap.txt",sep=""),row.names=F,quote=F,sep="\t")

firstnames = subset(namers,batchers == "firstbatch")
secondnames = subset(namers,batchers == "secondbatch")
first.clean = subset(names.clean,batchers == "firstbatch")
second.clean = subset(names.clean,batchers == "secondbatch")

results.heat.m = matrix(NA,length(namers),length(namers))
for(i in 1:length(namers)){
	for(j in i:length(namers)){
		if(i == j){
			#first.heat.m[i,j] = log10(min(as.numeric(results.m[,4])))
			results.heat.m[i,j] = -20
		}else{
			ier = grep(namers[i],results.m[,1])
			jer = grep(namers[j],results.m[,2])
			currval = as.numeric(results.m[intersect(ier,jer),4])
			results.heat.m[i,j] = -log10(currval)
			results.heat.m[j,i] = -log10(currval)
		}
	}
}
colnames(results.heat.m) = names.clean
rownames(results.heat.m) = names.clean

first.heat.m = matrix(NA,length(firstnames),length(firstnames))
for(i in 1:length(firstnames)){
	for(j in i:length(firstnames)){
		if(i == j){
			#first.heat.m[i,j] = log10(min(as.numeric(results.m[,4])))
			first.heat.m[i,j] = -20
		}else{
			ier = grep(firstnames[i],results.m[,1])
			jer = grep(firstnames[j],results.m[,2])
			currval = as.numeric(results.m[intersect(ier,jer),4])
			first.heat.m[i,j] = -log10(currval)
			first.heat.m[j,i] = -log10(currval)
		}
	}
}
colnames(first.heat.m) = first.clean
rownames(first.heat.m) = first.clean


second.heat.m = matrix(NA,length(secondnames),length(secondnames))
for(i in 1:length(secondnames)){
	for(j in i:length(secondnames)){
		if(i == j){
			#second.heat.m[i,j] = log10(min(as.numeric(results.m[,4])))
			second.heat.m[i,j] = -20
		}else{
			ier = grep(secondnames[i],results.m[,1])
			jer = grep(secondnames[j],results.m[,2])
			currval = as.numeric(results.m[intersect(ier,jer),4])
			second.heat.m[i,j] = -log10(currval)
			second.heat.m[j,i] = -log10(currval)
		}
	}
}
colnames(second.heat.m) = second.clean
rownames(second.heat.m) = second.clean

results.ma = matrix(NA,1,4)
for(i in 1:(length(names(pvals))-1)){
	for(j in (i+1):length(names(pvals))){
		ind1 = match(lessthanfive,rownames(pvals[[i]]))
		ind1 = ind1[which(!is.na(ind1))]
		ind2 = match(lessthanfive,rownames(pvals[[j]]))
		ind2 = ind2[which(!is.na(ind2))]
		genes1 = pvals[[i]][ind1,]
		genes2 = pvals[[j]][ind2,]
		siggenes1 = genes1[which(genes1$Qvalue < 0.05),5]
		siggenes2 = genes2[which(genes2$Qvalue < 0.05),5]
		q = length(intersect(siggenes1,siggenes2)) - 1
		m = length(siggenes1)
		n = length(intersect(genes1$ProbeID,genes2$ProbeID)) - m
		k = length(siggenes2)
		p = phyper(q,m,n,k, lower.tail=FALSE)
		newline = c(names(pvals)[i],names(pvals)[j],(q+1),p)
		results.ma = rbind(results.ma,newline)
	}
}
names(results.ma) = c("Factor1","Factor2","Overlap","Pvalue")
results.ma = results.ma[-1,]

results.heat.ma = matrix(NA,length(namers),length(namers))
for(i in 1:length(namers)){
	for(j in i:length(namers)){
		if(i == j){
			results.heat.ma[i,j] = -log10(min(as.numeric(results.ma[,4])))
			#results.heat.ma[i,j] = -20
		}else{
			ier = grep(namers[i],results.ma[,1])
			jer = grep(namers[j],results.ma[,2])
			currval = as.numeric(results.ma[intersect(ier,jer),4])
			results.heat.ma[i,j] = -log10(currval)
			results.heat.ma[j,i] = -log10(currval)
		}
	}
}
colnames(results.heat.ma) = names.clean
rownames(results.heat.ma) = names.clean

resultsmatrix = master[[2]][,4:dim(master[[2]])[2]]
resultsbinary = resultsmatrix<0.05
resultsbinary = resultsbinary + 0
keepers = resultsmatrix < 2
keepers = keepers + 0
batchcol2 = batchcol[which(colSums(resultsbinary) > 0)]
resultsbinary = resultsbinary[,which(colSums(resultsbinary) > 0)]
keepers = keepers[,which(colSums(resultsbinary) > 0)]

binary.phi = matrix(1,dim(resultsbinary)[2],dim(resultsbinary)[2])
for(i in 1:(dim(resultsbinary)[2])-1){
    for(j in (i+1):dim(resultsbinary)[2]){
        currkeepers = intersect(which(keepers[,i]>0),which(keepers[,j]>0))
        currresults = resultsbinary[currkeepers,c(i,j)]
        BC = length(which(resultsbinary[,i] == 1 & resultsbinary[,j] == 1))*length(which(resultsbinary[,i] == 0 & resultsbinary[,j] == 0))
        AD = length(which(resultsbinary[,i] == 1 & resultsbinary[,j] == 0))*length(which(resultsbinary[,i] == 0 & resultsbinary[,j] == 1))
        AB = sum(resultsbinary[,i])
        CD = length(resultsbinary[,i]) - sum(resultsbinary[,i])
        AC = length(resultsbinary[,j]) - sum(resultsbinary[,j])
        BD = sum(resultsbinary[,j])
        phi = (BC - AD)/(sqrt(AB)*sqrt(CD)*sqrt(AC)*sqrt(BD))
        binary.phi[i,j] = phi
        binary.phi[j,i] = phi
    }
}

colnames(binary.phi) = colnames(resultsbinary)
rownames(binary.phi) = colnames(resultsbinary)


ColorRamp <- rgb( seq(0,1,length=256),
	seq(0,1,length=256),
	seq(1,0,length=256))
mypalette = colorRampPalette(c("cornflowerblue", "white", "darkblue"), space = "Lab")

pdf(paste(resultsbin,"OverlapMaps.pdf",sep=""),height=11,width=8.5)
heatmap.2(results.heat.m,trace="none",col=ColorRamp,ColSideColors=batchcol,distfun=function(x) as.dist(1-abs(x)))
#heatmap.2(first.heat.m,trace="none",col=ColorRamp,ColSideColors=subset(batchcol,batchers == "firstbatch"))
#heatmap.2(second.heat.m,trace="none",col=ColorRamp,ColSideColors=subset(batchcol,batchers == "secondbatch"))
#heatmap.2(results.heat.ma,trace="none",col=ColorRamp,ColSideColors=batchcol,main="FET O/L only <5 Factors")
heatmap.2(binary.phi,trace="none",col=mypalette(20),ColSideColors=batchcol2,main="Phi Coefficients (Shared DE)",
          distfun=function(x) as.dist(1-abs(x)))
dev.off()

pdf(paste(resultsbin,"FigS6.pdf",sep=""),height=5.5,width=8.5)
barplot(as.numeric(summar[,2]),names.arg=names.clean,las=2,cex.names=0.8,col=batchcol,ylim=c(0,(max(as.numeric(summar[,2]))+3000)),
        main="Summary of No. Genes Expressed",ylab="No. Genes Expressed",legend.text=c("Batch 1","Batch 2","Batch 3"),
        args.legend = list(fill=c("indianred","dodgerblue2","goldenrod")))
dev.off()