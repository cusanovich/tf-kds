library('gplots')
library('plyr')
source('./config.R')

outbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Results/"
resultsmatrix = as.matrix(read.table(bindingmatrix,sep="\t"))
binary.phi = as.matrix(read.table(binarytable))
factors = as.matrix(read.table(factored))
nobfactors = as.matrix(read.table(nobfactored))
censusfactors = as.matrix(read.table(censused,sep="\t",fill=NA,header=T))
#resultsbin = "~/home/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/"
#outbin = "~/home/Kd_Arrays/CombinedBinding/Results/"
#resultsmatrix = as.matrix(read.table("~/home/Kd_Arrays/CombinedBinding/Binding/allbindingresults10kb.txt",sep="\t"))
#binary.phi = as.matrix(read.table("~/home/Kd_Arrays/CombinedBinding/PhiTables/AllFactorBinding10kbPhis.txt"))
#factors = as.matrix(read.table("~/home/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt"))
#nobfactors = as.matrix(read.table("~/home/Kd_Arrays/CombinedBinding/Annotations/nobinding_list.txt"))
#censusfactors = as.matrix(read.table("~/home/Kd_Arrays/CombinedBinding/Annotations/TFcensus.txt",sep="\t",fill=NA,header=T))
outplots = paste0("union_factorresultsplots_RUV2_NSAveraged_",windowname,".pdf")
outtables = paste0("union_RUV2_NSAveraged_FACTOR_",windowname,"BindingWindow_")
olapplots = paste0("union_overlapplots_RUV2_NSAveraged_",windowname,".pdf")
de_threshold = 0.05
relaxed_threshold = 0.2
resultsbinary = resultsmatrix
censusfactors = censusfactors[,c(1,2,6)]
censusfactors = censusfactors[censusfactors[,1] != 'x',]
censusfactors = censusfactors[censusfactors[,1] != 'c',]
censusensg = setdiff(censusfactors[,2],nobfactors[,2])
censusensg = setdiff(censusensg,factors[,3])
censusfactors = censusfactors[match(censusensg,censusfactors[,2]),]
nobbers = rbind(censusfactors[,c(3,2)],nobfactors)

#Create a list of all knockdown experiments from teh resultsbin
all.pvals <- list.files(path = resultsbin,pattern="Pvalues.txt")
pvals <- llply(paste(resultsbin,all.pvals,sep=""), read.table)
namers = c()
for(i in 1:length(all.pvals)){
  namers[i] = strsplit(all.pvals[i],"batch")[[1]][1]
}
names(pvals) = namers

#Build a master list of matrices that contains all Fold-changes, P-values,
#and Q-values for each knockdown.  Genes DE are coded with a 5 instead of 'NA'.
master = list()
for(j in 1:3){
  master[[j]] = pvals[[1]][,c(4:6,j)]
  names(master)[j] = colnames(pvals[[1]])[j]
  names(master[[j]])[4] = namers[1]
  for(i in 2:length(namers)){
    #if(i == 22 | i == 46){next}
    newbies = setdiff(as.character(pvals[[i]][,4]),master[[j]][,1])
    master[[j]] = merge(master[[j]],pvals[[i]][,c(4,j)],by="ProbeID",all=T)
    master[[j]][,1] = as.character(master[[j]][,1])
    master[[j]][,2] = as.character(master[[j]][,2])
    master[[j]][,3] = as.character(master[[j]][,3])
    names(master[[j]])[dim(master[[j]])[2]] = namers[i]
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

nobbers = nobbers[nobbers[,2] %in% master[[2]][,2],]
winning = unique(c(factors[,3],nobbers[,2]))


#Plot the gene binding correlations
pdf(paste0(outbin,outplots),height=11,width=8.5)
heatmap.2(binary.phi,trace="none",main="Binding Phi Correlations",
          distfun=function(x) as.dist(1-abs(x)),cexRow=0.6)
par(mfrow=c(3,2))

#The meat: this does all the overlap calculations and binding measures
bindingmatrix = matrix(NA,1,14)
colnames(bindingmatrix) = c("TF","Motif Set","Expressed Genes","Bound Genes",
                            "DE Genes","Relaxed DE Genes","Bound DE Genes",
                            "Relaxed Bound DE Genes","FET","Relaxed FET",
                            "Bound+ Genes","Bound+ DE Genes",
                            "Relaxed Bound+ DE Genes","FET+")
bindinglist = list()
winnertfs = list()
winnernobbers = list()
mastermatch = c()
for(i in 4:length(names(master[[2]]))){
  mastermatch[i-3] = strsplit(names(master[[2]])[i],"_")[[1]][1]
}
dtfs = matrix(0,59,3)
dtf.names = c()
boundps = list()
sampleps = list()
for(i in 4:dim(master[[2]])[2]){
  currgene = names(master[[2]])[i]
  matchgene = strsplit(currgene,"_")[[1]][1]
  commongenes = intersect(master[[2]][which(master[[2]][,i]< 5),2],rownames(resultsmatrix))
  commonind = match(commongenes,master[[2]][,2])
  currmaster = master[[2]][commonind,c(2,i)]
  winners = unique(currmaster[which(currmaster[,2]<= de_threshold),1])
  winnerfactors = intersect(winners,factors[,3])
  dtfs[i-3,1] = length(winners)
  dtfs[i-3,2] = length(winnerfactors)
  dtfs[i-3,3] = length(intersect(winners,winning))
  dtf.names[i-3] = matchgene
  winnerfactors = setdiff(winnerfactors,factors[match(matchgene,factors[,1]),3])
  winnerfactornames = as.character(factors[match(winnerfactors,factors[,3]),1])
  winnertfs[[i-3]] = c(matchgene,winnerfactornames)
  names(winnertfs)[i-3] = currgene
  winnernobbers[[i-3]] = intersect(winners,winning)
  names(winnernobbers)[i-3] = currgene
  currresultsmatrix = resultsmatrix[match(commongenes,rownames(resultsmatrix)),]
  matrixcol = grep(matchgene,colnames(currresultsmatrix))
  downstreamcol = matrixcol
  for(k in 1:length(winnerfactornames)){
    downstreamcol = c(downstreamcol,match(winnerfactornames[k],colnames(currresultsmatrix)))
  }
  if(length(winnerfactornames) == 0){
    downstreamcol = downstreamcol[1]
  }
  if(length(matrixcol) < 1 & length(downstreamcol) < 1){
    next
  }
  if(length(matrixcol) > 1){
    print(paste0("Whoa! ",matchgene," matches too many columns in the binding matrix!!!"))
    next
  }
  boundgenes = rownames(currresultsmatrix)[which(currresultsmatrix[,matrixcol] > 0)]
  downstreambound = rownames(currresultsmatrix)[which(rowSums(as.matrix(currresultsmatrix[,downstreamcol])) > 0)]
  if(length(boundgenes) == 0 & length(downstreambound) == 0){
    print(paste0("Whoa! ",matchgene," doesn't bind anything!!!"))
    next
  }
  q = length(intersect(winners,boundgenes)) - 1
  if(length(boundgenes) == 0){
    q = 0
  }
  m = length(winners)
  q.down = length(intersect(winners,downstreambound)) - 1
  relaxedwinners = currmaster[which(currmaster[,2]<=relaxed_threshold),1]
  m.relaxed = length(relaxedwinners)
  n = length(unique(currmaster[,1])) - m
  n.relaxed = length(unique(currmaster[,1])) - m.relaxed
  k = length(boundgenes)
  k.down = length(downstreambound)
  q.relaxed = length(intersect(relaxedwinners,boundgenes)) - 1
  p = phyper(q,m,n,k, lower.tail=FALSE)
  if(q == 0){
    p = 1
  }
  p.relaxed = phyper(q.relaxed,m.relaxed,n.relaxed,k, lower.tail=FALSE)
  p.down = phyper(q.down,m,n,k.down,lower.tail=FALSE)
  newline = c(currgene,matchgene,(n+m),k,m,m.relaxed,ifelse(q>0,q+1,q),
              (q.relaxed+1),p,p.relaxed,k.down,(q.down+1),
              length(intersect(relaxedwinners,downstreambound)),p.down)
  bindingmatrix = rbind(bindingmatrix,newline)
  listspot = match(currgene,names(master[[2]]))
  j = c((listspot*2-1),(listspot*2))
  currps = master[[1]][commonind,c(2,i)]
  bindind = match(boundgenes,currps[,1])
  downind = match(downstreambound,currps[,1])
  ssize = min(length(downind),(length(commonind) - length(downind)))
  boundps[[i-3]] = -log10(currps[sample(downind,ssize),2])
  sampleps[[i-3]] = -log10(currps[sample(c(1:dim(currps)[1])[-downind],ssize),2])
  top = max(max(boundps[[i-3]],sampleps[[i-3]]))
  dense = density(currps[,2])
  plot(dense,lwd=3,main=paste0(matchgene," P-values"),xlab="P-value")
  legend("topright",legend=c("FDR 0.05","FDR 0.20"),fill=c("black","dodgerblue2"))
  abline(v=max(currps[which(currmaster[,2] <= de_threshold),2]),lty=2,lwd=2)
  abline(v=max(currps[which(currmaster[,2] <= relaxed_threshold),2]),lty=2,lwd=2,col="dodgerblue2")
  plot(sort(-log10(ppoints(ssize))),sort(boundps[[i-3]]),
       main=paste0(matchgene," QQ Plot"),pch=20,col="indianred",ylim=c(0,top),
       xlab="Expected -log10(P-value)",ylab="Observed -log10(P-value)")
  legend("topleft",c("Bound Genes","Unbound Genes"),fill=c("indianred","black"))
  abline(a=0,b=1,lwd=2,lty="dashed",col="dodgerblue2")
  points(sort(-log10(ppoints(ssize))),sort(sampleps[[i-3]]),pch=20)
  points(sort(-log10(ppoints(ssize))),sort(boundps[[i-3]]),pch=20,col="indianred")
}

bindingmatrix = bindingmatrix[-1,]
percs = cbind(as.numeric(bindingmatrix[,4])/as.numeric(bindingmatrix[,3]),
              as.numeric(bindingmatrix[,7])/as.numeric(bindingmatrix[,4]),
              as.numeric(bindingmatrix[,8])/as.numeric(bindingmatrix[,4]),
              as.numeric(bindingmatrix[,11])/as.numeric(bindingmatrix[,3]),
              as.numeric(bindingmatrix[,12])/as.numeric(bindingmatrix[,11]),
              as.numeric(bindingmatrix[,13])/as.numeric(bindingmatrix[,11]),
              as.numeric(bindingmatrix[,5])/as.numeric(bindingmatrix[,3]))
colnames(percs) = c("Bound/Total","DE/Bound","RelaxedDE/Bound","BoundPlus/Total",
                    "DEPlus/Bound","RelaxedDEPlus/Bound","DE/Total")

write.table(bindingmatrix,paste0(outbin,outtables,"overlaptable.txt"),row.names=F,quote=F,sep="\t")

uwins = unique(unlist(winnertfs))
kds = names(winnertfs)[1:length(winnertfs)]
tfmatrix = matrix(0,length(uwins),length(kds))
colnames(tfmatrix) = kds
rownames(tfmatrix) = uwins
for(i in 1:length(kds)){
  for(j in 1:length(uwins)){
    currkd = kds[i]
    curruwins = uwins[j]
    currwintfs = match(currkd,names(winnertfs))
    if(curruwins %in% winnertfs[[currwintfs]]){
      tfmatrix[j,i] = 1
    }
  }
}

tfs.phi = matrix(1,dim(tfmatrix)[2],dim(tfmatrix)[2])
for(i in 1:(dim(tfmatrix)[2])-1){
  for(j in (i+1):dim(tfmatrix)[2]){
    BC = length(which(tfmatrix[,i] == 1 & tfmatrix[,j] == 1))*length(which(tfmatrix[,i] == 0 & tfmatrix[,j] == 0))
    AD = length(which(tfmatrix[,i] == 1 & tfmatrix[,j] == 0))*length(which(tfmatrix[,i] == 0 & tfmatrix[,j] == 1))
    AB = sum(tfmatrix[,i])
    CD = length(tfmatrix[,i]) - sum(tfmatrix[,i])
    AC = length(tfmatrix[,j]) - sum(tfmatrix[,j])
    BD = sum(tfmatrix[,j])
    phi = (BC - AD)/(sqrt(AB)*sqrt(CD)*sqrt(AC)*sqrt(BD))
    tfs.phi[i,j] = phi
    tfs.phi[j,i] = phi
  }
}

colnames(tfs.phi) = colnames(tfmatrix)
rownames(tfs.phi) = colnames(tfmatrix)
heatmap.2(tfs.phi,trace="none",main="TF DE Phi Correlations",
          distfun=function(x) as.dist(1-abs(x)))
dev.off()

commoners = intersect(master[[2]][,2],rownames(resultsmatrix))
masterind = match(commoners,master[[2]][,2])
comaster = master[[2]][masterind,]

#causalmatrix = matrix(NA,1830,10)
causalmatrix = matrix(NA,1711,12)
colnames(causalmatrix) = c("Gene1","Gene2","CommonTFsNOB","DisparateTFsNOB",
                           "CommonTFs","DisparateTFs","CommonDE","DisparateDE",
                           "CommonFET","DisparateFET","LeftFET","RightFET")
z=1
for(i in 1:(length(winnertfs)-1)){
  for(j in (i+1):length(winnertfs)){
    newline = c()
    newline[1] = names(winnertfs[i])
    newline[2] = names(winnertfs[j])
    de1 = match(names(winnertfs)[i],colnames(comaster))
    de2 = match(names(winnertfs)[j],colnames(comaster))
    commonexpr = comaster[which(comaster[,de1] < 5 & comaster[,de2] < 5),2]
    deind = match(commonexpr,comaster[,2])
    resultsind = match(commonexpr,rownames(resultsmatrix))
    commaster = comaster[deind,]
    comresults = resultsmatrix[resultsind,]
    tfs1 = winnertfs[[i]]
    tfs1 = intersect(tfs1,colnames(comresults))
    tfs2 = winnertfs[[j]]
    tfs2 = intersect(tfs2,colnames(comresults))
    alltfs1 = winnernobbers[[i]]
    alltfs2 = winnernobbers[[j]]
    commons = intersect(tfs1,tfs2)
    newline[3] = length(intersect(alltfs1,alltfs2))
    newline[4] = length(c(setdiff(alltfs1,alltfs2),setdiff(alltfs2,alltfs1)))
    newline[5] = length(commons)
    differs = c(setdiff(tfs1,tfs2),setdiff(tfs2,tfs1))
    newline[6] = length(differs)
    diffs1 = setdiff(tfs1,tfs2)
    diffs2 = setdiff(tfs2,tfs1)
    commonbind = c()
    bind1 = c()
    bind2 = c()
    for(a in 1:length(commons)){
      commonbind = c(commonbind,match(commons[a],colnames(comresults)))
    }
    for(b in 1:length(diffs1)){
      bind1 = c(bind1,match(diffs1[b],colnames(comresults)))
    }
    for(c in 1:length(diffs2)){
      bind2 = c(bind2,match(diffs2[c],colnames(comresults)))
    }
    winners1 = which(commaster[,de1] < 0.05)
    winners2 = which(commaster[,de2] < 0.05)
    cowinners = intersect(winners1,winners2)
    newline[7] = length(cowinners)
    colosers = c(setdiff(winners1,winners2),setdiff(winners2,winners1))
    newline[8] = length(colosers)
    losers1 = setdiff(winners1,winners2)
    losers2 = setdiff(winners2,winners1)
    commonbound = rowSums(as.matrix(comresults[,commonbind]))
    commonbound = commonbound>0
    commonbound = commonbound + 0
    bound1 = rowSums(as.matrix(comresults[,bind1]))
    bound1 = bound1>0
    bound1 = bound1 + 0
    bound2 = rowSums(as.matrix(comresults[,bind2]))
    bound2 = bound2>0
    bound2 = bound2 + 0
    diffbound = rowSums(as.matrix(comresults[,c(bind1,bind2)]))
    diffbound = diffbound>0
    diffbound = diffbound + 0
    qol = ifelse(is.na(sum(commonbound[cowinners])),0,ifelse(sum(commonbound[cowinners]) > 0,sum(commonbound[cowinners]) - 1,0))
    mol = ifelse(is.na(sum(commonbound)),0,sum(commonbound))
    nol = length(commonbound) - mol
    kol = length(cowinners)
    newline[9] = ifelse(phyper(qol,mol,nol,kol, lower.tail=FALSE)==0,1,phyper(qol,mol,nol,kol, lower.tail=FALSE))
    qdiff = ifelse(sum(diffbound[colosers]) > 0,sum(diffbound[colosers]) - 1,0)
    mdiff = sum(diffbound)
    ndiff = length(diffbound) - mdiff
    kdiff = length(colosers)
    newline[10] = phyper(qdiff,mdiff,ndiff,kdiff, lower.tail=FALSE)
    q1 = ifelse(sum(bound1[losers1]) > 0,sum(bound1[losers1]) - 1,0)
    m1 = sum(bound1)
    n1 = length(bound1) - m1
    k1 = length(losers1)
    newline[11] = phyper(q1,m1,n1,k1, lower.tail=FALSE)
    q2 = ifelse(sum(bound2[losers2]) > 0,sum(bound2[losers2]) - 1,0)
    m2 = sum(bound2)
    n2 = length(bound2) - m2
    k2 = length(losers2)
    newline[12] = phyper(q2,m2,n2,k2, lower.tail=FALSE)
    causalmatrix[z,] = newline
    z = z+1
  }
}
tfperccomm = as.numeric(causalmatrix[,5])/(as.numeric(causalmatrix[,5])+as.numeric(causalmatrix[,6]))
alltfperccomm = as.numeric(causalmatrix[,3])/(as.numeric(causalmatrix[,3])+as.numeric(causalmatrix[,4]))
deperccomm = as.numeric(causalmatrix[,7])/(as.numeric(causalmatrix[,7])+as.numeric(causalmatrix[,8]))
liner = lm(deperccomm ~ tfperccomm)
allliner = lm(deperccomm ~ alltfperccomm)
tfing = as.numeric(causalmatrix[,5])
tfing[which(tfing > 9)] = 10
ys = c()
yi = c()
meds = c()
xs = c(1:10,10:1)
for(i in 1:10){
  ys[i] = quantile(-log10(as.numeric(causalmatrix[tfing == i,9])),probs=0.75)
  yi[i] = quantile(-log10(as.numeric(causalmatrix[tfing == i,9])),probs=0.25)
  meds[i] = median(-log10(as.numeric(causalmatrix[tfing == i,9])))
}
ys = c(ys,yi[10:1])
weirdspots = which(tfperccomm>0.2 & deperccomm<0.15)

pdf(paste0(outbin,olapplots),width=10,height=6)
par(mar=c(5,5,2,2)+0.1)
#par(mfrow=c(2,1))
#par(mar=c(8, 8, 2, 4) + 0.1)
#par(mgp=c(3,1.5,0))
boxplot(percs[,c(4,7,5,6)],ylim=c(0,1),notch=T,col=c("indianred","goldenrod","dodgerblue2","mediumseagreen"),
        outpch=20,outcol=c("indianred","goldenrod","dodgerblue2","mediumseagreen"),
        ylab="Fraction",cex.lab=2,las=1,cex=2,boxlwd=3,medlwd=4,
        names=c("Bound Genes\nAll Genes","DE Genes\nAll Genes","DE Genes (FDR 0.05)\nBound Genes","DE Genes (FDR 0.20)\nBound Genes"))
abline(v=2.5,lty="dashed")
plot(dtfs[,2],dtfs[,1],las=1,cex.lab=2,cex=3,
     xlab="Differentially Expressed Transcription Factors",
     ylab="Differentially Expressed Genes",pch=20,col="dodgerblue2")
abline(lm(dtfs[,1] ~ dtfs[,2]),lwd=4,lty="dashed")
points(dtfs[,2],dtfs[,1],cex=3,pch=20,col="dodgerblue2")
text(10,3500,paste("R^2 = ",round(cor(dtfs[,2],dtfs[,1],
                                      use="pairwise.complete.obs")^2,2)),cex=2)
plot(dtfs[,3],dtfs[,1],las=1,cex.lab=2,cex=3,main="NoBindingTFs Incl.",
     xlab="Differentially Expressed Transcription Factors",
     ylab="Differentially Expressed Genes",pch=20,col="dodgerblue2")
abline(lm(dtfs[,1] ~ dtfs[,3]),lwd=4,lty="dashed")
points(dtfs[,3],dtfs[,1],cex=3,pch=20,col="dodgerblue2")
text(40,3500,paste0("R^2 = ",round(cor(dtfs[,3],dtfs[,1],
                                      use="pairwise.complete.obs")^2,3)),cex=2)

plot(tfperccomm,deperccomm,pch=20,cex.lab=2,las=1,
     xlab="Fraction of TFs DE in Common",
     ylab="Fraction of Genes DE in Common",col="indianred")
abline(liner,lwd=4,lty="dashed")
text(0.1,0.3,
     paste0("R^2 = ",round(cor(tfperccomm,deperccomm,use="pairwise.complete.obs")^2,2)),
     cex=2)
plot(alltfperccomm,deperccomm,pch=20,cex.lab=2,las=1,main="NoBindingTFs Incl.",
     xlab="Fraction of TFs DE in Common",
     ylab="Fraction of Genes DE in Common",col="indianred")
abline(allliner,lwd=4,lty="dashed")
text(0.05,0.3,
     paste0("R^2 = ",round(cor(alltfperccomm,deperccomm,use="pairwise.complete.obs")^2,2)),
     cex=2)

plot(1:10,1:10,type="n",ylim=c(min(ys),max(ys)),las=1,cex.lab=2,xaxt="n",
     xlab = "No. DE TFs in Common",ylab = "Degree of Binding")
axis(side=1,at=1:10,labels=c(1:9,"10+"))
polygon(xs,ys,col="mediumseagreen",border="mediumseagreen")
lines(meds,type="b",lwd=3)

inder = rep(10,length(alltfperccomm[causalmatrix[,3] != 0]))
reorg = order(alltfperccomm[causalmatrix[,3] != 0])
inder[1:(round(length(alltfperccomm[causalmatrix[,3] != 0])/10,0)*9)] = rep(c(1:9),
                                                                            each=round(length(alltfperccomm[causalmatrix[,3] != 0])/10,0))
boxer = cbind(as.factor(inder),
              -log10(as.numeric(causalmatrix[causalmatrix[,3] != 0,9][reorg])))
boxplot(boxer[,2] ~ boxer[,1],notch=T,outline=F,las=1,cex=2,
        xlab="Fraction Common TFs Decile",ylab="Degree of Co-occupancy",cex.lab=2,
        col="mediumseagreen",boxlwd=3,medlwd=4)
par(mfrow=c(1,2))
#par(mar=c(8, 8, 2, 4) + 0.1)
#par(mgp=c(3,1.5,0))
boxplot(percs[,4],ylim=c(0,1),notch=T,col="indianred",outpch=20,outcol="indianred",
        ylab="Fraction",cex.lab=2,las=1,cex=2,boxlwd=3,medlwd=4,names="Bound Genes\nAll Genes")
boxplot(percs[,c(5,6)],ylim=c(0,1),notch=T,col=c("dodgerblue2","mediumseagreen"),
        outpch=20,outcol=c("dodgerblue2","mediumseagreen"),ylab="Fraction",
        cex.lab=2,las=1,cex=2,boxlwd=3,medlwd=4,
        names=c("DE Genes (FDR 0.05)\nBound Genes","DE Genes (FDR 0.20)\nBound Genes"))
#permer = c()
#for(i in 4:dim(master[[2]])[2]){
#  permer[i-3] = length(which(master[[2]][,i] < 5))
#}
#perms = cbind(permer,dtfs)

#commoners = list()
#for(i in 4:dim(master[[2]])[2]){
#  commoners[[i-3]] = intersect(unique(master[[2]][which(master[[2]][,i]< 5),2]),rownames(resultsmatrix))
#}

#winning = unique(c(factors[,3],nobbers[,2]))

#rounds = 1000
#rpermingtf = matrix(NA,1711,rounds)
#rpermingde = matrix(NA,1711,rounds)
#for(t in 1:100){
#for(t in 1:rounds){
#   if(t%%100 == 0){print(t)}
#   holder = matrix(NA,1711,2)
#   s=1
#   for(i in 1:(dim(perms)[1]-1)){
#     for(j in (i+1):dim(perms)[1]){
#       deperms1 = sample(commoners[[i]],perms[i,2])
#       deperms2 = sample(commoners[[j]],perms[j,2])
#       tfs1 = intersect(deperms1,winning)
#       tfs2 = intersect(deperms2,winning)
#       holder[s,1] = length(intersect(deperms1,deperms2))/length(union(deperms1,deperms2))
#       holder[s,2] = length(intersect(tfs1,tfs2))/length(union(tfs1,tfs2))
#       if(length(union(tfs1,tfs2)) == 0){
#         holder[s,2] = 0
#       }
#       s = s+1
#     }
#   }
#   rpermingtf[,t] = holder[,2]
#   rpermingde[,t] = holder[,1]
# }
# 
# meanerstf = apply(rpermingtf,1,mean)
# sderstf = apply(rpermingtf,1,sd)
# differstf = abs(rpermingtf - meanerstf)
# compartf = differstf - abs(alltfperccomm - meanerstf)
# tfps = rowSums(compartf > 0)/dim(rpermingtf)[2]
# tfts = (alltfperccomm - meanerstf)/ifelse(sderstf > 0,sderstf,1)
# tfpermts = apply(rpermingtf,2,function(x){(x - meanerstf)/ifelse(sderstf > 0,sderstf,1)})
# tfpermsorts = apply(tfpermts,2,sort)
# conferlowtfts = apply(tfpermsorts,1,function(x){quantile(x,0.05)})
# conferhitfts = apply(tfpermsorts,1,function(x){quantile(x,0.95)})
# nullts = qt(ppoints(length(alltfperccomm)),df=length(alltfperccomm)-2)
# plot(sort(nullts),sort(tfts),main="All TFs",xlab="Expected",ylab="Observed",
#      pch=20,las=1,cex.lab=2,type="n")
# polygon(c(sort(nullts),sort(nullts,decreasing=T)),c(sort(conferlowtfts),sort(conferhitfts,decreasing=T)),col="gray",border="gray")
# abline(a=0,b=1,lwd=2,col="dodgerblue2",lty="dashed")
# points(sort(nullts),sort(tfts),col="indianred",pch=20)
# 
# meanersde = apply(rpermingde,1,mean)
# sdersde = apply(rpermingde,1,sd)
# differsde = abs(rpermingde - meanersde)
# comparde = differsde - abs(deperccomm - meanersde)
# deps = rowSums(comparde > 0)/dim(rpermingde)[2]
# dets = (deperccomm - meanersde)/sdersde
# depermts = apply(rpermingde,2,function(x){(x - meanersde)/ifelse(sdersde > 0,sdersde,1)})
# depermsorts = apply(depermts,2,sort)
# conferlowdets = apply(depermsorts,1,function(x){quantile(x,0.05)})
# conferhidets = apply(depermsorts,1,function(x){quantile(x,0.95)})
# nulldets = qt(ppoints(length(deperccomm)),df=length(deperccomm)-2)
# plot(sort(nulldets),sort(dets),xlim=c(min(nulldets,conferlowdets),max(nulldets,conferhidets)),
#      ylim=c(min(dets,conferlowdets),max(dets,conferhidets)),main="All DE",xlab="Expected",ylab="Observed",
#      pch=20,las=1,cex.lab=2,type="n")
# polygon(c(sort(nulldets),sort(nulldets,decreasing=T)),c(sort(conferlowdets),sort(conferhidets,decreasing=T)),col="gray",border="gray")
# abline(a=0,b=1,lwd=2,col="dodgerblue2",lty="dashed")
# points(sort(nulldets),sort(dets),col="indianred",pch=20)
# 
# medstf = apply(rpermingtf,2,median)
# diffingtf = abs(medstf - mean(medstf))
# comparatf = diffingtf - abs(median(alltfperccomm) - mean(medstf))
# tfp = length(which(comparatf > 0))/dim(rpermingtf)[2]
# 
# medsde = apply(rpermingde,2,median)
# diffingde = abs(medsde - mean(medsde))
# comparade = diffingde - abs(median(deperccomm) - mean(medsde))
# dep = length(which(comparade > 0))/dim(rpermingde)[2]
# 
# bottom = min(median(alltfperccomm),median(deperccomm),medsde,medstf)
# top = max(median(alltfperccomm),median(deperccomm),medsde,medstf)
# par(mfrow=c(1,1))
# par(mar=c(2, 8, 2, 6) + 0.1)
# par(mgp=c(4,1,0))
# boxplot(medstf,medsde,ylim=c(bottom,top),names=c("All TFs","All DE"),notch=T,
#         outline=F,las=1,cex=2,ylab="Median Degree Sharing",cex.lab=2,
#         col="dodgerblue2",boxlwd=3,medlwd=4)
# points(c(1,2),c(median(alltfperccomm),median(deperccomm)),pch=8,cex=3,lwd=3,
#        col="indianred")
dev.off()