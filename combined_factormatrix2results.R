library('gplots')
library('plyr')
resultsbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_Results/"
outbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Results/"
outplots = "3factorresultsplots_RUV2_NSAveraged_10kb.pdf"
outtables = "3RUV2_NSAveraged_FACTOR_10kbBindingWindow_"
olapplots = "3overlapplots_RUV2_NSAveraged_10kb.pdf"
de_threshold = 0.05
relaxed_threshold = 0.2
resultsmatrix = as.matrix(read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/allbindingresults10kb.txt",sep="\t"))
resultsbinary = resultsmatrix
binary.phi = as.matrix(read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/PhiTables/AllFactorBinding10kbPhis.txt"))
factors = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt")

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
rbindinglist = list()
winnertfs = list()
mastermatch = c()
for(i in 1:length(names(master[[2]]))){
  mastermatch[i] = strsplit(names(master[[2]])[i],"_")[[1]][1]
}
listgene = intersect(colnames(resultsmatrix),mastermatch)
sp1.ind = match("SP1",listgene)
listgene = c(listgene[1:sp1.ind],"SP1",listgene[(sp1.ind+1):length(listgene)])
sp1 = 0
dtfs = matrix(0,61,4)
dtf.names = c()
for(i in 4:dim(master[[2]])[2]){
  currgene = names(master[[2]])[i]
  matchgene = strsplit(currgene,"_")[[1]][1]
  commongenes = intersect(unique(master[[2]][which(master[[2]][,i]< 5),2]),rownames(resultsmatrix))
  commonind = match(commongenes,master[[2]][,2])
  currmaster = master[[2]][commonind,c(2,i)]
  winners = unique(currmaster[which(currmaster[,2]<= de_threshold),1])
  winnerfactors = factors$V3[factors$V3 %in% winners]
  winnerfactors = setdiff(winnerfactors,factors[match(matchgene,factors[,1]),3])
  winnerfactornames = as.character(factors[match(winnerfactors,factors$V3),1])
  dtfs[i-3,1] = length(winners)
  dtfs[i-3,2] = length(winnerfactors)
  dtf.names[i-3] = matchgene
  winnertfs[[i-3]] = c(matchgene,winnerfactornames)
  names(winnertfs)[i-3] = currgene
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
  newline = c(currgene,matchgene,(n+m),k,m,m.relaxed,ifelse(q>0,q+1,q),(q.relaxed+1),p,p.relaxed,k.down,(q.down+1),length(intersect(relaxedwinners,downstreambound)),p.down)
  bindingmatrix = rbind(bindingmatrix,newline)
  listspot = match(currgene,names(master[[2]]))
  j = c((listspot*2-1),(listspot*2))
  currps = master[[1]][commonind,c(2,i)]
  bindind = match(boundgenes,currps[,1])
  dense = density(currps[,2])
  plot(dense,lwd=3,main=paste0(matchgene," P-values"),xlab="P-value")
  legend("topright",legend=c("FDR 0.05","FDR 0.20"),fill=c("black","dodgerblue2"))
  abline(v=max(currps[which(currmaster[,2] <= de_threshold),2]),lty=2,lwd=2)
  abline(v=max(currps[which(currmaster[,2] <= relaxed_threshold),2]),lty=2,lwd=2,col="dodgerblue2")
  bindinglist[[j[1]]] = abs(master[[3]][bindind,i])
  bindinglist[[j[1]]] = bindinglist[[j[1]]][bindinglist[[j[1]]] != 5]
  bindinglist[[j[2]]] = abs(master[[3]][-bindind,i])
  bindinglist[[j[2]]] = bindinglist[[j[2]]][bindinglist[[j[2]]] != 5]
  signifier = "-"
  if(length(bindind) != 0){
    signifier = ifelse(t.test(bindinglist[[j[1]]],bindinglist[[j[2]]])$p.value < 0.05,"*","-")
  }
  names(bindinglist)[j] = signifier
}
bindingmatrix = bindingmatrix[-1,]
percs = cbind(as.numeric(bindingmatrix[,4])/as.numeric(bindingmatrix[,3]),
              as.numeric(bindingmatrix[,7])/as.numeric(bindingmatrix[,4]),
              as.numeric(bindingmatrix[,8])/as.numeric(bindingmatrix[,4]),
              as.numeric(bindingmatrix[,11])/as.numeric(bindingmatrix[,3]),
              as.numeric(bindingmatrix[,12])/as.numeric(bindingmatrix[,11]),
              as.numeric(bindingmatrix[,13])/as.numeric(bindingmatrix[,11]))
colnames(percs) = c("Bound/Total","DE/Bound","RelaxedDE/Bound","BoundPlus/Total",
                    "DEPlus/Bound","RelaxedDEPlus/Bound")
par(mfrow=c(2,1))
# boxplot(percs,ylim=c(0,1),notch=T,col=c("tomato","dodgerblue2","mediumseagreen"),
#         outpch=20,outcol=c("tomato","dodgerblue2","mediumseagreen"),ylab="Percent",
#         names=c("Genes\nBound","Bound Genes DE\n(FDR 0.05)",
#                 "Bound Genes DE\n(FDR 0.20)","Genes\nBound",
#                 "Bound Genes DE\n(FDR 0.05)","Bound Genes DE\n(FDR 0.20)"),las=2)
# abline(v=3.5,lty="dashed")

par(mar=c(8, 8, 2, 4) + 0.1)
par(mgp=c(3,1.5,0))
boxplot(percs[,4:6],ylim=c(0,1),notch=T,col=c("tomato","dodgerblue2","mediumseagreen"),
        outpch=20,outcol=c("tomato","dodgerblue2","mediumseagreen"),ylab="Fraction",
        cex.lab=2,las=1,names=c("Bound Genes\nAll Genes","DE Genes (FDR 0.05)\nBound Genes","DE Genes (FDR 0.20)\nBound Genes"))
abline(v=1.5,lty="dashed")
#mtext(adj=0.5,side=1,text=c("Bound Genes\nAll Genes","DE Genes (FDR 0.05)\nBound Genes","DE Genes (FDR 0.20)\nBound Genes"),at=1:3)

write.table(bindingmatrix,paste0(outtables,"overlaptable.txt"),row.names=F,quote=F,sep="\t")

uwins = unique(unlist(winnertfs))
kds = names(winnertfs)[4:length(winnertfs)]
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

causalmatrix = matrix(NA,1830,10)
colnames(causalmatrix) = c("Gene1","Gene2","CommonTFs","DisparateTFs","CommonDE","DisparateDE","CommonFET","DisparateFET","LeftFET","RightFET")
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
    commons = intersect(tfs1,tfs2)
    newline[3] = length(commons)
    differs = c(setdiff(tfs1,tfs2),setdiff(tfs2,tfs1))
    newline[4] = length(differs)
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
    newline[5] = length(cowinners)
    colosers = c(setdiff(winners1,winners2),setdiff(winners2,winners1))
    newline[6] = length(colosers)
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
    qol = ifelse(sum(commonbound[cowinners]) > 0,sum(commonbound[cowinners]) - 1,0)
    mol = sum(commonbound)
    nol = length(commonbound) - mol
    kol = length(cowinners)
    newline[7] = phyper(qol,mol,nol,kol, lower.tail=FALSE)
    qdiff = ifelse(sum(diffbound[colosers]) > 0,sum(diffbound[colosers]) - 1,0)
    mdiff = sum(diffbound)
    ndiff = length(diffbound) - mdiff
    kdiff = length(colosers)
    newline[8] = phyper(qdiff,mdiff,ndiff,kdiff, lower.tail=FALSE)
    q1 = ifelse(sum(bound1[losers1]) > 0,sum(bound1[losers1]) - 1,0)
    m1 = sum(bound1)
    n1 = length(bound1) - m1
    k1 = length(losers1)
    newline[9] = phyper(q1,m1,n1,k1, lower.tail=FALSE)
    q2 = ifelse(sum(bound2[losers2]) > 0,sum(bound2[losers2]) - 1,0)
    m2 = sum(bound2)
    n2 = length(bound2) - m2
    k2 = length(losers2)
    newline[10] = phyper(q2,m2,n2,k2, lower.tail=FALSE)
    causalmatrix[z,] = newline
    z = z+1
  }
}
tfperccomm = as.numeric(causalmatrix[,3])/(as.numeric(causalmatrix[,3])+as.numeric(causalmatrix[,4]))
deperccomm = as.numeric(causalmatrix[,5])/(as.numeric(causalmatrix[,5])+as.numeric(causalmatrix[,6]))
liner = lm(deperccomm ~ tfperccomm)
tfing = as.numeric(causalmatrix[,3])
tfing[which(tfing > 9)] = 10
ys = c()
yi = c()
meds = c()
xs = c(1:10,10:1)
for(i in 1:10){
  ys[i] = quantile(-log10(as.numeric(causalmatrix[tfing == i,7])),probs=0.75)
  yi[i] = quantile(-log10(as.numeric(causalmatrix[tfing == i,7])),probs=0.25)
  meds[i] = median(-log10(as.numeric(causalmatrix[tfing == i,7])))
}
ys = c(ys,yi[10:1])
weirdspots = which(tfperccomm>0.2 & deperccomm<0.15)

pdf(paste0(outbin,olapplots),width=10,height=6)
par(mar=c(5,5,2,2)+0.1)
plot(dtfs[,2],dtfs[,1],las=1,cex.lab=2,cex=2,xlab="Differentially Expressed Transcription Factors",ylab="Differentially Expressed Genes",pch=20,col="dodgerblue2")
abline(lm(dtfs[,1] ~ dtfs[,2]),lwd=4,lty="dashed")
text(5,3500,paste("R^2 = ",round(cor(dtfs[,2],dtfs[,1],use="pairwise.complete.obs")^2,2)))
plot(tfperccomm,deperccomm,pch=20,cex.lab=2,las=1,xlab="Fraction of DE TFs in Common",ylab="Fraction of DE Genes in Common")
abline(liner,col="tomato",lwd=4)
text(0.05,0.35,paste("R^2 = ",round(cor(tfperccomm,deperccomm,use="pairwise.complete.obs")^2,2)))

plot(1:10,1:10,type="n",ylim=c(min(ys),max(ys)),las=1,cex.lab=2,xaxt="n",
     xlab = "No. DE TFs in Common",ylab = "Degree of Binding")
axis(side=1,at=1:10,labels=c(1:9,"10+"))
polygon(xs,ys,col="mediumseagreen",border="mediumseagreen")
lines(meds,type="b",lwd=3)
dev.off()