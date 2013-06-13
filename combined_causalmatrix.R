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
tfperccomm = causalmatrix[,1]/(causalmatrix[,1]+causalmatrix[,2])
deperccomm = causalmatrix[,3]/(causalmatrix[,3]+causalmatrix[,4])
liner = lm(deperccomm ~ tfperccomm)
tfing = causalmatrix[,1]
tfing[which(tfing > 9)] = 10
ys = c()
yi = c()
meds = c()
xs = c(1:10,10:1)
for(i in 1:10){
  ys[i] = quantile(-log10(causalmatrix[tfing == i,5]),probs=0.75)
  yi[i] = quantile(-log10(causalmatrix[tfing == i,5]),probs=0.25)
  meds[i] = median(-log10(causalmatrix[tfing == i,5]))
}
ys = c(ys,yi[10:1])
weirdspots = which(tfperccomm>0.2 & deperccomm<0.15)

pdf("More TFs equals More DE.pdf",width=9,height=6)
plot(dtfs[,2],dtfs[,1],las=1,xlab="Differentially Expressed Transcription Factors",ylab="Differentially Expressed Genes",pch=20,col="dodgerblue2")
plot(tfperccomm,deperccomm,pch=20,las=1,xlab="Percent of DE TFs in Common",ylab="Percent of DE Genes in Common")
abline(liner,col="tomato",lwd=4)
text(0.05,0.35,paste("R^2 = ",round(cor(tfperccomm,deperccomm,use="pairwise.complete.obs")^2,2)))

plot(1:10,1:10,type="n",ylim=c(min(ys),max(ys)),las=1,xaxt="n",
     xlab = "No. DE TFs in Common",ylab = "Degree of Binding")
axis(side=1,at=1:10,labels=c(1:9,"10+"))
polygon(xs,ys,col="mediumseagreen",border="mediumseagreen")
lines(meds,type="b",lwd=3)
dev.off()