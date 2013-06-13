commoners = intersect(master[[2]][,2],rownames(resultsmatrix))
masterind = match(commoners,master[[2]][,2])
#resultsind = match(commoners,rownames(resultsmatrix))
comaster = master[[2]][masterind,]
#comresults = resultsmatrix[resultsind,]
de1 = grep("IRF3",colnames(commaster))
de2 = grep("RXRA",names(commaster))
commonexpr = comaster[which(comaster[,de1] < 5 & comaster[,de2] < 5),2]
deind = match(commonexpr,comaster[,2])
resultsind = match(commonexpr,rownames(resultsmatrix))
commaster = comaster[deind,]
comresults = resultsmatrix[resultsind,]

#bindind = match(commonexpr,rownames(comresults))
bind1 = grep("IRF3",colnames(comresults))
bind2 = grep("RXRA",colnames(comresults))

winners1 = which(commaster[,de1] < 0.05)
winners2 = which(commaster[,de2] < 0.05)
cowinners = intersect(winners1,winners2)
colosers = c(setdiff(winners1,winners2),setdiff(winners2,winners1))
irf3losers = setdiff(winners1,winners2)
rxralosers = setdiff(winners2,winners1)
irf3bound = sum(comresults[,bind1])
rxrabound = sum(comresults[,bind2])
eitherbound = rowSums(comresults[,c(bind1,bind2)])

#7714 genes expressed
#97 IRF3-DE
#219 RXRA-DE
#7 O/L
#2600 bound by either factor
#106 colosers are bound

#q = number of genes in o/l - 1
#m = number of winners
#n = number of non-winners
#k = number of genes picked
q=106-1
m=2600
n=5114
k=302
phyper(q,m,n,k, lower.tail=FALSE)

q=34-1
m=228
n=7486
k=90
phyper(q,m,n,k, lower.tail=FALSE)

q=72-1
m=2482
n=5232
k=212
phyper(q,m,n,k, lower.tail=FALSE)

q=7-1
m=219
n=7495
k=97
phyper(q,m,n,k, lower.tail=FALSE)

diffcheckups = list()
olapcheckups = list()
x=1
y=1
z=0
for(i in 1:(length(winnertfs)-1)){
  for(j in (i+1):length(winnertfs)){
    z = z+1
    tfs1 = winnertfs[[i]]
    tfs2 = winnertfs[[j]]
    commons = intersect(tfs1,tfs2)
    differs = c(setdiff(tfs1,tfs2),setdiff(tfs2,tfs1))
    overall = union(tfs1,tfs2)
    if(length(commons) > 0 & length(commons) < 5){
      olapcheckups[[(2*x-1)]] = commons
      olapcheckups[[(2*x)]] = differs
      names(olapcheckups)[c((2*x-1),(2*x))] = paste(names(winnertfs[i]),names(winnertfs[j]),sep="+")
      x = x+1
    }
    if(length(overall) > 0 & length(differs) < 6){
      diffcheckups[[(2*y-1)]] = commons
      diffcheckups[[(2*y)]] = differs
      names(diffcheckups)[c((2*y-1),(2*y))] = paste(names(winnertfs[i]),names(winnertfs[j]),sep="+")
      y = y+1
    }
  }
}