library('plyr')
library('stringr')
library('qvalue')

windowname = "10kb"
resultsbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/"
outdir = paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Results/",windowname,"/")
de_threshold = 0.05
resultsmatrices = list.files(path = paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Matrices/",windowname,"/"),
                             pattern = "_resultsmatrix.txt")
factors = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt")
masterresults = as.matrix(read.table(paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/allbindingresults",windowname,".txt")))

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
  
  for(k in 4:dim(master[[j]])[2]){
    curcol = master[[j]][,k]
    master[[j]][which(is.na(curcol)),k] = 5
  }
}
master[[2]] = master[[2]][match(master[[1]][,2],master[[2]][,2]),]
master[[3]] = master[[3]][match(master[[1]][,2],master[[3]][,2]),]

for(x in 1:length(resultsmatrices)){
  resultsmatrix = as.matrix(read.table(paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Matrices/",windowname,"/",resultsmatrices[x]),sep="\t"))
  chromstate = str_split(resultsmatrices[x],"_resultsmatrix")[[1]][1]
  resultsbinary = resultsmatrix>0
  resultsbinary = resultsbinary + 0
  bindingmatrix = matrix(NA,1,15)
  colnames(bindingmatrix) = c("TF","MotifSet","ChromatinState","BoundGenes",
                              "StateBoundGenes","BoundDEGenes","StateDEBoundGenes","FET",
                              "FETDepletion","BoundPlus Genes","State BoundPlus Genes",
                              "BoundPlus DE Genes","State DE BoundPlus Genes","FETPlus",
                              "FETPlusDepletion")
  bindinglist = list()
  rbindinglist = list()
  winnertfs = list()
  mastermatch = c()
  for(i in 1:length(names(master[[2]]))){
    mastermatch[i] = strsplit(names(master[[2]])[i],"_")[[1]][1]
  }
  listgene = intersect(colnames(resultsmatrix),mastermatch)
  #sp1.ind = match("SP1",listgene)
  #listgene = c(listgene[1:sp1.ind],"SP1",listgene[(sp1.ind+1):length(listgene)])
  #sp1 = 0
  dtfs = matrix(0,61,4)
  dtf.names = c()
  print(chromstate)
  for(i in 4:dim(master[[2]])[2]){
    currgene = names(master[[2]])[i]
    matchgene = strsplit(currgene,"_")[[1]][1]
    commongenes = intersect(unique(master[[2]][which(master[[2]][,i]< 5),2]),
                          rownames(resultsmatrix))
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
    currmasterresults = masterresults[match(commongenes,rownames(masterresults)),]
    matrixcol = grep(matchgene,colnames(currresultsmatrix))
    mastercol = grep(matchgene,colnames(currmasterresults))
    downstreamcol = matrixcol
    downstreammastercol = mastercol
    for(k in 1:length(winnerfactornames)){
      downstreamcol = c(downstreamcol,match(winnerfactornames[k],colnames(currresultsmatrix)))
      downstreammastercol = c(downstreammastercol,match(winnerfactornames[k],colnames(currmasterresults)))
    }
    if(length(winnerfactornames) == 0){
      downstreamcol = downstreamcol[1]
      downstreammastercol = downstreammastercol[1]
    }
    if(length(matrixcol) < 1 & length(downstreamcol) < 1){
      next
    }
    if(length(matrixcol) > 1){
      print(paste0("Whoa! ",matchgene,
                 " matches too many columns in the binding matrix!!!"))
      next
    }
    masterbound = rownames(currmasterresults)[which(currmasterresults[,mastercol] > 0)]
    downstreammasterbound = rownames(currmasterresults)[which(rowSums(as.matrix(currmasterresults[,downstreammastercol])) > 0)]
    boundgenes = rownames(currresultsmatrix)[which(currresultsmatrix[,matrixcol] > 0)]
    downstreambound = rownames(currresultsmatrix)[which(rowSums(as.matrix(currresultsmatrix[,downstreamcol])) > 0)]
    statebound = intersect(masterbound,boundgenes)
    downstreamstatebound = intersect(downstreammasterbound,downstreambound)
    winnerbound = intersect(winners,masterbound)
    winnerdownstreambound = intersect(winners,downstreammasterbound)
    if(length(masterbound) == 0 & length(downstreammasterbound) == 0){
      print(paste0("Whoa! ",matchgene," doesn't bind anything!!!"))
      next
    }
    q = length(intersect(statebound,winnerbound)) - 1
    if(length(statebound) == 0 | q < 0){
      q = 0
    }
    q.down = length(intersect(downstreamstatebound,winnerdownstreambound)) - 1
    if(length(downstreamstatebound) == 0 | q.down < 0){
      q.down = 0
    }
    m = length(winnerbound)
    n = length(masterbound) - m
    k = length(statebound)
    p = phyper(q,m,n,k, lower.tail=FALSE)
    p.dep = phyper(q+1,m,n,k)
    m.down = length(winnerdownstreambound)
    n.down = length(downstreammasterbound) - m.down
    k.down = length(downstreamstatebound)
    p.down = phyper(q.down,m.down,n.down,k.down, lower.tail=FALSE)
    p.down.dep = phyper(q.down+1,m.down,n.down,k.down)
    if(q == 0){
      p = 1
    }
    if(q.down == 0){
      p.down = 1
    }
    newline = c(currgene,matchgene,chromstate,(n+m),k,m,ifelse(q>0,q+1,q),p,p.dep,
                (n.down+m.down),k.down,m.down,ifelse(q.down>0,q.down+1,q.down),p.down,p.down.dep)
    bindingmatrix = rbind(bindingmatrix,newline)
  }
  bindingmatrix = bindingmatrix[-1,]

  write.table(bindingmatrix,paste0(outdir,chromstate,".overlaptable.txt"),
              row.names=F,quote=F,sep="\t")
}
overlaps = list.files(path=outdir,pattern=".overlaptable.txt")
for(j in 1:length(overlaps)){
  if(j == 1){
    masterol = read.table(paste0(outdir,overlaps[j]),sep="\t",header=T)
  }else{
    currol = read.table(paste0(outdir,overlaps[j]),sep="\t",header=T)
    masterol = rbind(masterol,currol)
  }
}
write.table(masterol,paste0(outdir,"overlaptable.txt"),row.names=F,quote=F,sep="\t")

#Cleaning up my tables a bit - want to remove any comparisons that just don't have binding
clean_and_q = function(x,colnums=c(1:7),filternums=4,thresh=0,pcol=8){
  cleaner = x[,colnums]
  if(length(filternums) > 1){
    cleaner = cleaner[which(rowSums(cleaner[,filternums]) > thresh),]
  }else{
    cleaner = cleaner[which(cleaner[,filternums] > thresh),]
  }
  fishing = apply(cleaner,1,fisherian)
  y <- sapply(fishing, rbind)
  fishing <- matrix(unlist(t(y)), ncol=2)
  fishing[,2] = ifelse(fishing[,2] > 1,1,fishing[,2])
  colnames(fishing) = c("Expected","FET")
  Qvalues = qvalue(fishing[,2])$qvalues
  cleaner = cbind(cleaner,fishing)
  cleaner = cbind(cleaner,Qvalues)
  cleaner = cleaner[order(cleaner$FET),]
  cleaner$UporDown = "Enriched"
  cleaner$Enriched = ifelse(cleaner[,7] > cleaner[,8],"Enriched","Depleted")
  return(cleaner)
}

fisherian = function(x,universe=4,overlap=7,group1=5,group2=6){
  fishers = as.numeric(c(x[overlap],(as.numeric(x[group1])-as.numeric(x[overlap])),
                         (as.numeric(x[group2])-as.numeric(x[overlap])),
                         (as.numeric(x[universe])+as.numeric(x[overlap])-as.numeric(x[group1])-as.numeric(x[group2]))))
  testers = matrix(fishers,2,2)
  group1prob = as.numeric(x[group1])/as.numeric(x[universe])
  group2prob = as.numeric(x[group2])/as.numeric(x[universe])
  return(list(expectation=as.numeric(group1prob*group2prob*as.numeric(x[universe])),
              fishers=as.numeric(fisher.test(testers)$p.value)))
}
enrichment = clean_and_q(masterol)
write.table(enrichment,paste0(outdir,"enrichmenttable.txt"),row.names=F,quote=F,sep="\t")
enrichmentplus = clean_and_q(masterol,colnums=c(1:3,10:13))
write.table(enrichmentplus,paste0(outdir,"enrichmentplustable.txt"),row.names=F,quote=F,sep="\t")