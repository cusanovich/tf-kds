library('plyr')
library('stringr')
library('qvalue')

windowname = "10kb"
resultsbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/"
outdir = paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Results/",windowname,"/")
de_threshold = 0.05
resultsmatrix = read.table(paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/ChromatinStates/Matrices/",windowname,"/Strong_Enhancer_resultsmatrix.txt"),sep="\t")
factors = read.table("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt")
masterresults = as.matrix(read.table(paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Binding/allunionbindingresults",windowname,".txt")))

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


resultsbinary = resultsmatrix>0
resultsbinary = resultsbinary + 0

bindinglist = list()
rbindinglist = list()
winnertfs = list()
mastermatch = c()
for(i in 1:length(names(master[[2]]))){
  mastermatch[i] = strsplit(names(master[[2]])[i],"_")[[1]][1]
}
listgene = intersect(colnames(resultsmatrix),mastermatch)

bindingmatrix = matrix(NA,1,15)
colnames(bindingmatrix) = c("TF","MotifSet","ChromatinState","BoundGenes",
                            "StateBoundGenes","BoundDEGenes","StateDEBoundGenes","FET",
                            "FETDepletion","BoundPlus Genes","State BoundPlus Genes",
                            "BoundPlus DE Genes","State DE BoundPlus Genes","FETPlus",
                            "FETPlusDepletion")
dtfs = matrix(0,61,4)
dtf.names = c()
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
    downstreamcol = downstreamcol[!is.na(downstreamcol)]
    downstreammastercol = downstreammastercol[!is.na(downstreammastercol)]
  }
  if(length(matrixcol) < 1 & length(downstreamcol) < 1){
    print(paste0("Whoa! ",matchgene," doesn't bind anything!!!"))
    next
  }
  if(length(matrixcol) > 1){
    print(paste0("Whoa! ",matchgene,
               " matches too many columns in the binding matrix!!!"))
    next
  }
  masterbound = rownames(currmasterresults)[which(currmasterresults[,mastercol] > 0)]
  downstreammasterbound = rownames(currmasterresults)[which(rowSums(as.matrix(currmasterresults[,downstreammastercol])) > 0)]
  if(length(matrixcol) > 0){
    boundgenes = rownames(currresultsmatrix[currresultsmatrix[,matrixcol] > 0,])
  }else{
    boundgenes = NA
  }
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
  newline = c(currgene,matchgene,"Strong_Enhancer",(n+m),k,m,ifelse(q>0,q+1,q),p,p.dep,
              (n.down+m.down),k.down,m.down,ifelse(q.down>0,q.down+1,q.down),p.down,p.down.dep)
  bindingmatrix = rbind(bindingmatrix,newline)
}
bindingmatrix = bindingmatrix[-1,]
write.table(bindingmatrix,"./strongenhancer_overlap.txt",row.names=F,sep="\t",quote=F)