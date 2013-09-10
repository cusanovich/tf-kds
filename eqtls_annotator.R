library('gplots')
library('plyr')
source('./config.R')

outbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/Bindings/"
resultsmatrix = as.matrix(read.table(bindingmatrix,sep="\t"))
factors = as.matrix(read.table(factored))
nobfactors = as.matrix(read.table(nobfactored))
#resultsbin = "~/home/Kd_Arrays/Analysis/Results/RUV2_NSAveraged_alt_Results/"
#outbin = "~/home/Kd_Arrays/CombinedBinding/Results/"
#resultsmatrix = as.matrix(read.table("~/home/Kd_Arrays/CombinedBinding/Binding/allunionbindingresults10kb.txt",sep="\t"))
#factors = as.matrix(read.table("~/home/Kd_Arrays/CombinedBinding/Annotations/allbinding_list.txt"))
#nobfactors = as.matrix(read.table("~/home/Kd_Arrays/CombinedBinding/Annotations/nobinding_list.txt"))
de_threshold = 0.05

#Create a list of all knockdown experiments from teh resultsbin
all.pvals <- list.files(path = resultsbin,pattern="Pvalues.txt")
pvals <- llply(paste(resultsbin,all.pvals,sep=""), read.table)
namers = c()
namers.clean = c()
for(i in 1:length(all.pvals)){
  namers[i] = strsplit(all.pvals[i],"batch")[[1]][1]
  namers.clean[i] = strsplit(namers[i],"_")[[1]][1]
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
}
master[[2]] = master[[2]][match(master[[1]][,2],master[[2]][,2]),]
master[[3]] = master[[3]][match(master[[1]][,2],master[[3]][,2]),]

qvals = master[[2]]
inder = match(master[[1]][,2],rownames(resultsmatrix))
newmatrix = resultsmatrix[inder,]
for(i in 1:dim(newmatrix)[2]){
  print(colnames(newmatrix)[i])
  currens = factors[match(colnames(newmatrix)[i],factors[,1]),3]
  if(is.na(currens)){next}
  decols = c(match(colnames(newmatrix)[i],namers.clean)+3,which(qvals[match(currens,qvals[,2]),] < 0.05))
  decols = decols[!is.na(decols)]
  if(length(decols) == 0){next}
  print(paste0("DE in ",length(decols)," experiments..."))
  deandbound = c()
  bound = c()
  unbound = c()
  for(j in 1:length(decols)){
    currcol = decols[j]
    newqvals = qvals[!is.na(qvals[,currcol]),c(2,currcol)]
    newestmatrix = newmatrix[!is.na(qvals[,currcol]),i]
    currdeb = newqvals[newqvals[,2] < 0.05 & newestmatrix > 0,1]
    currb = newqvals[newqvals[,2] >= 0.05 & newestmatrix > 0,1]
    curru = newqvals[newqvals[,2] >= 0.05 & newestmatrix == 0,1]
    deandbound = union(deandbound, currdeb)
    bound = union(bound, currb)
    unbound = union(unbound, curru)
  }
  write.table(paste(deandbound,collapse="|"),paste0(outbin,colnames(newmatrix)[i],".DEandBound.txt"),row.names=F,col.names=F,quote=F)
  write.table(paste(bound,collapse="|"),paste0(outbin,colnames(newmatrix)[i],".Bound.txt"),row.names=F,col.names=F,quote=F)
  write.table(paste(unbound,collapse="|"),paste0(outbin,colnames(newmatrix)[i],".Unbound.txt"),row.names=F,col.names=F,quote=F)
}