library('gplots')
library('plyr')
source('./config.R')

outbin = "/mnt/lustre/home/cusanovich/Kd_Arrays/Supp/"
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
}
master[[2]] = master[[2]][match(master[[1]][,2],master[[2]][,2]),]
master[[3]] = master[[3]][match(master[[1]][,2],master[[3]][,2]),]

newmastercolnames = c()
for(i in 1:dim(master[[1]])[2]){
  currfac = strsplit(colnames(master[[1]])[i],"_")[[1]][1]
  if(i > 3){
    currens = factors[match(currfac,factors[,1]),3]
    if(is.na(currens)){
      currens = nobfactors[match(currfac,nobfactors[,1]),2]
    }
    currfac = paste0(currfac,"_",currens,"_DE_p")
  }
  newmastercolnames[i] = currfac
}
newbindingcolnames = c()
for(i in 1:dim(resultsmatrix)[2]){
  currfac = colnames(resultsmatrix)[i]
  if(currfac == "NKX2.5"){
    currfac = "NKX2-5"
  }
  currens = factors[match(currfac,factors[,1]),3]
  currfac = paste0(currfac,"_",currens,"_binding")
  newbindingcolnames[i] = currfac
}
combops = master[[1]]
colnames(combops) = newmastercolnames
combobs = resultsmatrix
colnames(combobs) = newbindingcolnames
inder = match(combops[,2],rownames(combobs))
tables2 = cbind(combops,combobs[inder,])
write.table(tables2,paste0(outbin,"TableS2.txt"),row.names=F,quote=F,sep="\t")