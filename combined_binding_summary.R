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
degenes = c()
egenes = c()
for(i in 1:dim(newmatrix)[2]){
  print(colnames(newmatrix)[i])
  currens = factors[match(colnames(newmatrix)[i],factors[,1]),3]
  if(is.na(currens)){next}
  ecols = which(qvals[match(currens,qvals[,2]),] < 5)
  decols = which(qvals[match(currens,qvals[,2]),] < 0.05)
  ecols = ecols[!is.na(ecols)]
  decols = decols[!is.na(decols)]
  if(length(ecols) == 0){next}
  print(paste0("Expressed in ",length(ecols)," experiments..."))
  print(paste0("DE in ",length(decols)," experiments..."))
  egenes = c(egenes,colnames(newmatrix)[i])
  if(length(decols) > 0){
    degenes = c(degenes,colnames(newmatrix)[i])
  }
}

degenes = union(namers.clean,degenes)
indy = match(degenes,colnames(newmatrix))
indy = indy[!is.na(indy)]
newestermatrix = newmatrix[,indy]

range(colSums(newestermatrix>0))
median(colSums(newestermatrix>0))
pdf(paste0(outbin,"FigS7.pdf"),height=10,width=5)
par(mfrow=c(2,1))
hist(colSums(newestermatrix>0),col="dodgerblue2",las=1,main="A. Number of Genes Bound by Each Factor",xlab="Number of Genes")
hist(rowSums(newestermatrix>0),col="dodgerblue2",las=1,main="B. Number of Factors Bound by Each Gene",xlab="Number of Factors")
dev.off()
length(which(rowSums(newestermatrix)>0)) - dim(newmatrix)[1]
range(rowSums(newestermatrix[rowSums(newestermatrix)>0,]>0))
median(rowSums(newestermatrix[rowSums(newestermatrix)>0,]>0))