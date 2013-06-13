library('plyr')
resultsbin = "../Results/RUV2_NSAveraged_Results/"
all.pvals <- list.files(path = resultsbin,pattern="Pvalues.txt")
pvals <- llply(paste(resultsbin,all.pvals,sep=""), read.table)

genes = list()
for(i in 1:length(pvals)){
  genes[[i]] = pvals[[i]][,5]
}

exprgenes = sort(unique(as.character(unlist(genes))))

write.table(exprgenes,"../Annotations/uniqueHT12exprgenes.txt",row.names=F,col.names=F,quote=F)