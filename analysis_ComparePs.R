library('plyr')
resultsbin1 = "../Results/All3RUV2K8_NSAveraged_Results2/"
resultsbin2 = "../Results/RUV2_NSAveraged_Results/"
pvals1 <- list.files(path = resultsbin1,pattern="Pvalues.txt")
first <- llply(paste(resultsbin1,pvals1,sep=""), read.table)
namers1 = c()
for(i in 1:length(pvals1)){
    namers1[i] = strsplit(pvals1[i],"_")[[1]][1]
}
names(first) = namers1

pvals2 <- list.files(path = resultsbin2,pattern="Pvalues.txt")
second <- llply(paste(resultsbin2,pvals2,sep=""), read.table)
namers2 = c()
for(i in 1:length(pvals2)){
    namers2[i] = strsplit(pvals2[i],"batch")[[1]][1]
}
names(second) = namers2

#all.pdfs <- list.files(path = resultsbin1,pattern="_Results.pdf")
#batchers = c()
#for(i in 1:length(all.pdfs)){
#batchers[i] = strsplit(all.pdfs[i],"_")[[1]][2]
#}
#batchers=replace(batchers,grep("firstbatch",batchers),"black")
#batchers=replace(batchers,grep("secondbatch",batchers),"gray")


pdf("ComparingPvaluesOldandNew.pdf")
par(mfrow=c(2,2))
for(i in 1:length(pvals1)){
    curr.one = first[[i]]
    twos = grep(namers1[i],namers2)
    if(length(twos) == 0){
        next
    }
    for(j in 1:length(twos)){
        curr.two = second[[twos[j]]]
        two.ind = match(rownames(curr.one),rownames(curr.two))
        curr.cor = round(cor(-log10(curr.two[two.ind,1]),-log10(curr.one[,1]),use="pairwise.complete.obs")^2,2)
        plot(-log10(curr.two[two.ind,1]),-log10(curr.one[,1]),xlab=paste(namers2[twos[j]]," (All3)RUV2+SVA",sep=""),ylab=paste(namers1[i]," RUV2+SVA",sep=""),main=paste(namers1[i],"\nR^2 = ",curr.cor,sep=""),pch=20)
    }
}
dev.off()