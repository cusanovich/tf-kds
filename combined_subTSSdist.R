#read in bound genes
#read in bound and DE genes
#determine subsample size
#read in perm file
#subsample 1000 times
#plot polygon vs. real data
#sample bound genes, grep genes from bound genes, pull distances, unlist, hist, make matrix with bin counts,take 5th and 95th quantile for each bin

library('plyr')
library('beanplot')
source('./config.R')

currpath = paste0(grandpath,"GrepDistance")
altpath = "/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Grepers/Union/10kb"
#altpath = "~/home/Kd_Arrays/GenomeAnnotations/Grepers/Union/10kb"
currdownfiles = list.files(path=currpath,pattern = ".down.permadistance")
currdefiles = gsub("down.permadistance","DEandBoundPlus.txt",currdownfiles)
currboundfiles = gsub("down.permadistance","BoundPlus.txt",currdownfiles)
currdowns = llply(paste(currpath,currdownfiles,sep="/"), read.table)
currdes = llply(paste(altpath,currdefiles,sep="/"), read.table,sep="|",stringsAsFactors=F)
currdes = llply(currdes,unlist)
currbounds = llply(paste(altpath,currboundfiles,sep="/"), read.table,sep="|",stringsAsFactors=F)
currbounds = llply(currbounds,unlist)
delen = unlist(llply(currdes,length))
boundlen = unlist(llply(currbounds,length))
subing = apply(cbind(delen,boundlen),1,min)
breaking = seq(-10000,10000,250)
holdem = matrix(NA,1000,length(breaking)-1)
#holdem = matrix(NA,10,length(breaking)-1)
for(j in 1:1000){
  holder = c()
  for(i in 1:length(subing)){
    currgenes = sample(currbounds[[i]],subing[i])
    holder = c(holder,currdowns[[i]][currdowns[[i]][,1] %in% currgenes, 2])
    holder[which(holder > 10000)] = 10000
    holder[which(holder < -10000)] = -10000
  }
  holdem[j,] = hist(holder,plot=F,breaks=breaking)$counts
}

tops = apply(holdem,2,quantile,probs=0.95)
bottoms = apply(holdem,2,quantile,probs=0.05)
xs = c(hist(holder,breaks=breaking,plot=FALSE)$mids,rev(hist(holder,breaks=breaking,plot=FALSE)$mids))
ys = c(tops,rev(bottoms))

holder = c()
for(i in 1:length(subing)){
  currgenes = sample(currdes[[i]],subing[i])
  holder = c(holder,currdowns[[i]][currdowns[[i]][,1] %in% currgenes, 2])
  holder[which(holder > 10000)] = 10000
  holder[which(holder < -10000)] = -10000
}
obs = hist(holder,breaks=breaking,plot=F)$counts
breakers = hist(holder,breaks=breaking,plot=F)$mids
pdf("TSS Dist Not Integrated.pdf")
plot(xs,ys,type="n")
polygon(xs,ys,col="dodgerblue2")
lines(breakers,obs,col="indianred",lwd=2)
dev.off()