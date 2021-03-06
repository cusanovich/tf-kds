library('plyr')
library('stringr')
library('qvalue')
source('./config.R')
bounddir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/eQTLs/Bound/'
dedir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/eQTLs/DEandBound/'
unbounddir = '/mnt/lustre/home/cusanovich/Kd_Arrays/Trans_eQTLs/eQTLs/Unbound/'

qqer = function(de,bound,name){
  obsp <- -log10(de)
  N <- length(obsp)
  null <- -log10(ppoints(N))
  c95 <- rep(0,N)
  c05 <- rep(0,N)
  for(i in 1:N){
    c95[i] <- qbeta(0.95,i,N-i+1)
    c05[i] <- qbeta(0.05,i,N-i+1)
  }
  xx <- c(sort(null),sort(null,decreasing=T))
  yy <- c(sort(-log10(c95)),sort(-log10(c05),decreasing=T))
  plot(1,1, type="n", xlim=c(0,max(null)+0.5),ylim=c(0,max(obsp,-log10(bound))+1.5),
       xlab="Expected", ylab="Observed",main=name,las=1)
  polygon(xx, yy, col="gray",border="gray")
  lines(null,null,lwd=2)
  points(sort(null),sort(-log10(bound)),pch=20,col="mediumorchid3")
  points(sort(null),sort(obsp),pch=20,col="indianred")
  points(sort(null),sort(-log10(unbound)),pch=20,col="dodgerblue2")
  legend("topleft",legend=c("Candidate Genes","Bound\nControl Genes","Unbound\nControl Genes"),
         col=c("indianred","mediumorchid3","dodgerblue2"),pch=20)
}

all.boundeqtls = list.files(path = bounddir,pattern="_100kb_sidney_results.txt")
all.boundeqtls = all.boundeqtls[file.info(paste0(bounddir,all.boundeqtls))$size > 0]
all.boundeqtls = all.boundeqtls[file.info(paste0(dedir,all.boundeqtls))$size > 0]
bound = llply(paste0(bounddir,all.boundeqtls), read.table)
de = llply(paste0(dedir,all.boundeqtls), read.table)
unbound = llply(paste0(unbounddir,all.boundeqtls), read.table)
namers = c()
gemmacisers = c()
permcisers = c()
for(i in 1:length(all.boundeqtls)){
  namers[i] = strsplit(all.boundeqtls[i],"_")[[1]][1]
  gemmacisers[i] = de[[i]]$V4[1]
  permcisers[i] = de[[i]]$V5[1]
  bound[[i]] = bound[[i]][-1,]
  de[[i]] = de[[i]][-1,]
  unbound[[i]] = unbound[[i]][-1,]
}
names(bound) = namers
names(de) = namers
names(unbound) = namers

pdf("qqs.pdf")
par(mfrow=c(2,2))
gemmaallbound = c()
gemmaallde = c()
gemmaallunbound = c()
permallbound = c()
permallde = c()
permallunbound = c()
for(i in 1:length(bound)){
  print(namers[i])
  desize = length(de[[i]]$V4)
  boundsize = length(bound[[i]]$V4)
  unboundsize = length(unbound[[i]]$V4)
  samplesize = min(desize,boundsize,unboundsize)
  gemmacurrde = de[[i]]$V4[sample(desize,samplesize)]
  permcurrde = de[[i]]$V5[sample(desize,samplesize)]
  permcurrde[permcurrde == 0] = 0.0001
  gemmacurrbound = bound[[i]]$V4[sample(boundsize,samplesize)]
  permcurrbound = bound[[i]]$V5[sample(boundsize,samplesize)]
  permcurrbound[permcurrbound == 0] = 0.0001
  gemmacurrunbound = unbound[[i]]$V4[sample(unboundsize,samplesize)]
  permcurrunbound = unbound[[i]]$V5[sample(unboundsize,samplesize)]
  permcurrunbound[permcurrunbound == 0] = 0.0001
  gemmaallde = c(gemmaallde,gemmacurrde)
  gemmaallbound = c(gemmaallbound,gemmacurrbound)
  gemmaallunbound = c(gemmaallunbound,gemmacurrunbound)
  permallde = c(permallde,permcurrde)
  permallbound = c(permallbound,permcurrbound)
  permallunbound = c(permallunbound,permcurrunbound)
  qqer(gemmacurrde,gemmacurrbound,gemmacurrunbound,name=paste0("GEMMA ",namers[i],"\ncis P-value = ",gemmacisers[i]))
  qqer(permcurrde,permcurrbound,permcurrunbound,name=paste0("Permutation ",namers[i],"\ncis P-value = ",permcisers[i]))
}

gemmaall.q = qvalue(gemmaallde)$qvalue
qqer(gemmaallde,gemmaallbound,gemmaallunbound,
     name=paste0("GEMMA Trans QQ-Plot\nMin q-value = ",round(min(gemmaall.q),3)))
permall.q = qvalue(permallde)$qvalue
qqer(permallde,permallbound,permallunbound,
     name=paste0("Permutation Trans QQ-Plot\nMin q-value = ",round(min(permall.q),3)))
hist(gemmaallde,main="GEMMA P-values",las=1,xlab="P-value",col="indianred")
hist(gemmaallbound,col="mediumorchid3",add=T)
hist(gemmaallunbound,col="dodgerblue2",add=T)
hist(permallde,main="Permutation P-values",las=1,xlab="P-value",col="indianred")
hist(permallbound,col="mediumorchid3",add=T)
hist(permallbound,col="dodgerblue2",add=T)
plot(density(gemmaallde),main="GEMMA P-values",las=1,xlab="P-value",col="indianred",
     type="l",lwd=2)
lines(density(gemmaallbound),lwd=2,col="mediumorchid3")
lines(density(gemmaallunbound),lwd=2,col="dodgerblue2")
plot(density(permallde),main="Permutation P-values",las=1,xlab="P-value",col="indianred",
     type="l",lwd=2)
lines(density(permallbound),lwd=2,col="mediumorchid3")
lines(density(permallunbound),lwd=2,col="dodgerblue2")
dev.off()