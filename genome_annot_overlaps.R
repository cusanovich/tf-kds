library('plyr')
library('beanplot')
source('./config.R')

centispec = c("GrepPhastCons","GrepPosterior","GrepPWM")
coloring = list(c("mediumorchid3","goldenrod"),c("indianred","dodgerblue2"))
#grandpath = paste0("~/home/Kd_Arrays/GenomeAnnotations/GrepBeds/",windowname,"/")
#masterpath = "~/home/Kd_Arrays/GenomeAnnotations/"
tests = list.files(grandpath)

pdf(paste0(masterpath,windowname,"_union_Annotations.pdf"))
par(mfrow=c(2,2))
par(mar=c(4,4,6,4)+0.1)
for(i in 1:length(tests)){
  currtest = tests[i]
  print(currtest)
  if(currtest == "GrepChipvCenti"){next}
  currpath = paste0(grandpath,currtest)
#  currmaster = read.table(paste0(masterpath,"Perms/",windowname,"/",currtest,".master.perms"))
#  mastermedstf = apply(currmaster[grep("TF",currmaster[,2]),3:dim(currmaster)[2]],2,median)
#  mastermedsdown = apply(currmaster[grep("Downstream",currmaster[,2]),3:dim(currmaster)[2]],2,median)
#  altmastertf = read.table(paste0(masterpath,"Perms/",windowname,"/",currtest,".alt.tf.perms"))
#  altmasterdown = read.table(paste0(masterpath,"Perms/",windowname,"/",currtest,".alt.down.perms"))
  currtfboundfiles = list.files(path=currpath,pattern = ".tf.bound.")
  currtfdefiles = list.files(path=currpath,pattern = ".tf.debound.")
  currdownboundfiles = list.files(path=currpath,pattern = ".down.boundplus.")
  currdowndefiles = list.files(path=currpath,pattern = ".down.deboundplus.")
  currtfbounds = llply(paste(currpath,currtfboundfiles,sep="/"), read.table)
  currtfdes = llply(paste(currpath,currtfdefiles,sep="/"), read.table)
  currdownbounds = llply(paste(currpath,currdownboundfiles,sep="/"), read.table)
  currdowndes = llply(paste(currpath,currdowndefiles,sep="/"), read.table)
  if(currtest == "GrepChipvCenti"){
    currtfbounds = llply(currtfbounds, function(x){length(which(x == "centi"))/length(x[,1])})
    currtfdes = llply(currtfdes, function(x){length(which(x == "centi"))/length(x[,1])})
    currdownbounds = llply(currdownbounds, function(x){length(which(x == "centi"))/length(x[,1])})
    currdowndes = llply(currdowndes, function(x){length(which(x == "centi"))/length(x[,1])})
    for(i in 1:length(currtfbounds)){
      tfchecker = currtfbounds[[i]] + currtfdes[[i]]
      downchecker = currdownbounds[[i]] + currdowndes[[i]]
      if(tfchecker == 0 | tfchecker == 2){
        currtfbounds[[i]] = NA
        currtfdes[[i]] = NA
      }
      if(downchecker == 0 | downchecker == 2){
        currdownbounds[[i]] = NA
        currdowndes[[i]] = NA
      }
    }
  }
  if(currtest != "GrepChipvCenti"){
    tfmed = median(unlist(llply(currtfdes,function(x){median(x[,1],na.rm=T)})))
    downmed = median(unlist(llply(currdowndes,function(x){median(x[,1],na.rm=T)})))
  }else{
    tfmed = median(unlist(currtfdes),na.rm = T)
    downmed = median(unlist(currdowndes),na.rm = T)
  }
  tfmeds = median(unlist(currtfdes),na.rm = T) - median(unlist(currtfbounds),na.rm = T)
  downmeds = median(unlist(currdowndes),na.rm = T) - median(unlist(currdownbounds),na.rm = T)
  tfdir = ifelse(tfmeds > 0,"DE > Bound","DE < Bound")
  downdir = ifelse(downmeds > 0,"DE > Bound","DE < Bound")
#  tfempiricalp = length(which(abs(mastermedstf-mean(mastermedstf)) > abs(tfmed - mean(mastermedstf))))/1000
#  downempiricalp = length(which(abs(mastermedsdown-mean(mastermedsdown)) > abs(downmed - mean(mastermedsdown))))/1000
#  alttfempiricalp = length(which(abs(altmastertf[,1]-mean(altmastertf[,1])) > abs(tfmed - mean(altmastertf[,1]))))/1000
#  altdownempiricalp = length(which(abs(altmasterdown[,1]-mean(altmasterdown[,1])) > abs(downmed - mean(altmasterdown[,1]))))/1000
  if(currtest == "GrepDistance"){
    tfmed = median(unlist(llply(currtfdes,function(x){median(x[,1],na.rm=T)})))
    downmed = median(unlist(llply(currdowndes,function(x){median(x[,1],na.rm=T)})))
    alttfmed = median(unlist(currtfdes),na.rm = T)
    altdownmed = median(unlist(currdowndes),na.rm = T)
    tfmeds = median(abs(unlist(currtfdes)),na.rm = T) - median(abs(unlist(currtfbounds)),na.rm = T)
    downmeds = median(abs(unlist(currdowndes)),na.rm = T) - median(abs(unlist(currdownbounds)),na.rm = T)
#    tfempiricalp = length(which(abs(mastermedstf-mean(mastermedstf)) > abs(tfmed - mean(mastermedstf))))/1000
#    downempiricalp = length(which(abs(mastermedsdown-mean(mastermedsdown)) > abs(downmed - mean(mastermedsdown))))/1000
#    alttfempiricalp = length(which(abs(altmastertf[,1]-mean(altmastertf[,1])) > abs(alttfmed - mean(altmastertf[,1]))))/1000
#    altdownempiricalp = length(which(abs(altmasterdown[,1]-mean(altmasterdown[,1])) > abs(altdownmed - mean(altmasterdown[,1]))))/1000
    tfdir = ifelse(tfmeds > 0,"DE > Bound","DE < Bound")
    downdir = ifelse(downmeds > 0,"DE > Bound","DE < Bound")
    dense1 = density(unlist(currtfdes))
    dense2 = density(unlist(currtfbounds))
    plot(range(dense1$x, dense2$x), range(dense1$y, dense2$y), type = "n", xlab = "Distance",
         ylab = "Density", main="Distance from TSS (TFs)")
    lines(dense1, col = "indianred",lwd=3)
    lines(dense2, col = "dodgerblue2",lwd=3)
    legend("topright",legend=c("DE","Bound"),fill=c("indianred","dodgerblue2"))
    plot(range(-1000, 1000), range(dense1$y, dense2$y), type = "n", xlab = "Distance",
         ylab = "Density", main="Distance from TSS (TFs) - Zoom")
    abline(v=0,lty="dashed",lwd=3)
    lines(dense1, col = "indianred",lwd=3)
    lines(dense2, col = "dodgerblue2",lwd=3)
    legend("topright",legend=c("DE","Bound"),fill=c("indianred","dodgerblue2"))
    dense1 = density(unlist(currdowndes))
    dense2 = density(unlist(currdownbounds))
    plot(range(dense1$x, dense2$x), range(dense1$y, dense2$y), type = "n", xlab = "Distance",
         ylab = "Density", main="Distance from TSS (Downstream)")
    lines(dense1, col = "indianred",lwd=3)
    lines(dense2, col = "dodgerblue2",lwd=3)
    legend("topright",legend=c("DE","Bound"),fill=c("indianred","dodgerblue2"))
    plot(range(-1000, 1000), range(dense1$y, dense2$y), type = "n", xlab = "Distance",
         ylab = "Density", main="Distance from TSS (Downstream) - Zoom")
    abline(v=0,lty="dashed",lwd=3)
    lines(dense1, col = "indianred",lwd=3)
    lines(dense2, col = "dodgerblue2",lwd=3)
    legend("topright",legend=c("DE","Bound"),fill=c("indianred","dodgerblue2"))
    dense1 = density(unlist(currdowndes)[abs(unlist(currdowndes)) > 1000])
    dense2 = density(unlist(currdownbounds)[abs(unlist(currdownbounds)) > 1000])
    plot(range(dense1$x, dense2$x), range(dense1$y, dense2$y), type = "n", xlab = "Distance",
         ylab = "Density", main="Distance from TSS (Downstream)")
    lines(dense1, col = "indianred",lwd=3)
    lines(dense2, col = "dodgerblue2",lwd=3)
    legend("topright",legend=c("DE","Bound"),fill=c("indianred","dodgerblue2"))
    plot(range(-5000, 5000), range(dense1$y, dense2$y), type = "n", xlab = "Distance",
         ylab = "Density", main="Distance from TSS (Downstream) - Zoom")
    abline(v=0,lty="dashed",lwd=3)
    lines(dense1, col = "indianred",lwd=3)
    lines(dense2, col = "dodgerblue2",lwd=3)
    boxplot(log10(abs(unlist(currtfdes))+0.001),log10(abs(unlist(currtfbounds))+0.001),col=c("indianred","dodgerblue2"),outline=F,
            names=c("DE","Bound"),
            main=paste0(currtest," TF\nP-value = ",signif(wilcox.test(abs(unlist(currtfdes)),abs(unlist(currtfbounds)),na.action = na.omit)$p.value,4),"\n",tfdir),
            notch=T)
    boxplot(log10(abs(unlist(currdowndes))+0.001),log10(abs(unlist(currdownbounds))+0.001),col=c("indianred","dodgerblue2"),outline=F,
            names=c("DE","Bound"),
            main=paste0(currtest," Downstream\nP-value = ",signif(wilcox.test(abs(unlist(currdowndes)),abs(unlist(currdownbounds)),na.action = na.omit)$p.value,4),"\n",downdir),
            notch=T)
#    hist(mastermedstf,xlab="Absolute Distance to TSS",
#         main=paste0(currtest," TF\nEmpirical P-value = ",tfempiricalp),
#                     xlim=c(min(min(mastermedstf),tfmed),max(max(mastermedstf),tfmed)))
#    abline(v=tfmed,lwd=3,col="indianred")
#    hist(mastermedsdown,xlab="Absolute Distance to TSS",
#         main=paste0(currtest," Downstream\nEmpirical P-value = ",downempiricalp),
#                     xlim=c(min(min(mastermedsdown),downmed),max(max(mastermedsdown),downmed)))
#    abline(v=downmed,lwd=3,col="indianred")
#    hist(altmastertf[,1],xlab="Absolute Distance to TSS",
#         main=paste0(currtest," TF\nAlt Empirical P-value = ",alttfempiricalp),
#                     xlim=c(min(min(altmastertf),alttfmed),max(max(altmastertf),alttfmed)))
#    abline(v=alttfmed,lwd=3,col="indianred")
#    hist(altmasterdown[,1],xlab="Absolute Distance to TSS",
#         main=paste0(currtest," Downstream\nAlt Empirical P-value = ",altdownempiricalp),
#                     xlim=c(min(min(altmasterdown),altdownmed),max(max(altmasterdown),altdownmed)))
#    abline(v=altdownmed,lwd=3,col="indianred")
    currtfboundsb = llply(currtfbounds, function(x){length(which(x > 0))/length(x[,1])})
    currtfdesb = llply(currtfdes, function(x){length(which(x > 0))/length(x[,1])})
    currdownboundsb = llply(currdownbounds, function(x){length(which(x > 0))/length(x[,1])})
    currdowndesb = llply(currdowndes, function(x){length(which(x > 0))/length(x[,1])})
    boxplot(unlist(currtfdesb),unlist(currtfboundsb),col=c("indianred","dodgerblue2"),outline=F,
            names=c("DE","Bound"),
            main=paste0(currtest," TF Frac(+)\nP-value = ",signif(wilcox.test(abs(unlist(currtfdesb)),abs(unlist(currtfboundsb)),na.action = na.omit)$p.value,4),"\n",tfdir),
            notch=T)
    boxplot(unlist(currdowndesb),unlist(currdownboundsb),col=c("indianred","dodgerblue2"),outline=F,
            names=c("DE","Bound"),
            main=paste0(currtest," Downstream Frac(+)\nP-value = ",signif(wilcox.test(abs(unlist(currdowndesb)),abs(unlist(currdownboundsb)),na.action = na.omit)$p.value,4),"\n",downdir),
            notch=T)
  }else{
    if(currtest == "GrepCage"){
      currtfbounds = llply(currtfbounds, function(x){log10(x+0.001)})
      currtfdes = llply(currtfdes, function(x){log10(x+0.001)})
      currdownbounds = llply(currdownbounds, function(x){log10(x+0.001)})
      currdowndes = llply(currdowndes, function(x){log10(x+0.001)})
    }
    colorno = ifelse(currtest %in% centispec,1,2)
    boxplot(unlist(currtfdes),unlist(currtfbounds),
            col=coloring[[colorno]],names=c("DE","Bound"),outline=F,notch=T,
            main=paste0(currtest," TF\nP-value = ",signif(wilcox.test(unlist(currtfdes),unlist(currtfbounds),na.action = na.omit)$p.value,4),"\n",tfdir))
    boxplot(unlist(currdowndes),unlist(currdownbounds),
            col=coloring[[colorno]],names=c("DE","Bound"),outline=F,notch=T,
            main=paste0(currtest," Downstream\nP-value = ",signif(wilcox.test(unlist(currdowndes),unlist(currdownbounds),na.action = na.omit)$p.value,4),"\n",downdir))
#    hist(mastermedstf,main=paste0(currtest," TF\nEmpirical P-value = ",tfempiricalp),
#         xlim=c(min(min(mastermedstf),tfmed),max(max(mastermedstf),tfmed)))
#    abline(v=tfmed,lwd=3,col="indianred")
#    hist(mastermedsdown,main=paste0(currtest," Downstream\nEmpirical P-value = ",downempiricalp),
#         xlim=c(min(min(mastermedsdown),downmed),max(max(mastermedsdown),downmed)))
#    abline(v=downmed,lwd=3,col="indianred")
#    hist(altmastertf[,1],main=paste0(currtest," TF\nAlt Empirical P-value = ",alttfempiricalp),
#         xlim=c(min(min(mastermedstf),tfmed),max(max(mastermedstf),tfmed)))
#    abline(v=tfmed,lwd=3,col="indianred")
#    hist(altmasterdown[,1],main=paste0(currtest," Downstream\nAlt Empirical P-value = ",altdownempiricalp),
#         xlim=c(min(min(altmasterdown),downmed),max(max(altmasterdown),downmed)))
#    abline(v=downmed,lwd=3,col="indianred")
  }
}