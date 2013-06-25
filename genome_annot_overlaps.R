library('plyr')
grandpath = "/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/GrepBeds/"
tests = list.files("/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/GrepBeds/")
pdf("/mnt/lustre/home/cusanovich/Kd_Arrays/GenomeAnnotations/Annotations.pdf")
#grandpath = "~/home/Kd_Arrays/GenomeAnnotations/GrepBeds/"
#tests = list.files("~/home/Kd_Arrays/GenomeAnnotations/GrepBeds/")
par(mfrow=c(2,2))
par(mar=c(4,4,6,4)+0.1)
for(i in 1:length(tests)){
  currtest = tests[i]
  print(currtest)
  currpath = paste0(grandpath,currtest)
  currtfboundfiles = list.files(path=currpath,pattern = ".tf.bound.")
  currtfdefiles = list.files(path=currpath,pattern = ".tf.debound.")
  currdownboundfiles = list.files(path=currpath,pattern = ".down.bound.")
  currdowndefiles = list.files(path=currpath,pattern = ".down.debound.")
  currtfbounds = llply(paste(currpath,currtfboundfiles,sep="/"), read.table)
  currtfdes = llply(paste(currpath,currtfdefiles,sep="/"), read.table)
  currdownbounds = llply(paste(currpath,currdownboundfiles,sep="/"), read.table)
  currdowndes = llply(paste(currpath,currdowndefiles,sep="/"), read.table)
  outliers = T
  if(currtest == "GrepPosterior" | currtest == "GrepPWM"){
    outliers = F
  }
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
#  if(currtest == "GrepCage"){
#    currtfbounds = llply(currtfbounds, log10)
#    currtfdes = llply(currtfdes, log10)
#    currdownbounds = llply(currdownbounds, log10)
#    currdowndes = llply(currdowndes, log10)
#  }
#  if(currtest == "GrepPosterior"){
#    currtfbounds = llply(currtfbounds, function(x){1-x})
#    currtfdes = llply(currtfdes, function(x){1-x})
#    currdownbounds = llply(currdownbounds, function(x){1-x})
#    currdowndes = llply(currdowndes, function(x){1-x})
#  }
  tfmeds = median(unlist(currtfdes)) - median(unlist(currtfbounds))
  downmeds = median(unlist(currdowndes)) - median(unlist(currdownbounds))
  tfdir = ifelse(tfmeds > 0,"DE > Bound","DE < Bound")
  downdir = ifelse(downmeds > 0,"DE > Bound","DE < Bound")
  if(currtest == "GrepDistance"){
    dense1 = density(unlist(currtfdes))
    dense2 = density(unlist(currtfbounds))
    plot(range(dense1$x, dense2$x), range(dense1$y, dense2$y), type = "n", xlab = "Distance",
         ylab = "Density", main="Distance from TSS (TFs)")
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
    hist(unlist(currtfdes),col="indianred")
    hist(unlist(currtfbounds),col="dodgerblue2")
    hist(unlist(currdowndes),col="indianred")
    hist(unlist(currdownbounds),col="dodgerblue2")
    boxplot(log10(abs(unlist(currtfdes))+0.001),log10(abs(unlist(currtfbounds))+0.001),col=c("indianred","dodgerblue2"),
            names=c("DE","Bound"),main=paste0(currtest," TF\nP-value = ",signif(wilcox.test(abs(unlist(currtfdes)),abs(unlist(currtfbounds)),na.action = na.omit)$p.value,4),"\n",tfdir),notch=T)
    boxplot(log10(abs(unlist(currdowndes))+0.001),log10(abs(unlist(currdownbounds))+0.001),col=c("indianred","dodgerblue2"),
            names=c("DE","Bound"),main=paste0(currtest," Downstream\nP-value = ",signif(wilcox.test(abs(unlist(currdowndes)),abs(unlist(currdownbounds)),na.action = na.omit)$p.value,4),"\n",downdir),notch=T)
  }else{
    boxplot(log10(unlist(currtfdes)+0.001),log10(unlist(currtfbounds)+0.001),
            col=c("indianred","dodgerblue2"),names=c("DE","Bound"),outline=outliers,
            main=paste0(currtest," TF\nP-value = ",signif(wilcox.test(unlist(currtfdes),unlist(currtfbounds),na.action = na.omit)$p.value,4),"\n",tfdir),
            notch=T)
    boxplot(log10(unlist(currdowndes)+0.001),log10(unlist(currdownbounds)+0.001),
            col=c("indianred","dodgerblue2"),names=c("DE","Bound"),outline=outliers,
            main=paste0(currtest," Downstream\nP-value = ",signif(wilcox.test(unlist(currdowndes),unlist(currdownbounds),na.action = na.omit)$p.value,4),"\n",downdir),
            notch=T)
    hist(unlist(currtfdes),col="indianred",main=paste0(currtest," TF"))
    hist(unlist(currtfbounds),col="dodgerblue2",main=paste0(currtest," Downstream"))
    hist(unlist(currdowndes),col="indianred",main=paste0(currtest," TF"))
    hist(unlist(currdownbounds),col="dodgerblue2",main=paste0(currtest," Downstream"))
  }
}