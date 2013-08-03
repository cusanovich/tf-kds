library('VennDiagram')
windowsize="1kb"
bindingmatrix = read.table(paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Results/RUV2_NSAveraged_FACTOR_",windowsize,"BindingWindow_overlaptable.txt"),sep="\t",header=T)
pdf(paste0("/mnt/lustre/home/cusanovich/Kd_Arrays/CombinedBinding/Results/",windowsize,"_venns.pdf"))
for(i in 1:dim(bindingmatrix)[1]){
  grid.newpage()
  grid.draw(draw.pairwise.venn(as.numeric(bindingmatrix[i,11]),
                               as.numeric(bindingmatrix[i,5]),
                               as.numeric(bindingmatrix[i,12]),
                               category=c(paste(bindingmatrix[i,2],"Bound",sep=" "),paste(bindingmatrix[i,2],"DE",sep=" ")),
                               fill=c("indianred","dodgerblue2"),
                               col=NA,ext.percent = rep(0.005, 3),
                               alpha=c(.8,.6),cex=2,cat.cex=2,cat.dist=c(0.1,0.1),
                               cat.just=list(c(-1,-1),c(1,1))))
}
dev.off()