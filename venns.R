library('VennDiagram')
bindingmatrix = read.table("3RUV2_NSAveraged_FACTOR_10kbBindingWindow_overlaptable.txt",sep="\t",header=T)
pdf("venns.pdf")
for(i in 1:dim(bindingmatrix)[1]){
  grid.newpage()
  grid.draw(draw.pairwise.venn(as.numeric(bindingmatrix[i,11]),
                               as.numeric(bindingmatrix[i,5]),
                               as.numeric(bindingmatrix[i,12]),
                               category=c(paste(bindingmatrix[i,1],"Bound",sep="\n"),paste(bindingmatrix[i,1],"DE",sep="\n")),
                               fill=c("tomato","dodgerblue2"),
                               col=NA,ext.percent = rep(0.005, 3),
                               alpha=c(.8,.8)))
}
dev.off()