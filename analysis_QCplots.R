png("EP300 Sample SVA QC.png")
pairs(svatable)
pairs(cleansed.irw)
heatmap(cor(svatable,method="spearman"))
heatmap(cor(cleansed.irw,method="spearman"))
rlemaker = function(expr){
    exprmeds = apply(expr,1,median)
    exprrle = expr - exprmeds
    return(exprrle)
}
svatable.rle = rlemaker(svatable)
boxplot(svatable.rle)
svatable.rle = rlemaker(cleansed.irw)
boxplot(svatable.rle)
dev.off()

detection = read.table("Detection_Scores_All3.txt")