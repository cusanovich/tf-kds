bindingmatrix = read.table("~/home/Kd_Arrays/CombinedBinding/Results/union_RUV2_NSAveraged_FACTOR_10kbBindingWindow_overlaptable.txt",sep="\t",header=T)
mycor = cor(bindingmatrix[bindingmatrix[,4] > 0,4],bindingmatrix[bindingmatrix[,4] > 0,5],method="spearman")
permp = 0
xer = c()
for(i in 1:10000){
  currbound = sample(bindingmatrix[bindingmatrix[,4] > 0,4])
  currde = sample(bindingmatrix[bindingmatrix[,4] > 0,5])
  if(abs(cor(currbound,currde,method="spearman")) > mycor){
    permp = permp + 1
  }
  xer[i] = cor(currbound,currde,method="spearman")
}
permp/10000