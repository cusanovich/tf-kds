source('./config.R')

resultsmatrix = as.matrix(read.table(bindingmatrix),sep="\t")
resultsmatrix = resultsmatrix[which(rowSums(resultsmatrix) > 0),]
resultsmatrix = resultsmatrix[,which(colSums(resultsmatrix) > 0)]
resultsbinary = resultsmatrix>0
resultsbinary = resultsbinary + 0

binary.phi = matrix(1,dim(resultsbinary)[2],dim(resultsbinary)[2])
for(i in 1:(dim(resultsbinary)[2])-1){
    for(j in (i+1):dim(resultsbinary)[2]){
        BC = length(which(resultsbinary[,i] == 1 & resultsbinary[,j] == 1))*length(which(resultsbinary[,i] == 0 & resultsbinary[,j] == 0))
        AD = length(which(resultsbinary[,i] == 1 & resultsbinary[,j] == 0))*length(which(resultsbinary[,i] == 0 & resultsbinary[,j] == 1))
        AB = sum(resultsbinary[,i])
        CD = length(resultsbinary[,i]) - sum(resultsbinary[,i])
        AC = length(resultsbinary[,j]) - sum(resultsbinary[,j])
        BD = sum(resultsbinary[,j])
        phi = (BC - AD)/(sqrt(AB)*sqrt(CD)*sqrt(AC)*sqrt(BD))
        binary.phi[i,j] = phi
        binary.phi[j,i] = phi
    }
}

colnames(binary.phi) = colnames(resultsmatrix)
rownames(binary.phi) = colnames(resultsmatrix)

write.table(binary.phi,binarytable)