sigmoid = function(x){
  return(2/(1+(exp(-x)))-1)
  }

sig = function(mat){
  mns = colMeans(mat) ###Moyenne par gène sur les x iterations/
  DS = (t(mat) - mns) ^ 2 / (4 * nrow(mat))
  DS = colSums(DS)
  return(mean(DS))
  }


W = as.matrix(read.csv("W.txt", header = F))
Si = as.numeric(read.csv("Si.txt", header = F, as.is = T))

sums = NULL
for(i in 1:10){
  Si = Si %*% W
  Si = sigm(Si)  
  sums = rbind(sums, Si)
  }  
  

sig(sums[8:10,])
