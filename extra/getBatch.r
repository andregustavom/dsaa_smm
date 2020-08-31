getBatch <- function(test, vClass, vPerc, testSize){
  
  dts <- NULL
  
  vsizes <- trunc(testSize*vPerc)
  resto  <- testSize - sum(vsizes)
  vsizes[length(vsizes)] <- vsizes[length(vsizes)]+resto
  indices <- NULL
  
  for(i in 1:length(vClass)){
    aux <- which(test$class==vClass[i])
    
    ind <- sample(aux, min(vsizes[i],length(aux)) )
    dts <- rbind(dts, test[ind,])
    indices <- c(indices, ind)
    
  }
  return(list(dts, indices))
}
