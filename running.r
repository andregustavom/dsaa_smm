source("./discret_data.r")
source("./distance_DYS.r")
source("./getBatch.r")
source("./getTPRandFPRbyThreshold.r")
source("./TernarySearch.r")

library(microbenchmark)
library(philentropy)
library(CORElearn)


library("batch")
parseCommandArgs()


EMQ_method <- function(train, test){


  pTr      <- c(0.5, 0.5)#table(train$class)/nrow(train)
  predTe_s <- test#CC_method(test[,1])[[1]]#as.data.frame(cbind(test,1-test))
  p   <- list(predTe_s) # probabilities of the test samples for each iteration
  pTe <- list(pTr)    # distributions for each iteration
  nE  <- nrow(test) #number of test samples
  nC  <- 2#length(pTr)      #number of classes
  s   <- 1


  repeat{
    aux <- matrix(ncol=nC, nrow = nE)
    auxC <- c(1:nC)
    for (ic in 1:nC){
      for(ie in 1:nE){
        numerator   <- (pTe[[s]][ic]/pTr[ic])*predTe_s[ie,ic]
        denominator <- c(1:nC)
        for(ic2 in 1:nC)  denominator[ic2] <- (pTe[[s]][ic2]/pTr[ic2])*predTe_s[ie,ic2]
        aux[ie,ic] <- numerator/sum(denominator)
      }
      auxC[ic] <- sum(p[[s]][,ic])/nrow(test)
    }
    pTe <- c(pTe, list(auxC))
    p   <- c(p, list(aux))
    s <- s + 1
    if(s > 4 ){break}
  }

  re <- round(pTe[[s]],2)
  return(re)


}



SMM_method_perf <- function(p.score, n.score, test){
  
  #start
  Sty_1 <- sum(p.score)/length(p.score)
  Sty_2 <- sum(n.score)/length(n.score)
  Uy    <- sum(test)/length(test)
  
  result <- (Uy - Sty_2)/(Sty_1 - Sty_2)
  
  ifelse(result < 0, result <- 0, result <- result)
  ifelse(result > 1, result <- 1, result <- result)
  
  result <- c(result, 1 - result)
  
  return(result)
  
}

#----------------------------------------------------------------------------------------

HDy_method_perf <- function(p.score, n.score, test){
  
  #start
  
  alpha <- seq(0,1,by=0.01)
  b_sizes <- seq(10,110, length.out = 11)
  
  #result <- NULL
  result <- c(1:11)
  
  for(hi in 1:length(b_sizes)){
    
    v_prob <- seq(0,1,length.out = b_sizes[hi])
    v_prob <- c(v_prob[-length(v_prob)], v_prob[length(v_prob)]+0.1)
    
    Sty_1 <- discret_data(p.score,v_prob)
    Sty_2 <- discret_data(n.score,v_prob)
    Uy    <- discret_data(test,v_prob)
    
    #vDist <- NULL
    vDist <- c(1:length(alpha))
    for(k in 1:length(alpha)){
      aux <- (Sty_1*alpha[k])+ (Sty_2*(1-alpha[k]))
      #vDist <- c(vDist, as.numeric(distance(rbind(aux, Uy), method = "hellinger", p=0.5)))
      vDist[k] <- as.numeric(distance(rbind(aux, Uy), method = "hellinger", p=0.5))
    }
    #result <- c(result, alpha[which.min(vDist)])
    result[hi] <- alpha[which.min(vDist)]
  }  
  
  result <- median(result)
  result <- c(result, 1 - result)
  
  #stop
  
  #names(result) <- c("1", "2")
  
  return(result)
  
}


#----------------------------------------------------------------------------------------

SORD_perf_1 <- function(p.score, n.score, test){
  


  f <- function(x){
    return(PNTDiff_perf(p.score, n.score, test, x))
  }
  
  result <- TernarySearch(0, 1, f, 1e-2)
  result <- c(result, 1 - result)
  
  return(result)
}


PNTDiff_perf_1 <- function(pos, neg, test, pos_prop){
  p_w <- pos_prop / length(pos)
  n_w <- (1 - pos_prop) / length(neg)
  t_w <- -1 / length(test) # repare no -1 (menos um) !!
  
  # Cria listas com (valor da observacao), peso [fixo pra cada lista])
  p <- cbind(pos, rep(p_w, length(pos)))
  n <- cbind(neg, rep(n_w, length(neg)))
  t <- cbind(test, rep(t_w, length(test)))
  
  # Concatena numa lista soh e ordena pelos valores de obervacao
  #v <- sorted(p + n + t, key = lambda x: x[0])
  v <- rbind(p, n, t)
  v <- v[order(v[,1]),]
  acc <- v[1,2]
  total_cost <- 0
  
  # comeca a iterar no SEGUNDO elemento (python eh 0-indexado)
  for( i in 2:nrow(v)){
    cost_mul <- v[i,1] - v[i - 1, 1] # custo da movimentacao
    total_cost <- total_cost + abs(cost_mul * acc) #movimenta o que esta acumulado
    acc <- acc + v[i,2]
  }
  return(total_cost)
}

SORD_perf <- function(p.score, n.score, test){  


  f <- function(x){
    p_w <- x / length(p.score)
    n_w <- (1 - x) / length(n.score)
    t_w <- -1 / length(test) # repare no -1 (menos um) !!
  
    # Cria listas com (valor da observacao), peso [fixo pra cada lista])
    p <- cbind(p.score, rep(p_w, length(p.score)))
    n <- cbind(n.score, rep(n_w, length(n.score)))
    t <- cbind(test, rep(t_w, length(test)))
  
    # Concatena numa lista soh e ordena pelos valores de obervacao
    #v <- sorted(p + n + t, key = lambda x: x[0])
    v <- rbind(p, n, t)
    v <- v[order(v[,1]),]
    acc <- v[1,2]
    total_cost <- 0
  
    # comeca a iterar no SEGUNDO elemento (python eh 0-indexado)
    for( i in 2:nrow(v)){
      cost_mul <- v[i,1] - v[i - 1, 1] # custo da movimentacao
      total_cost <- total_cost + abs(cost_mul * acc) #movimenta o que esta acumulado
      acc <- acc + v[i,2]
    }
    return(total_cost)

  }
  
  result <- TernarySearch(0, 1, f, 1e-2)
  result <- c(result, 1 - result)
  
  return(result)
}




#----------------------------------------------------------------------------------------

DyS_method_perf <- function(p.score, n.score, test, measure="hellinger"){
  
  #start
  alpha <- seq(0,1,by=0.01)
  b_sizes <- seq(2,20,2)
  
  
  
  #result <- NULL
  result <- c(1:10)
  
  #for(hi in 1:length(b_sizes)){
  for(hi in 1:10){
    
    breaks <- seq(0,1,length.out = b_sizes[hi]+1)
    breaks <- c(breaks[-length(breaks)], 1.1)
    
    Sty_1 <- discret_data_normalized(p.score,breaks)
    Sty_2 <- discret_data_normalized(n.score,breaks)
    Uy    <- discret_data_normalized(test,breaks)
    
    #vDistAll <- NULL
    f <- function(x){
      return(distance_DYS(rbind((Sty_1*x)+ (Sty_2*(1-x)), Uy), method = measure, p=0.5))
    }
    
    #result <- c(result, TernarySearch(0, 1, f))
    result[hi] <- TernarySearch(0, 1, f)
    #vDistAll <- c(vDistAll, f(result))
  } 
  
  result <- median(result)
  result <- c(result, 1 - result)
  #stop
  
  #names(result) <- c("1", "2")
  return(result)
  
}

#----------------------------------------------------------------------------------------
X_method <- function(test, TprFpr){
  
  #x <- apply(TprFpr,2,as.numeric)
  
  x_sim <- cbind(round(1-TprFpr[,2],2), round(TprFpr[,3],2))
  
  idx_thr <- order(abs(x_sim[,1]-x_sim[,2]))[1]
  
  thr <- TprFpr[idx_thr,1]
  
  x <- as.data.frame(t(as.numeric(TprFpr[idx_thr,c(2,3)])))
  
  colnames(x) <- c("tpr", "fpr")
  
  TprFpr <- x
  
  #dC <- CC_method(test, thr)
  dC <- sum(test >= thr)/length(test)
  
  aux <- (dC[[1]] - TprFpr$fpr) / (TprFpr$tpr - TprFpr$fpr)
  
  dC <- aux  
  
  ifelse(dC < 0, dC <- 0, dC <- dC)
  ifelse(dC > 1, dC <- 1, dC <- dC)

  dC <- c(dC, 1 - dC)
  
  return(dC)
}

#----------------------------------------------------------------------------------------
PCC_method <- function(test){
  
  result <- mean(test)
  
  result <- c(result, 1 - result)
  #names(result) <- c("1", "2")
  
  return(result)
}
#----------------------------------------------------------------------------------------
PACC_method <- function(test, TprFpr, thr){
  #x <- apply(TprFpr,2,as.numeric)
  x <- as.data.frame(t(as.numeric(TprFpr[which(TprFpr[,1] == round(thr,2)),c(2,3)])))
  colnames(x) <- c("tpr", "fpr")
  TprFpr <- x
  
  #dC <- PCC_method(test)
  dC <- mean(test)
  
  dC <- c(dC, 1 - dC)
  
  
  aux <- (dC[1] - TprFpr$fpr) / (TprFpr$tpr - TprFpr$fpr)
  
  dC <- aux  
  
  ifelse(dC < 0, dC <- 0, dC <- dC)
  ifelse(dC > 1, dC <- 1, dC <- dC)

  dC <- c(dC, 1 - dC)
  
  return(dC)

}

#----------------------------------------------------------------------------------------
MS_method_old <- function(test, TprFpr, getThr = F){
  
  unique_scores <- TprFpr[,1]
  
  prevalances_array = c(1:length(unique_scores))

  for(i in 1:length(unique_scores)){
    pos <- which(TprFpr[,'tr'] == unique_scores[i])    
   
    tpr <- TprFpr[pos,'tpr']
    fpr <- TprFpr[pos,'fpr']
 
    estimated_positive_ratio <- sum(test>=unique_scores[i])/length(test)

    prevalances_array[i] <-  (abs(estimated_positive_ratio - fpr))/abs(tpr-fpr)
  }
    
  final_positive_prevalence <- median(prevalances_array)
  
  result <- c(final_positive_prevalence, 1 - final_positive_prevalence)


  
  return(result)
}

MS_method <- function(test, TprFpr, getThr = F){
  
  prevalances_array = c(1:nrow(TprFpr))
  
  for(i in 1:nrow(TprFpr)){
    
    estimated_positive_ratio <- sum(test>=TprFpr[i,'tr'])/length(test)
    
    prevalances_array[i] <-  (abs(estimated_positive_ratio - TprFpr[i,'fpr']))/abs(TprFpr[i,'tpr']-TprFpr[i,'fpr'])
  }

final_positive_prevalence <- median(prevalances_array)

result <- c(final_positive_prevalence, 1 - final_positive_prevalence)

return(result)
}


MS2_method_old <- function(test, TprFpr){
 
  TprFpr <- TprFpr[which(TprFpr[,2] - TprFpr[,3] >= 0.25),]

  unique_scores <- TprFpr[,1]

  prevalances_array = c(1:length(unique_scores))

  for(i in 1:length(unique_scores)){

    pos <- which(TprFpr[,'tr'] == unique_scores[i])    
   
    tpr <- TprFpr[pos,'tpr']
    fpr <- TprFpr[pos,'fpr']
 
    estimated_positive_ratio <- sum(test>=unique_scores[i])/length(test)

    prevalances_array[i] <-  (abs(estimated_positive_ratio - fpr))/abs(tpr-fpr)
  }

  final_positive_prevalence <- median(prevalances_array)

  result <- c(final_positive_prevalence, 1 - final_positive_prevalence)
  
  return(result)
}

MS2_method <- function(test, TprFpr){
 
  TprFpr <- TprFpr[which(TprFpr[,2] - TprFpr[,3] >= 0.25),]

  prevalances_array = c(1:nrow(TprFpr))
  
  for(i in 1:nrow(TprFpr)){
    
    estimated_positive_ratio <- sum(test>=TprFpr[i,'tr'])/length(test)
    
    prevalances_array[i] <-  (abs(estimated_positive_ratio - TprFpr[i,'fpr']))/abs(TprFpr[i,'tpr']-TprFpr[i,'fpr'])
  }

final_positive_prevalence <- median(prevalances_array)

result <- c(final_positive_prevalence, 1 - final_positive_prevalence)

return(result)
}


#----------------------------------------------------------------------------------------
CC_method <- function(test, thr=0.5){

  #size <- length(test)
  result <- sum(test >= thr)/length(test)
  return(result)    
  #return(c(result, 1 - result))    
  
  
}

#----------------------------------------------------------------------------------------
ACC_method <- function(test, TprFpr, thr=NULL){
  #dC <- CC_method(test, thr)

  dC <- sum(test >= thr)/length(test)
  
  #x <- apply(TprFpr,2,as.numeric)
  x <- as.data.frame(t(as.numeric(TprFpr[which(TprFpr[,1] == round(thr,2)),c(2,3)])))
  colnames(x) <- c("tpr", "fpr")
  TprFpr <- x
  
  aux <- (dC[1] - TprFpr$fpr) / (TprFpr$tpr - TprFpr$fpr)
  dC <- aux  
  
  ifelse(dC < 0, dC <- 0, dC <- dC)
  ifelse(dC > 1, dC <- 1, dC <- dC)

  dC <- c(dC, 1 - dC)
  
  return(dC)
}

#----------------------------------------------------------------------------------------
T50_method <- function(test, TprFpr){
  #thr <- round(quantile(test)[3], 2)
  
  #idx_thr <- order(abs(as.numeric(TprFpr[, 1]) - thr))[1]
  #x <- apply(TprFpr, 2, as.numeric)
  #x <- as.data.frame(t(as.numeric(TprFpr[idx_thr, c(2, 3)])))
  #colnames(x) <- c("tpr", "fpr")
  #TprFpr <- x
  TprFpr <- as.numeric(TprFpr[order(abs(as.numeric(TprFpr[,2])-0.5))[1],])
  dC <- CC_method(test, TprFpr[1])
  dC <- (dC[1] - TprFpr[3])/(TprFpr[2] - TprFpr[3])
  ifelse(dC < 0, dC <- 0, dC <- dC)
  ifelse(dC > 1, dC <- 1, dC <- dC)
  result <- c(dC, 1 - dC)
  return(result)
}

T50_method_old <- function(test, TprFpr){
  #thr <- quantile(test)[3]
  thr <- sort(test)[length(test)/2]
  #dC <- CC_method(test, thr)
  dC <- sum(test >= thr)/length(test)

  #TprFpr <- TprFpr[which.min(abs(TprFpr[,1]-thr)),]
  TprFpr <- TprFpr[TprFpr[,1]==thr,]
  
  result <- (dC[[1]] - TprFpr[3]) / (TprFpr[2] - TprFpr[3])

  result <- c(result, 1 - result)
  #result <- clamp(result)
  #names(result) <- c("+", "-")

  return(result)
}


#----------------------------------------------------------------------------------------
MAX_method <- function(test, TprFpr){
    
  ma <- abs(TprFpr[,2]-TprFpr[,3])
  id <- order(ma, decreasing = T)[1]
  TprFpr <- TprFpr[id,]

  dC <- (sum(test >= TprFpr[1])/length(test))

  result <- (dC[1] - TprFpr[3]) / (TprFpr[2] - TprFpr[3])

  result <- c(result, 1 - result)
  
  return(result)
}


#----------------------------------------------------------------------------------------

exec_eval_performance <- function(ni, co){

  var_perc <- 0.5#seq(0,1,0.05)
  #var_size <- c(10^2, 10^3, 10^4, 10^6, 10^6, 10^7, 10^8, 10^9, 10^10)[ni]#c(seq(10,100,10), seq(200,500,100)) #Variacoes no tamanho dos batchs de teste
  var_size <- c(100,1000,seq(10000, 100000, 10000))[ni]
  n_tests  <- 1#100 #repeticoes para cada vairiacao
  MF       <- 0.5#seq(0.05,1,0.05)
  MFtr     <- 0.5#seq(0.25,1,0.05)
  
  re <- NULL  
  x_times <- 100
    
  p.score <- runif(1000)**MFtr
  n.score <- 1 - p.score
  scores  <- cbind(c(p.score, n.score),c(rep(1,1000), rep(2,1000)))
  scores  <- cbind(scores[,1], scores[,1], scores[,2])

  if(co%in%c("ACC", "PACC", "MS", "MS2", "T50", "MAX", "X")){
      TprFpr  <- apply(getTPRandFPRbyThreshold(scores), 2, as.numeric)
      #reTPRFPR <- microbenchmark(TprFpr  <- getTPRandFPRbyThreshold(scores), times = x_times)$time
      #xq <- quantile(reTPRFPR, seq(0,1,0.05))
      #reTPRFPR <- mean(reTPRFPR[reTPRFPR <= xq[18] & reTPRFPR >= xq[4]])
      #saveRDS(reTPRFPR, paste0("reTPRFPR.rds"))
      #return(0)
      #reTPRFPR <- readRDS("reTPRFPR.rds")
  }
  
  i <- 10

  print(paste0("Class distr. ", var_size))
  aux <- NULL
        
  vclass     <- c(rep(1,round(var_size*var_perc)), rep(2, round(var_size*(1-var_perc))))
  test.p     <- runif(var_size*var_perc)**MF
  test.n     <- 1 - runif( round(var_size*(1-var_perc), digits=0) )**MF
  test_set   <- c(test.p, test.n)

  thr <- 0.5

  print(paste0("Running ", co, " ..."))
  
  if(co == "SMM"){
    re <- microbenchmark(SMM_method_perf(p.score = p.score, n.score = n.score, test = test_set), times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  if(co == "HDy"){
    re <- microbenchmark(HDy_method_perf(p.score = p.score, n.score = n.score, test = test_set), times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  if(co == "SORD"){
    re <- microbenchmark(SORD_perf(p.score = p.score, n.score = n.score, test = test_set), times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  if(co == "DyS-TS"){
    re <- microbenchmark(DyS_method_perf(p.score = p.score, n.score = n.score, test = test_set, measure = "topsoe"), times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  if(co == "DyS-ORD"){
    re <- microbenchmark(DyS_method_perf(p.score = p.score, n.score = n.score, test = test_set, measure = "ord"), times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  if(co == "CC"){
    test_set <- as.numeric(test_set)
    #re <- microbenchmark(CC_method(test = test_set), times = x_times)$time
    re <- microbenchmark(sum(test_set >= 0.5), times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  

  if(co%in%c("PACC", "PCC")){
    reCalibA <- microbenchmark(calib <- calibrate(as.factor(scores[,3]), as.numeric(scores[,1]), class1=1, method="isoReg",assumeProbabilities=TRUE), times = x_times)$time
    xq <- quantile(reCalibA, seq(0,1,0.05))
    reCalibA <- mean(reCalibA[reCalibA <= xq[18] & reCalibA >= xq[4]])
    
    reCalibB <- microbenchmark(calibProbs<- applyCalibration(as.numeric(test_set), calib),times = x_times)$time
    xq <- quantile(reCalibB, seq(0,1,0.05))
    reCalibB <- mean(reCalibB[reCalibB <= xq[18] & reCalibB >= xq[4]])
  }
  
  if(co == "PCC"){
    re <- microbenchmark(PCC_method(calibProbs), times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]]) + reCalibB
  }
  if(co == "PACC"){
    re <- microbenchmark(PACC_method(test = calibProbs, TprFpr = TprFpr, thr=thr) , times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]]) + reCalibB #+ reTPRFPR
  }


  if(co == "ACC"){
    re <- microbenchmark(ACC_method(test = test_set, TprFpr = TprFpr, thr = thr) , times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]]) #+ reTPRFPR
  }
  if(co == "T50"){
    re <- microbenchmark(T50_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time 
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]]) #+ reTPRFPR
  }
  if(co == "X"){
    re <- microbenchmark(X_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time 
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]]) #+ reTPRFPR
  }
  if(co == "MAX"){
    re <- microbenchmark(MAX_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time 
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]]) #+ reTPRFPR
  }
  if(co == "MS"){
    re <- microbenchmark(MS_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]]) #+ reTPRFPR
  }
  if(co == "MS2"){
    re <- microbenchmark(MS2_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time 
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]]) #+ reTPRFPR
  }

  if(co == "EMQ"){
    testEMQ <- as.data.frame(cbind(test_set,1-test_set))
    names(testEMQ) <- c("1", "2")
    re <- microbenchmark(EMQ_method(train = 0, test = testEMQ) , times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  
   

  
  #aux <- rbind(aux, c("SMM", mean(reSMM$time)))
  #aux <- rbind(aux, c("HDy",mean(reHDy$time)))
  #aux <- rbind(aux, c("SORD", mean(reSORD$time)))
  #aux <- rbind(aux, c("DyS-TS", mean(reDyS$time)))
  #aux <- rbind(aux, c("DyS-ORD", mean(reDyS_ORD$time)))

  #aux <- rbind(aux, c("CC", mean(reCC$time)))
  #aux <- rbind(aux, c("PCC", mean(rePCC$time) + mean(reCalibA$time) + mean(reCalibB$time)))
  #aux <- rbind(aux, c("ACC", mean(reACC$time) + mean(reTPRFPR$time)))
  #aux <- rbind(aux, c("T50", mean(reT50$time) + mean(reTPRFPR$time)))
  #aux <- rbind(aux, c("X", mean(reX$time) + mean(reTPRFPR$time)))
  #aux <- rbind(aux, c("MAX", mean(reMAX$time) + mean(reTPRFPR$time)))
  #aux <- rbind(aux, c("PACC", mean(rePACC$time) + mean(reTPRFPR$time) + mean(reCalibA$time) + mean(reCalibB$time) ))
  #aux <- rbind(aux, c("MS", mean(reMS$time) + mean(reTPRFPR$time)))
  #aux <- rbind(aux, c("MS2", mean(reMS2$time) + mean(reTPRFPR$time)))
  #aux <- rbind(aux, c("EMQ", mean(reEMQ$time)))
  
  print(re)
  aux <- as.data.frame(t(c(co, re, var_size)))

  names(aux) <- c("counter", "time", "size")

  saveRDS(aux, paste0(co,"_",var_size, ".rds"))
  return(1)
}

re <- exec_eval_performance(ni, co)

