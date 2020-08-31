library(microbenchmark)
library(philentropy)
library(CORElearn)
library("batch")
parseCommandArgs()

#----------------------------------------------------------------------------------------
discret_data <- function(x, inter){
  re <- NULL
  for(i in 2:length(inter)){
    re <- c(re, length(which(x >= inter[i-1] & x < inter[i]))/length(x) )
  }
  return(re)
}
#----------------------------------------------------------------------------------------
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
#----------------------------------------------------------------------------------------
SMM_method_perf <- function(p.score, n.score, test){
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
  alpha <- seq(0,1,by=0.01)
  b_sizes <- seq(10,110, length.out = 11)
  result <- c(1:11)
  for(hi in 1:length(b_sizes)){
    v_prob <- seq(0,1,length.out = b_sizes[hi])
    v_prob <- c(v_prob[-length(v_prob)], v_prob[length(v_prob)]+0.1)
    Sty_1 <- discret_data(p.score,v_prob)
    Sty_2 <- discret_data(n.score,v_prob)
    Uy    <- discret_data(test,v_prob)
  
    vDist <- c(1:length(alpha))
    for(k in 1:length(alpha)){
      aux <- (Sty_1*alpha[k])+ (Sty_2*(1-alpha[k]))
      vDist[k] <- as.numeric(distance(rbind(aux, Uy), method = "hellinger", p=0.5))
    }
    result[hi] <- alpha[which.min(vDist)]
  }  
  result <- median(result)
  result <- c(result, 1 - result)
  return(result)
}
#----------------------------------------------------------------------------------------
SORD_perf <- function(p.score, n.score, test){  
  f <- function(x){
    p_w <- x / length(p.score)
    n_w <- (1 - x) / length(n.score)
    t_w <- -1 / length(test)
    p <- cbind(p.score, rep(p_w, length(p.score)))
    n <- cbind(n.score, rep(n_w, length(n.score)))
    t <- cbind(test, rep(t_w, length(test)))

    v <- rbind(p, n, t)
    v <- v[order(v[,1]),]
    acc <- v[1,2]
    total_cost <- 0
  
    for( i in 2:nrow(v)){
      cost_mul <- v[i,1] - v[i - 1, 1] 
      total_cost <- total_cost + abs(cost_mul * acc)
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
  alpha <- seq(0,1,by=0.01)
  b_sizes <- seq(2,20,2)
  result <- c(1:10)
  for(hi in 1:10){
    breaks <- seq(0,1,length.out = b_sizes[hi]+1)
    breaks <- c(breaks[-length(breaks)], 1.1)
    
    Sty_1 <- discret_data_normalized(p.score,breaks)
    Sty_2 <- discret_data_normalized(n.score,breaks)
    Uy    <- discret_data_normalized(test,breaks)
  
    f <- function(x){
      return(distance_DYS(rbind((Sty_1*x)+ (Sty_2*(1-x)), Uy), method = measure, p=0.5))
    }
    result[hi] <- TernarySearch(0, 1, f)
  } 
  result <- median(result)
  result <- c(result, 1 - result)
  return(result)
}

#----------------------------------------------------------------------------------------
X_method <- function(test, TprFpr){
  x_sim <- cbind(round(1-TprFpr[,2],2), round(TprFpr[,3],2))
  idx_thr <- order(abs(x_sim[,1]-x_sim[,2]))[1]
  thr <- TprFpr[idx_thr,1]
  x <- as.data.frame(t(as.numeric(TprFpr[idx_thr,c(2,3)])))
  colnames(x) <- c("tpr", "fpr")
  TprFpr <- x
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
  return(result)
}
#----------------------------------------------------------------------------------------
PACC_method <- function(test, TprFpr, thr){
  x <- as.data.frame(t(as.numeric(TprFpr[which(TprFpr[,1] == round(thr,2)),c(2,3)])))
  colnames(x) <- c("tpr", "fpr")
  TprFpr <- x
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
#----------------------------------------------------------------------------------------
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
  result <- sum(test >= thr)/length(test)
  return(result)  
}
#----------------------------------------------------------------------------------------
ACC_method <- function(test, TprFpr, thr=NULL){
  dC <- sum(test >= thr)/length(test)
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
  TprFpr <- as.numeric(TprFpr[order(abs(as.numeric(TprFpr[,2])-0.5))[1],])
  dC <- CC_method(test, TprFpr[1])
  dC <- (dC[1] - TprFpr[3])/(TprFpr[2] - TprFpr[3])
  ifelse(dC < 0, dC <- 0, dC <- dC)
  ifelse(dC > 1, dC <- 1, dC <- dC)
  result <- c(dC, 1 - dC)
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
#########################################################################################
# MAIN PROCEDURE

exec_eval_performance <- function(ni, co){

  var_perc <- 0.5#seq(0,1,0.05)
  var_size <- ni
  x_times  <- 1000 #Increase this value for more confident results (repetitions)
  MF       <- 0.5
  MFtr     <- 0.5
  re <- NULL  
  # Generating artificial scores (1000 scores per class)
  p.score <- runif(1000)**MFtr
  n.score <- 1 - p.score
  scores  <- cbind(c(p.score, n.score),c(rep(1,1000), rep(2,1000)))
  scores  <- cbind(scores[,1], scores[,1], scores[,2])

  if(co%in%c("ACC", "PACC", "MS", "MS2", "T50", "MAX", "X"))
      TprFpr  <- apply(mlquantify::getTPRandFPRbyThreshold(scores), 2, as.numeric)

  print(paste0("Class distr. ", var_size))
        
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
  if(co %in%c("CC", "QT")){
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
    re <- mean(re[re <= xq[18] & re >= xq[4]]) 
  }
  if(co == "T50"){
    re <- microbenchmark(T50_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time 
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  if(co == "X"){
    re <- microbenchmark(X_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time 
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  if(co == "MAX"){
    re <- microbenchmark(MAX_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time 
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  if(co == "MS"){
    re <- microbenchmark(MS_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  if(co == "MS2"){
    re <- microbenchmark(MS2_method(test = test_set, TprFpr = TprFpr) , times = x_times)$time 
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }

  if(co == "EMQ"){
    testEMQ <- as.data.frame(cbind(test_set,1-test_set))
    names(testEMQ) <- c("1", "2")
    re <- microbenchmark(EMQ_method(train = 0, test = testEMQ) , times = x_times)$time
    xq <- quantile(re, seq(0,1,0.05))
    re <- mean(re[re <= xq[18] & re >= xq[4]])
  }
  
  aux <- as.data.frame(t(c(co, re, var_size)))
  names(aux) <- c("counter", "time", "size")
  #Saving the results for each test set size and quantifier separately
  saveRDS(aux, paste0("./performance/",co,"_",var_size, ".rds"))
  return(1)
}
#ni <- 1000 #Size of test set
#co <-  "QT" #Quantifier name
#re <- exec_eval_performance(ni, co)


counters <- c("T50","CC", "SORD","EMQ", "PACC",  "MAX","ACC","SMM","PCC","X","MS", "MS2","HDy","DyS-TS", "QT") 
test_set_sizes <- c(100,1000,seq(10000, 100000, 10000))

# uncommnet the next lines for getting the timing cost for all quantifiers
#for(ni in test_set_sizes)
#  for(co in counters)
#    re <- exec_eval_performance(ni, co)
    




