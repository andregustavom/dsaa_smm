#' Quantify a set of events by HDy method
#' 
#' Quantify a set of events given based on a fitted Weka classifier applying the Expectation Maximization for Quantification (EMQ) method
#' @param classifier A Weka classifier built on a training set
#' @param train A labeled set used to built the classifier
#' @param test A unlabeled set
#' @return The class distribution in the test set
#' @export

HDy_method <- function(p.score, n.score, test){
  
  alpha <- seq(0,1,by=0.01)
  
  sc_1 <- p.score
  sc_2 <- n.score
  
  Usc <- cbind(test, test)
  
  b_sizes <- seq(10,110, length.out = 11)
  
  result <- NULL
  
  for(hi in 1:length(b_sizes)){
    
    v_prob <- seq(0,1,length.out = b_sizes[hi])
    
    v_prob <- c(v_prob[-length(v_prob)], v_prob[length(v_prob)]+0.1)
      
    Sty_1 <- discret_data(sc_1,v_prob)
    Sty_2 <- discret_data(sc_2,v_prob)
    
    Uy <- discret_data(Usc[,1],v_prob)
    
    vDist <- NULL
    vDistAll <- NULL
    
    for(k in 1:length(alpha)){
      
      aux <- (Sty_1*alpha[k])+ (Sty_2*(1-alpha[k]))
      vDist <- c(vDist, as.numeric(distance(rbind(aux, Uy), method = "hellinger", p=0.5)))
    }
    
    vDistAll <- c(vDistAll, vDist)
    result <- c(result, alpha[which.min(vDist)])
    
  }  
  
  result <- median(result)
  result <- c(result, 1 - result)
  names(result) <- c("1", "2")
  
  return(list(round(result,2),min(vDistAll)))
    
}

HellingerDistance <-function(pd, qd){
  #takes two equal sized vectors and calculates the hellinger distance between the vectors
  
  # hellinger distance function
  vsum <- 0
  for(i in 1:length(pd))
    vsum <- vsum + (sqrt(pd[i]/sum(pd))- sqrt(qd[i]/sum(qd)))^2
    
  return(sqrt(vsum))
  
  
  return(sqrt(sum(((sqrt(pd) - sqrt(qd))^2)))/sqrt(2))
  
}

HDy_method_TS <- function(p.score, n.score, test){
  
  alpha <- seq(0,1,by=0.01)
  
  sc_1 <- p.score
  sc_2 <- n.score
  
  Usc <- cbind(test, test)
  
  b_sizes <- seq(10,110, length.out = 11)
  
  result <- NULL
  
  for(hi in 1:length(b_sizes)){
    
    v_prob <- seq(0,1,length.out = b_sizes[hi])
    
    v_prob <- c(v_prob[-length(v_prob)], v_prob[length(v_prob)]+0.1)
    
    Sty_1 <- discret_data(sc_1,v_prob)
    Sty_2 <- discret_data(sc_2,v_prob)
    
    Uy <- discret_data(Usc[,1],v_prob)
    
    vDist <- NULL
    vDistAll <- NULL
    
    f <- function(x){
      return(as.numeric(distance(rbind((Sty_1*x)+ (Sty_2*(1-x)), Uy), method = "hellinger", p=0.5)))
    }
    
    result <- c(result, TernarySearch(0, 1, f, 1e-2))
    vDistAll <- c(vDistAll, f(result))
    
  }  
  
  result <- median(result)
  result <- c(result, 1 - result)
  names(result) <- c("1", "2")
  
  return(list(round(result,2),min(vDistAll)))
  
}

