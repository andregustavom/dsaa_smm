#' Function used to perform experiments, given a quantifier method, varying the test size and the class distribution
#' 
#' Given a quantifier method (qntMethod), a scorer, and a test set, that function extract from the test set several 
#' subsets varying the size and for each size the class distribution is also varied
#' @author Andre G Maletzke - andregustavom at gmail.com
#' @param train data.frame containing the training set
#' @param test data.frame containing the test set
#' @param qntMethod the quantifier algorithm name
#' @param scorer A scorer built on a training set
#' @param TprFpr data.frame composed by True Positive and False Positive rates estimated on training set
#' @param scores data.frame contaning the scores for each instance of the training set, estimated by 10-folds
#' @param ml_alg_cl machine learning algorithm name
#' @param thr threshold used for some quantifier methos such as CC
#' @param measure parameter of the DyS framework
#' @return The class distribution for each simulated subset extracted from the test set
#' @export

run_qntMethod_over_dataset <- function(train, test, qntMethod, scorer, TprFpr, scores, thr=0.5, dts_id){

  
  var_perc <- seq(0,1,0.01)                       # range of values used to evaluate the impact of class distribution on error
  var_size <- c(seq(10,100,10), seq(200,500,100)) # range of test set sizes
  n_tests  <- 10                                  # for each combination of var_perc and var_size we build n_tests aleatory scenarios
  results  <- NULL
  
  idx_pos <- which(test$class==1) # positive class (1)
  idx_neg <- which(test$class==2) # negative class (2)
  
  DyS_conf <- c("topsoe", "jensen_difference", "prob_symm", "ord")
  names(DyS_conf) <- c("TS", "JD", "PS", "ORD")
  
  measure <- NA
  qi <- qntMethod
  aux <- strsplit(qntMethod, "DyS-")[[1]]
  if(length(aux)>1 ){
    qi <- "DyS"
    measure   <- DyS_conf[aux[2]]
  }
  
  if(!dir.exists(paste0("./tmp_", qntMethod)))  
    dir.create(paste0("./tmp_", qntMethod))
  
  if(!dir.exists(paste0("./tmp_",qntMethod,"/", dts_id)))
    dir.create(paste0("./tmp_",qntMethod,"/", dts_id))
  
  for(k in 1:length(var_size)){ 
    if(file.exists(paste0("./tmp_",qntMethod,"/", dts_id, "/", var_size[k], ".rds"))){
      print(paste0("Load data of size", var_size[k]))
      auxre <- readRDS(paste0("./tmp_",qntMethod,"/", dts_id, "/", var_size[k], ".rds"))
      results <- rbind(results, auxre)
    }else{
      print(paste0("Test size ", var_size[k]))
      for(i in 1:length(var_perc)){      
        for(j in 1:n_tests){
              n_pos <- round(i*k)
              n_neg <- k - n_pos
              test_set <- rbind(test[sample(idx_pos,n_pos),], test[sample(idx_neg,n_neg),])
              predTest <- predict(scorer, test_set, type = c("prob"))
              
              # PCC and PACC require calibrated scores
              if(qi%in%c("PCC", "PACC")){
                calib     <- calibrate(as.factor(scores[,3]), as.numeric(scores[,1]), class1=1, method="isoReg",assumeProbabilities=TRUE)
                calibProbs<- applyCalibration(as.numeric(predTest[,1]), calib)
              }else{
                calibProbs <- as.numeric(predTest[,1])
              }
              
              qnt_re <- apply.qntMethod(qntMethod = qi, 
                                        p.score = scores[scores[,3]==1,1],
                                        n.score = scores[scores[,3]==2,1], 
                                        test = calibProbs, 
                                        TprFpr = TprFpr, 
                                        thr = thr,
                                        measure = measure,
                                        train = train)
              
              freq_PRE <- round(qnt_re,2)
              freq_REAL<- table(test_set$class)/sum(table(test_set$class))
              
              results  <- rbind(results, unlist(c(freq_REAL,
                                                  freq_PRE,
                                                  round(abs(freq_REAL[1]-freq_PRE[1]),2),
                                                  nrow(test_set), 
                                                  qntMethod)))
              
        }      
      }
      saveRDS(results[which(results[,6]==var_size[k]),], paste0("./tmp_",qntMethod,"/", dts_id, "/", var_size[k], ".rds"))
    }
  }  
  results <- as.data.frame(results)
  names(results) <- c(paste("R", as.character(1:2), sep="_"), 
                      paste("P", as.character(1:2), sep="_"), 
                      "MAE-C1",
                      "Test_Size",
                      "Qnt")
  return(results)
}


