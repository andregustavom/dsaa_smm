#' Function used to perform experiments, given a quantifier method, varyng the test size and the class distribution
#' 
#' Given a quantifier method (qntMethod), a classifier, and a test set, that function extract from the test set several 
#' subsets varying the size and for each size the class distribution is also varied
#' @author Andre G Maletzke - andregustavom at gmail.com
#' @param test data.frame containing the test set
#' @param qntMethod the quantifier algorithm name
#' @param classifier A classifier built on a training set
#' @return The class distribution for each simulated subset extracted from the test set
#' @export

run_qntMethod_over_dataset_QT <- function(test, qntMethod, classifier, dts_id){
  
  
  var_perc <- seq(0,1,0.01)                       # range of values used to evaluate the impact of class distribution on error
  var_size <- c(seq(10,100,10), seq(200,500,100)) # range of test set sizes
  n_tests  <- 10                                  # for each combination of var_perc and var_size we build n_tests aleatory scenarios
  
  results     <- NULL
  
  idx_pos <- which(test$class==1) # positive class (1)
  idx_neg <- which(test$class==2) # negative class (2)

  if(!dir.exists(paste0("./tmp_", qntMethod)))  
    dir.create(paste0("./tmp_", qntMethod))
  
  if(!dir.exists(paste0("./tmp_",qntMethod,"/", dts_id)))
    dir.create(paste0("./tmp_",qntMethod,"/", dts_id))

  if(!dir.exists(paste0("./test_QT/",dts_id)))  
    dir.create(paste0("./test_QT/",dts_id))

  
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
          
          pp <- which(test_set$class==1)
          nca <- 1:nrow(test_set)

          if(length(pp)==0){
            nca[] <- 1
          }else{
            nca[pp] <- 0
            nca[-pp] <- 1
          }
          
          nn <- names(test_set)
          test_set <- test_set[,-1]
          test_set <- cbind(nca,test_set)
          names(test_set) <- nn    
          
          test_set$class <- factor(test_set$class, levels=c("0", "1"))
          nn <- names(test_set)
          nn <- c(nn[-1], nn[1])
          test_set <- cbind(test_set[,-1], test_set[,1])
          names(test_set) <- nn
           
          write.arff(test_set, paste0("./test_QT/",dts_id,"/test_set",".arff"))                  
          command <- paste0("java  -Xmx5G -cp quantify.jar:weka.jar:. weka.classifiers.trees.RandomForest -l ./models_train_test/",dts_id, "/", classifier,  
                            " -T ", "./test_QT/",dts_id,"/test_set",".arff")          
          system(paste0(command," > ", "./test_QT/",dts_id,"/","re.txt"))          
          x <- read.delim(paste0("./test_QT/",dts_id,"/","re.txt"))      
          
          pos <- as.vector(na.omit(as.numeric(strsplit(as.character(x[4,])," ")[[1]])))[1]          
          pos <- pos + as.vector(na.omit(as.numeric(strsplit(as.character(x[5,])," ")[[1]])))[1]          
          pos <- pos/nrow(test_set)          
          freq_PRE <- as.data.frame(round(cbind(pos, 1-pos),2))         
          names(freq_PRE) <- c("1", "2")
          freq_REAL<- table(test_set$class)/sum(table(test_set$class))
          results  <- rbind(results, unlist(c(freq_REAL,
                                              freq_PRE,
                                              round(abs(freq_REAL[1]-freq_PRE[1]),2),
                                              nrow(test_set), 
                                              qntMethod)))
          browser()
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


