source("./load_libraries.r")
library("batch")
parseCommandArgs()

run_experiment_QT <- function(dts_id, counter){
  
  ml_alg_cl    <- "RF" # machine learning algorithm used to induce the classifier for each dataset
  dts          <- list.files("./datasets/")[dts_id]
  thr          <- 0.5 # classifier threshold
  
  print(paste0("Quantifying ----->>>> ",dts))
  db_files <- paste0("./models_train_test/", dts, "/train_arff.arff")

  train     <- readRDS(paste0("./models_train_test/",dts,"/train.rds"))
  test      <- readRDS(paste0("./models_train_test/",dts,"/test.rds"))

  classifier<- paste0("classifier_RF_QT_", dts)
  
  if(!file.exists(db_files)){
    
    pp <- which(train$class==1)
    nca <- 1:nrow(train)
    nca[pp] <- 0
    nca[-pp] <- 1
    
    nn <- names(train)
    train <- train[,-1]
    train <- cbind(nca,train)
    names(train) <- nn
    
    train$class <- as.factor(train$class)
    nn <- names(train)
    nn <- c(nn[-1], nn[1])
    train <- cbind(train[,-1], train[,1])
    names(train) <- nn
    
    write.arff(train, paste0("./models_train_test/", dts, "/train_arff.arff"))
    
    command <- paste0("java -Xmx6G -cp quantify.jar:weka.jar:. weka.classifiers.trees.RandomForest -I 500 -P 1 -t ", "./models_train_test/",dts,"/train_arff", 
                      ".arff", " -d ./models_train_test/",dts, "/", classifier)
    
    system(command)
  }
  
  if(!dir.exists(paste0("./results_", counter)))  
    dir.create(paste0("./results_", counter))
  
  print(paste0("------>>>> Starting ", counter, " quantifier"))
  ran_f <- paste0("./results_", counter,"/QE_results_",counter, "_",dts,".rds")   
  print(ran_f) 
  if(!file.exists(ran_f)){
      print(paste0("Running ", counter))                
    re <- run_qntMethod_over_dataset_QT(test = test, 
                                     qntMethod = counter, 
                                     classifier = classifier, 
                                     dts_id=dts)
    re$dataset <- dts
    saveRDS(re, paste0("./results_", counter,"/QE_results_",counter, "_",dts,".rds")) 
    
  }else{  print(ran_f) }
  
  print(paste0(dts, " ----->>>> Concluded!"))
  
  return(1)
  
}

#results <- run_experiment_QT(dts_id=ndts, counter=counter, measure=measure)
#results <- run_experiment_QT(dts_id=ndts, counter=counter)
results <- run_experiment_QT(dts_id = 2, counter="QT")

print("Finished!")