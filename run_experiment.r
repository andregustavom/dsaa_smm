source("./load_libraries.r")
library("batch")
parseCommandArgs()

run_experiment <- function(dts_id, counter){
  
  dts <- list.files("./datasets/")[dts_id]
  thr <- 0.5 # classifier threshold
  
  print(paste0("Quantifying ----->>>> ",dts))
  db_dir <- paste0("./models_train_test/",dts)
    
  if(!file.exists(db_dir)){
    dir.create(db_dir)
    # Training - 50% and Test - 50%
    db    <- as.data.frame(fread(paste0("./datasets/",dts)))
    db    <- db %>% mutate_if(is.character, as.factor)  
    db    <- as.data.frame(db)
    db$class <- as.factor(db$class)
    cv   <- createFolds(y = db$class, k = 2)
    train <- db[cv[[1]],]
    test  <- db[cv[[2]],]

     if(dts == "click-prediction.data"){
      train <- apply(train, 2, as.numeric)
      train <- as.data.frame(train)
      train$class <- as.factor(train$class)
      test <- apply(test, 2, as.numeric)
      test <- as.data.frame(test)
      test$class <- as.factor(test$class)
    }

    tryCatch(scorer <- randomForest::randomForest(class~., data=train, ntree=200), error = function(e) { print("ERROR 1 - scorer error!")})
    tryCatch(scores     <- getScore_using_K_folds(train, 10), error = function(e) { print("ERROR 2 - Scores error!")})
    tryCatch(TprFpr     <- mlquantify::getTPRandFPRbyThreshold(scores), error = function(e) { print("ERROR 3 - TprFpr error!")})
    
    # save in .RDS files scorer, training set, test set, scores and TprFpr for each dataset
    saveRDS(scorer, paste0("./models_train_test/",dts,"/classifier.rds"))
    saveRDS(train     , paste0("./models_train_test/",dts,"/train.rds"))
    saveRDS(test      , paste0("./models_train_test/",dts,"/test.rds"))
    saveRDS(scores    , paste0("./models_train_test/",dts,"/scores.rds"))
    saveRDS(TprFpr    , paste0("./models_train_test/",dts,"/TprFpr.rds"))
  }else{
      # load the .RDS file of each dataset
      scorer<- readRDS(paste0("./models_train_test/",dts,"/classifier.rds"))
      train     <- readRDS(paste0("./models_train_test/",dts,"/train.rds"))
      test      <- readRDS(paste0("./models_train_test/",dts,"/test.rds"))
      scores    <- readRDS(paste0("./models_train_test/",dts,"/scores.rds"))
      TprFpr    <- readRDS(paste0("./models_train_test/",dts,"/TprFpr.rds"))
      TprFpr    <- apply(TprFpr, 2, as.numeric)
    }

    db_files <- paste0("./models_train_test/",dts,"/sampled_scores.rds")
    # To reduce the time of experiments we use only a sample of the positive and negative scores. This step has only efect over DyS and HDy 
    if(length(which(counter%in%c("DyS-ORD", "MLQ", "RND", "SORD", "TOP", "RTQ")))>0 ){
      if(!file.exists(db_files)){
        ifelse(length(which(scores[,3]==1))<1000 , p_scores <- scores[which(scores[,3]==1),], p_scores <- scores[sample(which(scores[,3]==1),1000),])
        ifelse(length(which(scores[,3]==2))<1000 , n_scores <- scores[which(scores[,3]==2),], n_scores <- scores[sample(which(scores[,3]==2),1000),])
        scores   <- rbind(p_scores, n_scores)
        print(nrow(scores))
        saveRDS(scores    , paste0("./models_train_test/",dts,"/sampled_scores.rds"))
      }else{
          scores    <- readRDS(paste0("./models_train_test/",dts,"/sampled_scores.rds"))
          print("2000 scores loaded!")
        }
        }else{
          print("All scores used!")
        }

if(!dir.exists(paste0("./results_", counter)))  
    dir.create(paste0("./results_", counter))

  print(paste0("------>>>> Starting ", counter, " quantifier"))
  ran_f <- paste0("./results_", counter,"/QE_results_",counter, "_",dts,".rds")   
  print(ran_f) 
  if(!file.exists(ran_f)){
      print(paste0("Running ", counter))                
      re <- run_qntMethod_over_dataset(train = train, 
                                        test = test, 
                                   qntMethod = counter, 
                                      scorer = scorer, 
                                      TprFpr = TprFpr, 
                                      scores = scores, 
                                         thr = thr,
                                      dts_id = dts)
      re$dataset <- dts
      saveRDS(re, paste0("./results_", counter,"/QE_results_",counter, "_",dts,".rds")) 
  
  }else{  print(ran_f) }
  
  print(paste0(dts, " ----->>>> Concluded!"))
  
  return(1)
  
}

results <- run_experiment(dts_id = 10, counter="DyS-PS")

print("Finished!")