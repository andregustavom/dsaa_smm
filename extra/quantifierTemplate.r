apply.qntMethod <- function(qntMethod, 
                            p.score,
                            n.score,
                            test, 
                            TprFpr = NULL, 
                            thr = NULL,
                            measure="hellinger",
                            train = NULL
                            ){
  if(qntMethod == "CC")
    return(mlquantify::CC(test = test, thr = thr))
  
  if(qntMethod == "PCC")
    return(mlquantify::PCC(test))

  if(qntMethod == "ACC")
    return(mlquantify::ACC(test = test, TprFpr = TprFpr, thr = thr))
  
  if(qntMethod == "T50")
    return(T50_new(test = test, TprFpr = TprFpr))
    #return(mlquantify::T50(test = test, TprFpr = TprFpr))
  
  if(qntMethod == "X")
    return(mlquantify::X(test = test, TprFpr = TprFpr))
  
  if(qntMethod == "MAX")
    return(mlquantify::MAX(test = test, TprFpr = TprFpr))
  
  if(qntMethod == "PACC")
    return(mlquantify::PACC(test = test, TprFpr = TprFpr, thr=thr))
  
  if(qntMethod == "HDy")
    return(HDy_method(p.score = p.score, n.score = n.score, test = test))

  if(qntMethod == "EMQ"){
    test <- as.data.frame(cbind(test,1-test))
    names(test) <- c("1", "2")
    return(mlquantify::EMQ(train = train, test = test))
  }
  
  if(qntMethod == "SMM")
    return(SMM(p.score = p.score, n.score = n.score, test = test))
  
  if(qntMethod == "HDy-LP")
    return(mlquantify::HDy_LP(p.score = p.score, n.score = n.score, test = test))
  
  if(qntMethod == "DyS")
    return(mlquantify::DyS(p.score = p.score, n.score = n.score, test = test, measure = measure))

  if(qntMethod == "SORD")
    return(mlquantify::SORD(p.score = p.score, n.score = n.score, test = test))
  
  if(qntMethod == "MS")
    return(mlquantify::MS(test = test, TprFpr = TprFpr))
    
  if(qntMethod == "MS2")
    return(mlquantify::MS2(test = test, TprFpr = TprFpr))
    
  print("ERROR - Quantification method was not applied!") 
  return(NULL)
}

#----------------------------------------------------------------------------------------------------------------------
