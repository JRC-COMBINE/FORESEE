#' Apply the trained ForeseeModel on a TestObject
#'
#' The Foreseer applies the ForeseeModel that was trained on a FORESEE TrainObject to the test data to gain a prediction of the TestObject's response.
#'
#' @param TestObject Object that contains all data that the model is to be tested on, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type) and drug response data
#' @param ForeseeModel Model that has been trained on a TrainObject with the function ForeseeTrain.
#'
#' @return \item{Foreseen}{Predicted drug response of the samples listed in the TestObject obtained by applying the ForeseeModel to the molecular data.}
#' @export

Foreseer <- function(TestObject, ForeseeModel, BlackBox){
  if (is.function(BlackBox)){
    TestObject_test<- as.data.frame(as.matrix(t(TestObject$Features)))
    # For some weird reason the object still contains duplicates? Check duplication handler
    # Just take the first occuring gene name (here: in columns!) for now
    TestObject_test <- TestObject_test[,!duplicated(colnames(TestObject_test))]

    Foreseen <- BlackBox(ForeseeModel, TestObject_test)
  }
  else if (any(BlackBox==c("lasso","elasticnet"))){

    TestObject_test<- t(TestObject$Features)
    Foreseen <- predict(object=ForeseeModel, newx=TestObject_test)
    #Fixing some warnings in Validator: Some functions in validation, e.g. pROC::roc expect a
    #numeric vector, while predict of lasso and elasticnet return a matrix, so:
    if(is.matrix(Foreseen) & ncol(Foreseen) == 1){
      Foreseen <- as.numeric(Foreseen)
    } else {
      stop("Unexpected Foreseen class and/or dimensions!")
    }

  }
  else if(BlackBox=="tandem"){

    TestObject_test<- t(TestObject$Features)
    ForeseenTmp <- predict(object=ForeseeModel, newx=TestObject_test)
    Foreseen <- ForeseenTmp@x
    names(Foreseen) <- ForeseenTmp@Dimnames[[1]]

  }
  else if(BlackBox=="rf"){
    TestObject_test<- as.data.frame(as.matrix(t(TestObject$Features)))
    # For some weird reason the object still contains duplicates? Check duplication handler
    # Just take the first occuring gene name (here: in columns!) for now
    TestObject_test <- TestObject_test[,!duplicated(colnames(TestObject_test))]

    ##randomForest and ranger have a problem with features that are named starting with an integer, like Entrez IDs that we used
    #Have to fix the names:
    names(TestObject_test) <- make.names(names(TestObject_test))

    Foreseen <- predict(ForeseeModel, TestObject_test)
  } else if(BlackBox=="rf_ranger"){
    TestObject_test<- as.data.frame(as.matrix(t(TestObject$Features)))
    # For some weird reason the object still contains duplicates? Check duplication handler
    # Just take the first occuring gene name (here: in columns!) for now
    TestObject_test <- TestObject_test[,!duplicated(colnames(TestObject_test))]

    ##randomForest and ranger have a problem with features that are named starting with an integer, like Entrez IDs that we used
    #Have to fix the names:
    names(TestObject_test) <- make.names(names(TestObject_test))

    Foreseen <- predict(ForeseeModel, TestObject_test)$predictions
  } else{
  TestObject_test<- as.data.frame(as.matrix(t(TestObject$Features)))
  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TestObject_test <- TestObject_test[,!duplicated(colnames(TestObject_test))]

  Foreseen <- predict(ForeseeModel, TestObject_test)
  }

  # Update Object in the Environment
  return(Foreseen)
}

