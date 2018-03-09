#' Apply the trained ForeseeModel on a TestObject
#'
#' The Foreseer applies the ForeseeModel that was trained on a FORESEE TrainObject to the test data to gain a prediction of the TestObject's response.
#'
#' @param TestObject Object that contains all data that the model is to be tested on, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param ForeseeModel Model that has been trained on a TrainObject with ForeseeTrain.
#'
#' @return \item{Foreseen}{Predicted drug response of the TestObject obtained by applying the ForeseeModel.}
#' @export

Foreseer <- function(TestObject, ForeseeModel, BlackBox){

  if (BlackBox=="lasso"){

    TestObject_test<- t(TestObject$Features)
    Foreseen <- predict(object=ForeseeModel, newx=TestObject_test)

  }
  else if (BlackBox=="elasticnet"){

    TestObject_test<- t(TestObject$Features)
    Foreseen <- predict(object=ForeseeModel, newx=TestObject_test)

  } else if(BlackBox=="rf"){
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
  # Plots

  # Update Object in the Environment
  return(Foreseen)
}

