#' Apply the trained ForeseeModel on a TestObject
#'
#' The Foreseer applies the ForeseeModel that was trained on the features of a FORESEE TrainObject to the test data to gain a prediction of the TestObject's drug response.
#'
#' @param TestObject Object that contains all data that the model is to be tested on, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type, etc. ) and drug response data
#' @param ForeseeModel Model that has been trained on a TrainObject with the function ForeseeTrain().
#'
#' @return \item{Foreseen}{Predicted drug response of the samples listed in the TestObject obtained by applying the ForeseeModel to the features of the TestObject.}
#' @export

#########################
# This file is part of the FORESEE R-package
# File authors: Lisa-Katrin Turnhoff <turnhoff@combine.rwth-aachen.de> and Ali Hadizadeh Esfahani <hadizadeh@combine.rwth-aachen.de>
# Distributed under the GNU General Public License v3.0.(http://www.gnu.org/licenses/gpl-3.0.html)
#########################

Foreseer <- function(TestObject, ForeseeModel, BlackBox){
  if (is.function(BlackBox)){
    TestObject_test<- as.data.frame(as.matrix(t(TestObject$Features)))
    # Check if there are still duplications (to avoid that the pipeline breaks)
    # Just take the first occuring gene name (here: in columns)
    TestObject_test <- TestObject_test[,!duplicated(colnames(TestObject_test))]

    Foreseen <- BlackBox(ForeseeModel, TestObject_test)
  }
  else if (any(BlackBox==c("lasso","elasticnet"))){

    TestObject_test<- t(TestObject$Features)
    Foreseen <- predict(object=ForeseeModel, newx=TestObject_test)
    # Fixing some warnings in Validator: Some functions in validation, e.g. pROC::roc expect a
    # numeric vector, while predict of lasso and elasticnet return a matrix, so:
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
    # Check if there are still duplications (to avoid that the pipeline breaks)
    # Just take the first occuring gene name (here: in columns)
    TestObject_test <- TestObject_test[,!duplicated(colnames(TestObject_test))]

    # randomForest and ranger have a problem with features that are named starting with an integer, like Entrez IDs that we used
    # Have to fix the names:
    names(TestObject_test) <- make.names(names(TestObject_test))

    Foreseen <- predict(ForeseeModel, TestObject_test)
  } else if(BlackBox=="rf_ranger"){
    TestObject_test<- as.data.frame(as.matrix(t(TestObject$Features)))
    # Check if there are still duplications (to avoid that the pipeline breaks)
    # Just take the first occuring gene name (here: in columns)
    TestObject_test <- TestObject_test[,!duplicated(colnames(TestObject_test))]

    # randomForest and ranger have a problem with features that are named starting with an integer, like Entrez IDs that we used
    # Have to fix the names:
    names(TestObject_test) <- make.names(names(TestObject_test))

    Foreseen <- predict(ForeseeModel, TestObject_test)$predictions
  } else{
  TestObject_test<- as.data.frame(as.matrix(t(TestObject$Features)))
  # Check if there are still duplications (to avoid that the pipeline breaks)
  # Just take the first occuring gene name (here: in columns)
  TestObject_test <- TestObject_test[,!duplicated(colnames(TestObject_test))]

  Foreseen <- predict(ForeseeModel, TestObject_test)
  }

  # Update Object in the Environment
  return(Foreseen)
}

