#' Train a Black Box Model for Drug Efficacy Prediction
#'
#' The BlackBoxFilter applies a machine learning algorithm to the data of the TrainObject to create a model that is predictive of the drug response.
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param BlackBox Modeling algorithm for training:
#' The function 'linear' fits a linear regression model to the training data,
#' The function 'ridge' fits a linear ridge regression model by Cule et al. (2012) to the training data,
#' The function 'lasso' fits a lasso regression model from the glmnet package by Friedman et al. (2008) to the training data,
#' The function 'elasticnet' fits an elastic net regression model from the glmnet package by Friedman et al. (2008) to the training data,
#' The function 'svm' fits a support vector regression model from the e1071 package by Meyer and Chih-Chung (2017) to the training data,
#' The function 'rf' fits a random forest regression model by Breiman (2001) to the training data
#' The function 'rf_ranger' fits a fast random forest regression model by Marvin N. Wright (2018) to the training data
#' The function 'tandem' fits a two-stage regression model by Nanne Aben (2017) to the training data
#'
#'
#' @param nfoldCrossvalidation # folds to use for crossvalidation while training the model. If put to zero, the complete data of the TrainObject is used for training.

#' @return \item{ForeseeModel}{A black box model trained on the TrainObject data that can be applied to new test data.}
#'         \item{TrainObject}{The TrainObject that was used to train the model.}
#' @export


# ToDo: User-defined function


BlackBoxFilter <- function(TrainObject, BlackBox = "ridge", nfoldCrossvalidation = 1, ...){
  if (nfoldCrossvalidation==1 && class(nfoldCrossvalidation)=="numeric"){
    UseMethod("BlackBoxFilter", object = BlackBox)
  }

  else if (nfoldCrossvalidation>1 && class(nfoldCrossvalidation)=="numeric"){

    # Establish backup in environment
    TrainObject_backup<-TrainObject
    assign("TrainObject_backup", value = TrainObject, envir = parent.frame())

    train <- TrainObject
    test  <- TrainObject

    # Shuffle the data
    sample <- sample.int(n = ncol(TrainObject_backup$Features))
    TrainObject$Features<-TrainObject_backup$Features[,sample]

    # Create n equally size folds
    folds <- cut(seq(1,ncol(TrainObject_backup$Features)),breaks=nfoldCrossvalidation,labels=FALSE)

    Performance_bestfold <- 0
    ForeseeModel_bestfold <- list()
    test_ids_bestfold <- integer()
    # Perform 10 fold cross validation
    for (i in 1:nfoldCrossvalidation){

      #Segement data
      test_ids <-  which(folds==i,arr.ind=TRUE)
      test$Features <- TrainObject_backup$Features[,test_ids]
      test$DrugResponse <- TrainObject_backup$DrugResponse[test_ids]
      train$Features <- TrainObject_backup$Features[,-test_ids]
      train$DrugResponse <- TrainObject_backup$DrugResponse[-test_ids]

      # Update TrainObject in environment
      assign("TrainObject", value = train, envir = parent.frame())

      #Train on subset of data (as the TrainObject is now updated in the environment)
      UseMethod("BlackBoxFilter", object = BlackBox)

      # Apply Model to cross-validation test set
      ForeseeTest(test, ForeseeModel, BlackBox, Evaluation=Crossvalidation_Criterion)
      # Compare current performance to previous once to just keep the best one in the environment
      if (Performance>Performance_bestfold){
        # Update objects in the environment
        assign("Performance_bestfold", value = Performance, envir = parent.frame())
        assign("ForeseeModel_bestfold", value = ForeseeModel, envir = parent.frame())
        assign("test_ids_bestfold", value = test_ids, envir = parent.frame())
      }

      TrainObject$Features <- TrainObject_backup$Features[,!test_ids_bestfold]
      TrainObject$DrugResponse <- TrainObject_backup$DrugResponse[!test_ids_bestfold]

      assign("ForeseeModel", value = ForeseeModel_bestfold, envir = parent.frame())
      assign("TrainObject", value = TrainObject, envir = parent.frame())

    }

  }

  else {
    stop("nfoldCrossvalidation needs to be a positive integer")
  }
}

BlackBoxFilter.character <- function(TrainObject, BlackBox, nfoldCrossvalidation){
  class(BlackBox) <- BlackBox
  UseMethod("BlackBoxFilter", object = BlackBox)
}

# BlackBoxFilter.function <- function(TrainObject, BlackBox, nfoldCrossvalidation){
#   UseMethod("BlackBoxFilter", object = BlackBox)
# }

################################################################################
### Function "linear" to train a linear regression model
BlackBoxFilter.linear <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$Features),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TrainObject_train <- as.data.frame(TrainObject_train)
  TrainObject_train <- TrainObject_train[,!duplicated(colnames(TrainObject_train))]
  TrainObject_train <- TrainObject_train[!duplicated(rownames(TrainObject_train)),]
  lm_fit <- lm(formula = DrugResponse~., TrainObject_train)

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <- TrainObject_train
  assign("TrainObject", value = TrainObject, envir = parent.frame())
  assign("ForeseeModel", value = lm_fit, envir = parent.frame())
}


################################################################################
### Function "ridge" to train a linear ridge regression model
BlackBoxFilter.ridge <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$Features),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TrainObject_train <- as.data.frame(TrainObject_train)
  TrainObject_train <- TrainObject_train[,!duplicated(colnames(TrainObject_train))]
  TrainObject_train <- TrainObject_train[!duplicated(rownames(TrainObject_train)),]

  # Ridge package by Cule, E. and De Iorio, M., A semi-automatic method to guide the choice of ridge parameter in ridge regression. (2012) arXiv:1205.0686v1
  require(ridge)
  ridge_fit <- linearRidge(formula = DrugResponse~., TrainObject_train)

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <- TrainObject_train
  assign("TrainObject", value = TrainObject, envir = parent.frame())
  assign("ForeseeModel", value = ridge_fit, envir = parent.frame())
}


################################################################################
### Function "lasso" to train a lasso regression model
BlackBoxFilter.lasso <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$Features),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # Package glmnet by Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, https://web.stanford.edu/~hastie/Papers/glmnet.pdf
  require(glmnet)
  lasso_fit <- glmnet(x = t(TrainObject$Features), y=TrainObject$DrugResponse, alpha = 1,lambda=cv.glmnet(x = t(TrainObject$Features), y=TrainObject$DrugResponse, alpha = 1)$lambda.min)

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <- TrainObject_train
  assign("TrainObject", value = TrainObject, envir = parent.frame())
  assign("ForeseeModel", value = lasso_fit, envir = parent.frame())
}


################################################################################
### Function "elasticnet" to train a elasticnet regression model
BlackBoxFilter.elasticnet <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$Features),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TrainObject_train <- as.data.frame(TrainObject_train)
  TrainObject_train <- TrainObject_train[,!duplicated(colnames(TrainObject_train))]
  TrainObject_train <- TrainObject_train[!duplicated(rownames(TrainObject_train)),]

  # Package glmnet by Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, https://web.stanford.edu/~hastie/Papers/glmnet.pdf
  require(glmnet)
  require(glmnetUtils)
  elasticnet_fit <- glmnet(x = t(TrainObject$Features), y=TrainObject$DrugResponse,  alpha = 0.5 ,lambda=cv.glmnet(x = t(TrainObject$Features), y=TrainObject$DrugResponse, alpha = 0.5)$lambda.min)

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <- TrainObject_train
  assign("TrainObject", value = TrainObject, envir = parent.frame())
  assign("ForeseeModel", value = elasticnet_fit, envir = parent.frame())
}


################################################################################
### Function "svm" to train a support vector regression model
BlackBoxFilter.svm <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$Features),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TrainObject_train <- as.data.frame(TrainObject_train)
  TrainObject_train <- TrainObject_train[,!duplicated(colnames(TrainObject_train))]
  TrainObject_train <- TrainObject_train[!duplicated(rownames(TrainObject_train)),]

  require(e1071)
  svm_fit <- svm(formula = DrugResponse~., data=data.frame(TrainObject_train))

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <- TrainObject_train
  assign("TrainObject", value = TrainObject, envir = parent.frame())
  assign("ForeseeModel", value = svm_fit, envir = parent.frame())
}


################################################################################
### Function "rf" to train a random forest regression model
BlackBoxFilter.rf <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$Features),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TrainObject_train <- as.data.frame(TrainObject_train)
  TrainObject_train <- TrainObject_train[,!duplicated(colnames(TrainObject_train))]
  TrainObject_train <- TrainObject_train[!duplicated(rownames(TrainObject_train)),]

  # Random Forest Package by Breiman, L. (2001), Random Forests, Machine Learning 45(1), 5-32.
  require(randomForest)
  rf_fit <- randomForest(formula = DrugResponse~., data.frame(TrainObject_train))

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <- TrainObject_train
  assign("TrainObject", value = TrainObject, envir = parent.frame())
  assign("ForeseeModel", value = rf_fit, envir = parent.frame())
}

################################################################################
### Function "rf" to train a random forest regression model
BlackBoxFilter.rf_ranger <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$Features),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TrainObject_train <- as.data.frame(TrainObject_train)
  TrainObject_train <- TrainObject_train[,!duplicated(colnames(TrainObject_train))]
  TrainObject_train <- TrainObject_train[!duplicated(rownames(TrainObject_train)),]

  # Random Forest Package by Marvin N. Wright (2018)
  require(ranger)
  rf_ranger_fit <- ranger(formula = DrugResponse~., data.frame(TrainObject_train))

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <- TrainObject_train
  assign("TrainObject", value = TrainObject, envir = parent.frame())
  assign("ForeseeModel", value = rd_ranger_fit, envir = parent.frame())

}


################################################################################
### Function "tandem" to train a lasso regression model
BlackBoxFilter.tandem <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(t(TrainObject$Features))
  upstream_index <- (TrainObject$FeatureTypes[1,])!="GeneExpression"

  # Package glmnet by Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, https://web.stanford.edu/~hastie/Papers/glmnet.pdf
  require(TANDEM)
  require(glmnet)

  tandem_fit <- tandem(x = TrainObject_train, y=TrainObject$DrugResponse, upstream=upstream_index)

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <- TrainObject_train
  TrainObject[["UpstreamIndex"]] <- upstream_index
  assign("TrainObject", value = TrainObject, envir = parent.frame())
  assign("ForeseeModel", value = tandem_fit, envir = parent.frame())
}
