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
#' The function 'rf' fits a random forest regression model by by Breiman (2001) to the training data

#' @param nfoldCrossvalidation # folds to use for crossvalidation while training the model. If put to zero, the complete data of the TrainObject is used for training.

#' @return \item{ForeseeModel}{A black box model trained on the TrainObject data that can be applied to new test data.}
#'         \item{TrainObject}{The TrainObject that was used to train the model.}
#' @export

# ToDo: nfoldCrossvalidation!
# ToDo: Check why objects still contain duplicates!
# ToDo: User-defined function
# ToDo: Tandem


BlackBoxFilter <- function(TrainObject, BlackBox, nfoldCrossvalidation){
  UseMethod("BlackBoxFilter", object = BlackBox)
  # Update Objects in the Environment
  TrainObject <<- TrainObject
  ForeseeModel <<- ForeseeModel

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

  TrainObject_train<- as.matrix(cbind(t(TrainObject$GeneExpression),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TrainObject_train <- as.data.frame(TrainObject_train)
  TrainObject_train <- TrainObject_train[,!duplicated(colnames(TrainObject_train))]
  TrainObject_train <- TrainObject_train[!duplicated(rownames(TrainObject_train)),]
  lm_fit <- lm(formula = DrugResponse~., TrainObject_train)

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <<- TrainObject_train
  ForeseeModel <<- lm_fit
}


################################################################################
### Function "ridge" to train a linear ridge regression model
BlackBoxFilter.ridge <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$GeneExpression),TrainObject$DrugResponse))
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
  TrainObject[["TrainFrame"]] <<- TrainObject_train
  ForeseeModel <<- ridge_fit
}


################################################################################
### Function "lasso" to train a lasso regression model
BlackBoxFilter.lasso <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$GeneExpression),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # Package glmnet by Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, https://web.stanford.edu/~hastie/Papers/glmnet.pdf
  require(glmnet)
  lasso_fit <- glmnet(x = t(TrainObject$GeneExpression), y=TrainObject$DrugResponse, alpha = 1)

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <<- TrainObject_train
  ForeseeModel <<- lasso_fit
}


################################################################################
### Function "elasticnet" to train a elasticnet regression model
BlackBoxFilter.elasticnet <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$GeneExpression),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TrainObject_train <- as.data.frame(TrainObject_train)
  TrainObject_train <- TrainObject_train[,!duplicated(colnames(TrainObject_train))]
  TrainObject_train <- TrainObject_train[!duplicated(rownames(TrainObject_train)),]

  # Package glmnet by Friedman, J., Hastie, T. and Tibshirani, R. (2008) Regularization Paths for Generalized Linear Models via Coordinate Descent, https://web.stanford.edu/~hastie/Papers/glmnet.pdf
  require(glmnet)
  elasticnet_fit <- glmnet(x = t(TrainObject$GeneExpression), y=TrainObject$DrugResponse, alpha = 0.5)

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <<- TrainObject_train
  ForeseeModel <<- elasticnet_fit
}


################################################################################
### Function "svm" to train a support vector regression model
BlackBoxFilter.svm <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$GeneExpression),TrainObject$DrugResponse))
  colnames(TrainObject_train)[dim(TrainObject_train)[2]]<-"DrugResponse"

  # For some weird reason the object still contains duplicates? Check duplication handler
  # Just take the first occuring gene name (here: in columns!) for now
  TrainObject_train <- as.data.frame(TrainObject_train)
  TrainObject_train <- TrainObject_train[,!duplicated(colnames(TrainObject_train))]
  TrainObject_train <- TrainObject_train[!duplicated(rownames(TrainObject_train)),]

  require(e1071)
  svm_fit <- svm(formula = DrugResponse~., data=data.frame(TrainObject_train))

  # Update Objects in the Environment
  TrainObject[["TrainFrame"]] <<- TrainObject_train
  ForeseeModel <<- svm_fit
}


################################################################################
### Function "rf" to train a random forest regression model
BlackBoxFilter.rf <- function(TrainObject, BlackBox, nfoldCrossvalidation){

  TrainObject_train<- as.matrix(cbind(t(TrainObject$GeneExpression),TrainObject$DrugResponse))
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
  TrainObject[["TrainFrame"]] <<- TrainObject_train
  ForeseeModel <<- rf_fit
}


