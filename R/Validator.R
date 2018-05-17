#' Validate the Performance of a ForeseeModel on a new TestObject
#'
#' The Validator evaluates the performance of the model by comparing the predicted response Foreseen to the reported, true response from the TestObject's Annotation.
#'
#' @param TestObject Object that contains all data that the model is to be tested on, especially the true, measured drug response.
#' @param Foreseen Predicted drug response of the TestObject obtained by applying the ForeseeModel.
#' @param Evaluation Measure for evaluating the model performance, such as ROC-Curve, AUC or p-value of ROC-Curve, Rsquared, MSE, Correlation, F-Test, etc.
#' The function 'listInputOptions("Validator")' returns a list of the possible options.
#' Instead of choosing one of the implemented options, a user-defined function can be used as an input.
#' @return \item{Performance}{Evaluation Measure of the Predictability of the ForeseeModel trained on the TrainObject and tested on the TestObject.}
#' @examples
#' Validator(rep(1,length(GSE6434$Annotation)),GSE6434,"rocauc")
#' @export

#########################
# This file is part of the FORESEE R-package
# File authors: Lisa-Katrin Turnhoff <turnhoff@combine.rwth-aachen.de> and Ali Hadizadeh Esfahani <hadizadeh@combine.rwth-aachen.de>
# Distributed under the GNU General Public License v3.0.(http://www.gnu.org/licenses/gpl-3.0.html)
#########################

Validator <- function(Foreseen, TestObject, Evaluation){
  UseMethod("Validator", object = Evaluation)
}

#' @export
Validator.character <- function(Foreseen, TestObject, Evaluation){
  class(Evaluation) <- Evaluation;
  UseMethod("Validator", object = Evaluation)
  return(Performance)
}

################################################################################
### Function "function" applies the function in "Evaluation" to predicted values
#' @export
Validator.function <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  return(Evaluation(Foreseen, ANNOTATIONS))
}

################################################################################
### Function "rocauc" calculates the AUC of the ROC curve (binarizes TestObj annotation on median if they have continues values)
#' @export
Validator.rocauc <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  if(is.numeric(ANNOTATIONS)){
    message("Annotation of the test set is binarized for calculating ROC")
    ANNOTATIONS <- ifelse(ANNOTATIONS-median(ANNOTATIONS) > 0, TRUE, FALSE)
  }
  requireForesee(pROC) ## glmnet masks auc so had to use functions of pROC as pROC::function (pROC::roc and pROC::auc)
  ROCObj <- pROC::roc(-(as.numeric(ANNOTATIONS)), Foreseen, direction="<")
  AUCofROC <- pROC::auc(pROC::roc(ANNOTATIONS, Foreseen))[[1]]
  return(AUCofROC)

  ### Plot ROC Curve
  #jpeg(filename = filename_roc, width=10, height=10, units="in", res=600)
  #plot(ROCObj, main = "Prediction of Patient Response after Training on Cell Line Data",legacy.axes = TRUE, print.auc = TRUE)
  #dev.off()


}

################################################################################
### Function "rocpvalue" calculates the t.test p value of the ROC curve AUC (binarizes TestObj annotation on median if they have continues values)
## versus the 10000 permutated annotation values
#' @export
Validator.rocpvalue <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  if(is.numeric(ANNOTATIONS)){
    message("Annotation of the test set is binarized for calculating ROC")
    ANNOTATIONS <- ifelse(ANNOTATIONS-median(ANNOTATIONS) > 0, TRUE, FALSE)
  }
  requireForesee(pROC) ## glmnet masks auc so had to use functions of pROC as pROC::function (pROC::roc and pROC::auc)
  AUCofROC <- pROC::auc(pROC::roc(ANNOTATIONS, Foreseen))[[1]]
  AUCofROCRandom <- numeric(length = 10000)
  for(K in 1:10000){
    set.seed(0)
    AUCofROCRandom <- pROC::auc(pROC::roc(sample(ANNOTATIONS), Foreseen))[[1]]
  }
  return(sum(AUCofROC<=AUCofROCRandom)/length(AUCofROCRandom))
}


################################################################################
### Function "prauc" calculates the AUC of the precision-recall curve
#' @export
Validator.prauc <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  if(is.numeric(ANNOTATIONS)){
    message("Annotation of the test set is binarized for calculating ROC")
    ANNOTATIONS <- ifelse(ANNOTATIONS-median(ANNOTATIONS) > 0, TRUE, FALSE)
  }

  requireForesee(PRROC)
  PRObj <- PRROC::pr.curve(-(as.numeric(ANNOTATIONS)), Foreseen)
  AUCofPR <- PRObj$auc.integral
  return(AUCofPR)
}

################################################################################
### Function "rsquared" calculates the fraction of variance explained by a linear model between predictions and
## actual annotations, !! only for when there is continuous numeric (not binary) annotation availble on TestObj!!
#' @export
Validator.rsquared <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  if(is.logical(ANNOTATIONS)){
    warning(paste("Annotation of the test set is binary! Is", Evaluation,"the correct validation method?"))
    ANNOTATIONS <- as.numeric(ANNOTATIONS)
  }
  return(summary(lm(ANNOTATIONS~Foreseen))$r.squared)
}

################################################################################
### Function "rsquared_adjusted" calculates the fraction of variance explained by a linear model between predictions and
## actual annotations, corrected by the p-value of F-test,
##!! only for when there is continuous numeric (not binary) annotation availble on TestObj!!
#' @export
Validator.rsquared_adjusted <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  if(is.logical(ANNOTATIONS)){
    warning(paste("Annotation of the test set is binary! Is", Evaluation,"the correct validation method?"))
    ANNOTATIONS <- as.numeric(ANNOTATIONS)
  }
  return(summary(lm(ANNOTATIONS~Foreseen))$adj.r.squared)
}

################################################################################
### Function "fpvalue" calculates the p value of an F test on a linear model between predictions and
## actual annotations, !! only for when there is continuous numeric (not binary) annotation availble on TestObj!!
#' @export
Validator.fpvalue <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  if(is.logical(ANNOTATIONS)){
    warning(paste("Annotation of the test set is binary! Is", Evaluation,"the correct validation method?"))
    ANNOTATIONS <- as.numeric(ANNOTATIONS)
  }
  FStats <- summary(lm(ANNOTATIONS~Foreseen))$fstatistic

  if (is.null(FStats)==TRUE){
  return(NA)
  }
  else{
  return(pf(FStats[1],FStats[2],FStats[3],lower.tail=F))
  }
}

################################################################################
### Function "mse" calculates the mean square error of a linear model between predictions and
## actual annotations, !! only for when there is continuous numeric (not binary) annotation availble on TestObj!!
#' @export
Validator.mse <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  if(is.logical(ANNOTATIONS)){
    warning(paste("Annotation of the test set is binary! Is", Evaluation,"the correct validation method?"))
    ANNOTATIONS <- as.numeric(ANNOTATIONS)
  }
  WholeSummary <- summary(lm(ANNOTATIONS~Foreseen))
  return(mean(WholeSummary$residuals^2))
}

################################################################################
### Function "spearman" calculates the spearman correlation between predictions and
## actual annotations, !! only for when there is continuous numeric (not binary) annotation availble on TestObj!!
#' @export
Validator.spearman <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  if(is.logical(ANNOTATIONS)){
    warning(paste("Annotation of the test set is binary! Is", Evaluation,"the correct validation method?"))
    ANNOTATIONS <- as.numeric(ANNOTATIONS)
  }
  return(cor(ANNOTATIONS,Foreseen,method = Evaluation))
}

################################################################################
### Function "pearson" calculates the spearman correlation between predictions and
## actual annotations, !! only for when there is continuous numeric (not binary) annotation availble on TestObj!!
#' @export
Validator.pearson <- function(Foreseen, TestObject, Evaluation) {
  ANNOTATIONS <- if(CellorPatient(TestObject)) TestObject$DrugResponse else TestObject$Annotation
  if(is.logical(ANNOTATIONS)){
    warning(paste("Annotation of the test set is binary! Is", Evaluation,"the correct validation method?"))
    ANNOTATIONS <- as.numeric(ANNOTATIONS)
  }
  return(cor(ANNOTATIONS,Foreseen,method = Evaluation))
}

################################################################################
### Function "default" is called in case method in "Evaluation" is
# unknown to Validator
#' @export
Validator.default <- function(Foreseen, TestObject, Evaluation){
  stop(paste("Method",Evaluation,"is not defined as an evaluation method!"))
}
