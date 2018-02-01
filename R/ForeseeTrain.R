#' Train a Drug Efficacy Prediction Model
#'
#' ForeseeTrain uses the data of the TrainObject to train a black box model that can later applied to new data in order to predict drug efficacy.
#' First, duplicates in the gene names are removed using the Foresee DuplicationHandler.
#' Second, the Homogenizer function removes batch effects and homogenizes train and test data.
#' Third, the FeatureSelector filters the input and selects features to be used for the prediction.
#' Fourth, the FeaturePreprocessor converts the original features into predictive features.
#' Fifth, the CellResponseProcessor prepares the response data of the TrainObject for prediction.
#' Last, a machine learning algorithm is applied to the data to create a model that is predictive of the drug response.
#'
#'
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param TestObject Object that contains all data that the model is to be tested on, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param DrugName Name of the drug whose efficacy is supposed to be predicted with the model
#' @param CellResponseType Format of the drug response data of the TrainObject, such as IC50, AUC, GI50, etc.
#' @param CellResponseTransformation Method that is to be used to transform the drug response data of the TrainObject, such as power transform, logarithm, binarization, user defined functions, etc.
#' @param InputDataTypes Data types of the TrainObject that are to be used to train the model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param DuplicationHandling Method for handling duplicates of gene names, such as taking none, the average, the first hit, etc.
#' @param HomogenizationMethod Method for homogenizing data of the TrainObject and TestObject, such as ComBat, quantile normalization, limma, RUV, etc.
#' @param GeneFilter Set of genes to be considered for training the model, such as all, a certain percantage based on variance or p-value, specific gene sets like landmark genes, gene ontologies or pathways, etc.
#' @param FeaturePreprocessing Method for preprocessing the inputs of the model, such as z-score, principal component analysis, PhysioSpace similarity, etc.
#' @param BlackBox Modeling algorithm for training, such as linear regression, elastic net, lasso regression, ridge regression, tandem, support vector machines, random forests, user defined functions, etc.
#' @param nfoldCrossvalidation # folds to use for crossvalidation while training the model. If put to zero, the complete data of the TrainObject is used for training.

#' @return \item{ForeseeModel}{A black box model trained on the TrainObject data that can be applied to new test data.}
#'         \item{TrainObject}{The TrainObject with preprocessed and filtered features.}
#'         \item{TestObject}{The TestObject with preprocessed and filtered features.}
#' @export

ForeseeTrain <- function(TrainObject, TestObject, DrugName, CellResponseType, CellResponseTransformation, InputDataTypes,
                         DuplicationHandling, HomogenizationMethod, GeneFilter, FeaturePreprocessing, BlackBox, nfoldCrossvalidation){
  #Handling duplicated genes:
  TrainObject <- DuplicationHandler(Object = TrainObject, DuplicationHandling = DuplicationHandling)
  TestObject <- DuplicationHandler(Object = TestObject, DuplicationHandling = DuplicationHandling)

  #Homogenizing Train and Test Objects:
  Homogenizer(TrainObject, TestObject, HomogenizationMethod)

  return(list("TrainObject"=TrainObject, "TestObject"=TestObject))
}
