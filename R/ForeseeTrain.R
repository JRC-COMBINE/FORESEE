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

  TrainObject <- TrainObject
  TestObject <- TestObject
  DrugName <- DrugName
  CellResponseType <- CellResponseType
  CellResponseTransformation <- CellResponseTransformation
  InputDataTypes <- InputDataTypes
  DuplicationHandling <- DuplicationHandling
  HomogenizationMethod <- HomogenizationMethod
  GeneFilter <- GeneFilter
  FeaturePreprocessing <- FeaturePreprocessing
  BlackBox <- BlackBox
  nfoldCrossvalidation <- nfoldCrossvalidation
  ForeseeModel <- 0

#################################################################################################################################
# 1. Cell Response Processing

  # Do CellResponseProcessor in the beginning to avoid doing the rest on duplicated cell line names etc.
  # the function 'powertransform' power transforms the drug response data,
  # the function 'logarithm' returns the natural logarithm of the drug response data,
  # the function 'binarization_kmeans' returns a binarized drug response vector based on 2 kmeans clusters,
  # the function 'binarization_cutoff' returns a binarized drug response vector based on a cutoff at the median,
  # the function 'none' returns the unchanged drug response data,
  # the function 'user-defined function' is determined by the function in the input

  # CellResponseType_options <- c("IC50", "AUC")
  # CellResponseTransformation_options <- c("powertransform", "logarithm", "binarization_kmeans", "binarization_cutoff", "none")

  # Process Cell Response
  CellResponseProcessor(TrainObject, DrugName, CellResponseType, CellResponseTransformation)

#################################################################################################################################
# 2. Remove Duplicates in gene names
# Options
  # The function 'mean' calculates the mean of all rows that have the same gene name,
  # The function 'first' chooses the first occuring row of duplicated genes only,
  # The function 'none' removes all genes that occur more than once.

  # DuplicationHandling_options <- c("first", "mean", "none")

# Train matrix
  TrainObject <- DuplicationHandler(TrainObject, DuplicationHandling)
# Test matrix
  TestObject <- DuplicationHandler(TestObject, DuplicationHandling)

#################################################################################################################################
# 3. Homogenize training and test data set
  # Options
  # The function 'ComBat' uses the batch effect removal ComBat of the sva package,
  # The function 'quantile' uses the quantile normalization of the preprocessCore package,
  # The function 'limma' uses the removeBatchEffect function of the limma package,
  # The function 'RUV' uses the RUV normalization, requiring a set of negative control genes,
  # The function 'none' doesn't do any batch effect correction,
  # If the user wants to implement a user-defined function batch effect removal function, the input should be the function.

  # HomogenizationMethod_options <- c("ComBat", "quantile", "limma", "RUV", "none")

  # Homogenize
  Homogenizer(TrainObject, TestObject, HomogenizationMethod)

#################################################################################################################################
# 4. Feature Selection

  # Options
  # The option 'variance' removes the 20 % genes of lowest variance across samples in the TrainObject
  # The option 'pvalue' removes the 20 % genes of lowest p-value (ttest) across samples in the TrainObject
  # The option 'landmarkgenes' uses the L1000 gene set downloaded from CLUE Command App
  # The option 'ontology' uses a specific set of genes included in the chosen ontology? -> Ali
  # The option 'pathway' uses a specific set of genes included in the chosen pathway? -> Ali
  # The option 'all' keeps all genes as features
  # If the user inserts a list as an input, the list is considered as chosen features.

  # FeatureSelector_options <- c("variance", "pvalue", "landmarkgenes", "ontology", "pathway", "all")

  FeatureSelector(TrainObject,TestObject, GeneFilter)

#################################################################################################################################
# 5. Feature Preprocessing

  # Options
  # The function 'zscore_genewise' calculates the zscore normalizing each gene over all samples,
  # The function 'zscore_samplewise' calculates the zscore normalizing each sample over all genes,
  # The function 'pca' does principal component analysis,
  # The function 'physio' does physiospace analysis with the samples using cell line gene expression of the gdsc data base as physiological references,
  # The function 'none' keeps the gene expression values unchanged,
  # If the user wants to implement a user-defined function batch effect removal function, the input should be the function.

  # FeaturePreprocessing_options <- c("zscore_genewise", "zscore_samplewise", "pca", "physio", "none")

  FeaturePreprocessor(TrainObject, TestObject, FeaturePreprocessing)


#################################################################################################################################
# 6. Black Box Filter

  # Options
  # The function 'linear' fits a linear regression model to the training data,
  # The function 'ridge' fits a linear ridge regression model by Cule et al. (2012) to the training data,
  # The function 'lasso' fits a lasso regression model from the glmnet package by Friedman et al. (2008) to the training data,
  # The function 'elasticnet' fits an elastic net regression model from the glmnet package by Friedman et al. (2008) to the training data,
  # The function 'svm' fits a support vector regression model ffrom the e1071 package by Meyer and Chih-Chung (2017) to the training data,
  # The function 'rf' fits a random forest regression model by by Breiman (2001) to the training data

  # BlackBox_options <- c("linear", "ridge", "lasso", "elasticnet", "svm", "rf")

  BlackBoxFilter(TrainObject, BlackBox, nfoldCrossvalidation)

  ForeseeModel <<- ForeseeModel
  TrainObject <<- TrainObject
  TestObject <<- TestObject

}
