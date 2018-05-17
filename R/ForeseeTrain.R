#' Train a Drug Efficacy Prediction Model
#'
#' ForeseeTrain uses the data of the TrainObject to train a black box model that can later be applied to new data in order to predict drug efficacy.
#' Duplicates in the gene names are removed using the Foresee DuplicationHandler.
#' The Homogenizer function removes batch effects and homogenizes train and test data.
#' The FeatureSelector filters the input and selects features to be used for the prediction.
#' The FeaturePreprocessor converts the original features into predictive features.
#' The CellResponseProcessor prepares the response data of the TrainObject for prediction.
#' A machine learning algorithm is applied to the preprocessed data to create a model that is predictive of drug response.
#'
#'
#'
#' @param TrainObject Object that contains all data needed to train a model, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type) and drug response data
#' @param TestObject Object that contains all data that the model is to be tested on, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type) and drug response data
#' @param DrugName Name of the drug whose efficacy is supposed to be predicted with the model.
#' You can get all possible values with listDrugs(OBJ) or listInputOptions("DrugName", OBJ), where OBJ is the object you want to use as TrainObject.
#' @param CellResponseType Format of the drug response data of the TrainObject, such as IC50, AUC, GI50, etc.
#' You can get all possible values with listInputOptions("CellResponseType", OBJ), where OBJ is the object you want to use as TrainObject.
#' @param CellResponseTransformation Method that is to be used to transform the drug response data of the TrainObject, such as power transform, logarithm, binarization, user defined functions, etc.
#' Get all possible values with listInputOptions("CellResponseProcessor").
#' @param InputDataTypes Data types of the TrainObject that are to be used to train the model, such as GeneExpression, Mutation, CopyNumberVariation, Methylation, Cancertype, etc.
#' You can get all possible values with listInputOptions("InputDataTypes", OBJ), where OBJ is the object you want to use as TrainObject.
#' @param DuplicationHandling Method for handling duplicates of gene names, such as considering none, the mean, the first hit, etc.
#' Get all possible values with listInputOptions("DuplicationHandler").
#' @param HomogenizationMethod Method for homogenizing data of the TrainObject and TestObject, such as ComBat, quantile normalization, limma, RUV, etc.
#' Get all possible values with listInputOptions("Homogenizer").
#' @param TrainingTissue Tissue type that the cell lines of the TrainObject should be of, such as pancreas or lung. Default is "all" for pancancer analysis.
#' You can get all possible values with listInputOptions("TrainingTissue", OBJ), where OBJ is the object you want to use as TrainObject.
#' @param TestingTissue Tissue type that the cell lines or samples of the TestObject should be of, such as pancreas or lung. Default is "all" for analysis of all samples.
#' You can get all possible values with listInputOptions("TestingTissue", OBJ), where OBJ is the object you want to use as TestObject.
#' @param GeneFilter Set of genes to be considered for training the model, such as all, a certain percantage based on variance or p-value, specific gene sets like landmark genes, gene ontologies or pathways, etc.
#' Get all possible values with listInputOptions("FeatureSelector").
#' @param FeaturePreprocessing Method for preprocessing the inputs of the model, such as z-score, principal component analysis, PhysioSpace similarity, etc.
#' Get all possible values with listInputOptions("FeaturePreprocessor").
#' @param BlackBox Modeling algorithm for training, such as linear regression, elastic net, lasso regression, ridge regression, tandem, support vector machines, random forests, user defined functions, etc.
#' Get all possible values with listInputOptions("BlackBoxFilter").
#' @param nfoldCrossvalidation # folds to use for crossvalidation while training the model. If put to one, the complete data of the TrainObject is used for training.

#' @return \item{ForeseeModel}{A black box model trained on the TrainObject data that can be applied to new test data.}
#'         \item{TrainObject}{The TrainObject with preprocessed and filtered features.}
#'         \item{TestObject}{The TestObject with preprocessed and filtered features.}
#' @export

#########################
# This file is part of the FORESEE R-package
# File authors: Lisa-Katrin Turnhoff <turnhoff@combine.rwth-aachen.de> and Ali Hadizadeh Esfahani <hadizadeh@combine.rwth-aachen.de>
# Distributed under the GNU General Public License v3.0.(http://www.gnu.org/licenses/gpl-3.0.html)
#########################

ForeseeTrain <- function(TrainObject, TestObject, DrugName, CellResponseType, CellResponseTransformation = "powertransform",
                         InputDataTypes = "GeneExpression", TrainingTissue = "all", TestingTissue = "all",
                         DuplicationHandling = "first", HomogenizationMethod = "ComBat",
                         GeneFilter = "all", FeaturePreprocessing = "none", BlackBox = "ridge", nfoldCrossvalidation = 1,...){


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
  TrainObject <- CellResponseProcessor(TrainObject, DrugName, CellResponseType, CellResponseTransformation)


  # If FORESEE is used to do cell2cell prediction, the drug response data of the TestObject is processed in the same manner as the TrainObject
  if (class(TestObject)=="ForeseeCell"){
    TestObject <- CellResponseProcessor(TestObject, DrugName=DrugName, CellResponseType=CellResponseType, CellResponseTransformation=CellResponseTransformation)
  }

  #################################################################################################################################
  # 2. Selecting Samples

  # SampleSelector_options <- c("all",as.character(unique(TrainObject$TissueInfo$Site)))

  TrainObject <- SampleSelector(TrainObject=TrainObject,TrainingTissue=TrainingTissue, InputDataTypes=InputDataTypes)

  # If FORESEE is used to do cell2cell prediction, relevant samples of the TestObject is extracted in the same manner as the TrainObject
  if (class(TestObject)=="ForeseeCell"){
    TestObject <- SampleSelector(TrainObject=TestObject,TrainingTissue=TestingTissue, InputDataTypes=InputDataTypes)
  }

  #################################################################################################################################
  # 3. Preprocess Gene Expression Data

  if (("GeneExpression" %in% InputDataTypes)==TRUE){


    #################################################################################################################################
    # 3a. Remove Duplicates in gene names
    # Options
      # The function 'mean' calculates the mean of all rows that have the same gene name,
      # The function 'first' chooses the first occuring row of duplicated genes only,
      # The function 'none' removes all genes that occur more than once.

      # DuplicationHandling_options <- c("first", "mean", "none")

    # Train matrix
      TrainObject <- DuplicationHandler(Object=TrainObject, DuplicationHandling=DuplicationHandling)
    # Test matrix
      TestObject <- DuplicationHandler(Object=TestObject, DuplicationHandling=DuplicationHandling)

    #################################################################################################################################
    # 3b. Homogenize training and test data set
      # Options
      # The function 'ComBat' uses the batch effect removal ComBat of the sva package,
      # The function 'quantile' uses the quantile normalization of the preprocessCore package,
      # The function 'limma' uses the removeBatchEffect function of the limma package,
      # The function 'RUV' uses the RUV normalization, requiring a set of negative control genes,
      # The function 'none' doesn't do any batch effect correction,
      # If the user wants to implement a user-defined function batch effect removal function, the input should be the function.

      #   HomogenizationMethod_options <- c("ComBat", "quantile", "limma", "YuGene", "RUV", "RUV4", "none")

      # Homogenize
      Homogenizer(TrainObject=TrainObject, TestObject=TestObject, HomogenizationMethod=HomogenizationMethod)

    #################################################################################################################################
    # 3c. Feature Selection

      # Options
      # The option 'variance' removes the 20 % genes of lowest variance across samples in the TrainObject
      # The option 'pvalue' removes the 20 % genes of lowest p-value (ttest) across samples in the TrainObject
      # The option 'landmarkgenes' uses the L1000 gene set downloaded from CLUE Command App
      # The option 'ontology' uses a specific set of genes included in the chosen ontology? -> Ali
      # The option 'pathway' uses a specific set of genes included in the chosen pathway? -> Ali
      # The option 'all' keeps all genes as features
      # If the user inserts a list as an input, the list is considered as chosen features.

      # FeatureSelector_options <- c("variance", "pvalue", "landmarkgenes", "ontology", "pathway", "all")

      FeatureSelector(TrainObject=TrainObject,TestObject=TestObject, GeneFilter=GeneFilter, DrugName=DrugName)

    #################################################################################################################################
    # 3d. Feature Preprocessing

      # Options
      # The function 'zscore_genewise' calculates the zscore normalizing each gene over all samples,
      # The function 'zscore_samplewise' calculates the zscore normalizing each sample over all genes,
      # The function 'pca' does principal component analysis,
      # The function 'physio' does physiospace analysis with the samples using cell line gene expression of the gdsc data base as physiological references,
      # The function 'none' keeps the gene expression values unchanged,
      # If the user wants to implement a user-defined function batch effect removal function, the input should be the function.

      # FeaturePreprocessing_options <- c("zscore_genewise", "zscore_samplewise", "pca", "physio", "none")

      FeaturePreprocessor(TrainObject=TrainObject, TestObject=TestObject, FeaturePreprocessing=FeaturePreprocessing)


  }


  #################################################################################################################################
  # 4. Creating the Feature Matrix

  FeatureCombiner(TrainObject=TrainObject, TestObject=TestObject, InputDataTypes=InputDataTypes)

#################################################################################################################################
# 5. Black Box Filter

  # Options
  # The function 'linear' fits a linear regression model to the training data,
  # The function 'ridge' fits a linear ridge regression model by Cule et al. (2012) to the training data,
  # The function 'lasso' fits a lasso regression model from the glmnet package by Friedman et al. (2008) to the training data,
  # The function 'elasticnet' fits an elastic net regression model from the glmnet package by Friedman et al. (2008) to the training data,
  # The function 'svm' fits a support vector regression model ffrom the e1071 package by Meyer and Chih-Chung (2017) to the training data,
  # The function 'rf' fits a random forest regression model by by Breiman (2001) to the training data
  # The function 'rf_ranger' fits a fast random forest regression model by Marvin N. Wright (2018) to the training data

  # BlackBox_options <- c("linear", "ridge", "lasso", "elasticnet", "svm", "rf", "rf_ranger")

  BlackBoxFilter(TrainObject=TrainObject, BlackBox=BlackBox, nfoldCrossvalidation=nfoldCrossvalidation)


  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject, envir = parent.frame())
  assign("TestObject", value = TestObject, envir = parent.frame())
  assign("ForeseeModel", value = ForeseeModel, envir = parent.frame())
}
