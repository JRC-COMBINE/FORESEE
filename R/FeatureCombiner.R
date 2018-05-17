#' FeatureCombiner
#'
#' The Feature Combiner combines all selected Input Data Types into one Feature Matrix
#'
#' @param TrainObject Object that contains all data needed to train a model, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type, etc. ) and drug response data
#' @param TestObject Object that contains all data that the model is to be tested on, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type, etc. ) and drug response data
#' @param InputDataTypes Data types of the TrainObject that are to be used to train the model, such as "GeneExpression", "Mutation", "CopyNumberVariation", "Methylation", "CancerType", etc.
#' @return \item{TrainObject}{The TrainObject with a new Feature matrix combining all specified input data types and a Featuretype Vector indicating the molecular data type of each feature}
#'         \item{TestObject}{The TestObject with a new Feature matrix combining all specified input data types and a Featuretype Vector indicating the molecular data type of each feature}
#' @examples
#' FeatureCombiner(GDSC,GSE6434,c("GeneExpression", "Mutation"))
#' @export

#########################
# This file is part of the FORESEE R-package
# File authors: Lisa-Katrin Turnhoff <turnhoff@combine.rwth-aachen.de> and Ali Hadizadeh Esfahani <hadizadeh@combine.rwth-aachen.de>
# Distributed under the GNU General Public License v3.0.(http://www.gnu.org/licenses/gpl-3.0.html)
#########################

FeatureCombiner<- function(TrainObject, TestObject, InputDataTypes){

  ################################################################################
  # Check if both objects include the specified InputDataTypes

  missing_data_train <- c()
  missing_data_index_train <- c()
  missing_data_test <- c()
  missing_data_index_test <- c()

  for (m in 1:length(InputDataTypes)){

    if(is.null(TrainObject[[InputDataTypes[m]]])==TRUE){
      warning(paste0("The ForeseeTrainObject does not comprise data of the type '",InputDataTypes[m],"'. Therefore, the training will proceed without that data type."))
      missing_data_index_train <- m

    }

    if(is.null(TestObject[[InputDataTypes[m]]])==TRUE){
    warning(paste0("The ForeseeTest Object does not comprise data of the type '",InputDataTypes[m],"'. Therefore, the training will proceed without that data type."))
    missing_data_index_test <- m
      }

    missing_data_train <- cbind(missing_data_train,missing_data_index_train)
    missing_data_test <- cbind(missing_data_test,missing_data_index_test)
  }

  missing_data<-Reduce(union,c(missing_data_train,missing_data_test))
  if(length(missing_data)>0) InputDataTypes<-InputDataTypes[-missing_data]

  ################################################################################
  # Reduce data types to the features that are common in TrainObject and TestObjecct

  JointExtractedFeatures_Train<-c()
  JointExtractedFeatures_Test<-c()

  for (m in 1:length(InputDataTypes)){

    FeatureType <- InputDataTypes[m]
    ExtractedFeatures_Train <- rownames(TrainObject[[FeatureType]])
    ExtractedFeatures_Test <- rownames(TestObject[[FeatureType]])
    JointExtractedFeatures <- Reduce(intersect, list(ExtractedFeatures_Train,ExtractedFeatures_Test))

    TrainObject[[FeatureType]]<-TrainObject[[FeatureType]][JointExtractedFeatures,]
    TestObject[[FeatureType]]<-TestObject[[FeatureType]][JointExtractedFeatures,]
  }


  ################################################################################
  # Create Feature Matrices
  Features_Train <- c()
  FeatureTypes_Train  <- c()

  for (i in 1:length(InputDataTypes)){
    FeatureType <- InputDataTypes[i]
    Features_i <- TrainObject[[FeatureType]]
    FeatureTypes_Train  <- cbind(FeatureTypes_Train, rbind(rep(x=FeatureType,times=dim(Features_i)[1])))
    Features_Train  <- rbind(Features_Train,Features_i)
  }

  Features_Test <- c()
  FeatureTypes_Test  <- c()

  for (i in 1:length(InputDataTypes)){
    FeatureType <- InputDataTypes[i]
    Features_i <- TestObject[[FeatureType]]
    FeatureTypes_Test  <- cbind(FeatureTypes_Test, rbind(rep(x=FeatureType,times=dim(Features_i)[1])))
    Features_Test  <- rbind(Features_Test,Features_i)
  }

################################################################################
# Create new feature matrices and feature type vectors in Foresee objects
TrainObject[["Features"]] <- Features_Train
TrainObject[["FeatureTypes"]] <- FeatureTypes_Train
TestObject[["Features"]] <- Features_Test
TestObject[["FeatureTypes"]] <- FeatureTypes_Test

# Update Objects in the Environment
assign("TrainObject", value = TrainObject, envir = parent.frame())
assign("TestObject", value = TestObject, envir = parent.frame())
}
