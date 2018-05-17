#' SampleSelector
#'
#' The Sample Selector returns features of only those samples that have a specified drug response and are available for all chosen input data types. If a specific tissue is selected, the samples are reduced to cell lines of that tissue only.
#'
#' @param TrainObject Object that contains all data needed to train a model, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type) and drug response data
#' @param TrainingTissue Tissue type that the cell lines of the TrainObject should be of, such as lung. Default should be "all" for pancancer analysis.
#' @param InputDataTypes Data types of the TrainObject that are to be used to train the model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @return \item{TrainObject}{The TrainObject with only those samples that have a specified drug response and are available for all chosen input data types.}
#' @examples
#' SampleSelector(GDSC, "pancreas", "GeneExpression")
#' SampleSelector(GDSC, "all", c("GeneExpression","Mutation"))
#' @export

#########################
# This file is part of the FORESEE R-package
# File authors: Lisa-Katrin Turnhoff <turnhoff@combine.rwth-aachen.de> and Ali Hadizadeh Esfahani <hadizadeh@combine.rwth-aachen.de>
# Distributed under the GNU General Public License v3.0.(http://www.gnu.org/licenses/gpl-3.0.html)
#########################

SampleSelector<- function(TrainObject,TrainingTissue, InputDataTypes){

  #################################################################################################################################
  # Check user's choice: Use all cell lines for training or only those of a specific tissue
  #################################################################################################################################
  # Check if the user wants to train on all cell lines (pan-cancer analysis)
  if (TrainingTissue=="all"){

    if (("GeneExpression" %in% InputDataTypes) == TRUE){
      CommonSamples <- colnames(TrainObject$GeneExpression)[colnames(TrainObject$GeneExpression) %in% names(TrainObject$DrugResponse)]
      TrainObject$GeneExpression <- TrainObject$GeneExpression[,CommonSamples]
      TrainObject$DrugResponse <- TrainObject$DrugResponse[CommonSamples]
    }
    if (("Mutation" %in% InputDataTypes) == TRUE){
      CommonSamples <- colnames(TrainObject$Mutation)[colnames(TrainObject$Mutation) %in% names(TrainObject$DrugResponse)]
      TrainObject$Mutation <- TrainObject$Mutation[,CommonSamples]
      TrainObject$DrugResponse <- TrainObject$DrugResponse[CommonSamples]
    }
    if (("Methylation" %in% InputDataTypes) == TRUE){
      CommonSamples <- colnames(TrainObject$Methylation)[colnames(TrainObject$Methylation) %in% names(TrainObject$DrugResponse)]
      TrainObject$Methylation <- TrainObject$Methylation[,CommonSamples]
      TrainObject$DrugResponse <- TrainObject$DrugResponse[CommonSamples]
    }
    if (("CopyNumberVaration" %in% InputDataTypes) == TRUE){
      CommonSamples <- colnames(TrainObject$CNVGain)[colnames(TrainObject$CNVGain) %in% names(TrainObject$DrugResponse)]
      TrainObject$CNVGain <- TrainObject$CNVGain[,CommonSamples]
      CommonSamples <- colnames(TrainObject$CNVLoss)[colnames(TrainObject$CNVLoss) %in% names(TrainObject$DrugResponse)]
      TrainObject$CNVLoss <- TrainObject$CNVLoss[,CommonSamples]
      TrainObject$DrugResponse <- TrainObject$DrugResponse[CommonSamples]
    }
    if (("CancerType" %in% InputDataTypes) == TRUE){

    CommonSamples <- colnames(TrainObject$CancerType)[colnames(TrainObject$CancerType) %in% names(TrainObject$DrugResponse)]
    TrainObject$CancerType <- TrainObject$CancerType[,CommonSamples]
    TrainObject$DrugResponse <- TrainObject$DrugResponse[CommonSamples]
    }
  }

  # Check if the user wants to train on cell linesof a specific tissue
  else if ((TrainingTissue %in% TrainObject$TissueInfo$Site)==TRUE){

    CellsofTissue <-TrainObject$TissueInfo$Line[TrainObject$TissueInfo$Site==TrainingTissue]
    RelevantSamples <- CellsofTissue[CellsofTissue %in% names(TrainObject$DrugResponse)]

    if (("GeneExpression" %in% InputDataTypes) == TRUE){
      TissueSamples <- colnames(TrainObject$GeneExpression)[colnames(TrainObject$GeneExpression) %in% RelevantSamples]
      TrainObject$GeneExpression <- TrainObject$GeneExpression[,TissueSamples]
      TrainObject$DrugResponse <- TrainObject$DrugResponse[TissueSamples]
    }
    if (("Mutation" %in% InputDataTypes) == TRUE){
      TissueSamples <- colnames(TrainObject$Mutation)[colnames(TrainObject$Mutation) %in% RelevantSamples]
      TrainObject$Mutation <- TrainObject$Mutation[,TissueSamples]
      TrainObject$DrugResponse <- TrainObject$DrugResponse[TissueSamples]
    }
    if (("Methylation" %in% InputDataTypes) == TRUE){
      TissueSamples <- colnames(TrainObject$Methylation)[colnames(TrainObject$Methylation) %in% RelevantSamples]
      TrainObject$Methylation <- TrainObject$Methylation[,TissueSamples]
      TrainObject$DrugResponse <- TrainObject$DrugResponse[TissueSamples]
    }
    if (("CopyNumberVaration" %in% InputDataTypes) == TRUE){
      TissueSamples <- colnames(TrainObject$CNVGain)[colnames(TrainObject$CNVGain) %in% RelevantSamples]
      TrainObject$CNVGain <- TrainObject$CNVGain[,TissueSamples]
      TissueSamples <- colnames(TrainObject$CNVLoss)[colnames(TrainObject$CNVLoss) %in% RelevantSamples]
      TrainObject$CNVLoss <- TrainObject$CNVLoss[,TissueSamples]
      TrainObject$DrugResponse <- TrainObject$DrugResponse[TissueSamples]
    }
    if (("CancerType" %in% InputDataTypes) == TRUE){

      TissueSamples <- colnames(TrainObject$CancerType)[colnames(TrainObject$CancerType) %in% RelevantSamples]
      TrainObject$CancerType <- TrainObject$CancerType[,TissueSamples]
      TrainObject$DrugResponse <- TrainObject$DrugResponse[TissueSamples]
    }
  }

  # Stop if specified tissue is not part of the tissue types listed in the Foresee TrainObject
  else{
    stop(paste0("The specified tissue '",TrainingTissue,"' is not included as a tissue type in the ForeseeObject."))
  }


  #################################################################################################################################
  # Identify those samples that are characterized for all chosen input data types
  #################################################################################################################################

  if (length(InputDataTypes)>1){

    JointSamples<-c()

    for (i in 1:length(InputDataTypes)){

      FeatureType <- InputDataTypes[i]
      JointSamples[[paste0("Samples_Features_",i)]] <- as.list(colnames(TrainObject[[FeatureType]]))
    }

    JointSamples <- as.character(Reduce(intersect, JointSamples))

    for (i in 1:length(InputDataTypes)){
      FeatureType <- InputDataTypes[i]
      TrainObject[[FeatureType]] <- TrainObject[[FeatureType]][,JointSamples]
    }

    TrainObject$DrugResponse <- TrainObject$DrugResponse[JointSamples]

    if (!length(JointSamples)>1){ # Stop if samples for the selected input data types do not overlap
      stop(paste0("The subsets of samples for the specified input data types '",InputDataTypes,"' do not overlap."))
    }
  }
  # Returning the new object
  return(TrainObject)
}



