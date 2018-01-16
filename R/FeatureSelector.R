#' Select the Genes that are used as Model Features
#'
#' The FeatureSelector selects features from all genes to be used for drug efficacy prediction.
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param DrugName Name of the drug whose efficacy is supposed to be predicted with the model
#' @param GeneFilter Set of genes to be considered for training the model, such as all, a certain percantage based on variance or p-value, specific gene sets like landmark genes, gene ontologies or pathways, etc.

#' @return \item{TrainObject}{The TrainObject with the selected gene set as features.}
#' @export

FeatureSelector <- function(TrainObject, DrugName, GeneFilter){
  return(TrainObject)
}
