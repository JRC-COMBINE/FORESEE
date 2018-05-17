#' Remove Duplicated Gene Names from a FORESEE Object
#'
#' DuplicationHandler finds duplicates in the gene names (features) from the FORESEE Object and summarizes or deletes them according to the user's preferences.
#' @param Object FORESEE Object (ForeseeCell or ForeseeTrain) that contains all data needed to train a model, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type, etc. ) and drug response data
#' @param DuplicationHandling Method for handling duplicates of gene names.
#' The function 'mean' calculates the mean of all rows that have the same gene name,
#' The function 'first' chooses the first hit of duplicated genes and discards the rest of genes with the same name,
#' The function 'none' removes all gene names that occur more than once.
#' The function 'listInputOptions("DuplicationHandler")' returns a list of the possible options.
#' Instead of chosing one of the implemented options, a user-defined function can be used as an input.
#' @return \item{Object}{The object with unique gene names}
#' @examples
#' DuplicationHandler(GDSC,"first")
#' @export

#########################
# This file is part of the FORESEE R-package
# File authors: Lisa-Katrin Turnhoff <turnhoff@combine.rwth-aachen.de> and Ali Hadizadeh Esfahani <hadizadeh@combine.rwth-aachen.de>
# Distributed under the GNU General Public License v3.0.(http://www.gnu.org/licenses/gpl-3.0.html)
#########################

DuplicationHandler <- function(Object, DuplicationHandling){
  UseMethod("DuplicationHandler", object = DuplicationHandling)
}

#' @export
DuplicationHandler.character <- function(Object, DuplicationHandling){
  class(DuplicationHandling) <- DuplicationHandling;
  UseMethod("DuplicationHandler", object = DuplicationHandling)
}

################################################################################
### Function "function" applies the function in "DuplicationHandling"
# to all rows that have the same gene name
#' @export
DuplicationHandler.function <- function(Object, DuplicationHandling){

  Object_withoutDuplicates <- Object

  # Counts the number of duplicated gene names in the gene expression matrix
  dim_before <- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Removes duplicates
  Object_withoutDuplicates$GeneExpression <- data.matrix(aggregate(x = Object_withoutDuplicates$GeneExpression, by = list(rownames(Object_withoutDuplicates$GeneExpression)), FUN = DuplicationHandling))
  rownames(Object_withoutDuplicates$GeneExpression) <- Object_withoutDuplicates$GeneExpression[,1]
  Object_withoutDuplicates$GeneExpression <- Object_withoutDuplicates$GeneExpression[,2:dim(Object_withoutDuplicates$GeneExpression)[2]]

  # Calculates the dimension of the gene expression matrix after removing duplicates
  dim_after<- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The removal of duplicates reduced the number of genes in the Foresee Object from ", dim_before, " to ", dim_after))

  # Update Object in the Environment
  Object <- Object_withoutDuplicates
  return(Object_withoutDuplicates)

}

################################################################################
### Function "first" to choose the first occuring row of duplicated genes only
#' @export
DuplicationHandler.first <- function(Object, DuplicationHandling){

  Object_withoutDuplicates <- Object

  # Counts the number of gene names in the gene expression matrix
  dim_before <- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Removes duplicates
  Object_withoutDuplicates$GeneExpression<-Object_withoutDuplicates$GeneExpression[!duplicated(rownames(Object_withoutDuplicates$GeneExpression)),]

  # Calculates the dimension of the gene expression matrix after removing duplicates
  dim_after <- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The removal of duplicates reduced the number of genes in the Foresee Object from ", dim_before, " to ", dim_after))

  # Update Object in the Environment
  Object <- Object_withoutDuplicates
  return(Object)

}

################################################################################
### Function none" removes all genes that occur more than once.
#' @export
DuplicationHandler.none <- function(Object, DuplicationHandling){

  Object_withoutDuplicates <- Object

  # Counts the number of duplicated gene names in the gene expression matrix
  dim_before <- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Removes duplicates
  Object_withoutDuplicates$GeneExpression<-Object_withoutDuplicates$GeneExpression[!(duplicated(rownames(Object_withoutDuplicates$GeneExpression)) | duplicated(rownames(Object_withoutDuplicates$GeneExpression), fromLast = TRUE)),]

  # Calculates the dimension of the gene expression matrix after removing duplicates
  dim_after<- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The removal of duplicates reduced the number of genes in the Foresee Object from ", dim_before, " to ", dim_after))

  # Update Object in the Environment

  Object <- Object_withoutDuplicates
  return(Object)

}


################################################################################
### Function "mean" calculates the mean of all rows that have the same gene name
#' @export
DuplicationHandler.mean <- function(Object, DuplicationHandling){

  Object_withoutDuplicates <- Object

  # Counts the number of duplicated gene names in the gene expression matrix
  dim_before <- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Removes duplicates
  Object_withoutDuplicates$GeneExpression <- aggregate(x = Object_withoutDuplicates$GeneExpression, by = list(rownames(Object_withoutDuplicates$GeneExpression)), FUN = mean)
  rownames(Object_withoutDuplicates$GeneExpression) <- Object_withoutDuplicates$GeneExpression[,1]
  Object_withoutDuplicates$GeneExpression <- data.matrix(Object_withoutDuplicates$GeneExpression[,-1])

  # Calculates the dimension of the gene expression matrix after removing duplicates
  dim_after<- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The removal of duplicates reduced the number of genes in the Foresee Object from ", dim_before, " to ", dim_after))

  # Update Object in the Environment
  Object <- Object_withoutDuplicates
  return(Object)

}

################################################################################
### Function "default" is called in case method in "DuplicationHandling" is
# unknown to DuplicationHandler
#' @export
DuplicationHandler.default <- function(Object, DuplicationHandling){
  stop(paste("Method",DuplicationHandling,"is not defined for handling duplicated genes!"))
}
