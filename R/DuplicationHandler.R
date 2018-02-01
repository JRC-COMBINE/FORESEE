#' Remove Duplicated Gene Names from a FORESEE Object
#'
#' DuplicationHandler removes duplicates in the gene names from the FORESEE Object.
#'
#' @param Object FORESEE Object (ForeseeCell or ForeseeTrain) that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param DuplicationHandling Method for handling duplicates of gene names.
#' The function 'mean' calculates the mean of all rows that have the same gene name,
#' The function 'first' chooses the first occuring row of duplicated genes only,
#' The function 'none' removes all genes that occur more than once.

#' @return \item{Object}{The object without duplicated gene names}
#' @export

#ToDo1: Should we include a line to remove NAs?
#ToDo2: Remove duplicates from molecular data other than gene expression as well!

DuplicationHandler <- function(Object, DuplicationHandling){
  UseMethod("DuplicationHandler", object = DuplicationHandling)
}

DuplicationHandler.character <- function(Object, DuplicationHandling){
  class(DuplicationHandling) <- DuplicationHandling;
  UseMethod("DuplicationHandler", object = DuplicationHandling)
}

################################################################################
### Function "function" applies the function in "DuplicationHandling"
# to all rows that have the same gene name
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
  print(paste0("The removal of duplicates reduced the number of genes in the Foresee Object from ", dim_before, " to ", dim_after))

  # Update Object in the Environment
  return(Object_withoutDuplicates)

}

################################################################################
### Function "first" to choose the first occuring row of duplicated genes only
DuplicationHandler.first <- function(Object, DuplicationHandling){

  Object_withoutDuplicates <- Object

  # Counts the number of gene names in the gene expression matrix
  dim_before <- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Removes duplicates
  Object_withoutDuplicates$GeneExpression<-Object_withoutDuplicates$GeneExpression[!duplicated(rownames(Object_withoutDuplicates$GeneExpression)),]

  # Calculates the dimension of the gene expression matrix after removing duplicates
  dim_after <- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Prints the reduction of gene names
  print(paste0("The removal of duplicates reduced the number of genes in the Foresee Object from ", dim_before, " to ", dim_after))

  # Update Object in the Environment
  return(Object_withoutDuplicates)
}

################################################################################
### Function none" removes all genes that occur more than once.
DuplicationHandler.none <- function(Object, DuplicationHandling){

  Object_withoutDuplicates <- Object

  # Counts the number of duplicated gene names in the gene expression matrix
  dim_before <- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Removes duplicates
  Object_withoutDuplicates$GeneExpression<-Object_withoutDuplicates$GeneExpression[!(duplicated(rownames(Object_withoutDuplicates$GeneExpression)) | duplicated(rownames(Object_withoutDuplicates$GeneExpression), fromLast = TRUE)),]

  # Calculates the dimension of the gene expression matrix after removing duplicates
  dim_after<- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Prints the reduction of gene names
  print(paste0("The removal of duplicates reduced the number of genes in the Foresee Object from ", dim_before, " to ", dim_after))

  # Update Object in the Environment
  return(Object_withoutDuplicates)
}


################################################################################
### Function "mean" calculates the mean of all rows that have the same gene name
DuplicationHandler.mean <- function(Object, DuplicationHandling){

  Object_withoutDuplicates <- Object

  # Counts the number of duplicated gene names in the gene expression matrix
  dim_before <- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Removes duplicates
  Object_withoutDuplicates$GeneExpression <- data.matrix(aggregate(x = Object_withoutDuplicates$GeneExpression, by = list(rownames(Object_withoutDuplicates$GeneExpression)), FUN = mean))
  rownames(Object_withoutDuplicates$GeneExpression) <- Object_withoutDuplicates$GeneExpression[,1]
  Object_withoutDuplicates$GeneExpression <- Object_withoutDuplicates$GeneExpression[,2:dim(Object_withoutDuplicates$GeneExpression)[2]]

  # Calculates the dimension of the gene expression matrix after removing duplicates
  dim_after<- dim(Object_withoutDuplicates$GeneExpression)[1]

  # Prints the reduction of gene names
  print(paste0("The removal of duplicates reduced the number of genes in the Foresee Object from ", dim_before, " to ", dim_after))

  # Update Object in the Environment
  return(Object_withoutDuplicates)

}

################################################################################
### Function "default" is called in case method in "DuplicationHandling" is
# unknown to DuplicationHandler
DuplicationHandler.default <- function(Object, DuplicationHandling){
  stop(paste("Method",DuplicationHandling,"is not defined for handling duplicated genes!"))
}
