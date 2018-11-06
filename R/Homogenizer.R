#' Homogenize a Pair of two FORESEE Objects
#'
#' The Homogenizer function reduces batch effects and homogenizes the data of a given TrainObject and TestObject.
#'
#' @param TrainObject Object that contains all data needed to train a model, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type, etc. ) and drug response data
#' @param TestObject Object that contains all data that the model is to be tested on, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type, etc. ) and drug response data
#' @param HomogenizationMethod Method for homogenizing data of the TrainObject and TestObject.
#' The function 'ComBat' uses the batch effect removal ComBat of the sva package,
#' The function 'quantile' uses the quantile normalization of the preprocessCore package,
#' The function 'limma' uses the removeBatchEffect function of the limma package,
#' The function 'YuGene' uses the function of the YuGene package by Le Cao,
#' The function 'RUV' regresses out unwanted variation using the 10 principal components of negative control genes (here: list of human housekeeping by Eisenberg and Levanon (2013))
#' The function 'RUV4' regresses out unwanted variation using the ruv package,
#' The function 'none' does not do any batch effect correction.
#' The function 'listInputOptions("Homogenizer")' returns a list of the possible options.
#' Instead of choosing one of the implemented options, a user-defined function can be used as an input.
#' @return \item{TrainObject}{The TrainObject with homogenized features according to the chosen TestObject.}
#'         \item{TestObject}{The TestObject with homogenized features according to the chosen TrainObject.}
#' @examples
#' Homogenizer(GDSC,GSE6434,"quantile")
#' @export

#########################
# This file is part of the FORESEE R-package
# File authors: Lisa-Katrin Turnhoff <turnhoff@combine.rwth-aachen.de> and Ali Hadizadeh Esfahani <hadizadeh@combine.rwth-aachen.de>
# Distributed under the GNU General Public License v3.0.(http://www.gnu.org/licenses/gpl-3.0.html)
#########################

Homogenizer <- function(TrainObject, TestObject, HomogenizationMethod){
  UseMethod("Homogenizer", object = HomogenizationMethod)
}

#' @export
Homogenizer.character <- function(TrainObject, TestObject, HomogenizationMethod){
  class(HomogenizationMethod) <- HomogenizationMethod;
  UseMethod("Homogenizer", object = HomogenizationMethod)
}

################################################################################
### Function "function" applies the function in "HomogenizationMethod"
# to Train and Test objects
#' @export
Homogenizer.function <- function(TrainObject, TestObject, HomogenizationMethod) {

  TrainObject_homogenized <- TrainObject
  TestObject_homogenized <- TestObject

  # Removes gene expression matrices to common gene names only
  CommonGenes <- rownames(TrainObject_homogenized$GeneExpression)[rownames(TrainObject_homogenized$GeneExpression) %in% rownames(TestObject_homogenized$GeneExpression)]
  TrainObject_homogenized$GeneExpression <- TrainObject_homogenized$GeneExpression[CommonGenes,]
  TestObject_homogenized$GeneExpression <- TestObject_homogenized$GeneExpression[CommonGenes,]

  ### What to do with User-defined function goes here
  homogenizedTrainAndTestGex <- HomogenizationMethod(TrainObject_homogenized$GeneExpression, TestObject_homogenized$GeneExpression)
  TrainObject_homogenized$GeneExpression <- homogenizedTrainAndTestGex[[1]]
  TestObject_homogenized$GeneExpression <- homogenizedTrainAndTestGex[[2]]

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_homogenized, envir = parent.frame())
  assign("TestObject", value = TestObject_homogenized, envir = parent.frame())
}


################################################################################
### Function "ComBat"
#' @export
Homogenizer.ComBat <- function(TrainObject, TestObject, HomogenizationMethod){

  # Load package
  requireForesee(sva)

  TrainObject_homogenized <- TrainObject
  TestObject_homogenized <- TestObject

  # Counts the number of gene names in the objects' gene expression matrices
  dim_before_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_before_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Removes gene expression matrices to common gene names only
  CommonGenes <- rownames(TrainObject_homogenized$GeneExpression)[rownames(TrainObject_homogenized$GeneExpression) %in% rownames(TestObject_homogenized$GeneExpression)]
  TrainObject_homogenized$GeneExpression <- TrainObject_homogenized$GeneExpression[CommonGenes,]
  TestObject_homogenized$GeneExpression <- TestObject_homogenized$GeneExpression[CommonGenes,]

  # Homogenizes both gene expression Data sets
  BothBatches <- cbind(TrainObject_homogenized$GeneExpression,TestObject_homogenized$GeneExpression)
  BatchIndex <- as.factor(c(rep("trainobject", ncol(TrainObject_homogenized$GeneExpression)), rep("testobject", ncol(TestObject_homogenized$GeneExpression))))
  ComBatObject <- ComBat(BothBatches, BatchIndex)

  TrainObject_homogenized$GeneExpression <- ComBatObject[, BatchIndex=="trainobject"]
  TestObject_homogenized$GeneExpression <- ComBatObject[, BatchIndex=="testobject"]

  # Counts the number of gene names in the objects' gene expression matrices after homogenization
  dim_after_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_after_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The homogenization of both gene expression matrices reduced the number of common genes in the Foresee objects to ", dim_after_train))

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_homogenized, envir = parent.frame())
  assign("TestObject", value = TestObject_homogenized, envir = parent.frame())
}



################################################################################
### Function "quantile"
#' @export
Homogenizer.quantile <- function(TrainObject, TestObject, HomogenizationMethod){

  # Load package
  requireForesee(preprocessCore)

  TrainObject_homogenized <- TrainObject
  TestObject_homogenized <- TestObject

  # Counts the number of gene names in the objects' gene expression matrices
  dim_before_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_before_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Removes gene expression matrices to common gene names only
  CommonGenes <- rownames(TrainObject_homogenized$GeneExpression)[rownames(TrainObject_homogenized$GeneExpression) %in% rownames(TestObject_homogenized$GeneExpression)]
  TrainObject_homogenized$GeneExpression <- TrainObject_homogenized$GeneExpression[CommonGenes,]
  TestObject_homogenized$GeneExpression <- TestObject_homogenized$GeneExpression[CommonGenes,]

  # Homogenizes both gene expression Data sets
  BothBatches <- cbind(TrainObject_homogenized$GeneExpression,TestObject_homogenized$GeneExpression)
  BatchIndex <- as.factor(c(rep("trainobject", ncol(TrainObject_homogenized$GeneExpression)), rep("testobject", ncol(TestObject_homogenized$GeneExpression))))
  quantileObject <- normalize.quantiles(BothBatches)

  # "normalize.quantiles" doesn't keep dimnames, have to transfer row and column names manually:
  dimnames(quantileObject) <- dimnames(BothBatches)

  TrainObject_homogenized$GeneExpression <- quantileObject[, BatchIndex=="trainobject"]
  TestObject_homogenized$GeneExpression <- quantileObject[, BatchIndex=="testobject"]

  # Counts the number of gene names in the objects' gene expression matrices after homogenization
  dim_after_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_after_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The homogenization of both gene expression matrices reduced the number of common genes in the Foresee objects to ", dim_after_train))

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_homogenized, envir = parent.frame())
  assign("TestObject", value = TestObject_homogenized, envir = parent.frame())
}


################################################################################
### Function "limma"
#' @export
Homogenizer.limma <- function(TrainObject, TestObject, HomogenizationMethod){

  # Load package
  requireForesee(limma)

  TrainObject_homogenized <- TrainObject
  TestObject_homogenized <- TestObject

  # Counts the number of gene names in the objects' gene expression matrices
  dim_before_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_before_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Removes gene expression matrices to common gene names only
  CommonGenes <- rownames(TrainObject_homogenized$GeneExpression)[rownames(TrainObject_homogenized$GeneExpression) %in% rownames(TestObject_homogenized$GeneExpression)]
  TrainObject_homogenized$GeneExpression <- TrainObject_homogenized$GeneExpression[CommonGenes,]
  TestObject_homogenized$GeneExpression <- TestObject_homogenized$GeneExpression[CommonGenes,]

  # Homogenizes both gene expression Data sets
  BothBatches <- cbind(TrainObject_homogenized$GeneExpression,TestObject_homogenized$GeneExpression)
  BatchIndex <- as.factor(c(rep("trainobject", ncol(TrainObject_homogenized$GeneExpression)), rep("testobject", ncol(TestObject_homogenized$GeneExpression))))
  limmaObject <- removeBatchEffect(x=BothBatches, batch=BatchIndex)

  TrainObject_homogenized$GeneExpression <- limmaObject[, BatchIndex=="trainobject"]
  TestObject_homogenized$GeneExpression <- limmaObject[, BatchIndex=="testobject"]

  # Counts the number of gene names in the objects' gene expression matrices after homogenization
  dim_after_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_after_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The homogenization of both gene expression matrices reduced the number of common genes in the Foresee objects to ", dim_after_train))

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_homogenized, envir = parent.frame())
  assign("TestObject", value = TestObject_homogenized, envir = parent.frame())
}

################################################################################
### Function "YuGene"
#' @export
Homogenizer.YuGene <- function(TrainObject, TestObject, HomogenizationMethod){

  # Load package
  requireForesee(YuGene)

  TrainObject_homogenized <- TrainObject
  TestObject_homogenized <- TestObject

  # Counts the number of gene names in the objects' gene expression matrices
  dim_before_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_before_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Removes gene expression matrices to common gene names only
  CommonGenes <- rownames(TrainObject_homogenized$GeneExpression)[rownames(TrainObject_homogenized$GeneExpression) %in% rownames(TestObject_homogenized$GeneExpression)]
  TrainObject_homogenized$GeneExpression <- TrainObject_homogenized$GeneExpression[CommonGenes,]
  TestObject_homogenized$GeneExpression <- TestObject_homogenized$GeneExpression[CommonGenes,]

  # Homogenizes both gene expression Data sets
  BothBatches <- cbind(TrainObject_homogenized$GeneExpression,TestObject_homogenized$GeneExpression)
  BatchIndex <- as.factor(c(rep("trainobject", ncol(TrainObject_homogenized$GeneExpression)), rep("testobject", ncol(TestObject_homogenized$GeneExpression))))
  YuGeneObject <- YuGene(data.prop = BothBatches, progressBar = TRUE)

  TrainObject_homogenized$GeneExpression <-YuGeneObject[, BatchIndex=="trainobject"]
  TestObject_homogenized$GeneExpression <- YuGeneObject[, BatchIndex=="testobject"]

  # Counts the number of gene names in the objects' gene expression matrices after homogenization
  dim_after_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_after_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The homogenization of both gene expression matrices reduced the number of common genes in the Foresee objects to ", dim_after_train))

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_homogenized, envir = parent.frame())
  assign("TestObject", value = TestObject_homogenized, envir = parent.frame())
}

################################################################################
### Function "RUV"
#' @export
Homogenizer.RUV <- function(TrainObject, TestObject, HomogenizationMethod){

  # requireForesee(ruv)
  TrainObject_homogenized <- TrainObject
  TestObject_homogenized <- TestObject

  # Counts the number of gene names in the objects' gene expression matrices
  dim_before_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_before_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Removes gene expression matrices to common gene names only
  CommonGenes <- rownames(TrainObject_homogenized$GeneExpression)[rownames(TrainObject_homogenized$GeneExpression) %in% rownames(TestObject_homogenized$GeneExpression)]
  TrainObject_homogenized$GeneExpression <- TrainObject_homogenized$GeneExpression[CommonGenes,]
  TestObject_homogenized$GeneExpression <- TestObject_homogenized$GeneExpression[CommonGenes,]

  # Load list of housekeeping genes
  # "Human housekeeping genes revisited", E. Eisenberg and E.Y. Levanon, Trends in Genetics, 29 (2013)
  # saved as LM_genes_entrez in package internal data:
  NegativeControl <-  rownames(TestObject_homogenized$GeneExpression) %in% HK_genes_entrez

  # Create Data out of both batches
  BothBatches <- cbind(TrainObject_homogenized$GeneExpression,TestObject_homogenized$GeneExpression)
  BatchIndex <- as.factor(c(rep("trainobject", ncol(TrainObject_homogenized$GeneExpression)), rep("testobject", ncol(TestObject_homogenized$GeneExpression))))

  # Apply principal component analysis to gene expression matrix of only the control genes
  pca_NegativeControl <- princomp(BothBatches[NegativeControl,])
  loadings_10PCs<- pca_NegativeControl$loadings[,1:10]
  scores_10PCs<- pca_NegativeControl$scores[,1:10]
  LinearModel_RUV<-BothBatches
  for (i in 1:dim(LinearModel_RUV)[1]){
  LinearModel_RUV[i,] <- lm(BothBatches[i,]~loadings_10PCs)$residuals
  }

  TrainObject_homogenized$GeneExpression <- LinearModel_RUV[, BatchIndex=="trainobject"]
  TestObject_homogenized$GeneExpression <- LinearModel_RUV[, BatchIndex=="testobject"]

  # Counts the number of gene names in the objects' gene expression matrices after homogenization
  dim_after_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_after_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The homogenization of both gene expression matrices reduced the number of common genes in the Foresee objects to ", dim_after_train))

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_homogenized, envir = parent.frame())
  assign("TestObject", value = TestObject_homogenized, envir = parent.frame())
}

################################################################################
### Function "RUV", based on RUV4 in package "ruv"
#' @export
Homogenizer.RUV4 <- function(TrainObject, TestObject, HomogenizationMethod){

  requireForesee(ruv)
  TrainObject_homogenized <- TrainObject
  TestObject_homogenized <- TestObject

  # Counts the number of gene names in the objects' gene expression matrices
  dim_before_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_before_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Removes gene expression matrices to common gene names only
  CommonGenes <- rownames(TrainObject_homogenized$GeneExpression)[rownames(TrainObject_homogenized$GeneExpression) %in% rownames(TestObject_homogenized$GeneExpression)]
  TrainObject_homogenized$GeneExpression <- TrainObject_homogenized$GeneExpression[CommonGenes,]
  TestObject_homogenized$GeneExpression <- TestObject_homogenized$GeneExpression[CommonGenes,]


  # Load list of housekeeping genes
  # "Human housekeeping genes revisited", E. Eisenberg and E.Y. Levanon, Trends in Genetics, 29 (2013)
  # load("./data/HK_genes_entrez.rda")
  NegativeControl <-  rownames(TestObject_homogenized$GeneExpression) %in% HK_genes_entrez

  # Create Data out of both batches
  BothBatches <- cbind(TrainObject_homogenized$GeneExpression,TestObject_homogenized$GeneExpression)
  BatchIndex <- as.factor(c(rep("trainobject", ncol(TrainObject_homogenized$GeneExpression)), rep("testobject", ncol(TestObject_homogenized$GeneExpression))))

  RUVObject <- RUV4(Y=t(BothBatches), X=matrix(c(rep(1, ncol(TrainObject_homogenized$GeneExpression)),
                                          rep(2, ncol(TestObject_homogenized$GeneExpression))),
                                          ncol = 1),
                    ctl = NegativeControl, k = 10)

  # BothBatchesHomogenized <- t(as.matrix(c(rep(1, ncol(TrainObject_homogenized$GeneExpression)),
  #             rep(2, ncol(TestObject_homogenized$GeneExpression)))) %*% RUVObject$betahat)
  BothBatchesHomogenized <- t(matrix(1, nrow = ncol(TrainObject_homogenized$GeneExpression)+
                                       ncol(TestObject_homogenized$GeneExpression), ncol = 1) %*% RUVObject$betahat) +
                            t(RUVObject$W %*% RUVObject$alpha)
  TrainObject_homogenized$GeneExpression <- BothBatchesHomogenized[, BatchIndex=="trainobject"]
  TestObject_homogenized$GeneExpression <- BothBatchesHomogenized[, BatchIndex=="testobject"]

  # Counts the number of gene names in the objects' gene expression matrices after homogenization
  dim_after_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_after_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Recover the lost 'colnames' in the process:
  colnames(TrainObject_homogenized$GeneExpression) <- colnames(TrainObject$GeneExpression)
  colnames(TestObject_homogenized$GeneExpression) <- colnames(TestObject$GeneExpression)

  # Prints the reduction of gene names
  message(paste0("The homogenization of both gene expression matrices reduced the number of common genes in the Foresee objects to ", dim_after_train))

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_homogenized, envir = parent.frame())
  assign("TestObject", value = TestObject_homogenized, envir = parent.frame())
}

################################################################################
### Function "none"
#' @export
Homogenizer.none <- function(TrainObject, TestObject, HomogenizationMethod){

  TrainObject_homogenized <- TrainObject
  TestObject_homogenized <- TestObject

  # Counts the number of gene names in the objects' gene expression matrices
  dim_before_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_before_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Removes gene expression matrices to common gene names only
  CommonGenes <- rownames(TrainObject_homogenized$GeneExpression)[rownames(TrainObject_homogenized$GeneExpression) %in% rownames(TestObject_homogenized$GeneExpression)]
  TrainObject_homogenized$GeneExpression <- TrainObject_homogenized$GeneExpression[CommonGenes,]
  TestObject_homogenized$GeneExpression <- TestObject_homogenized$GeneExpression[CommonGenes,]

  # Don't do any batch correction

  # Counts the number of gene names in the objects' gene expression matrices after homogenization
  dim_after_train <- dim(TrainObject_homogenized$GeneExpression)[1]
  dim_after_test <- dim(TestObject_homogenized$GeneExpression)[1]

  # Prints the reduction of gene names
  message(paste0("The homogenization of both gene expression matrices reduced the number of common genes in the Foresee objects to ", dim_after_train))

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_homogenized, envir = parent.frame())
  assign("TestObject", value = TestObject_homogenized, envir = parent.frame())
}

################################################################################
### Function "default" is called in case method in "HomogenizationMethod" is
# unknown to Homogenizer
#' @export
Homogenizer.default <- function(TrainObject, TestObject, HomogenizationMethod){
  stop(paste("Method",HomogenizationMethod,"is not defined as a homogenization method!"))
}

