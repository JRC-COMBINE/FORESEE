#' Select the Genes that are used as Model Features
#'
#' The FeatureSelector selects features from all genes to be used for drug efficacy prediction.
#'
#' @param TrainObject Object that contains all data needed to train a model, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param TestObject Object that contains all data that the model is to be tested on, such as gene expression, mutation, copy number variation, methylation, cancer type, drug response data, etc.
#' @param GeneFilter Set of genes to be considered for training the model, such as all, a certain percantage based on variance or p-value, specific gene sets like landmark genes, gene ontologies or pathways, etc.
#' The option 'variance' removes the 20 % genes of lowest variance across samples in the TrainObject
#' The option 'pvalue' removes the 20 % genes of lowest p-value (ttest) across samples in the TrainObject
#' The option 'landmarkgenes' uses the L1000 gene set downloaded from CLUE Command App
#' The option 'ontology' uses a specific set of genes included in the chosen ontology? -> Ali
#' The option 'pathway' uses a specific set of genes included in the chosen pathway? -> Ali
#' The option 'all' keeps all genes as features
#' If the user inserts a list as an input, the list is considered as chosen features.
#' @return \item{TrainObject}{The TrainObject with the selected gene set as features.}
#'         \item{TestObject}{The TestObject with homogenized features according to the chosen TrainObject.}
#' @export


#ToDo3: Ontologies
#ToDo4: Pathways
#ToDo5: User-defined gene lists


FeatureSelector <- function(TrainObject,TestObject, GeneFilter){
  UseMethod("FeatureSelector", object = GeneFilter)
}

FeatureSelector.character <- function(TrainObject, TestObject, GeneFilter){
  class(GeneFilter) <- GeneFilter;
  UseMethod("FeatureSelector", object = GeneFilter)
}

# FeatureSelector.list <- function(TrainObject, TestObject, GeneFilter){
#   UseMethod("FeatureSelector", object = GeneFilter)
# }

################################################################################
### Function "function" applies the function in "GeneFilter"
# to Train and Test objects
FeatureSelector.function <- function(TrainObject, TestObject, GeneFilter){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Calculate a measure by GeneFilter  of each gene across all samples in TrainObject
  measure_TrainObject<-c()
  for (i in 1:dim(TrainObject_selectedfeatures$GeneExpression)[1]){
    measure_TrainObject[i] <- GeneFilter(TrainObject_selectedfeatures$GeneExpression[i,])
  }
  names(measure_TrainObject) <- rownames(TrainObject_selectedfeatures$GeneExpression)
  highest_measured_genes <- sort(measure_TrainObject,decreasing = TRUE)
  features <- names(highest_measured_genes[1:(0.8*length(highest_measured_genes))])

  TrainObject_selectedfeatures$GeneExpression<-TrainObject_selectedfeatures$GeneExpression[features,]
  TestObject_selectedfeatures$GeneExpression<-TestObject_selectedfeatures$GeneExpression[features,]

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
  assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
}

################################################################################
### Option "variance" to keep variance genes as features in the TrainObject
FeatureSelector.variance <- function(TrainObject, TestObject, GeneFilter){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Calculate variance of each gene across all samples in TrainObject
  variance_TrainObject<-c()
  for (i in 1:dim(TrainObject_selectedfeatures$GeneExpression)[1]){
    variance_TrainObject[i] <- var(TrainObject_selectedfeatures$GeneExpression[i,])
  }
  names(variance_TrainObject) <- rownames(TrainObject_selectedfeatures$GeneExpression)
  highest_variant_genes <- sort(variance_TrainObject,decreasing = TRUE)
  features <- names(highest_variant_genes[1:(0.8*length(highest_variant_genes))])

  TrainObject_selectedfeatures$GeneExpression<-TrainObject_selectedfeatures$GeneExpression[features,]
  TestObject_selectedfeatures$GeneExpression<-TestObject_selectedfeatures$GeneExpression[features,]

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
  assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
}

################################################################################
### Option "pvalue" to keep pvalue genes as features in the TrainObject
FeatureSelector.pvalue <- function(TrainObject, TestObject, GeneFilter){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Calculate ttest-p-value of each gene between sensitive and resistant samples in TrainObject

  # Order samples in TrainObject according to their DrugResponse
  num <- 50
  sensInd <- order(TrainObject$DrugResponse)[1:num]
  resInd <- order(TrainObject$DrugResponse)[(length(TrainObject$DrugResponse)-num):length(TrainObject$DrugResponse)]

  # Calculate t-test between most sensitive and most resistant samples to find significant genes
  pvalue_TrainObject<-c()
  # for each gene:
  for (i in 1:dim(TrainObject_selectedfeatures$GeneExpression)[1]){
    pvalue_TrainObject[i] <- t.test(TrainObject$GeneExpression[i,sensInd],
                                    TrainObject$GeneExpression[i,resInd])$p.value
  }

  # Get names of genes with lowest p-value:
  names(pvalue_TrainObject) <- rownames(TrainObject$GeneExpression)
  lowest_pvalued_genes <- sort(pvalue_TrainObject,decreasing = FALSE)
  features <- names(lowest_pvalued_genes[1:(0.8*length(lowest_pvalued_genes))])

  TrainObject_selectedfeatures$GeneExpression<-TrainObject_selectedfeatures$GeneExpression[features,]
  TestObject_selectedfeatures$GeneExpression<-TestObject_selectedfeatures$GeneExpression[features,]

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
  assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
}

################################################################################
### Option "landmarkgenes" to keep landmarkgenes genes as features in the TrainObject
FeatureSelector.landmarkgenes <- function(TrainObject, TestObject, GeneFilter){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Load list of landmark genes
  load("./data/LM_genes_entrez.rda")
  TrainObject_selectedfeatures$GeneExpression<-TrainObject_selectedfeatures$GeneExpression[rownames(TrainObject_selectedfeatures$GeneExpression) %in% LM_genes_entrez,]
  TestObject_selectedfeatures$GeneExpression<-TestObject_selectedfeatures$GeneExpression[rownames(TestObject_selectedfeatures$GeneExpression) %in% LM_genes_entrez,]

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
  assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
}

################################################################################
### Option "ontology" to keep ontology genes as features in the TrainObject
FeatureSelector.ontology <- function(TrainObject, TestObject, GeneFilter){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Ali's implementation

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
  assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
}

################################################################################
### Option "pathway" to keep pathway genes as features in the TrainObject
FeatureSelector.pathway <- function(TrainObject, TestObject, GeneFilter){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Ali's implementation

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
  assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
}

################################################################################
### Option "all" to keep all genes as features in the TrainObject
FeatureSelector.all <- function(TrainObject, TestObject, GeneFilter){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Don't remove any genes

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
  assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
}

