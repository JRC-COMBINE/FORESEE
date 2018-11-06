#' Select the Genes that are used as Model Features
#'
#' The FeatureSelector selects a subset of features from all genes to be used for drug efficacy prediction.
#'
#' @param TrainObject Object that contains all data needed to train a model, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type, etc.) and drug response data
#' @param TestObject Object that contains all data that the model is to be tested on, including molecular data (such as gene expression, mutation, copy number variation, methylation, cancer type, etc.) and drug response data
#' @param GeneFilter Set of genes to be considered for training the model, such as all, a certain percantage based on variance or p-value, specific gene sets like landmark genes, gene ontologies or pathways, etc.
#' The option 'variance' removes the 20 % genes of lowest variance across samples in the TrainObject
#' The option 'pvalue' removes the 20 % genes of lowest p-value (ttest) across samples in the TrainObject
#' The option 'landmarkgenes' uses the L1000 gene set downloaded from CLUE Command App
#' The option 'ontology' uses a specific set of genes included in the ontology associated with the drug
#' The option 'pathway' uses a specific set of genes included in the pathway associated with the drug
#' The option 'all' keeps all genes as features
#' The function 'listInputOptions("FeatureSelector")' returns a list of the possible options.
#' Instead of choosing one of the implemented options, a user-defined function can be used as an input.
#' If the user inserts a list as an input, the list is considered as chosen features.
#' @return \item{TrainObject}{The TrainObject with the selected gene set as features.}
#'         \item{TestObject}{The TestObject with homogenized features according to the chosen TrainObject.}
#' @examples
#' FeatureSelector(GDSC,GSE6434,"variance","Docetaxel")
#' FeatureSelector(GDSC,GSE6434,"pathway","Tamoxifen")
#' @export

FeatureSelector <- function(TrainObject, TestObject, GeneFilter, DrugName){
  UseMethod("FeatureSelector", object = GeneFilter)
}

#' @export
FeatureSelector.character <- function(TrainObject, TestObject, GeneFilter, DrugName){
  class(GeneFilter) <- GeneFilter;
  UseMethod("FeatureSelector", object = GeneFilter)
}


################################################################################
### Function "function" applies the function in "GeneFilter"
# to Train and Test objects
#' @export
FeatureSelector.function <- function(TrainObject, TestObject, GeneFilter, DrugName){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Calculate a measure by GeneFilter of each gene across all samples in TrainObject
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
#' @export
FeatureSelector.variance <- function(TrainObject, TestObject, GeneFilter, DrugName){

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
#' @export
FeatureSelector.pvalue <- function(TrainObject, TestObject, GeneFilter, DrugName){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Calculate ttest-p-value of each gene between sensitive and resistant samples in TrainObject

  # Order samples in TrainObject according to their DrugResponse
  num <- min(50,floor(length(TrainObject$DrugResponse)/2)) # in case length(DrugResponse) < 100
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
#' @export
FeatureSelector.landmarkgenes <- function(TrainObject, TestObject, GeneFilter, DrugName){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Using list of landmark genes save as LM_genes_entrez in package internal data:
  TrainObject_selectedfeatures$GeneExpression<-TrainObject_selectedfeatures$GeneExpression[rownames(TrainObject_selectedfeatures$GeneExpression) %in% LM_genes_entrez,]
  TestObject_selectedfeatures$GeneExpression<-TestObject_selectedfeatures$GeneExpression[rownames(TestObject_selectedfeatures$GeneExpression) %in% LM_genes_entrez,]

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
  assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
}

################################################################################
### Option "ontology" to keep ontology genes as features in the TrainObject
#' @export
FeatureSelector.ontology <- function(TrainObject, TestObject, GeneFilter, DrugName){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject
  if(!any(names(TrainObject)=="DrugInfo") | !any(names(TrainObject$DrugInfo)=="DRUG_NAME") | !any(names(TrainObject$DrugInfo)=="TARGET")) {
    warning("TrainObject doesn't have the information needed for ontology feature selection, will continue with all features...")
    # Update Objects in the Environment
    assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
    assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
  } else {
    if(!any(TrainObject$DrugInfo$DRUG_NAME == DrugName)){
      stop("Matching DrugName failed!!")
    } else {
      TargetGene <- TrainObject$DrugInfo$TARGET[TrainObject$DrugInfo$DRUG_NAME == DrugName]
    }
    if(!any(ConvTableGo2Sym$hgnc_symbol == TargetGene)){
      warning(
        paste0(
          "Matching Target-Gene failed, probably because the drug doesn't have a single target gene, the target we tried to match was '",TargetGene,"'"
        )
      )
      GenesWithRelevantGos <- character(length = 0)
    } else {
      TargetGos <- ConvTableGo2Sym$go_id[ConvTableGo2Sym$hgnc_symbol == TargetGene]
      TargetGos <- TargetGos[TargetGos!=""] #removing empty match(es)
      GenesWithRelevantGos <- sort(table(ConvTableGo2Sym$hgnc_symbol[ConvTableGo2Sym$go_id %in% TargetGos]), decreasing = T)
    }
    if(length(GenesWithRelevantGos) < 1000){ ## ToDo: Have this 1000 as a variable user feeds
      warning("Not enough features matched, will continue with all features")
      # Update Objects in the Environment
      assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
      assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
    } else {
      FirstFewEntrezGenesWithRelevantGos <- ConvTableSym2Entrez$entrezgene[match(names(GenesWithRelevantGos)[1:1000],ConvTableSym2Entrez$hgnc_symbol,nomatch = 0)]## ToDo: Have this 1000 as a variable user feeds
      TrainObject_selectedfeatures$GeneExpression<-TrainObject_selectedfeatures$GeneExpression[rownames(TrainObject_selectedfeatures$GeneExpression) %in% FirstFewEntrezGenesWithRelevantGos,]
      TestObject_selectedfeatures$GeneExpression<-TestObject_selectedfeatures$GeneExpression[rownames(TestObject_selectedfeatures$GeneExpression) %in% FirstFewEntrezGenesWithRelevantGos,]
      # Update Objects in the Environment
      assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
      assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
    }
  }
}

################################################################################
### Option "pathway" to keep pathway genes as features in the TrainObject
#' @export
FeatureSelector.pathway <- function(TrainObject, TestObject, GeneFilter, DrugName){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject
  if(!any(names(TrainObject)=="DrugInfo") | !any(names(TrainObject$DrugInfo)=="DRUG_NAME") | !any(names(TrainObject$DrugInfo)=="TARGET")) {
    warning("TrainObject doesn't have the information needed for pathway feature selection, will continue with all features...")
    # Update Objects in the Environment
    assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
    assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
  } else {
    if(!any(TrainObject$DrugInfo$DRUG_NAME == DrugName)){
      stop("Matching DrugName failed!!")
    } else {
      TargetGene <- TrainObject$DrugInfo$TARGET[TrainObject$DrugInfo$DRUG_NAME == DrugName]
      TargetGeneEntrez <- ConvTableSym2Entrez$entrezgene[match(TargetGene,ConvTableSym2Entrez$hgnc_symbol)]
      #TargetGeneEntrez <- ConvTableSym2Entrez$entrezgene[match(TargetGene,ConvTableSym2Entrez$hgnc_symbol)]
    }
    if(is.na(TargetGeneEntrez) | (!any(names(Entrez2PathID) == TargetGeneEntrez))){
      warning(
        paste0(
          "Matching Target-Gene (to find the target pathway) failed, probably because the drug doesn't have a single target gene, the target we tried to match was '",TargetGene,"'"
        )
      )
      GenesWithRelevantPaths <- 0
    } else {
      TargetPaths <- Entrez2PathID[[as.character(TargetGeneEntrez)]] #TargetGeneEntrez is a factor, and indexing with a factor return weird stuff hence the 'as.character'
      GenesWithRelevantPaths <- sort(table(unlist(PathID2Entrez[TargetPaths])), decreasing = T)
    }
    if(length(GenesWithRelevantPaths) < 1000){ ## ToDo: Have this 1000 as a variable user feeds
      warning("Not enough features matched, will continue with all features")
      # Update Objects in the Environment
      assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
      assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
    } else {
      FirstFewEntrezGenesWithRelevantPaths <- names(GenesWithRelevantPaths)[1:1000]## ToDo: Have this 1000 as a variable user feeds
      TrainObject_selectedfeatures$GeneExpression<-TrainObject_selectedfeatures$GeneExpression[rownames(TrainObject_selectedfeatures$GeneExpression) %in% FirstFewEntrezGenesWithRelevantPaths,]
      TestObject_selectedfeatures$GeneExpression<-TestObject_selectedfeatures$GeneExpression[rownames(TestObject_selectedfeatures$GeneExpression) %in% FirstFewEntrezGenesWithRelevantPaths,]
      # Update Objects in the Environment
      assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
      assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
    }
  }
}

################################################################################
### Option "all" to keep all genes as features in the TrainObject
#' @export
FeatureSelector.all <- function(TrainObject, TestObject, GeneFilter, DrugName){

  TrainObject_selectedfeatures <- TrainObject
  TestObject_selectedfeatures <- TestObject

  # Don't remove any genes

  # Update Objects in the Environment
  assign("TrainObject", value = TrainObject_selectedfeatures, envir = parent.frame())
  assign("TestObject", value = TestObject_selectedfeatures, envir = parent.frame())
}

