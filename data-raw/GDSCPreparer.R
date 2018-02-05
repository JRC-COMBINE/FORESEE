## Script that will download, process and add GDSC to this package
# by esfahani@aices.rwth-aachen.de, 19.1.2018

##Downloading files directly from sanger ftp server seemed to be buggy and unstable, so for now we assume user downloaded the last versions of:
#1-
#from sanger and put them in data-raw

RespondInfo <- read.csv2("data-raw/v17_fitted_dose_response.csv",  header = T, sep = ";") #read.csv will have problem with comma in place of decimal point, Damn excel!!
DrugInfo <- read.csv("data-raw/Screened_Compounds.csv", header = T, sep = ";")
CelllineInfo <- read.csv("data-raw/Cell_Lines_Details_Sheet1.csv", header = T, sep = ";")
TissueInfo <- read.csv("data-raw/Cell_Lines_Details_Sheet2.csv", header = T, sep = ";")

##Before matching drug names: There are drugs that have two matching IDs,
# gonna have to do something about them for now (ask what to do later!!)
DrugInfo$Drug.Name <- as.character(DrugInfo$Drug.Name)
DrugInfo$Drug.Name[duplicated(DrugInfo$Drug.Name)] <- paste0(DrugInfo$Drug.Name[duplicated(DrugInfo$Drug.Name)],"(2)")
DrugInfo$Drug.Name <- as.factor(DrugInfo$Drug.Name)
# now name matching:
RespondInfo$DRUG_NAME <- DrugInfo$Drug.Name[match(RespondInfo$DRUG_ID, DrugInfo$Drug.ID)]

RespondInfo$CELLLINE_NAME <- CelllineInfo$Sample.Name[match(RespondInfo$COSMIC_ID, CelllineInfo$COSMIC.identifier)]
##There are cosmic ids that doesn't have corrisponding cell name??!!
IC50 <- AUC <- RMSE <- matrix(NA, nrow = nlevels(RespondInfo$CELLLINE_NAME),
               ncol = nlevels(RespondInfo$DRUG_NAME),
               dimnames = list(levels(RespondInfo$CELLLINE_NAME),levels(RespondInfo$DRUG_NAME)))

require(progress)
pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                       total = length(IC50), clear = FALSE)

for(CellLine in levels(RespondInfo$CELLLINE_NAME)){
  for(DrUg in levels(RespondInfo$DRUG_NAME)){
    bullseyeIndex <- RespondInfo$DRUG_NAME %in% DrUg & RespondInfo$CELLLINE_NAME %in% CellLine
    if(sum(bullseyeIndex)>1) stop("Something's wrong: couldn't match properly!!")
    if(sum(bullseyeIndex)==1){
      IC50[CellLine, DrUg] <- RespondInfo$LN_IC50[bullseyeIndex]
      AUC[CellLine, DrUg] <- RespondInfo$AUC[bullseyeIndex]
      RMSE[CellLine, DrUg] <- RespondInfo$RMSE[bullseyeIndex]
    }
    pb$tick()
  }
}

# #Let's do the above parallel:
# require(doParallel)
# CL <- makeCluster(4, outfile="")
# registerDoParallel(CL)
# require(progress)
# pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
#                        total = length(IC50), clear = FALSE)
#
# foreach(CellLine=levels(RespondInfo$CELLLINE_NAME)[1:nlevels(RespondInfo$CELLLINE_NAME)], .combine = "c") %dopar% {
#   for(DrUg in levels(RespondInfo$DRUG_NAME)){
#     bullseyeIndex <- RespondInfo$DRUG_NAME %in% DrUg & RespondInfo$CELLLINE_NAME %in% CellLine
#     if(sum(bullseyeIndex)>1) stop("Something's wrong: couldn't match properly!!")
#     if(sum(bullseyeIndex)==1){
#       IC50[CellLine, DrUg] <- RespondInfo$LN_IC50[bullseyeIndex]
#       AUC[CellLine, DrUg] <- RespondInfo$AUC[bullseyeIndex]
#       RMSE[CellLine, DrUg] <- RespondInfo$RMSE[bullseyeIndex]
#     }
#     pb$tick()
#   }
# }
# stopCluster(CL)

GeneExpression <- read.csv("data-raw/sanger1018_brainarray_ensemblgene_rma.txt",
                           sep = "\t", row.names = 1)
colnames(GeneExpression) <- CelllineInfo$Sample.Name[match(as.integer(substr(x = colnames(GeneExpression),
                                                                             start = 2, stop = nchar(colnames(GeneExpression)))),
                                                           CelllineInfo$COSMIC.identifier)]

require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTabelle <- getBM(attributes = c("ensembl_gene_id","entrezgene"),
      filters = "ensembl_gene_id", values = rownames(GeneExpression), mart = humaRt)
GeneExpression <- as.matrix(GeneExpression)
rownames(GeneExpression) <- ConvTabelle$entrezgene[match(rownames(GeneExpression),ConvTabelle$ensembl_gene_id)]
GeneExpression <- GeneExpression[!is.na(rownames(GeneExpression)),] #Losing 673 genes, is Entrez a good idea? doublecheck with others!

#Loading Iorio et. al. 2016 info for mut, cna, methylation:
IorioBinaryMat <- read.csv(file = "data-raw/mmc4S3B.csv", header = T, sep = ";", skip = 6, row.names = 1)
splittedRowNames <- strsplit(x = rownames(IorioBinaryMat), split = "_")
CNAsIndx <- !sapply(splittedRowNames,length) == 2
CNAs <- IorioBinaryMat[CNAsIndx,]
CNAsWithMatchedGenes <- CNAs[grep(rownames(CNAs), pattern = " ", fixed = T),]
CNAsWithMatchedGenesGain <- CNAsWithMatchedGenes[substr(x = rownames(CNAsWithMatchedGenes),
                                                        start = 1, stop = 4)=="gain",]
CNAsWithMatchedGenesLoss <- CNAsWithMatchedGenes[substr(x = rownames(CNAsWithMatchedGenes),
                                                        start = 1, stop = 4)=="loss",]
#Have to convert cause there're duplications in rownames:
CNAsWithMatchedGenesGain <- as.matrix(CNAsWithMatchedGenesGain)
CNAsWithMatchedGenesLoss <- as.matrix(CNAsWithMatchedGenesLoss)

rownames(CNAsWithMatchedGenesGain) <- sapply(strsplit(rownames(CNAsWithMatchedGenesGain),
                                                      split = " ", fixed = T),
                                             function(x) x[2])
rownames(CNAsWithMatchedGenesLoss) <- sapply(strsplit(rownames(CNAsWithMatchedGenesLoss),
                                                      split = " ", fixed = T),
                                             function(x) x[2])
rownames(CNAsWithMatchedGenesGain) <- substr(x = rownames(CNAsWithMatchedGenesGain),
                                             start = 2,
                                             stop = nchar(rownames(CNAsWithMatchedGenesGain))-1)
rownames(CNAsWithMatchedGenesLoss) <- substr(x = rownames(CNAsWithMatchedGenesLoss),
                                             start = 2,
                                             stop = nchar(rownames(CNAsWithMatchedGenesLoss))-1)
CNAsWithMatchedGenesGainSignleGenes <- matrix(NA, nrow = 0, ncol = ncol(CNAsWithMatchedGenesGain))
for(K in 1:nrow(CNAsWithMatchedGenesGain)){
  GenesInRowN <- strsplit(rownames(CNAsWithMatchedGenesGain)[K], split = ",", fixed = T)[[1]]
  for(GENE in GenesInRowN){
    tmpMat <- CNAsWithMatchedGenesGain[K,,drop=F]
    rownames(tmpMat) <- GENE
    CNAsWithMatchedGenesGainSignleGenes <- rbind(CNAsWithMatchedGenesGainSignleGenes,
                                                 tmpMat)
  }
}
CNAsWithMatchedGenesLossSignleGenes <- matrix(NA, nrow = 0, ncol = ncol(CNAsWithMatchedGenesLoss))
for(K in 1:nrow(CNAsWithMatchedGenesLoss)){
  GenesInRowN <- strsplit(rownames(CNAsWithMatchedGenesLoss)[K], split = ",", fixed = T)[[1]]
  for(GENE in GenesInRowN){
    tmpMat <- CNAsWithMatchedGenesLoss[K,,drop=F]
    rownames(tmpMat) <- GENE
    CNAsWithMatchedGenesLossSignleGenes <- rbind(CNAsWithMatchedGenesLossSignleGenes,
                                                 tmpMat)
  }
}

ConvTabelle <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                     filters = "hgnc_symbol",
                     values = rownames(CNAsWithMatchedGenesGainSignleGenes), mart = humaRt)
rownames(CNAsWithMatchedGenesGainSignleGenes) <-
  ConvTabelle$entrezgene[match(rownames(CNAsWithMatchedGenesGainSignleGenes),
                               ConvTabelle$hgnc_symbol)]
CNAsWithMatchedGenesGainSignleGenes <-
  CNAsWithMatchedGenesGainSignleGenes[!is.na(rownames(CNAsWithMatchedGenesGainSignleGenes)),]
ConvTabelle <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                     filters = "hgnc_symbol",
                     values = rownames(CNAsWithMatchedGenesLossSignleGenes), mart = humaRt)
rownames(CNAsWithMatchedGenesLossSignleGenes) <-
  ConvTabelle$entrezgene[match(rownames(CNAsWithMatchedGenesLossSignleGenes),
                               ConvTabelle$hgnc_symbol)]
CNAsWithMatchedGenesLossSignleGenes <-
  CNAsWithMatchedGenesLossSignleGenes[!is.na(rownames(CNAsWithMatchedGenesLossSignleGenes)),]

CellLineNames <- read.csv(file = "data-raw/mmc4S3B.csv", header = F, sep = ";", skip = 3, nrows = 1, row.names = 1)
colnames(CNAsWithMatchedGenesGainSignleGenes) <- as.character(as.matrix(CellLineNames))
colnames(CNAsWithMatchedGenesLossSignleGenes) <- as.character(as.matrix(CellLineNames))
##Have to match the cell-line names better, right now not good!!! see: match(colnames(CNAsWithMatchedGenesGainSignleGenes), colnames(GDSC$GeneExpression))

#Now to mutation and methylation:
MutsAndMethyls <- IorioBinaryMat[!CNAsIndx,]
MutsAndMethylsRowDecomposed <- strsplit(x = rownames(MutsAndMethyls),
                                               split = "_", fixed = T)
Muts <- MutsAndMethyls[sapply(MutsAndMethylsRowDecomposed, function(x) x[2]=="mut"),]
rownames(Muts) <- sapply(strsplit(x = rownames(Muts),split = "_", fixed = T), function(x) x[1])

ConvTabelle <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                     filters = "hgnc_symbol",
                     values = rownames(Muts), mart = humaRt)
#Duplication in row names so:
Muts <- as.matrix(Muts)

rownames(Muts) <- ConvTabelle$entrezgene[match(rownames(Muts), ConvTabelle$hgnc_symbol)]
Muts <- Muts[!is.na(rownames(Muts)),]
colnames(Muts) <- as.character(as.matrix(CellLineNames))


Methyls <- MutsAndMethyls[sapply(MutsAndMethylsRowDecomposed, function(x) x[2]=="HypMET"),]
rownames(Methyls) <- sapply(strsplit(x = rownames(Methyls),split = "_", fixed = T),
                            function(x) x[1])
anyMatchingGeneNames <- sapply(strsplit(rownames(Methyls),
                                    split = "[(]"),function(x) substr(x = x[2],
                                                                      start = 1,
                                                                      stop = nchar(x[2])-1))
#Duplication in row names so:
Methyls <- as.matrix(Methyls)

rownames(Methyls) <- anyMatchingGeneNames
MethylsWithMatchedGenes <- Methyls[rownames(Methyls)!="", ]

ConvTabelle <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                     filters = "hgnc_symbol",
                     values = rownames(MethylsWithMatchedGenes), mart = humaRt)
rownames(MethylsWithMatchedGenes) <- ConvTabelle$entrezgene[match(rownames(MethylsWithMatchedGenes), ConvTabelle$hgnc_symbol)]
MethylsWithMatchedGenes <- MethylsWithMatchedGenes[!is.na(rownames(MethylsWithMatchedGenes)),]
colnames(MethylsWithMatchedGenes) <- as.character(as.matrix(CellLineNames))


#Making the Foresee object:
GDSC <- list()
class(GDSC) <- "ForeseeCell"
GDSC[["GeneExpression"]] <- GeneExpression
GDSC[["CNVGain"]] <- CNAsWithMatchedGenesGainSignleGenes
GDSC[["CNVLoss"]] <- CNAsWithMatchedGenesLossSignleGenes
GDSC[["Mutation"]] <- Muts
GDSC[["Methylation"]] <- MethylsWithMatchedGenes
GDSC[["IC50"]] <- IC50
GDSC[["AUC"]] <- AUC
GDSC[["RMSE"]] <- RMSE
GDSC[["DrugInfo"]] <- DrugInfo
GDSC[["CelllineInfo"]] <- CelllineInfo
GDSC[["TissueInfo"]] <- TissueInfo


devtools::use_data(GDSC, overwrite = T)




