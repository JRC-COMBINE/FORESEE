## Script that will (download,) process and add GDSC to this package
# by esfahani@aices.rwth-aachen.de, 19.1.2018

##Downloading files directly from sanger ftp server seemed to be buggy and unstable,
#so for now we assume user downloaded the relevant versions of data files
#from sanger and put them in GDSC-raw.

## List and link of downloaded data, all downloaded on 25.4.2018 from https://www.cancerrxgene.org/downloads:
#1- log(IC50) and AUC values, named v17.3_fitted_dose_response.xlsx, last updated on March 27th 2018,
#with the link ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/v17.3_fitted_dose_response.xlsx
#2- Screened compounds, named Screened_Compounds.xlsx, last updated on March 27th 2018,
#with the link ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Screened_Compounds.xlsx
#3- Annotated list of Cell lines, named Cell_Lines_Details.xlsx, last updated on July 4th 2016,
# with the link ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/Cell_Lines_Details.xlsx
#4- RMA normalised gene expression data for Cell lines, named sanger1018_brainarray_ensemblgene_rma.txt.gz, last updated on March 2nd 2017,
#with the link ftp://ftp.sanger.ac.uk/pub/project/cancerrxgene/releases/current_release/sanger1018_brainarray_ensemblgene_rma.txt.gz
#5*- RACSs CNV BEMs for cell lines, named CellLines_CNV_BEMs.zip, last updated on July 4th 2016,
#with the link http://www.cancerrxgene.org/gdsc1000/GDSC1000_WebResources/Data/BEMs/CellLines/CellLines_CNV_BEMs.zip
## * -> As an alternative to 5, we used the original binary matrix published in Iorio F, et al. Cell. 2016. rather than the files provided by
# cancerrxgene.org since in the original publication the binary CNV, Mutation and Methylation of all cell lines are in one big binary matrix,
# making it easier to work with for our case. Here's the info of the file we worked with:
#*- Global binary matrix summarising the status of the Selected functional events across cell lines and tumours both, related to figure 3, named mmc4.xlsx,
#with the link https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4967469/bin/mmc4.xlsx

# We use readxl package for reading xlsx files, so readxl must be installed beforehand.

###Preparing Meta information first:
##Importing Meta Info into R:
RespondInfo <- readxl::read_excel(path = "./data-raw/GDSC/GDSC-raw/v17.3_fitted_dose_response.xlsx", sheet = 1)
DrugInfo <- as.data.frame(readxl::read_excel(path = "./data-raw/GDSC/GDSC-raw/Screened_Compounds.xlsx", sheet = 1))
CelllineInfo <- as.data.frame(readxl::read_excel(path = "./data-raw/GDSC/GDSC-raw/Cell_Lines_Details.xlsx", sheet = "Cell line details"))
TissueInfo <- as.data.frame(readxl::read_excel(path = "./data-raw/GDSC/GDSC-raw/Cell_Lines_Details.xlsx", sheet = "COSMIC tissue classification"))

##There are drugs that have two matching IDs, for now we add a '(2)' suffix to the drug name that already had another ID:
DrugInfo$DRUG_NAME[duplicated(DrugInfo$DRUG_NAME)] <- paste0(DrugInfo$DRUG_NAME[duplicated(DrugInfo$DRUG_NAME)],"(2)")
RespondInfo$DRUG_NAME <- DrugInfo$DRUG_NAME[match(RespondInfo$DRUG_ID,DrugInfo$DRUG_ID)]

## We rearrange drug responses (IC50, AUC, RMSE, ...) into a matrix, with cell lines as rows and drugs as columns:
#Initializing drug response matrices:
ALL_CELL_LINE_NAME <- unique(RespondInfo$CELL_LINE_NAME)
ALL_DRUG_NAME <- unique(RespondInfo$DRUG_NAME)
LN_IC50 <- AUC <- RMSE <- Z_SCORE <-
  MAX_CONC_MICROMOLAR <- MIN_CONC_MICROMOLAR <- matrix(NA, nrow = length(ALL_CELL_LINE_NAME),
                                                       ncol = length(ALL_DRUG_NAME),
                                                       dimnames = list(ALL_CELL_LINE_NAME,ALL_DRUG_NAME))

#We use progress package for making a more sophisticated progress bar for our loops:
require(progress)
pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                       total = length(LN_IC50), clear = FALSE)

#These loops fill drug response matrices we initialized above:
for(CellLine in ALL_CELL_LINE_NAME){
  for(DrUg in ALL_DRUG_NAME){
    bullseyeIndex <- RespondInfo$DRUG_NAME %in% DrUg & RespondInfo$CELL_LINE_NAME %in% CellLine
    if(sum(bullseyeIndex)>1) stop("Something's wrong: couldn't match properly!!")
    if(sum(bullseyeIndex)==1){
      LN_IC50[CellLine, DrUg] <- RespondInfo$LN_IC50[bullseyeIndex]
      AUC[CellLine, DrUg] <- RespondInfo$AUC[bullseyeIndex]
      RMSE[CellLine, DrUg] <- RespondInfo$RMSE[bullseyeIndex]
      Z_SCORE[CellLine, DrUg] <- RespondInfo$Z_SCORE[bullseyeIndex]
      MAX_CONC_MICROMOLAR[CellLine, DrUg] <- RespondInfo$MAX_CONC_MICROMOLAR[bullseyeIndex]
      MIN_CONC_MICROMOLAR[CellLine, DrUg] <- RespondInfo$MIN_CONC_MICROMOLAR[bullseyeIndex]
    }
    pb$tick()
  }
}



###Preparing Gene Expression Data:
#Importing:
GeneExpression <- read.delim("./data-raw/GDSC/GDSC-raw/sanger1018_brainarray_ensemblgene_rma.txt.gz", row.names = 1)

#Converting colnames from COSMIC IDs to Cell line names:
colnames(GeneExpression) <-
  CelllineInfo$`Sample Name`[match(as.integer(substr(x = colnames(GeneExpression),
    start = 2,
    stop = nchar(colnames(GeneExpression))
  )), CelllineInfo$`COSMIC identifier`)]

#Converting rownames from ensemble IDs to Entrez:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTabelle <- getBM(attributes = c("ensembl_gene_id","entrezgene"),
      filters = "ensembl_gene_id", values = rownames(GeneExpression), mart = humaRt)
GeneExpression <- as.matrix(GeneExpression) #Cause it's a dataframe and dataframes have problem with duplicated row names!
rownames(GeneExpression) <- ConvTabelle$entrezgene[match(rownames(GeneExpression),ConvTabelle$ensembl_gene_id)]
GeneExpression <- GeneExpression[!is.na(rownames(GeneExpression)),]

###Preparing CNA, Mutation and Methylation:
#Importing Iorio et. al. 2016 table:
IorioBinary <- readxl::read_excel(path = "./data-raw/GDSC/GDSC-raw/mmc4.xlsx", sheet = "TableS3B-GlobalBEM", skip = 6)
colnames(IorioBinary) <- as.character(readxl::read_excel(path = "./data-raw/GDSC/GDSC-raw/mmc4.xlsx",
                                                         sheet = "TableS3B-GlobalBEM", skip = 3,n_max = 1, col_names = F))
IorioBinaryMat <- as.matrix(IorioBinary[,2:1002])
rownames(IorioBinaryMat) <- IorioBinary$`Sample identifier`

##Based on rownames in IorioBinaryMat, we want to extract CNA, Mutation and Methylation into different matrices:
splittedRowNames <- strsplit(x = rownames(IorioBinaryMat), split = "_")
#starting with CNA:
CNAsIndx <- !sapply(splittedRowNames,length) == 2
CNAs <- IorioBinaryMat[CNAsIndx,]
CNAsWithMatchedGenes <- CNAs[grep(rownames(CNAs), pattern = " ", fixed = T),]
CNAsWithMatchedGenesGain <- CNAsWithMatchedGenes[substr(x = rownames(CNAsWithMatchedGenes),
                                                        start = 1, stop = 4)=="gain",]
CNAsWithMatchedGenesLoss <- CNAsWithMatchedGenes[substr(x = rownames(CNAsWithMatchedGenes),
                                                        start = 1, stop = 4)=="loss",]

#Have to convert to matrix cause there're duplications in rownames:
CNAsWithMatchedGenesGain <- as.matrix(CNAsWithMatchedGenesGain)
CNAsWithMatchedGenesLoss <- as.matrix(CNAsWithMatchedGenesLoss)

#Extracting gene names that are in parenthesis:
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

#Some rows have more than one gene per row, making a new matrix with single genes in rows:
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

#Converting gene symbols to Entrez in row names:
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


#Now to mutation and methylation:
#Extracting Mutation matrix:
MutsAndMethyls <- IorioBinaryMat[!CNAsIndx,]
MutsAndMethylsRowDecomposed <- strsplit(x = rownames(MutsAndMethyls),
                                               split = "_", fixed = T)
Muts <- MutsAndMethyls[sapply(MutsAndMethylsRowDecomposed, function(x) x[2]=="mut"),]
rownames(Muts) <- sapply(strsplit(x = rownames(Muts),split = "_", fixed = T), function(x) x[1])

#Converting gene symbols to Entrez in row names:
ConvTabelle <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                     filters = "hgnc_symbol",
                     values = rownames(Muts), mart = humaRt)
#Duplication in row names so:
Muts <- as.matrix(Muts)

rownames(Muts) <- ConvTabelle$entrezgene[match(rownames(Muts), ConvTabelle$hgnc_symbol)]
Muts <- Muts[!is.na(rownames(Muts)),]

#Extracting Methylation matrix:
Methyls <- MutsAndMethyls[sapply(MutsAndMethylsRowDecomposed, function(x) x[2]=="HypMET"),]
#Extracting gene names hidden inside parenthesis in the rownames:
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

#Converting gene symbols to Entrez in row names:
ConvTabelle <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                     filters = "hgnc_symbol",
                     values = rownames(MethylsWithMatchedGenes), mart = humaRt)
rownames(MethylsWithMatchedGenes) <- ConvTabelle$entrezgene[match(rownames(MethylsWithMatchedGenes), ConvTabelle$hgnc_symbol)]
MethylsWithMatchedGenes <- MethylsWithMatchedGenes[!is.na(rownames(MethylsWithMatchedGenes)),]

###Preparing 'ResponseTypes' slot:
##ResponseTypes is a data frame containing variable names and description of availble ForeseeCell drug response variables:
Response <- c("LN_IC50","AUC","RMSE","Z_SCORE","MAX_CONC_MICROMOLAR","MIN_CONC_MICROMOLAR")
Response <- cbind(Response,c("Natural log of the fitted IC50",
                             "Area Under the Curve for the fitted model.
                             Presented as a fraction of the total area between
                             the highest and lowest screening concentration.",
                             "Root Mean Squared Error, a measurement of how well
                             the modelled curve fits the data points. Curves
                             with RMSE > 0.3 are excluded prior to release as
                             part of quality control.",
                             "Z score of the LN_IC50 (x) comparing it to the
                             mean (μ) and standard deviation ( σ 2 ) of
                             the LN_IC50 values for the drug in question
                             over all cell lines treated.",
                             "Maximum micromolar screening concentration of the drug",
                             "Minimum micromolar screening concentration of the drug"))
Response <- as.data.frame(Response);colnames(Response) <- c("Name","Description")

###Preparing 'InputTypes' slot:
##InputTypes is a data frame containing variable names (and description) of availble ForeseeCell input data variables:
InputTypes <- c("GeneExpression","CNVGain","CNVLoss","Mutation","Methylation")
InputTypes <- cbind(InputTypes,c("RMA-normalized gene expression measured by DNA array",
                             "Gaining of copy number variation, a binary matrix",
                             "Losing of copy number variation, a binary matrix",
                             "Mutation, a binary matrix",
                             "Methylation, a binary matrix"))
InputTypes <- as.data.frame(InputTypes);colnames(InputTypes) <- c("Name","Description")

###Making and saving the Foresee object:
GDSC <- list()
class(GDSC) <- "ForeseeCell"
GDSC[["GeneExpression"]] <- GeneExpression
GDSC[["CNVGain"]] <- CNAsWithMatchedGenesGainSignleGenes
GDSC[["CNVLoss"]] <- CNAsWithMatchedGenesLossSignleGenes
GDSC[["Mutation"]] <- Muts
GDSC[["Methylation"]] <- MethylsWithMatchedGenes
GDSC[["LN_IC50"]] <- LN_IC50
GDSC[["AUC"]] <- AUC
GDSC[["RMSE"]] <- RMSE
GDSC[["Z_SCORE"]] <- Z_SCORE
GDSC[["MAX_CONC_MICROMOLAR"]] <- MAX_CONC_MICROMOLAR
GDSC[["MIN_CONC_MICROMOLAR"]] <- MIN_CONC_MICROMOLAR
GDSC[["DrugInfo"]] <- DrugInfo
GDSC[["CelllineInfo"]] <- CelllineInfo
GDSC[["TissueInfo"]] <- TissueInfo
GDSC[["ResponseTypes"]] <- Response
GDSC[["InputTypes"]] <- InputTypes


devtools::use_data(GDSC, overwrite = T)




