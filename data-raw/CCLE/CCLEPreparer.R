## Script that will (download,) process and add CCLE to this package
# by esfahani@aices.rwth-aachen.de, XX.2.2018

## We downloaded CCLE from Broad institute's portal https://portals.broadinstitute.org/ccle, since the portal
#requires registration and authentication, we don't directly download the raw files here but assume user
#downloaded the relevant versions of data files and put them in CCLE-raw.

## List and link of downloaded data, all downloaded on 5.2018 from https://portals.broadinstitute.org/ccle/data:
#1- RMA normalised gene expression data from DNA array, named CCLE_Expression_Entrez_2012-10-18.res, last updated on October 18th 2012,
#with the link https://data.broadinstitute.org/ccle_legacy_data/mRNA_expression/CCLE_Expression_Entrez_2012-10-18.res
#2- CCLE RNAseq gene expression data (RPKM), named CCLE_DepMap_18Q1_RNAseq_RPKM_20180214.gct, last updated on February 14th 2018,
#with the link https://data.broadinstitute.org/ccle/CCLE_DepMap_18Q1_RNAseq_RPKM_20180214.gct
#3- Binary Calls for Copy Number and Mutation Data, named CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct, last updated on February 29th 2016,
#with the link https://data.broadinstitute.org/ccle_legacy_data/binary_calls_for_copy_number_and_mutation_data/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct
#(We used) Legacy Data because the latest release lacked Copy number variation data)
#4- Reverse Phase Protein Array data, named ccle2maf_081117.txt, last updated on January 23rd 2018,
#with the link https://data.broadinstitute.org/ccle/CCLE_RPPA_20180123.csv
#5- Reverse Phase Protein Array antibody information, named ccle2maf_081117-2.txt, last updated on January 24rd 2018,
#with the link https://data.broadinstitute.org/ccle/CCLE_RPPA_Ab_info_20180123.csv
#6- Cell Line Annotations, named CCLE_sample_info_file_2012-10-18.txt, last updated on October 18th 2012,
#with the link https://data.broadinstitute.org/ccle_legacy_data/cell_line_annotations/CCLE_sample_info_file_2012-10-18.txt
#7- List of the 24 drugs profiled across 504 CCLE lines, named CCLE_NP24.2009_profiling_2012.02.20.csv, last updated on April 17th 2012,
#with the link https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_profiling_2012.02.20.csv
#8- Pharmacologic profiles for 24 anticancer drugs across 504 CCLE lines, named CCLE_NP24.2009_Drug_data_2015.02.24.csv, last updated on February 24th 2015,
#with the link https://data.broadinstitute.org/ccle_legacy_data/pharmacological_profiling/CCLE_NP24.2009_Drug_data_2015.02.24.csv


###Preparing Gene Expression Data, starting with DNA array data:
#Importing and cleaning up:
GEX <- read.delim("./data-raw/CCLE/CCLE-raw/CCLE_Expression_Entrez_2012-10-18.res", header = T, row.names = NULL, check.names = F)
colnames(GEX)[seq(from = 3, to = 2073, by = 2)] <- colnames(GEX)[seq(from = 4, to = 2074, by = 2)]
GEX <- GEX[,-seq(from = 4, to = 2074, by = 2)]
GexRowNames <- GEX[-1,1]
GEX <- GEX[,-(1:2)]
GEX <- as.matrix(GEX)
GEX <- GEX[-1,]
rownames(GEX) <- GexRowNames

#Converting gene symbols to Entrez IDs:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"), filters = "hgnc_symbol",values = rownames(GEX), mart = humaRt)
rownames(GEX) <- ConvTableSym2Entrez$entrezgene[match(rownames(GEX),ConvTableSym2Entrez$hgnc_symbol)]
GEX <- GEX[!is.na(rownames(GEX)),]


###Preparing Gene Expression Data, continuing with RNA-seq data:
#Importing and cleaning up:
GEXRNAseq <- read.delim("./data-raw/CCLE/CCLE-raw/CCLE_DepMap_18Q2_RNAseq_RPKM_20180502.gct",
                        skip = 2, header = T, check.names = F)
GEXRNAseqMat <- as.matrix(GEXRNAseq[,-c(1:2)])
rownames(GEXRNAseqMat) <- GEXRNAseq$Description

#Converting gene symbols to Entrez IDs:
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"), filters = "hgnc_symbol",values = rownames(GEXRNAseqMat), mart = humaRt)
rownames(GEXRNAseqMat) <- ConvTableSym2Entrez$entrezgene[match(rownames(GEXRNAseqMat),ConvTableSym2Entrez$hgnc_symbol)]
GEXRNAseqMat <- GEXRNAseqMat[!is.na(rownames(GEXRNAseqMat)),]


###Preparing CNA and Mutation:
#Importing:
MutAndCNABinary <- read.delim("./data-raw/CCLE/CCLE-raw/CCLE_MUT_CNA_AMP_DEL_binary_Revealer.gct",
                              skip = 2, header = T)

#Based on rownames, we want to extract CNA and Mutation into different matrices:
splittedNames <- strsplit(as.character(MutAndCNABinary$Name), split = "_")
#Three splits fail cause there's _ in the gene names, have to correct the the ruined gene names:
failedOnes <- sapply(splittedNames,length)==3
splittedNames[failedOnes] <- lapply(splittedNames[failedOnes], function(x) c(paste(x[1],x[2], sep = "_"),x[3]))

#starting with Mutation matrix:
Muts <- as.matrix(MutAndCNABinary[sapply(splittedNames,function(x) x[2]) == "MUT",-c(1:2)])
rownames(Muts) <- as.character(MutAndCNABinary[sapply(splittedNames,function(x) x[2]) == "MUT",1])
rownames(Muts) <- sapply(strsplit(x = rownames(Muts), split = "_"), function(x) x[1])

#Converting gene symbols to Entrez IDs:
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                             filters = "hgnc_symbol",values = rownames(Muts), mart = humaRt)
rownames(Muts) <- ConvTableSym2Entrez$entrezgene[match(rownames(Muts),ConvTableSym2Entrez$hgnc_symbol)]
Muts <- Muts[!is.na(rownames(Muts)),]

#continuing with Amplification in copy number alteration matrix:
CNAAmp <- as.matrix(MutAndCNABinary[sapply(splittedNames,function(x) x[2]) == "AMP",-c(1:2)])
rownames(CNAAmp) <- as.character(MutAndCNABinary[sapply(splittedNames,function(x) x[2]) == "AMP",1])
rownames(CNAAmp) <- sapply(strsplit(x = rownames(CNAAmp), split = "_"), function(x) x[1])

#Converting gene symbols to Entrez IDs:
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                             filters = "hgnc_symbol",values = rownames(CNAAmp), mart = humaRt)
rownames(CNAAmp) <- ConvTableSym2Entrez$entrezgene[match(rownames(CNAAmp),ConvTableSym2Entrez$hgnc_symbol)]
CNAAmp <- CNAAmp[!is.na(rownames(CNAAmp)),]

#last but not least: Deletions in copy number alteration matrix:
CNADel <- as.matrix(MutAndCNABinary[sapply(splittedNames,function(x) x[2]) == "DEL",-c(1:2)])
rownames(CNADel) <- as.character(MutAndCNABinary[sapply(splittedNames,function(x) x[2]) == "DEL",1])
rownames(CNADel) <- sapply(strsplit(x = rownames(CNADel), split = "_"), function(x) x[1])

#Converting gene symbols to Entrez IDs:
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                             filters = "hgnc_symbol",values = rownames(CNADel), mart = humaRt)
rownames(CNADel) <- ConvTableSym2Entrez$entrezgene[match(rownames(CNADel),ConvTableSym2Entrez$hgnc_symbol)]
CNADel <- CNADel[!is.na(rownames(CNADel)),]


###Preparing Reverse Phase Protein Array data:
#Importing:
RPPA <- t(as.matrix(read.csv("./data-raw/CCLE/CCLE-raw/ccle2maf_081117.txt",
                           header = T, row.names = 1, check.names = F)))
RPPAAntiInfo <- read.csv("./data-raw/CCLE/CCLE-raw/ccle2maf_081117-2.txt", header = T, check.names = F)

#Replacing anibody with their corrisponding target gene name(s):
rownames(RPPA) <- RPPAAntiInfo$Target_Genes

#Some rows have more than one gene in the rownames, have to fix it into one gene per row:
ROW <- 1
while(ROW <= nrow(RPPA)){
  ROWinHAND <- strsplit(rownames(RPPA)[ROW], split = " ", fixed = T)[[1]]
  if(length(ROWinHAND) > 1){
    RPPA <- RPPA[c(1:ROW,rep(ROW,(length(ROWinHAND)-2)),ROW:nrow(RPPA)),]
    rownames(RPPA)[ROW:(ROW+length(ROWinHAND)-1)] <- ROWinHAND
  }
  ROW <- ROW + 1
}

#There are duplication in measured genes, some anibodies have same targets,
#we will get the mean of the duplicated genes:
RPPAAggregated <- aggregate(x = RPPA, by = list(rownames(RPPA)), FUN = mean)
RPPA <- data.matrix(RPPAAggregated[,-1])
rownames(RPPA) <- RPPAAggregated$Group.1

#Converting gene symbols to Entrez IDs:
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                             filters = "hgnc_symbol",values = rownames(RPPA), mart = humaRt)
rownames(RPPA) <- ConvTableSym2Entrez$entrezgene[match(rownames(RPPA),ConvTableSym2Entrez$hgnc_symbol)]
RPPA <- RPPA[!is.na(rownames(RPPA)),]


###Preparing Meta information:
#Importing sample information:
SInfo <- read.delim("./data-raw/CCLE/CCLE-raw/CCLE_sample_info_file_2012-10-18.txt")
names(SInfo)[5] <- "Site" #For being compatible with SampleSelector function in FORESEE
#Drug Info:
DInfo <- read.csv("./data-raw/CCLE/CCLE-raw/CCLE_NP24.2009_profiling_2012.02.20.csv")
names(DInfo)[c(1,3)] <- c("DRUG_NAME","TARGET") #For being compatible with FeatureSelector function in FORESEE

###Response data curation:
ResponseData <- read.csv("./data-raw/CCLE/CCLE-raw/CCLE_NP24.2009_Drug_data_2015.02.24.csv")

## We rearrange drug responses (IC50, AUC, RMSE, ...) into a matrix, with cell lines as rows and drugs as columns:
#Initializing drug response matrices:
IC50 <- EC50 <- ActArea <- Amax <- matrix(NA, nrow = nlevels(ResponseData$CCLE.Cell.Line.Name),
                              ncol = nlevels(ResponseData$Compound),
                              dimnames = list(levels(ResponseData$CCLE.Cell.Line.Name),levels(ResponseData$Compound)))

#We use progress package for making a more sophisticated progress bar for our loops:
require(progress)
pb <- progress_bar$new(format = "(:spin) [:bar] :percent eta: :eta",
                       total = length(IC50), clear = FALSE)

#These loops fill drug response matrices we initialized above:
for(CellLine in levels(ResponseData$CCLE.Cell.Line.Name)){
  for(DrUg in levels(ResponseData$Compound)){
    bullseyeIndex <- ResponseData$Compound %in% DrUg & ResponseData$CCLE.Cell.Line.Name %in% CellLine
    if(sum(bullseyeIndex)>1) stop("Something's wrong: couldn't match properly!!")
    if(sum(bullseyeIndex)==1){
      IC50[CellLine, DrUg] <- ResponseData$"IC50..uM."[bullseyeIndex]
      EC50[CellLine, DrUg] <- ResponseData$"EC50..uM."[bullseyeIndex]
      ActArea[CellLine, DrUg] <- ResponseData$ActArea[bullseyeIndex]
      Amax[CellLine, DrUg] <- ResponseData$Amax[bullseyeIndex]
    }
    pb$tick()
  }
}


###Preparing 'ResponseTypes' slot:
##ResponseTypes is a data frame containing variable names and description of availble ForeseeCell drug response variables:
Response <- c("IC50","EC50","ActArea","Amax")
Response <- cbind(Response,c("Concentration at which the drug response
                             reached an absolute inhibition of 50%",
                             "Concentration at half-maximal activity of the compound",
                             "Area between the pharmacologic dose-response
                             curve and the zero inhibition level",
                             "Inhibition observed at the maximal drug concentration"))
Response <- as.data.frame(Response);colnames(Response) <- c("Name","Description")

###Preparing 'InputTypes' slot:
##InputTypes is a data frame containing variable names (and description) of availble ForeseeCell input data variables:
InputTypes <- c("GeneExpression","GeneExpressionRNAseq","ProteinExpression","CNVGain","CNVLoss","Mutation")
InputTypes <- cbind(InputTypes,c("RMA-normalized gene expression measured by DNA array",
                                 "Gene expression measured by sequencing in RPKM",
                                 "Protein expression measured by Reverse Phase Protein Array",
                                 "Gaining of copy number variation, a binary matrix",
                                 "Losing of copy number variation, a binary matrix",
                                 "Mutation, a binary matrix"))
InputTypes <- as.data.frame(InputTypes);colnames(InputTypes) <- c("Name","Description")


###Making and saving the Foresee object:
CCLE <- list()
class(CCLE) <- "ForeseeCell"
CCLE[["GeneExpression"]] <- GEX
CCLE[["GeneExpressionRNAseq"]] <- GEXRNAseqMat
CCLE[["ProteinExpression"]]<- RPPA
CCLE[["CNVGain"]] <- CNAAmp
CCLE[["CNVLoss"]] <- CNADel
CCLE[["Mutation"]] <- Muts
CCLE[["IC50"]] <- IC50
CCLE[["EC50"]] <- EC50
CCLE[["ActArea"]] <- ActArea
CCLE[["Amax"]] <- Amax
CCLE[["DrugInfo"]] <- DInfo
CCLE[["TissueInfo"]] <- SInfo
CCLE[["ResponseTypes"]] <- Response
CCLE[["InputTypes"]] <- InputTypes

devtools::use_data(CCLE, overwrite = T)

