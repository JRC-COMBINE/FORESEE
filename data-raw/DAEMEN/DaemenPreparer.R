## Script that will (download,) process and add Daemen
# dataset (based on the paper "Modeling precision treatment of breast cancer" by Daemen et. al.)
# to this package.
# by esfahani@aices.rwth-aachen.de, 5.2018

## We aimed to replicate the results of Daemen et. al. paper, so
#we downloaded the data provided for "Modeling precision treatment of breast cancer."
#on the paper, with the link: https://www.synapse.org/#!Synapse:syn2179898
#Since the portal requires registration and authentication, we don't directly download the
# raw files here but assume user downloaded the relevant data files and put them in DAEMEN-raw.

## List and link of downloaded data, all downloaded on 5.2018 from https://www.synapse.org/#!Synapse:syn2184886 :
#1- RMA normalised gene expression data from DNA array, named Neve_AffyRMA_genelevel_maxvar_stringent.csv, last updated on August 30th 2013,
#with the link https://www.synapse.org/Portal/filehandleassociation?associatedObjectId=syn2184894&associatedObjectType=FileEntity&fileHandleId=150628
#2- Gene expression counts based on sequencing data, named breastRNAseq_genelevel_stringent.txt, last updated on August 30th 2013,
#with the link https://www.synapse.org/Portal/filehandleassociation?associatedObjectId=syn2184895&associatedObjectType=FileEntity&fileHandleId=150630
#3- Methylation data, named Methylation_stringent.csv, last updated on August 30th 2013,
#with the link https://www.synapse.org/Portal/filehandleassociation?associatedObjectId=syn2184893&associatedObjectType=FileEntity&fileHandleId=150626
#4- Methylation annotation data, named Methylation_annotation_stringent.csv, last updated on August 30th 2013,
#with the link https://www.synapse.org/Portal/filehandleassociation?associatedObjectId=syn2184892&associatedObjectType=FileEntity&fileHandleId=150624
#5- SNP data, named SNP6_genelevel_stringent_std0.7.csv, last updated on August 30th 2013,
#with the link https://www.synapse.org/Portal/filehandleassociation?associatedObjectId=syn2184911&associatedObjectType=FileEntity&fileHandleId=150660
#(The next file is downloaded as a supplement file of the original paper, accesible via https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3937590/#B41, Additional file 1)
#6- GI50 drug response, named gb-2013-14-10-r110-S1.xlsx, last updated on ---
#with the link https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3937590/bin/gb-2013-14-10-r110-S1.xlsx


###Preparing Gene Expression Data, starting with DNA array data:
#Importing and cleaning up:
GEX <- read.csv("./data-raw/DAEMEN/DAEMEN-raw/Neve_AffyRMA_genelevel_maxvar_stringent.csv", check.names = F)
GEXMat <- as.matrix(GEX[,-1])
rownames(GEXMat) <- GEX$SampleNames

#Converting gene symbols to Entrez IDs:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"), filters = "hgnc_symbol",
                             values = rownames(GEXMat), mart = humaRt)
rownames(GEXMat) <- ConvTableSym2Entrez$entrezgene[match(rownames(GEXMat),ConvTableSym2Entrez$hgnc_symbol)]
GEXMat <- GEXMat[!is.na(rownames(GEXMat)),]


###Preparing Gene Expression Data, Continuing with RNA-seq data:
#Importing and cleaning up:
GEXRNAseq <- read.delim("./data-raw/DAEMEN/DAEMEN-raw/breastRNAseq_genelevel_stringent.txt", check.names = F)
GEXRNAseqMat <- as.matrix(GEXRNAseq[,-c(1:4)])
rownames(GEXRNAseqMat) <- GEXRNAseq$Seq_Name

#Converting gene symbols to Entrez IDs:
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"), filters = "hgnc_symbol",
                             values = rownames(GEXRNAseqMat), mart = humaRt)
rownames(GEXRNAseqMat) <- ConvTableSym2Entrez$entrezgene[match(rownames(GEXRNAseqMat),
                                                               ConvTableSym2Entrez$hgnc_symbol)]
GEXRNAseqMat <- GEXRNAseqMat[!is.na(rownames(GEXRNAseqMat)),]

##We move the seq values into logarithmic (base 2) scale so that we'd be having semi-normaly distributed values, necessary for linear models
#We cut all lower that 1 values, and replace them with 1.
GEXRNAseqMat[GEXRNAseqMat < 1] <- 1
GEXRNAseqMat <- log2(GEXRNAseqMat)

###Preparing methylation data:
#Importing and cleaning up:
Meth <- read.csv("./data-raw/DAEMEN/DAEMEN-raw/Methylation_stringent.csv", check.names = F)
MethNameConversion <- read.csv("./data-raw/DAEMEN/DAEMEN-raw/Methylation_annotation_stringent.csv",
                               check.names = F)
MethMat <- as.matrix(Meth[,-1])

#Assigning gene symbols as rownames:
rownames(MethMat) <- MethNameConversion$Symbol #Because the row orders of Meth and MethNameConversion are the same (check all(Meth$SampleNames == MethNameConversion$Name))

#Converting gene symbols to Entrez IDs:
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"), filters = "hgnc_symbol",
                             values = rownames(MethMat), mart = humaRt)
rownames(MethMat) <- ConvTableSym2Entrez$entrezgene[match(rownames(MethMat),
                                                               ConvTableSym2Entrez$hgnc_symbol)]
MethMat <- MethMat[!is.na(rownames(MethMat)),]


###Preparing SNP data:
#Importing and cleaning up:
SNP <- read.csv("./data-raw/DAEMEN/DAEMEN-raw/SNP6_genelevel_stringent_std0.7.csv", check.names = F)
SNPMat <- as.matrix(SNP[,-(1:5)])
#Already comes with EntrezIDs so no conversion is required:
rownames(SNPMat) <- SNP$GeneID

###Response data curation:
#Importing and cleaning up(you need readxl package for importing the excel file):
ResponseData <- readxl::read_excel("./data-raw/DAEMEN/DAEMEN-raw/gb-2013-14-10-r110-S1.xlsx", skip = 2)

ResponseData <- read.csv2("./data-raw/DAEMEN/DAEMEN-raw/13059_2013_3164_MOESM1_ESM.csv",
                         skip = 2, check.names = F)
ResponseDataMat <- data.matrix(ResponseData[,-c(1:11)])
rownames(ResponseDataMat) <- ResponseData$"Cell line"

##Exporting (short) TissueInfo from Response Data:
TissueInfo <- as.data.frame(ResponseData[,1:3])
names(TissueInfo)[2] <- "Site" #For being compatible with SampleSelector function in FORESEE


###Preparing 'ResponseTypes' slot:
##ResponseTypes is a data frame containing variable names and description of availble ForeseeCell drug response variables:
Response <- c("GI50")
Response <- cbind(Response,c("-log10 of GI50 (Concentration required to inhibit growth by 50%)"))
Response <- as.data.frame(Response);colnames(Response) <- c("Name","Description")

###Preparing 'InputTypes' slot:
##InputTypes is a data frame containing variable names (and description) of availble ForeseeCell input data variables:
InputTypes <- c("GeneExpression","GeneExpressionRNAseq","Methylation","SNP6")
InputTypes <- cbind(InputTypes,c("RMA-normalized gene expression measured by DNA array",
                                 "Gene expression measured by sequencing in counts",
                                 "Methylation Data",
                                 "SNP Data"))
InputTypes <- as.data.frame(InputTypes);colnames(InputTypes) <- c("Name","Description")


###Making and saving the Foresee object:
DAEMEN <- list()
class(DAEMEN) <- "ForeseeCell"
DAEMEN[["GeneExpression"]] <- GEXMat
DAEMEN[["GeneExpressionRNAseq"]] <- GEXRNAseqMat
DAEMEN[["Methylation"]] <- MethMat
DAEMEN[["SNP6"]] <- SNPMat
DAEMEN[["GI50"]] <- ResponseDataMat
DAEMEN[["TissueInfo"]] <- TissueInfo
DAEMEN[["ResponseTypes"]] <- Response
DAEMEN[["InputTypes"]] <- InputTypes


devtools::use_data(DAEMEN, overwrite = T)



