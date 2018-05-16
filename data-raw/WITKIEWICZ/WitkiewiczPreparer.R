## Script that will (download,) process and add Witkiewicz Xenograft dataset to this package
## data published in https://www.ncbi.nlm.nih.gov/pubmed/27498862
# by esfahani@aices.rwth-aachen.de, 19.3.2018


## We downloaded Witkiewicz Xenograft dataset as supplement files of the paper
#"Integrated Patient-Derived Models Delineate Individualized Therapeutic Vulnerabilities of Pancreatic Cancer"
#by Witkiewicz et. al. published in 2016 in Cell Rep. accessible via the link https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5287055/
# Except the gene expression dataset, which was downloaded from GEO via https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84023

##Raw data can be downloaded using the following lines:
download.file(url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE84023&format=file&file=GSE84023%5Fseq%5FProcessedData%2Etxt%2Egz",
              destfile = "./data-raw/WITKIEWICZ/WITKIEWICZ-raw/GSE84023_seq_ProcessedData.txt.gz")
download.file(url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5287055/bin/NIHMS803574-supplement-2.xlsx",
              destfile = "./data-raw/WITKIEWICZ/WITKIEWICZ-raw/NIHMS803574-supplement-2.xlsx")
download.file(url = "https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5287055/bin/NIHMS803574-supplement-3.xlsx",
              destfile = "./data-raw/WITKIEWICZ/WITKIEWICZ-raw/NIHMS803574-supplement-3.xlsx")


###Preparing Gene Expression Data based on RNA-seq: (Pipeline based on info of GSE84023: Illumina Casava1.7 software used for basecalling.
#Sequenced reads were trimmed for adaptor sequence, mapped to hg19 genome using bowtie TopHatCounts per gene was obtained
#using HTseq counts and normalized using edgeR package in RGenome_build: hg19, files_format_and_content: tab-delimited text file
#include matrix of normalized log counts per million for each sample)
#Importing and cleaning up:
RNAseq <- read.table("./data-raw/WITKIEWICZ/WITKIEWICZ-raw/GSE84023_seq_ProcessedData.txt.gz", header = T)
RNAseqMat <- as.matrix(RNAseq[,-1])
rownames(RNAseqMat) <- RNAseq$Gene

#Feature identifiers are in gene symbols, we convert them to Entrez IDs:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTabelle <- getBM(attributes = c("hgnc_symbol","entrezgene"),
                     filters = "hgnc_symbol", values = rownames(RNAseqMat), mart = humaRt)
rownames(RNAseqMat) <- ConvTabelle$entrezgene[match(rownames(RNAseqMat),ConvTabelle$hgnc_symbol)]
RNAseqMat <- RNAseqMat[!is.na(rownames(RNAseqMat)),] # No need to log2; as said above, values are normalized log counts per million.


###Response data curation:
#Importing and cleaning up:
##There are single and double drug treatment in this dataset; we curate both into two different metrics, with model in rows and treatment in columns:
#Single response data:
require(openxlsx)
AUC_single <- read.xlsx(xlsxFile = "./data-raw/WITKIEWICZ/WITKIEWICZ-raw/NIHMS803574-supplement-3.xlsx", sheet = 3, rowNames = T)
AUC_single_Mat <- as.matrix(AUC_single)
AUC_single_Mat <- t(AUC_single_Mat) #We want the models to be in rows and drugs in columns.

#Double response data:
require(openxlsx)
AUC_combo <- read.xlsx(xlsxFile = "./data-raw/WITKIEWICZ/WITKIEWICZ-raw/NIHMS803574-supplement-2.xlsx", sheet = 2, rowNames = T)
AUC_combo_Mat <- as.matrix(AUC_combo)
AUC_combo_Mat <- t(AUC_combo_Mat) #We want the models to be in rows and drugs in columns.

###Preparing Meta information:
#Some (maybe?) drug info:
DrugInfo <- read.xlsx(xlsxFile = "./data-raw/WITKIEWICZ/WITKIEWICZ-raw/NIHMS803574-supplement-3.xlsx", sheet = 1)


###Preparing 'ResponseTypes' slot:
##ResponseTypes is a data frame containing variable names and description of availble ForeseeCell drug response variables:
Response <- c("AUC","AUCCombo")
Response <- cbind(Response,c("patient-derived cell line models were screened in total at
                             100 nM-1 μM dose range, and area under the curve (AUC)
                             was calculated per drug per cell line (range 0.08 – 4.95), For single treatment",
                             "patient-derived cell line models were screened in total at
                             100 nM-1 μM dose range, and area under the curve (AUC)
                             was calculated per drug per cell line (range 0.08 – 4.95), For combinatory treatment"))
Response <- as.data.frame(Response);colnames(Response) <- c("Name","Description")

###Preparing 'InputTypes' slot:
##InputTypes is a data frame containing variable names (and description) of availble ForeseeCell input data variables:
InputTypes <- c("GeneExpression")
InputTypes <- cbind(InputTypes,c("Matrix of normalized log counts per million for each sample,
                                 detailed pipeline explained in https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE84023"))
InputTypes <- as.data.frame(InputTypes);colnames(InputTypes) <- c("Name","Description")


###Making and saving the Foresee object:
WITKIEWICZ <- list()
class(WITKIEWICZ) <- "ForeseeCell" ##We assign Xenografts as Cell objects!
WITKIEWICZ[["GeneExpression"]] <- RNAseqMat
WITKIEWICZ[["AUC"]] <- AUC_single_Mat
WITKIEWICZ[["AUCCombo"]] <- AUC_combo_Mat

WITKIEWICZ[["DrugInfo"]] <- DrugInfo

WITKIEWICZ[["ResponseTypes"]] <- Response
WITKIEWICZ[["InputTypes"]] <- InputTypes

devtools::use_data(WITKIEWICZ, overwrite = T)
