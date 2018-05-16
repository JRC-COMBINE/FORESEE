##This script prepares and saves internal data required by FORESEE,
# written by esfahani@aices.rwth-aachen.de, 11.5.2018

###Preparing internal data for FeatureSelector.ontology
## In FeatureSelector.ontology we need a conversion table between GO (gene ontology) and gene symbols (because usually target of drugs
#in ForeseeCell objects are save as gene symbols), we also need conversion table between gene symbols and Entrez IDs (because
#genes in expression matrix in ForeseeCell objects are save as Entrez IDs)
##We used biomaRt (https://bioconductor.org/packages/release/bioc/html/biomaRt.html) to make these two conversion tables:
require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"),values = "*", mart = humaRt)
ConvTableGo2Sym <- getBM(attributes = c("go_id","hgnc_symbol"),values = "*", mart = humaRt)

###Preparing internal data for FeatureSelector.pathway:
## In FeatureSelector.pathway we need a conversion table from pathway to gene Entrez IDs and vice versa.
##We used reactome database (https://bioconductor.org/packages/release/data/annotation/html/reactome.db.html) to make these two conversion tables:
require(reactome.db)
PathID2Entrez <- as.list(reactomePATHID2EXTID)
Entrez2PathID <- as.list(reactomeEXTID2PATHID)

###Preparing internal data for FeatureSelector.landmarkgenes:
## In FeatureSelector.landmarkgenes we needed a list of landmark genes in Entrez ID, which we obtained as symbols (from https://clue.io/command?q=/gene-space%20lm)
#and converted ourselves to Entrez:
load("./data-raw/internalData/internalData-raw/LM_genes_entrez.rda")
#for more info and list of lankmark genes check https://clue.io/connectopedia//what_are_landmark_genes.

###Preparing internal data for Homogenizer.RUV and Homogenizer.RUV4:
## In Homogenizer.RUV and Homogenizer.RUV4, a list of house keeping genes are used. This list is based on "'Human housekeeping
#genes revisited', E. Eisenberg and E.Y. Levanon, Trends in Genetics, 29 (2013)' accessible by the link:
# https://www.tau.ac.il/~elieis/HKG/, we downloaded and converted the list from symbols to Entrez IDs:
load("./data-raw/internalData/internalData-raw/HK_genes_entrez.rda")



###At the end we save all the variables explained above as internal data for the package:
devtools::use_data(ConvTableSym2Entrez, ConvTableGo2Sym, PathID2Entrez, Entrez2PathID,
                   HK_genes_entrez, LM_genes_entrez,internal = T)
