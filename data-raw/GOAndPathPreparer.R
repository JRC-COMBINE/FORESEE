#

require(biomaRt)
humaRt <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
ConvTableSym2Entrez <- getBM(attributes = c("hgnc_symbol","entrezgene"),values = "*", mart = humaRt)
ConvTableGo2Sym <- getBM(attributes = c("go_id","hgnc_symbol"),values = "*", mart = humaRt)
ConvTableKEGG2Sym <- getBM(attributes = c("kegg_enzyme","hgnc_symbol"),values = "*", mart = humaRt)
## sum(ConvTableKEGG2Sym$kegg_enzyme == "") -> mostly empty! gonna switch to reactome:
require(reactome.db)
PathID2Entrez <- as.list(reactomePATHID2EXTID)
Entrez2PathID <- as.list(reactomeEXTID2PATHID)

load("./data-raw/HK_genes_entrez.rda") #Need the script that made these lists
load("./data-raw/LM_genes_entrez.rda") #Need the script that made these lists

devtools::use_data(ConvTableSym2Entrez, ConvTableGo2Sym, PathID2Entrez, Entrez2PathID,
                   HK_genes_entrez, LM_genes_entrez,internal = T)
