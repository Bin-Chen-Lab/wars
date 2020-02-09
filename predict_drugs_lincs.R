args <- commandArgs(trailingOnly=T)

disease <- args[1]
dz_sig_path <- args[2]

## ----setup, include=TRUE, echo=TRUE, message=FALSE, warning=FALSE----------------------------------------------------------------------------
start.time <- Sys.time()
# knitr::opts_knit$set(root.dir = '../..')
base.folder <-  getwd() #'~/Documents/stanford/tumor_cell_line/pipeline/octad_desktop//' # should be set to installation parent folder
setwd(base.folder)
library(dplyr)
library(ggplot2)
#library("RUVSeq")
library("GSVA")
library(compiler)
library(data.table)
library(devtools)
library(limma)
library(rhdf5)
# library(doParallel)
library(pheatmap)

enableJIT(3)

outputFolder = paste0(base.folder,'/', disease, "/")
if (!dir.exists(outputFolder)) {
  dir.create(outputFolder)
}
dataFolder = paste0(base.folder,'/raw/')
CodeFolder = paste0('../code/')

load(paste0(dataFolder,'metadata.RData')) # loads 'ensemble_info', 'merged_gene_info', 'breastpam50', 'tsne', and 'phenoDF'

source(paste0(CodeFolder,'core_functions.R'))

# avoiding logging errors
tsne.time    <- "Not run"
GEA.time     <- "Not run"
VDH.time     <- "Not run"
enrich.time  <- "Not run"
ranking.time <- "Not run"


## ----parameters------------------------------------------------------------------------------------------------------------------------------
setwd(base.folder)
# DE params:

normalize_samples = T
k = 1
n_topGenes = 10000
#this demo sample set does not have too much genes. The usual default for n_topGenes = 10000. Set this higher if you think there are more significant genes
DE_method = 'edgeR' # Other choice is 'limma' for DE_method.


# disease signature params:
log2fchange <- 2 # the cutoff for "significant gene impact". Default = 1, so any DE gene with log2foldchange < 1 will be omitted. Improved results with 2.
padjusted <- 0.001 # the cutoff for "significant genes". Default = 0.001, so any DE genes with padj > 0.001 will be omitted.


# sRGES params:
# Only mess with "max_gene_size", which is how many up-regulated genes and down-regulated genes are allowed. ("100" = 100up + 100down = 200 total)
landmark = 0
choose_fda_drugs = F
max_gene_size = 100
weight_cell_line = F

# logging parameters
write('x', file = paste0(outputFolder,"parameters.txt"))
fileConn<-file(paste0(outputFolder,"parameters.txt"))
writeLines(c("Normalization Parameters:","normalize_samples",normalize_samples,"","k",k,"","n_topGenes",n_topGenes,"","DE_method",DE_method,"","--------------------------","","Disease Signature Parameters:","log2fchange",log2fchange,"","padjusted",padjusted,"","--------------------------","","sRGES Parameters:","landmark",landmark,"","choose_fda_drugs",choose_fda_drugs,"","max_gene_size",max_gene_size,"","weight_cell_line",weight_cell_line), fileConn)
close(fileConn)



## ----dataFolder, results='hide'--------------------------------------------------------------------------------------------------------------
setwd(base.folder)
dir(dataFolder)
#check to make sure you have these files


## ----CodeFolder------------------------------------------------------------------------------------------------------------------------------
setwd(base.folder)
dir(CodeFolder)
#check to make sure you have these files


## --------------------------------------------------------------------------------------------------------------------------------------------
setwd(base.folder)
dz_signature = read.csv(dz_sig_path, sep = "\t", stringsAsFactors = F)
dz_signature$log2FoldChange = dz_signature$value #pseudo fc]
dz_signature = dz_signature[abs(dz_signature$log2FoldChange) > 1.2, ]
dz_signature$Symbol = sapply(dz_signature$Symbol, function(x) unlist(strsplit(x, " "))[1])

#lincs genes
if (landmark == 1) {
  load(paste0(base.folder, "/raw/lincs_signatures_cmpd_landmark_symbol.RData"))
} else {
  load(paste0(base.folder, "/raw/lincs_signatures_cmpd_landmark_GSE92742.RData"))
}
dz_signature_used = dz_signature[dz_signature$Symbol %in% rownames(lincs_signatures), ]
write.csv(dz_signature_used, paste0(outputFolder, "/dz_signature_lincs_used.csv"))

# #visualize mapped signature
# load(paste(disease, "/dz_expr.RData", sep=""))
# comparison_frame_subset <- dz_expr[rownames(dz_expr) %in% dz_signature_used$GeneID, ]
# pheatmap(t(scale(t(comparison_frame_subset))), col = my.cols, annotation = annotation,  annotation_colors = anno_colors,
#          show_colnames=F, legend=T, show_rownames=F, filename=paste(disease, "/dz_sig_validation_lincs.pdf", sep="")
# )

## ----sRGES-----------------------------------------------------------------------------------------------------------------------------------
setwd(base.folder)
sRGES.start <- Sys.time()
# source(paste0(CodeFolder,'core_functions.R'))
#functions needed to run RGES
# source(paste0(CodeFolder,'runRGES_dz_arrangeMax_compiler.R'))
runsRGES(
  dz_signature = dz_signature,
  choose_fda_drugs = choose_fda_drugs,
  parallel = F,                         # causes issues for some machines.
  max_gene_size = max_gene_size,
  landmark = landmark
  )
#,cells='HEPG2'
sRGES.end <- Sys.time()
sRGES.time <- sRGES.end - sRGES.start
print(sRGES.time)

## --------------------------------------------------------------------------------------------------------------------------------------------
setwd(base.folder)
#repurposing = read.csv( paste0(base.folder, "/raw/repurposing_drugs_20170327.csv"))
repurposing = read.csv( paste0(base.folder, "/raw/repurposing_drugs_20180907.txt"), skip = 9, sep = "\t")
repurposing$name = tolower(repurposing$pert_iname)

sRGES = read.csv(paste0(outputFolder, "/sRGES.csv"))
sRGES$name = tolower(sRGES$pert_iname)

sRGES_drug = merge(sRGES, repurposing, by = "name")
sRGES_drug = sRGES_drug[order(sRGES_drug$sRGES), ]
write.csv(sRGES_drug, paste0(outputFolder, "/sRGES_drugs.csv"))


## --------------------------------------------------------------------------------------------------------------------------------------------
setwd(base.folder)
enrich.time.start <- Sys.time()

source(paste0(CodeFolder,'core_functions.R'))
#targets = c('chembl_targets','mesh','sea_targets','ChemCluster') # include as many or few as you prefer. 
targets = c('mesh', 'chembl_targets')

enrichFolder <- paste0(outputFolder,'enrichment_analysis/')
if (!dir.exists(enrichFolder)) {
  dir.create(enrichFolder)
}

sRGES = read.csv(paste0(outputFolder,'/sRGES.csv'),stringsAsFactors = F)
load(paste0(dataFolder,"random_gsea_score.RData"))

for (target_type in targets){
  drug_enrichment(sRGES = sRGES, target_type = target_type)
}


enrich.time.stop <- Sys.time()
enrich.time <- enrich.time.stop - enrich.time.start
enrich.time

visualize_dz_sig_pathway(dz_signature)

visualize_drug_hits()

