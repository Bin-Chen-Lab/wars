setwd("~/Documents/stanford/wars/cmap/data//")
data_dir = "~/Documents/stanford/tumor_cell_line/pipeline/octad_desktop/results/"
GSE17400 = read.csv(paste0(data_dir, "GSE17400", "/sRGES_drugs.csv"), stringsAsFactors = F)
GSE30589 = read.csv(paste0(data_dir, "GSE30589", "/sRGES_drugs.csv"), stringsAsFactors = F)
GSE45042 = read.csv(paste0(data_dir, "GSE45042", "/sRGES_drugs.csv"), stringsAsFactors = F)
GSE47960 = read.csv(paste0(data_dir, "GSE47960", "/sRGES_drugs.csv"), stringsAsFactors = F)
GSE22581 = read.csv(paste0(data_dir, "GSE22581", "/sRGES_drugs.csv"), stringsAsFactors = F)
GSE68820 = read.csv(paste0(data_dir, "GSE68820", "/sRGES_drugs.csv"), stringsAsFactors = F)
GSE36016 = read.csv(paste0(data_dir, "GSE36016", "/sRGES_drugs.csv"), stringsAsFactors = F)

all = data.frame( GSE17400 = GSE17400$sRGES, GSE30589 = GSE30589$sRGES, GSE45042 = GSE45042$sRGES, GSE47960 = GSE47960$sRGES,
                 GSE22581 = GSE22581$sRGES, GSE68820 = GSE68820$sRGES, GSE36016 = GSE36016$sRGES)
rownames(all) = GSE17400$pert_iname
all = all[GSE17400$n > 1, ]

all_rank = all
for( i in 1:ncol(all)){
  all_rank[,i] = rank(all[,i])
}

all_rank["deferiprone", ]
drug_rank = apply(all_rank, 1, median)
tail(sort(drug_rank))
head(sort(drug_rank))
drug_rank["lopinavir"]
drug_rank["niclosamide"]

