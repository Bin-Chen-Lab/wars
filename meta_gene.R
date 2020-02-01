setwd("~/Documents/stanford/wars/cmap/data//")
GSE17400 = read.csv("sars_GSE17400/dz_signature_all.txt", sep = "\t", stringsAsFactors = F)
GSE22581 = read.csv("GSE22581/dz_signature_all_mapped.txt", sep = "\t", stringsAsFactors = F)
GSE22581 = GSE22581[, colnames(GSE17400)]
GSE30589 = read.csv("GSE30589/dz_signature_all.txt", sep = "\t", stringsAsFactors = F)
GSE36016 = read.csv("GSE36016/dz_signature_all_mapped.txt", sep = "\t", stringsAsFactors = F)
GSE36016 = GSE36016[, colnames(GSE17400)]
GSE45042 = read.csv("GSE45042/dz_signature_all.txt", sep = "\t", stringsAsFactors = F)
GSE47960 = read.csv("GSE47960/dz_signature_all.txt", sep = "\t", stringsAsFactors = F)
GSE68820 = read.csv("GSE68820/dz_signature_all_mapped.txt", sep = "\t", stringsAsFactors = F)
GSE68820 = GSE68820[, colnames(GSE17400)]


all = data.frame(GSE = "GSE17400", type="human", GSE17400)
all = rbind(all, data.frame(GSE = "GSE22581", type = "other", GSE22581))
all = rbind(all, data.frame(GSE = "GSE30589", type = "human", GSE30589))
all = rbind(all, data.frame(GSE = "GSE36016", type = "other", GSE36016))
all = rbind(all, data.frame(GSE = "GSE45042", type = "human", GSE45042))
all = rbind(all, data.frame(GSE = "GSE47960", type = "human", GSE47960))
all = rbind(all, data.frame(GSE = "GSE68820", type = "other", GSE68820))

all$Symbol = sapply(all$Symbol, function(x) unlist(strsplit(x, " "))[1])
all = unique(all[, -4])

tail(sort(table(all$GeneID)), 10)

genes = (table(all$GeneID))
top_genes = top_genes[top_genes > 3]

sum(genes > 4)
all[all$GeneID == "3454",]
#all[all$Symbol == "NPIPL3",]

all_select = all[all$GeneID %in% names(genes[genes > 4]), ]
all_select$direction = as.numeric(as.factor(all_select$up_down))
all_select_collapse = aggregate(cbind(value, q.value, direction) ~ GeneID,  all_select, mean)
all_select_collapse = all_select_collapse[all_select_collapse$direction == 2, ]
gene_mapping = read.csv("raw/cmap/gene_info_hs.csv")[, c("GeneID", "Symbol")]
all_select_collapse = merge(all_select_collapse, gene_mapping,by="GeneID")

write.csv(all_select_collapse, "top_gene.csv")


#H1NI
GSE47960_H1N1 = read.csv("GSE47960_H1N1/dz_signature_all.txt", sep = "\t", stringsAsFactors = F)
GSE17400_DOHV = read.csv("DOHV_GSE17400/dz_signature_all.txt", sep = "\t", stringsAsFactors = F)
GSE17400_DOHV$Symbol = sapply(GSE17400_DOHV$Symbol, function(x) unlist(strsplit(x, " "))[1])
GSE47960_H1N1$Symbol = sapply(GSE47960_H1N1$Symbol, function(x) unlist(strsplit(x, " "))[1])

all_select_collapse_specific = all_select_collapse[!all_select_collapse$GeneID%in% 
                                                     c(as.character(GSE17400_DOHV$GeneID), as.character(GSE47960_H1N1$GeneID)), ]

all_select_collapse_specific = merge(all_select_collapse_specific, gene_mapping,by="GeneID")

write.csv(all_select_collapse_specific, "top_gene_specific.csv")



