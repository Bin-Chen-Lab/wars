setwd("~/Documents/stanford/wars/cmap/data//")

library(stringr)
pathway_set_name = "Virus_Perturbations_from_GEO"
pathway_gene_set_up = read.delim(paste0( pathway_set_name, "_up.txt"),  sep ="$", header=F, stringsAsFactors = F)
pathway_gene_set_down = read.delim(paste0(pathway_set_name, "_down.txt"),  sep ="$", header=F, stringsAsFactors = F)

if (!file.exists(paste0(pathway_set_name))) {dir.create(paste0(pathway_set_name))}

for (pathway_id in 1:nrow(pathway_gene_set_up)){  #nrow(pathway_gene_set) nrow(pathway_gene_set_up
  pathways_up = unlist(strsplit(pathway_gene_set_up[pathway_id,1], "\t\t"))
  pathway_name_up = pathways_up[1]
  
  pathways_down = unlist(strsplit(pathway_gene_set_down[pathway_id,1], "\t\t"))
  pathway_name_down = pathways_down[1]
  
  pathway_genes_up = as.character(sapply(unlist(strsplit(pathways_up[2], "\t")), function(x) {unlist(strsplit(x, ","))[1]}))
  pathway_genes_up_value = as.character(sapply(unlist(strsplit(pathways_up[2], "\t")), function(x) {unlist(strsplit(x, ","))[2]}))
  
  pathway_genes_down = as.character(sapply(unlist(strsplit(pathways_down[2], "\t")), function(x) {unlist(strsplit(x, ","))[1]}))
  pathway_genes_down_value = as.character(sapply(unlist(strsplit(pathways_down[2], "\t")), function(x) {unlist(strsplit(x, ","))[2]}))
  
  if (pathway_name_up != pathway_name_down){
    stop(paste(pathway_id, " name inconsitent"))
  }else{
    disease = pathway_name_up
  }
  disease = str_replace_all(disease, " ", "_")
  
  dz_signature = rbind(data.frame(Symbol = pathway_genes_up, up_down = "up", value = as.numeric(pathway_genes_up_value), GeneID = NA, p.value = 0.01, q.value = 0.01),
                       data.frame(Symbol = pathway_genes_down, up_down = "down", value = as.numeric(pathway_genes_down_value), GeneID = NA, p.value = 0.01, q.value = 0.01))
  
  if (!file.exists(paste0(pathway_set_name, "/", disease))) {dir.create(paste0( pathway_set_name, "/", disease))}

  write.table(dz_signature, paste0( pathway_set_name, "/", disease, "/dz_signature_all.txt"), col.names = T, row.names = F, sep = "\t", quote=F)
}

diseases = list.files(pathway_set_name)
for (disease in diseases){
  if (file.exists(paste0(pathway_set_name, "/", disease, "/sRGES.csv"))) next
  cmd <- paste("Rscript ../code/predict_drugs_lincs_mini.R", disease,  paste0(pathway_set_name, "/", disease, "/dz_signature_all.txt", sep=""))
  print(cmd)
  system(cmd)
}
        
