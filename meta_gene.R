setwd("~/Documents/stanford/wars/cmap/data//")

meta_files = list.files(pattern =  "meta_")
dz_signatures = NULL
meta_files_valid = NULL
for (meta_file  in meta_files){
  if (file.exists(paste0(getwd(), "/", meta_file, "/sRGES_drugs.csv"))) {
    dz_signature = read.table(paste0(getwd(), "/", meta_file, "/dz_signature_all.txt"), sep="\t", header=T)
    
    if (nrow(dz_signature) < 5) next
    dz_signatures = rbind(dz_signatures, data.frame(name = meta_file, dz_signature))
    meta_files_valid = c(meta_files_valid, meta_file)
  }
}

dz_signatures_mock = (dz_signatures[grep("MOCK", dz_signatures$name), ])
dz_signatures_cov = (dz_signatures[grep("SARS|MERS", dz_signatures$name), ])
dz_signatures_mock = unique(dz_signatures_mock[dz_signatures_mock$up_down == "up", c("Symbol", "name")])
dz_signatures_cov = unique(dz_signatures_cov[dz_signatures_cov$up_down == "up", c("Symbol", "name")])

dz_signatures_cov_genes = table(dz_signatures_cov$Symbol)/length(unique(dz_signatures_cov$name))
dz_signatures_mock_genes = table(dz_signatures_mock$Symbol)/length(unique(dz_signatures_mock$name))

tail(sort(dz_signatures_cov_genes))
tail(sort(dz_signatures_mock_genes))

cov_genes_up = intersect(names(dz_signatures_cov_genes[dz_signatures_cov_genes > 0.5]),
          names(dz_signatures_mock_genes[dz_signatures_mock_genes < 0.5]))


dz_signatures_mock = (dz_signatures[grep("MOCK", dz_signatures$name), ])
dz_signatures_cov = (dz_signatures[grep("SARS|MERS", dz_signatures$name), ])
dz_signatures_mock = unique(dz_signatures_mock[dz_signatures_mock$up_down == "down", c("Symbol", "name")])
dz_signatures_cov = unique(dz_signatures_cov[dz_signatures_cov$up_down == "down", c("Symbol", "name")])

dz_signatures_cov_genes = table(dz_signatures_cov$Symbol)/length(unique(dz_signatures_cov$name))
dz_signatures_mock_genes = table(dz_signatures_mock$Symbol)/length(unique(dz_signatures_mock$name))

tail(sort(dz_signatures_cov_genes))
tail(sort(dz_signatures_mock_genes))

cov_genes_down = intersect(names(dz_signatures_cov_genes[dz_signatures_cov_genes > 0.5]),
                      names(dz_signatures_mock_genes[dz_signatures_mock_genes < 0.5]))

cov_genes = rbind(data.frame(gene = cov_genes_up, up_down = "up"), data.frame(gene = cov_genes_down, up_down = "down"))

write.csv(cov_genes, "meta_cov_genes.csv")
