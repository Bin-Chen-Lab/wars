setwd("~/Documents/stanford/wars/cmap/data//")

meta_files = list.files(pattern =  "meta_")

load("raw/cmpd_sets_mesh.RData")

moas = cmpd_sets$cmpd.set.names

for (meta_file  in meta_files){
  if (file.exists(paste0(getwd(), "/", meta_file, "/sRGES.csv"))) {
    prediction = read.csv(paste0(getwd(), "/", meta_file, "/sRGES.csv"))

    moa_enriched = NULL
    for (moa in moas){
      test = wilcox.test(prediction$sRGES[prediction$pert_iname %in% tolower(cmpd_sets$cmpd.sets[[which(cmpd_sets$cmpd.set.names == moa)]])],
                         prediction$sRGES[!prediction$pert_iname %in% tolower(cmpd_sets$cmpd.sets[[which(cmpd_sets$cmpd.set.names == moa)]])])
      moa_enriched = c(moa_enriched, test$p.value)
    }
    moa_enriched = data.frame(moa = moas, p = p.adjust(moa_enriched))
    
    write.csv(moa_enriched, paste0(getwd(), "/", meta_file, "/sRGES_mesh.csv"))
  }
}

meta_moas = data.frame(moa = moas)
meta_files_valid = NULL
for (meta_file  in meta_files){
  if (file.exists(paste0(getwd(), "/", meta_file, "/sRGES_drugs.csv"))) {
    prediction = read.csv(paste0(getwd(), "/", meta_file, "/sRGES_mesh.csv"))
    meta_moas = cbind(meta_moas, prediction$p)
    meta_files_valid = c(meta_files_valid, meta_file)
  }
}

rownames(meta_moas) = meta_moas$moa
meta_moas = meta_moas[, -1]
colnames(meta_moas) = meta_files_valid

types = sapply(rownames(t(meta_moas)), function(x) unlist(strsplit(x, "_"))[3])
write.csv(cbind(t(meta_moas), types), "meta_mesh.csv")


