setwd("~/Documents/stanford/wars/cmap/data//")

meta_files = list.files(pattern =  "meta_")

prediction = read.csv(paste0(getwd(), "/", "meta_GSE17400_DOHV_Calu-3_12_24", "/sRGES_drugs.csv"), stringsAsFactors = F)
prediction = prediction[order(prediction$name), ]

prediction = prediction[prediction$moa != "", ]
drug_moa = NULL
for (i in 1:nrow(prediction)){
  moas = unlist(strsplit(prediction$moa[i], "\\|"))
  drug_moa = rbind(drug_moa, data.frame(drug = prediction$name[i], moa = moas))
}

moas = table(drug_moa$moa)
moas = names(moas[moas > 3])

for (meta_file  in meta_files){
  if (file.exists(paste0(getwd(), "/", meta_file, "/sRGES_drugs.csv"))) {
    prediction = read.csv(paste0(getwd(), "/", meta_file, "/sRGES_drugs.csv"))

    moa_enriched = NULL
    for (moa in moas){
      test = wilcox.test(prediction$sRGES[prediction$name %in% drug_moa$drug[drug_moa$moa == moa]],
                         prediction$sRGES[!prediction$name %in% drug_moa$drug[drug_moa$moa == moa]])
      moa_enriched = c(moa_enriched, test$p.value)
    }
    moa_enriched = data.frame(moa = moas, p = p.adjust(moa_enriched))
    
    write.csv(moa_enriched, paste0(getwd(), "/", meta_file, "/sRGES_drug_MOA.csv"))
  }
}

meta_moas = data.frame(moa = moas)
meta_files_valid = NULL
for (meta_file  in meta_files){
  if (file.exists(paste0(getwd(), "/", meta_file, "/sRGES_drugs.csv"))) {
    prediction = read.csv(paste0(getwd(), "/", meta_file, "/sRGES_drug_MOA.csv"))
    meta_moas = cbind(meta_moas, prediction$p)
    meta_files_valid = c(meta_files_valid, meta_file)
  }
}

rownames(meta_moas) = meta_moas$moa
meta_moas = meta_moas[, -1]
colnames(meta_moas) = meta_files_valid

write.csv(t(meta_moas), "meta_moas.csv")
