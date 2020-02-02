setwd("~/Documents/stanford/wars/cmap/data//")

meta_files = list.files(pattern =  "GSE")

prediction = read.csv(paste0(getwd(), "/", "meta_GSE17400_DOHV_Calu-3_12_24", "/sRGES_drugs.csv"))
prediction = prediction[order(prediction$name), ]

scores = data.frame(name = prediction$name)
meta_files_valid = NULL
for (meta_file  in meta_files){
  if (file.exists(paste0(getwd(), "/", meta_file, "/sRGES_drugs.csv"))) {
    prediction = read.csv(paste0(getwd(), "/", meta_file, "/sRGES_drugs.csv"))
    prediction = prediction[order(prediction$name), ]
    scores = cbind(scores, prediction$sRGES)
    meta_files_valid = c(meta_files_valid, meta_file)
  }
}
rownames(scores) = scores$name
scores = scores[, -1]
colnames(scores) = meta_files_valid
#remove n = 1 profiles
scores = scores[prediction$n > 1, ]
scores_rank = scores
for (i in 1:ncol(scores)){
  scores_rank[,i] = rank(scores[,i])
}

scores_sars = scores[, grep("SARS", colnames(scores))]
scores_sars_rank = scores_rank[, grep("SARS", colnames(scores_rank))]

scores_sars["ritonavir",]
scores_sars_rank["ritonavir",]

scores_rank["ritonavir",]
scores_rank["niclosamide",]

scores_sars_rank[grep("vir", rownames(scores_sars_rank)), ]

#choose comparisons where any positive hit is on the top 100 ritonavir/niclosamide/lopinavir
scores_rank_pos = scores_rank[rownames(scores_rank) %in% c("ritonavir", "niclosamide", "lopinavir"), ]
scores_rank_pos_min = apply(scores_rank_pos, 2, min)
scores_rank_pos_selected = names(scores_rank_pos_min[scores_rank_pos_min < 100])
scores_rank_pos_selected = scores_rank_pos_selected[!scores_rank_pos_selected %in% scores_rank_pos_selected[grep("MOCK", scores_rank_pos_selected)]]

#among those selected comparisons, choose drugs that are better
scores_rank_selected = scores_rank[, scores_rank_pos_selected]
drug_rank_selected = apply(scores_rank_selected, 1, median)
head(sort(drug_rank_selected), 10)

scores_rank["dabrafenib",scores_rank_pos_selected]


