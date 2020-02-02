#this pipeline is a simple version of predicting drug hits for a given disease
#to get started, you need to select one disease and find one dataset from GEO that includes both case (disease) and control (healthy/non-tumor) samples
#depending the question you are asking, you can also compare responders vs non-responders, mutation vs wild type, etc.

setwd("~/Documents/stanford/wars/cmap/data")

#insall packages before running the code
#source("http://www.bioconductor.org/biocLite.R")
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocInstaller")
#biocLite("impute")
#biocLite("siggenes")
#biocLite("RankProd")
#biocLite("preprocessCore")
#biocLite("GEOquery")
#biocLite("qvalue")

#BiocManager::install(c("impute", "siggenes", "RankProd", "preprocessCore", "GEOquery", "qvalue"))
#installed.packages(c("pheatmap", "gplots", "ggplot2", "RColorBrewer", "plyr"))

library("pheatmap")
library("gplots")
library("ggplot2")
library("RColorBrewer")
library("plyr")
library("stringr")
###############################
#parameters

###############################

meta = read.csv("../code/meta.csv", stringsAsFactors = F)
GSEs = unique(meta$GSE)
#GSEs = c("GSE47960")
for (GSE in GSEs){
  diseases = unique(meta$disease[meta$GSE == GSE])
  for (dz in diseases){
    models = unique(meta$model[meta$GSE == GSE & meta$disease == dz])
    for (model in models){
      times = sort(unique(meta$time[meta$GSE == GSE & meta$disease == dz & meta$model == model]))
      if (length(times) <= 1) next
       for (i in 1:(length(times) - 1)){
         control_time = times[i]
         case_time = times[i + 1]
         disease = paste("meta", GSE, dz, model, control_time, case_time, sep = "_")
         disease = str_replace_all(disease, " ", "_")
         dz_sig.file <- paste(disease, "/dz_signature_cmap.txt", sep="") #paste(GSE, "_dz_sig_",method_id,".txt",sep="")
         dz_sig_all.file <- paste(disease, "/dz_signature_all.txt", sep="") #paste(GSE, "_dz_sig_",method_id,".txt",sep="")
         
         if (!dir.exists(disease)) dir.create(disease)
         if (file.exists(paste0(disease, "/cmap_drug_predictions.csv"))) next
         if (disease %in% c("meta_GSE45042_EMC_Calu-3_3_7", "meta_GSE45042_EMC_Calu-3_7_12", "meta_GSE17400_SARS_Calu-3_12_24",
                            "meta_GSE17400_MOCK_Calu-3_12_24", "meta_GSE17400_MOCK_Calu-3_24_48")) next #skip problematic disease
         
         print(disease)
         
         case_reg <- meta$GSM[meta$GSE == GSE & meta$disease == dz & meta$model == model & meta$time == case_time]
         control_reg <- meta$GSM[meta$GSE == GSE & meta$disease == dz & meta$model == model & meta$time == control_time]
         
         method_id <- 3 #2: rankprod, 3: siggenes sam, 
         q_thresh <- 0.05 #fdr; may loose the threshold if there are no differentially expressed genes
         fold_change <- 1 #fold change cutoff
         
         source("../code/create_dz_signature_from_GEO.R")
         
         #if method 3 does not give enough signature genes, try method 2, which is more sensitive
         dz_signature = read.csv(dz_sig.file)
         if (nrow(dz_signature) < 30){
           #method_id <- 2 #2: rankprod, 3: siggenes sam, 
           #source("../code/create_dz_signature_from_GEO.R")
           next
         }
         
         #map gene to Homo sapiens
         if (meta$organism[meta$GSE == GSE][1] != "Homo sapiens"){
           gene_mapping = read.csv("raw/cmap/gene_info_hs.csv")[, c("GeneID", "Symbol")]
           
           dz_signature = read.csv(dz_sig_all.file, sep = "\t", stringsAsFactors = F)
           dz_signature$Symbol = toupper(sapply(dz_signature$Symbol, function(x) unlist(strsplit(x, " "))[1]))
           dz_signature = merge(dz_signature, gene_mapping, by = "Symbol")
           dz_signature$GeneID = dz_signature$GeneID.y
           write.table(dz_signature,dz_sig_all.file,sep="\t",quote=F,row.names=F,col.names=T)
           
           dz_signature = dz_signature[order(dz_signature$value), ]
           dz_signature = rbind(head(dz_signature, 150), tail(dz_signature, 150))
           write.table(dz_signature, dz_sig.file,sep="\t",quote=F,row.names=F,col.names=T)
         }
         
         dz_signature = read.csv(dz_sig.file)
         if (sum(dz_signature$up_down == "up") < 3 | sum(dz_signature$up_down == "down") < 3) next
         #predict  drugs using connectivity map data
         cmd <- paste("Rscript ../code/predict_drugs_cmap.R", disease, "cmap", paste0(disease, "/dz_signature_cmap.txt", sep=""))
         system(cmd)
         source("../code/drug_repurpose.R")
         
         dz_signature = read.csv(dz_sig_all.file)
         if (nrow(dz_signature) < 5) next
         #predict  drugs using lincs data
         cmd <- paste("Rscript ../code/predict_drugs_lincs.R", disease,  paste0(disease, "/dz_signature_all.txt", sep=""))
         system(cmd)
         
         gc()
       }
      
    }
  }
}






