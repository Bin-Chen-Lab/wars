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
  diseases = unique(meta$disease[meta$GSE == GSE & meta$disease != "MOCK"])
  for (dz in diseases){
    models = unique(meta$model[meta$GSE == GSE & meta$disease == dz])
    for (model in models){
      times = sort(unique(meta$time[meta$GSE == GSE & meta$disease == dz & meta$model == model]))
       for (control_time in 1:(length(times) -1 )){
         for (case_time in (control_time + 1):length(times)){
           print(i)
           #first compute signature between time points and then filter those appeared in the mock group
            disease = paste("all_virus_meta", GSE, dz, model, times[control_time], times[case_time], sep = "_")
            disease = str_replace_all(disease, " ", "_")
            
            case_reg <- meta$GSM[meta$GSE == GSE & meta$disease == dz & meta$model == model & meta$time == times[case_time]]
            control_reg <- meta$GSM[meta$GSE == GSE & meta$disease == dz & meta$model == model & meta$time == times[control_time]]
            
            if (length(case_reg) < 3 | length(control_reg) < 3) next
            
            dz_sig.file <- paste(disease, "/dz_signature_cmap.txt", sep="") #paste(GSE, "_dz_sig_",method_id,".txt",sep="")
            dz_sig_all.file <- paste(disease, "/dz_signature_all.txt", sep="") #paste(GSE, "_dz_sig_",method_id,".txt",sep="")
         
            if (!dir.exists(disease)) dir.create(disease)
            if (file.exists(paste0(disease, "/sRGES.csv"))) next
    
            method_id <- 3 #2: rankprod, 3: siggenes sam, 
            q_thresh <- 0.05 #fdr; may loose the threshold if there are no differentially expressed genes
            fold_change <- 1 #fold change cutoff
            
            source("../code/create_dz_signature_from_GEO.R")
         
            dz_signature = read.table(dz_sig.file, sep = "\t")
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
            
            
            ########
            #MOCK group
            disease = paste("all_mock_meta", GSE, dz, model, times[control_time], times[case_time], sep = "_")
            disease = str_replace_all(disease, " ", "_")
            
            case_reg <- meta$GSM[meta$GSE == GSE & meta$disease == "MOCK" & meta$model == model & meta$time == times[case_time]]
            control_reg <- meta$GSM[meta$GSE == GSE & meta$disease == "MOCK" & meta$model == model & meta$time == times[control_time]]
            
            if (length(case_reg) < 3 | length(control_reg) < 3) next
            
            dz_sig.file <- paste(disease, "/dz_signature_cmap.txt", sep="") #paste(GSE, "_dz_sig_",method_id,".txt",sep="")
            dz_sig_all.file <- paste(disease, "/dz_signature_all.txt", sep="") #paste(GSE, "_dz_sig_",method_id,".txt",sep="")
            
            if (!dir.exists(disease)) dir.create(disease)
            if (file.exists(paste0(disease, "/sRGES.csv"))) next
            
            method_id <- 3 #2: rankprod, 3: siggenes sam, 
            q_thresh <- 0.05 #fdr; may loose the threshold if there are no differentially expressed genes
            fold_change <- 1 #fold change cutoff
            
            source("../code/create_dz_signature_from_GEO.R")
            
            dz_signature = read.table(dz_sig.file, sep = "\t")
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
            
            
            ########
            #virus vs MOCK group
            disease = paste("all_virus_mock_meta", GSE, dz, model, times[control_time], times[case_time], sep = "_")
            disease = str_replace_all(disease, " ", "_")
            
            case_dz = str_replace_all(paste0(paste("all_virus_meta", GSE, dz, model, times[control_time], times[case_time], sep = "_"), "/dz_signature_all.txt"), " ", "_")
            control_dz = str_replace_all(paste0(paste("all_mock_meta", GSE, dz, model, times[control_time], times[case_time], sep = "_"), "/dz_signature_all.txt"), " ", "_")
            if (!file.exists(case_dz) | !file.exists(control_dz)) next
            
            dz_signature = read.table(case_dz, sep = "\t", header=T)
            mock_signature = read.table(control_dz, sep = "\t", header=T)
            
            #print(paste("dz", nrow(dz_signature)))
            #print(paste("mock", nrow(mock_signature)))
            if (!dir.exists(disease)) dir.create(disease)
            
            dz_signature = dz_signature[!dz_signature$Symbol %in% mock_signature$Symbol,]
            print(nrow(dz_signature))
            if (nrow(dz_signature) > 30){
              print(disease)
              write.table(dz_signature, paste0(disease, "/dz_signature_all.txt"), quote=F, col.names = T, row.names = F, sep = "\t")
            }
             gc()
       }
      
    }
    }
  }
}


#prediction
diseases = list.files(pattern =  "all_")
diseases = sample(diseases, length(diseases))
for (disease in diseases){
  if (file.exists(paste0(getwd(), "/", disease, "/sRGES.csv"))) next
  if (file.exists(paste0(getwd(), "/", disease, "/dz_signature_all.txt"))) {
    cmd <- paste("Rscript ../code/predict_drugs_lincs_mini.R", disease,  paste0(getwd(), "/", disease, "/dz_signature_all.txt", sep=""))
    print(cmd)
    system(cmd)
    gc()
  }
}





