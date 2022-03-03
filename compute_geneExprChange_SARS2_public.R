library('octad')
setwd('~/Documents/OneDrive - Michigan State University/CoV/WARS/data/RNA_seq_expression_data/GEC_210818/')

load('../COVID1783.RData')
RDCT = COVID1783_read.count.matrix
rm(COVID1783_fpkm.matrix)
rm(COVID1783_read.count.matrix)
rm(COVID1783_tpm.matrix)

comparison_table = read.csv('../SARS2_public_comparisons.csv', stringsAsFactors=F)
for (i in 1:dim(comparison_table)[1]) {
  comp = comparison_table$Comparison_Name[i]
  print(comp)
  #if (comp != 'SRP297170__Caco2__SARS-CoV-2') next
  dir.create(comp)
  setwd(comp)
  
  ctrl_ids = unlist(strsplit(comparison_table[i, 'Control_IDs'], '__'))
  case_ids = unlist(strsplit(comparison_table[i, 'Case_IDs'], '__'))
  tryCatch({
    DE = diffExp(case_ids, ctrl_ids, source='side', output=T, n_topGenes=10000, expSet=as.matrix(RDCT), annotate=T)
    write.csv(DE, file=paste0('log2fc__', comp, '.csv'))
  
    #DE = subset(DE, abs(log2FoldChange) > 0.58 & padj < 0.2)
    DE$FDR = p.adjust(DE$pvalue, method='fdr')
    DE = subset(DE, abs(log2FoldChange) >= 0.585 & pvalue < 0.05 & FDR < 0.05)
    write.csv(DE, file=paste0('DZSIG__', comp, '.csv'))
  
    srges = runsRGES(DE, max_gene_size=100, permutations=10000, choose_fda_drugs=T)
    file.rename('sRGES.csv', paste0('sRGES__', comp, '.csv'))
  }, error=function(err){print(paste('ERROR in ', comp, ': ', err))})
  
  setwd('../')
}






