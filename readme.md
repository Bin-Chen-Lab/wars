# Reversal of Infected Host Gene Expression Identifies Repurposed Drug Candidates for COVID-19

## Abstract
Repurposing existing drugs is a timely option to cope with COVID-19. We sought to target infection-induced genes in the host cells, hoping to mitigate disease progression and alleviate symptoms. Based on our previous experience that reversal of gene expression (namely **sRGES**) correlates to drug efficacy, we utilized **a systems-based approach that employs gene expression profiles of SARS-/MERS-CoV infected samples** (as COVID-19 samples not available at the time of writing) **and drug-induced gene expression profiles from cell lines to discover new thereapeutic candidates for COVID-19**. 430 samples from 12 studies of SARS-/MERS-CoV infection were collected from [GEO](https://www.ncbi.nlm.nih.gov/geo/), 215 comparisons across all the infection time duration were employed to create different infection signatures and predict drug candidates. Drugs effective to SARS or MERS served as positive control to select valid signatures that captured virus induced biology change in the host cells. Finally we computed the consensus score of drug predictions derived from all the valid infection signatures to generate the final prediction.

![](Figure_1.pdf)

## Instructions

### Before Running
Create data folder and code folder. Download data from Chen lab server (/home/ubuntu/chenlab_v2/chenlab_data/raw.zip), and unzip to data folder. Clone code and unzip to the code folder.

### Enumerate All Comparisons & Create Disease Signatures
Three ways of comparisons:
1. Between virus and mock at the same time point




#choose disease name (more specifically the comparison of name, I use GSE number + comparison method), here we use GSE17400_SARS_48_vs_12
#method_id = 3: by default we use SAM (3), but if we could not find enough signature genes, use rankProd (2)
#case_reg: can be keywords (or patten) or a list of GSM ids.

#create_dz_signature_from_GEO.R: create disase signature
#predict_drugs_cmap.R: predict drugs using cmap data
#predict_drugs_lincs:R: predict drugs using lincs data

