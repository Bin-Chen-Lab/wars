### Development of biomarkers
## Started on Date: March 26, 2020
## By Rama
#############################
##################
## data for GSE47960
##################
library(readxl)
meta_data = read_xlsx("Validation_data/Validation_metadata.xlsx", sheet = 1)

TPM_GSE47960 <- read.csv("Validation_data/GSE47960/GSE47960_normalized_exp.csv", row.names = 1, stringsAsFactors = F)
TPM_GSE47960_1 <- data.frame(t(TPM_GSE47960))
TPM_GSE47960_2 = TPM_GSE47960_1[, colnames(TPM_GSE47960_1) %in% RF_importance_biomarker$gene]
TPM_GSE47960_3 = merge(TPM_GSE47960_2, meta_data %>% dplyr::select(GSM, infection), by.x = 0, by.y = "GSM")
TPM_GSE47960_3$treatment <- ifelse(TPM_GSE47960_3$infection=="MOCK",0,1)
rownames(TPM_GSE47960_3) <- TPM_GSE47960_3[,1]
TPM_GSE47960_3 = TPM_GSE47960_3[,-c(1,13)]




TPM_GSE17400 <- read.csv("Validation_data/GSE17400/GSE17400_normalized_exp.csv", row.names = 1, stringsAsFactors = F)
TPM_GSE17400_1 <- data.frame(t(TPM_GSE17400))
TPM_GSE17400_2 = TPM_GSE17400_1[, colnames(TPM_GSE17400_1) %in% RF_importance_biomarker$gene]
library(readxl)
meta_data = read_xlsx("Validation_data/Validation_metadata.xlsx", sheet = 1)
library(dplyr)
TPM_GSE17400_3 = merge(TPM_GSE17400_2, meta_data %>% dplyr::select(GSM, infection), by.x = 0, by.y = "GSM")

TPM_GSE17400_3$treatment <- ifelse(TPM_GSE17400_3$infection=="MOCK",0,1)
rownames(TPM_GSE17400_3) <- TPM_GSE17400_3[,1]
TPM_GSE17400_3 = TPM_GSE17400_3[,-c(1,13)]
TPM_GSE17400_3$treatment = as.factor(TPM_GSE17400_3$treatment)
########
## prediction for model
#########
pred_GSE17400 <- predict(model, TPM_GSE17400_3, type = "class")

#plot(model, main = "Error rate of random forest")

#y_pred_num <- ifelse(predValid > 0.5, 1,0)
y_pred <- factor(pred_GSE17400, levels = c(0,1))
y_act <- TPM_GSE17400_3$treatment
## Accuracy
mean(y_pred==y_act)
confusionMatrix(pred_GSE17400, TPM_GSE17400_3$treatment)

## plot AUC
GSE17400_test.forest <- predict(model, newdata = TPM_GSE17400_3, type = "prob")
GSE17400_forestprediction <- prediction(GSE17400_test.forest[,2], TPM_GSE17400_3$treatment)
GSE17400_forestperf = performance(GSE17400_forestprediction, "tpr", "fpr")
plot(GSE17400_forestperf, main="ROC", colorize=F)

GSE17400_perf_AUC<- performance(GSE17400_forestprediction,"auc")
GSE17400_auc <- GSE17400_perf_AUC@y.values[[1]]


##########
## data for GSE30589
##########
TPM_GSE30589 <- read.csv("Validation_data/GSE30589/GSE30589_normalized_exp.csv", row.names = 1, stringsAsFactors = F)
TPM_GSE30589_1 <- data.frame(t(TPM_GSE30589))
TPM_GSE30589_2 = TPM_GSE30589_1[, colnames(TPM_GSE30589_1) %in% RF_importance_biomarker$gene]
TPM_GSE30589_3 = merge(TPM_GSE30589_2, meta_data %>% dplyr::select(GSM, infection), by.x = 0, by.y = "GSM")
TPM_GSE30589_3$treatment <- ifelse(TPM_GSE30589_3$infection=="MOCK",0,1)
rownames(TPM_GSE30589_3) <- TPM_GSE30589_3[,1]
TPM_GSE30589_3 = TPM_GSE30589_3[,-c(1,13)]
TPM_GSE30589_3$treatment <- as.factor(TPM_GSE30589_3$treatment)
# prediction

pred_GSE30589 <- predict(model, TPM_GSE30589_3, type = "class")

#plot(model, main = "Error rate of random forest")

#y_pred_num <- ifelse(predValid > 0.5, 1,0)
y_pred <- factor(pred_GSE30589, levels = c(0,1))
y_act <- TPM_GSE30589_3$treatment
## Accuracy
mean(y_pred==y_act)
confusionMatrix(pred_GSE30589, TPM_GSE30589_3$treatment)

## plot AUC
GSE30589_test.forest <- predict(model, newdata = TPM_GSE30589_3, type = "prob")
GSE30589_forestprediction <- prediction(GSE30589_test.forest[,2], TPM_GSE30589_3$treatment)
GSE30589_forestperf = performance(GSE30589_forestprediction, "tpr", "fpr")
plot(GSE30589_forestperf, main="ROC", colorize=F)

GSE30589_perf_AUC<- performance(GSE30589_forestprediction,"auc")
GSE30589_auc <- GSE30589_perf_AUC@y.values[[1]]


###############
## data for GSE45042
##############
TPM_GSE45042 <- read.csv("Validation_data/GSE45042/GSE45042_normalized_exp.csv", row.names = 1, stringsAsFactors = F)
TPM_GSE45042_1 <- data.frame(t(TPM_GSE45042))
TPM_GSE45042_2 = TPM_GSE45042_1[, colnames(TPM_GSE45042_1) %in% RF_importance_biomarker$gene]
TPM_GSE45042_3 = merge(TPM_GSE45042_2, meta_data %>% dplyr::select(GSM, infection), by.x = 0, by.y = "GSM")
TPM_GSE45042_3$treatment <- ifelse(TPM_GSE45042_3$infection=="MOCK",0,1)
rownames(TPM_GSE45042_3) <- TPM_GSE45042_3[,1]
TPM_GSE45042_3 = TPM_GSE45042_3[,-c(1,13)]




##################
## data for GSE79172
##################
TPM_GSE79172 <- read.csv("Validation_data/GSE79172/GSE79172_normalized_exp.csv", row.names = 1, stringsAsFactors = F)
TPM_GSE79172_1 <- data.frame(t(TPM_GSE79172))
TPM_GSE79172_2 = TPM_GSE79172_1[, colnames(TPM_GSE79172_1) %in% RF_importance_biomarker$gene]
TPM_GSE79172_3 = merge(TPM_GSE79172_2, meta_data %>% dplyr::select(GSM, infection), by.x = 0, by.y = "GSM")
TPM_GSE79172_3$treatment <- ifelse(TPM_GSE79172_3$infection=="MOCK",0,1)
rownames(TPM_GSE79172_3) <- TPM_GSE79172_3[,1]
TPM_GSE79172_3 = TPM_GSE79172_3[,-c(1,13)]

##################
## data for GSE79218
##################
TPM_GSE79218 <- read.csv("Validation_data/GSE79218/GSE79218_normalized_exp.csv", row.names = 1, stringsAsFactors = F)
TPM_GSE79218_1 <- data.frame(t(TPM_GSE79218))
TPM_GSE79218_2 = TPM_GSE79218_1[, colnames(TPM_GSE79218_1) %in% RF_importance_biomarker$gene]
TPM_GSE79218_3 = merge(TPM_GSE79218_2, meta_data %>% dplyr::select(GSM, infection), by.x = 0, by.y = "GSM")
TPM_GSE79218_3$treatment <- ifelse(TPM_GSE79218_3$infection=="MOCK",0,1)
rownames(TPM_GSE79218_3) <- TPM_GSE79218_3[,1]
TPM_GSE79218_3 = TPM_GSE79218_3[,-c(1,13)]










predValid_GSE17400 <- predict(model, TPM_GSE17400_3, type = "class")

plot(model, main = "Error rate of random forest")

#y_pred_num <- ifelse(predValid > 0.5, 1,0)
y_pred <- factor(predValid, levels = c(0,1))
y_act <- ValidSet_for_prediction$treatment
## Accuracy
mean(y_pred==y_act)
confusionMatrix(predValid, ValidSet_for_prediction$treatment)

## plot AUC
test.forest <- predict(model, newdata = ValidSet_for_prediction, type = "prob")
forestprediction <- prediction(test.forest[,2], ValidSet_for_prediction$treatment)
forestperf = performance(forestprediction, "tpr", "fpr")
plot(forestperf, main="ROC", colorize=F)

perf_AUC<- performance(forestprediction,"auc")
auc <- perf_AUC@y.values[[1]]

