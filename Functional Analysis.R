library(readxl)
library(tidyverse)
library(maSigPro)
library(dplyr)
library(gbutils)

#Convert FBGN to valid ENTREZ Gene IDs = https://www.biotools.fr/drosophila/fbgn_converter
#GO Analysis on ENTREZ Gene IDs = https://maayanlab.cloud/Enrichr/

#Intersection w/ Modifier Dataset

modifier_data = read_excel('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/Modifier Dataset.xlsx')
modifier_data = modifier_data[complete.cases(modifier_data[ , 1]),]
intersectedModifierData = modifier_data[FALSE,]

full_intersection = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/full_intersection_list.txt')
full_intersection_list <- as.list(as.data.frame(full_intersection))
full_intersection_list <- full_intersection_list$x

#sigs_summary = sigs2[["summary"]]
for (d in full_intersection_list){
  x <- modifier_data[grep(d, modifier_data$Fly.Base.GeneID), ]
  intersectedModifierData <- rbind(intersectedModifierData, x)
}

intersectedModifierData$`Modifier Construct` <- paste(intersectedModifierData$Allele.Class," ",intersectedModifierData$Stock.Number)

intersectedModifierData$Strain.Source = NULL
intersectedModifierData$Drosophila.Gene = NULL
intersectedModifierData$`gmr>Tau`= NULL
intersectedModifierData$`repo>FL-HTT`= NULL
intersectedModifierData$`repo>NT-HTT`= NULL
intersectedModifierData$`Any Human Homolog`= NULL
intersectedModifierData$Gene.Symbol= NULL
intersectedModifierData$Stock.Number= NULL
intersectedModifierData$Allele.Class= NULL

#remember to make results object independently available to this script!
cluster_membership <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  cluster_membership <- rbind (cluster_membership, cluster_results[[i]])
}
intersectedModifierData$`Cluster Membership` <- as.numeric(cluster_membership)
intersectedModifierData = intersectedModifierData[,c(8, 1, 10, 9, 2, 3, 4, 5, 6, 7)]
intersectedModifierData$`elav>APP` = NULL

#Add Log2FC values from previously calculated files
AD_B42_D2_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_10_AD-model–Secretedβ42_vs_NegativeControl(driver)_Age_2_DESeq_all_genes.txt')
rownames(AD_B42_D2_file) <- AD_B42_D2_file[,1]
AD_B42_D2_file$Gene_id = NULL
Log2FC_AD_B42_D2 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_B42_D2 <- c(Log2FC_AD_B42_D2, AD_B42_D2_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_B42_D2` <- Log2FC_AD_B42_D2

AD_B42_D5_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_12_AD-model–Secretedβ42_vs_NegativeControl(driver)_Age_5_DESeq_all_genes.txt')
rownames(AD_B42_D5_file) <- AD_B42_D5_file[,1]
AD_B42_D5_file$Gene_id = NULL
Log2FC_AD_B42_D5 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_B42_D5 <- c(Log2FC_AD_B42_D5, AD_B42_D5_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_B42_D5` <- Log2FC_AD_B42_D5

AD_B42_D7_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_14_AD-model–Secretedβ42_vs_NegativeControl(driver)_Age_7_DESeq_all_genes.txt')
rownames(AD_B42_D7_file) <- AD_B42_D7_file[,1]
AD_B42_D7_file$Gene_id = NULL
Log2FC_AD_B42_D7 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_B42_D7 <- c(Log2FC_AD_B42_D7, AD_B42_D7_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_B42_D7` <- Log2FC_AD_B42_D7

AD_B42_D10_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_16_AD-model–Secretedβ42_vs_NegativeControl(driver)_Age_10_DESeq_all_genes.txt')
rownames(AD_B42_D10_file) <- AD_B42_D10_file[,1]
AD_B42_D10_file$Gene_id = NULL
Log2FC_AD_B42_D10 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_B42_D10 <- c(Log2FC_AD_B42_D10, AD_B42_D10_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_B42_D10` <- Log2FC_AD_B42_D10

AD_B42_D14_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_18_AD-model–Secretedβ42_vs_NegativeControl(driver)_Age_14_DESeq_all_genes.txt')
rownames(AD_B42_D14_file) <- AD_B42_D14_file[,1]
AD_B42_D14_file$Gene_id = NULL
Log2FC_AD_B42_D14 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_B42_D14 <- c(Log2FC_AD_B42_D14, AD_B42_D14_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_B42_D14` <- Log2FC_AD_B42_D14

AD_Tau_D2_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_22_AD-model-Tauopathy_vs_NegativeControl(driver)_Age_2_DESeq_all_genes.txt')
rownames(AD_Tau_D2_file) <- AD_Tau_D2_file[,1]
AD_Tau_D2_file$Gene_id = NULL
Log2FC_AD_Tau_D2 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_Tau_D2 <- c(Log2FC_AD_Tau_D2, AD_Tau_D2_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_Tau_D2` <- Log2FC_AD_Tau_D2

AD_Tau_D5_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_24_AD-model-Tauopathy_vs_NegativeControl(driver)_Age_5_DESeq_all_genes.txt')
rownames(AD_Tau_D5_file) <- AD_Tau_D5_file[,1]
AD_Tau_D5_file$Gene_id = NULL
Log2FC_AD_Tau_D5 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_Tau_D5 <- c(Log2FC_AD_Tau_D5, AD_Tau_D5_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_Tau_D5` <- Log2FC_AD_Tau_D5

AD_Tau_D7_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_26_AD-model-Tauopathy_vs_NegativeControl(driver)_Age_7_DESeq_all_genes.txt')
rownames(AD_Tau_D7_file) <- AD_Tau_D7_file[,1]
AD_Tau_D7_file$Gene_id = NULL
Log2FC_AD_Tau_D7 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_Tau_D7 <- c(Log2FC_AD_Tau_D7, AD_Tau_D7_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_Tau_D7` <- Log2FC_AD_Tau_D7

AD_Tau_D10_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_28_AD-model-Tauopathy_vs_NegativeControl(driver)_Age_10_DESeq_all_genes.txt')
rownames(AD_Tau_D10_file) <- AD_Tau_D10_file[,1]
AD_Tau_D10_file$Gene_id = NULL
Log2FC_AD_Tau_D10 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_Tau_D10 <- c(Log2FC_AD_Tau_D10, AD_Tau_D10_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_Tau_D10` <- Log2FC_AD_Tau_D10

AD_Tau_D14_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_30_AD-model-Tauopathy_vs_NegativeControl(driver)_Age_14_DESeq_all_genes.txt')
rownames(AD_Tau_D14_file) <- AD_Tau_D14_file[,1]
AD_Tau_D14_file$Gene_id = NULL
Log2FC_AD_Tau_D14 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_AD_Tau_D14 <- c(Log2FC_AD_Tau_D14, AD_Tau_D14_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_AD_Tau_D14` <- Log2FC_AD_Tau_D14

Parkinson_D2_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_38_Parkinson’smodel-Synuclein_vs_NegativeControl(driver)_Age_2_DESeq_all_genes.txt')
rownames(Parkinson_D2_file) <- Parkinson_D2_file[,1]
Parkinson_D2_file$Gene_id = NULL
Log2FC_Parkinson_D2 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_Parkinson_D2 <- c(Log2FC_Parkinson_D2, Parkinson_D2_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_Parkinson_D2` <- Log2FC_Parkinson_D2

Parkinson_D5_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_40_Parkinson’smodel-Synuclein_vs_NegativeControl(driver)_Age_5_DESeq_all_genes.txt')
rownames(Parkinson_D5_file) <- Parkinson_D5_file[,1]
Parkinson_D5_file$Gene_id = NULL
Log2FC_Parkinson_D5 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_Parkinson_D5 <- c(Log2FC_Parkinson_D5, Parkinson_D5_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_Parkinson_D5` <- Log2FC_Parkinson_D5

Parkinson_D7_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_42_Parkinson’smodel-Synuclein_vs_NegativeControl(driver)_Age_7_DESeq_all_genes.txt')
rownames(Parkinson_D7_file) <- Parkinson_D7_file[,1]
Parkinson_D7_file$Gene_id = NULL
Log2FC_Parkinson_D7 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_Parkinson_D7 <- c(Log2FC_Parkinson_D7, Parkinson_D7_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_Parkinson_D7` <- Log2FC_Parkinson_D7

Parkinson_D10_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_44_Parkinson’smodel-Synuclein_vs_NegativeControl(driver)_Age_10_DESeq_all_genes.txt')
rownames(Parkinson_D10_file) <- Parkinson_D10_file[,1]
Parkinson_D10_file$Gene_id = NULL
Log2FC_Parkinson_D10 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_Parkinson_D10 <- c(Log2FC_Parkinson_D10, Parkinson_D10_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_Parkinson_D10` <- Log2FC_Parkinson_D10

Parkinson_D14_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_46_Parkinson’smodel-Synuclein_vs_NegativeControl(driver)_Age_14_DESeq_all_genes.txt')
rownames(Parkinson_D14_file) <- Parkinson_D14_file[,1]
Parkinson_D14_file$Gene_id = NULL
Log2FC_Parkinson_D14 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_Parkinson_D14 <- c(Log2FC_Parkinson_D14, Parkinson_D14_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_Parkinson_D14` <- Log2FC_Parkinson_D14

NT_HTT_D2_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_50_Huntington’sModel–NT-HTT-128Q_vs_NegativeControl(driver)_Age_2_DESeq_all_genes.txt')
rownames(NT_HTT_D2_file) <- NT_HTT_D2_file[,1]
NT_HTT_D2_file$Gene_id = NULL
Log2FC_NT_HTT_D2 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_NT_HTT_D2 <- c(Log2FC_NT_HTT_D2, NT_HTT_D2_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_NT_HTT_D2` <- Log2FC_NT_HTT_D2

NT_HTT_D5_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_52_Huntington’sModel–NT-HTT-128Q_vs_NegativeControl(driver)_Age_5_DESeq_all_genes.txt')
rownames(NT_HTT_D5_file) <- NT_HTT_D5_file[,1]
NT_HTT_D5_file$Gene_id = NULL
Log2FC_NT_HTT_D5 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_NT_HTT_D5 <- c(Log2FC_NT_HTT_D5, NT_HTT_D5_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_NT_HTT_D5` <- Log2FC_NT_HTT_D5

NT_HTT_D7_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_54_Huntington’sModel–NT-HTT-128Q_vs_NegativeControl(driver)_Age_7_DESeq_all_genes.txt')
rownames(NT_HTT_D7_file) <- NT_HTT_D7_file[,1]
NT_HTT_D7_file$Gene_id = NULL
Log2FC_NT_HTT_D7 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_NT_HTT_D7 <- c(Log2FC_NT_HTT_D7, NT_HTT_D7_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_NT_HTT_D7` <- Log2FC_NT_HTT_D7

NT_HTT_D10_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_56_Huntington’sModel–NT-HTT-128Q_vs_NegativeControl(driver)_Age_10_DESeq_all_genes.txt')
rownames(NT_HTT_D10_file) <- NT_HTT_D10_file[,1]
NT_HTT_D10_file$Gene_id = NULL
Log2FC_NT_HTT_D10 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_NT_HTT_D10 <- c(Log2FC_NT_HTT_D10, NT_HTT_D10_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_NT_HTT_D10` <- Log2FC_NT_HTT_D10

NT_HTT_D14_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_58_Huntington’sModel–NT-HTT-128Q_vs_NegativeControl(driver)_Age_14_DESeq_all_genes.txt')
rownames(NT_HTT_D14_file) <- NT_HTT_D14_file[,1]
NT_HTT_D14_file$Gene_id = NULL
Log2FC_NT_HTT_D14 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_NT_HTT_D14 <- c(Log2FC_NT_HTT_D14, NT_HTT_D14_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_NT_HTT_D14` <- Log2FC_NT_HTT_D14

FL_HTT_D2_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_62_Huntington’sModel–FL-HTT-200Q_vs_NegativeControl(driver)_Age_2_DESeq_all_genes.txt')
rownames(FL_HTT_D2_file) <- FL_HTT_D2_file[,1]
FL_HTT_D2_file$Gene_id = NULL
Log2FC_FL_HTT_D2 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_FL_HTT_D2 <- c(Log2FC_FL_HTT_D2, FL_HTT_D2_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_FL_HTT_D2` <- Log2FC_FL_HTT_D2

FL_HTT_D5_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_64_Huntington’sModel–FL-HTT-200Q_vs_NegativeControl(driver)_Age_5_DESeq_all_genes.txt')
rownames(FL_HTT_D5_file) <- FL_HTT_D5_file[,1]
FL_HTT_D5_file$Gene_id = NULL
Log2FC_FL_HTT_D5 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_FL_HTT_D5 <- c(Log2FC_FL_HTT_D5, FL_HTT_D5_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_FL_HTT_D5` <- Log2FC_FL_HTT_D5

FL_HTT_D7_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_66_Huntington’sModel–FL-HTT-200Q_vs_NegativeControl(driver)_Age_7_DESeq_all_genes.txt')
rownames(FL_HTT_D7_file) <- FL_HTT_D7_file[,1]
FL_HTT_D7_file$Gene_id = NULL
Log2FC_FL_HTT_D7 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_FL_HTT_D7 <- c(Log2FC_FL_HTT_D7, FL_HTT_D7_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_FL_HTT_D7` <- Log2FC_FL_HTT_D7

FL_HTT_D10_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_68_Huntington’sModel–FL-HTT-200Q_vs_NegativeControl(driver)_Age_10_DESeq_all_genes.txt')
rownames(FL_HTT_D10_file) <- FL_HTT_D10_file[,1]
FL_HTT_D10_file$Gene_id = NULL
Log2FC_FL_HTT_D10 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_FL_HTT_D10 <- c(Log2FC_FL_HTT_D10, FL_HTT_D10_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_FL_HTT_D10` <- Log2FC_FL_HTT_D10

FL_HTT_D14_file = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/DESeq/contrast_70_Huntington’sModel–FL-HTT-200Q_vs_NegativeControl(driver)_Age_14_DESeq_all_genes.txt')
rownames(FL_HTT_D14_file) <- FL_HTT_D14_file[,1]
FL_HTT_D14_file$Gene_id = NULL
Log2FC_FL_HTT_D14 <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  Log2FC_FL_HTT_D14 <- c(Log2FC_FL_HTT_D14, FL_HTT_D14_file[i,'log2FoldChange'])
}
intersectedModifierData$`Log2FC_FL_HTT_D14` <- Log2FC_FL_HTT_D14

intersectedModifierData = intersectedModifierData %>%relocate(Interpro.Annotation)
intersectedModifierData = intersectedModifierData %>%relocate(`Modifier Construct`, `Cluster Membership`, .before = `elav>B42`)

#Calculate DEG Direction per Disease Model
DEG_Direction_AD_B42 <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  sum = as.double(intersectedModifierData[[i, 'Log2FC_AD_B42_D2']]) + as.double(intersectedModifierData[[i, 'Log2FC_AD_B42_D5']]) + 
    as.double(intersectedModifierData[[i, 'Log2FC_AD_B42_D7']]) + as.double(intersectedModifierData[[i, 'Log2FC_AD_B42_D10']]) +
    as.double(intersectedModifierData[[i, 'Log2FC_AD_B42_D14']])
  mean = sum / 5
  if (mean > 0 ){
    DEG_Direction_AD_B42 = c(DEG_Direction_AD_B42, "UP DEG")
  } else {
    DEG_Direction_AD_B42 = c(DEG_Direction_AD_B42, "DOWN DEG")
  }
  i = i + 1
}
intersectedModifierData$`DEG_Direction_AD_B42` <- DEG_Direction_AD_B42
intersectedModifierData = intersectedModifierData %>%relocate(`DEG_Direction_AD_B42`, .before = `elav>B42`)

DEG_Direction_AD_Tau <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  sum = as.double(intersectedModifierData[[i, 'Log2FC_AD_Tau_D2']]) + as.double(intersectedModifierData[[i, 'Log2FC_AD_Tau_D5']]) + 
    as.double(intersectedModifierData[[i, 'Log2FC_AD_Tau_D7']]) + as.double(intersectedModifierData[[i, 'Log2FC_AD_Tau_D10']]) +
    as.double(intersectedModifierData[[i, 'Log2FC_AD_Tau_D14']])
  mean = sum / 5
  if (mean > 0 ){
    DEG_Direction_AD_Tau = c(DEG_Direction_AD_Tau, "UP DEG")
  } else {
    DEG_Direction_AD_Tau = c(DEG_Direction_AD_Tau, "DOWN DEG")
  }
  i = i + 1
}
intersectedModifierData$`DEG_Direction_AD_Tau` <- DEG_Direction_AD_Tau
intersectedModifierData = intersectedModifierData %>%relocate(`DEG_Direction_AD_Tau`, .before = `elav>Tau`)

DEG_Direction_PD <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  sum = as.double(intersectedModifierData[[i, 'Log2FC_Parkinson_D2']]) + as.double(intersectedModifierData[[i, 'Log2FC_Parkinson_D5']]) + 
    as.double(intersectedModifierData[[i, 'Log2FC_Parkinson_D7']]) + as.double(intersectedModifierData[[i, 'Log2FC_Parkinson_D10']]) +
    as.double(intersectedModifierData[[i, 'Log2FC_Parkinson_D14']])
  mean = sum / 5
  if (mean > 0 ){
    DEG_Direction_PD = c(DEG_Direction_PD, "UP DEG")
  } else {
    DEG_Direction_PD = c(DEG_Direction_PD, "DOWN DEG")
  }
  i = i + 1
}
intersectedModifierData$`DEG_Direction_PD` <- DEG_Direction_PD
intersectedModifierData = intersectedModifierData %>%relocate(`DEG_Direction_PD`, .before = `elav>aSyn`)

DEG_Direction_FL_HTT <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  sum = as.double(intersectedModifierData[[i, 'Log2FC_FL_HTT_D2']]) + as.double(intersectedModifierData[[i, 'Log2FC_FL_HTT_D5']]) + 
    as.double(intersectedModifierData[[i, 'Log2FC_FL_HTT_D7']]) + as.double(intersectedModifierData[[i, 'Log2FC_FL_HTT_D10']]) +
    as.double(intersectedModifierData[[i, 'Log2FC_FL_HTT_D14']])
  mean = sum / 5
  if (mean > 0 ){
    DEG_Direction_FL_HTT = c(DEG_Direction_FL_HTT, "UP DEG")
  } else {
    DEG_Direction_FL_HTT = c(DEG_Direction_FL_HTT, "DOWN DEG")
  }
  i = i + 1
}
intersectedModifierData$`DEG_Direction_FL_HTT` <- DEG_Direction_FL_HTT
intersectedModifierData = intersectedModifierData %>%relocate(`DEG_Direction_FL_HTT`, .before = `elav>FL-HTT`)

DEG_Direction_NT_HTT <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  sum = as.double(intersectedModifierData[[i, 'Log2FC_NT_HTT_D2']]) + as.double(intersectedModifierData[[i, 'Log2FC_NT_HTT_D5']]) + 
    as.double(intersectedModifierData[[i, 'Log2FC_NT_HTT_D7']]) + as.double(intersectedModifierData[[i, 'Log2FC_NT_HTT_D10']]) +
    as.double(intersectedModifierData[[i, 'Log2FC_NT_HTT_D14']])
  mean = sum / 5
  if (mean > 0 ){
    DEG_Direction_NT_HTT = c(DEG_Direction_NT_HTT, "UP DEG")
  } else {
    DEG_Direction_NT_HTT = c(DEG_Direction_NT_HTT, "DOWN DEG")
  }
  i = i + 1
}
intersectedModifierData$`DEG_Direction_NT_HTT` <- DEG_Direction_NT_HTT
intersectedModifierData = intersectedModifierData %>%relocate(`DEG_Direction_NT_HTT`, .before = `elav>NT-HTT`)

#Predict whether gene is associated or causal (compensator or pathogenic) change for AD_B42
Gene_Change_NT_HTT <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  if(isNA(intersectedModifierData[[i, 'elav>NT-HTT']])){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Unknown ???")  
  } else if (grepl("UN", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)) {
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Unknown ???")
  } else if (grepl("NM", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Associated Change")
  }
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Compensatory") 
  } 
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } 
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_NT_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>NT-HTT']], fixed = TRUE)){
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "Causal Change - Pathogenic") 
  }
  
  else{
    Gene_Change_NT_HTT = c(Gene_Change_NT_HTT, "EXCEPTION!!!!!!!!") 
  }
  i = i+1
}
intersectedModifierData$`Gene_Change_NT_HTT` <- Gene_Change_NT_HTT
intersectedModifierData = intersectedModifierData %>%relocate(`Gene_Change_NT_HTT`, .after = `elav>NT-HTT`)

Gene_Change_FL_HTT <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  if(isNA(intersectedModifierData[[i, 'elav>FL-HTT']])){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Unknown ???")  
  } else if (grepl("UN", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)) {
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Unknown ???")
  } else if (grepl("NM", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Associated Change")
  }
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Compensatory") 
  } 
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Pathogenic") 
  } 
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_FL_HTT']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>FL-HTT']], fixed = TRUE)){
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "Causal Change - Pathogenic") 
  }
  
  else{
    Gene_Change_FL_HTT = c(Gene_Change_FL_HTT, "EXCEPTION!!!!!!!!") 
  }
  i = i+1
}
intersectedModifierData$`Gene_Change_FL_HTT` <- Gene_Change_FL_HTT
intersectedModifierData = intersectedModifierData %>%relocate(`Gene_Change_FL_HTT`, .after = `elav>FL-HTT`)

Gene_Change_PD <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  if(isNA(intersectedModifierData[[i, 'elav>aSyn']])){
    Gene_Change_PD = c(Gene_Change_PD, "Unknown ???")  
  } else if (grepl("UN", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)) {
    Gene_Change_PD = c(Gene_Change_PD, "Unknown ???")
  } else if (grepl("NM", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Associated Change")
  }
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Compensatory") 
  } 
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Pathogenic") 
  } 
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_PD']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>aSyn']], fixed = TRUE)){
    Gene_Change_PD = c(Gene_Change_PD, "Causal Change - Pathogenic") 
  }
  
  else{
    Gene_Change_PD = c(Gene_Change_PD, "EXCEPTION!!!!!!!!") 
  }
  i = i+1
}
intersectedModifierData$`Gene_Change_PD` <- Gene_Change_PD
intersectedModifierData = intersectedModifierData %>%relocate(`Gene_Change_PD`, .after = `elav>aSyn`)


Gene_Change_AD_Tau <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  if(isNA(intersectedModifierData[[i, 'elav>Tau']])){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Unknown ???")  
  } else if (grepl("UN", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)) {
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Unknown ???")
  } else if (grepl("NM", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Associated Change")
  }
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Compensatory") 
  } 
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Pathogenic") 
  } 
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_Tau']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>Tau']], fixed = TRUE)){
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "Causal Change - Pathogenic") 
  }
  
  else{
    Gene_Change_AD_Tau = c(Gene_Change_AD_Tau, "EXCEPTION!!!!!!!!") 
  }
  i = i+1
}
intersectedModifierData$`Gene_Change_AD_Tau` <- Gene_Change_AD_Tau
intersectedModifierData = intersectedModifierData %>%relocate(`Gene_Change_AD_Tau`, .after = `elav>Tau`)

Gene_Change_AD_B42 <- vector(mode = "list")
i = 1
while (i <= length(intersectedModifierData)+1)
{
  if(isNA(intersectedModifierData[[i, 'elav>B42']])){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Unknown ???")  
  } else if (grepl("UN", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)) {
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Unknown ???")
  } else if (grepl("NM", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Associated Change")
  }
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("S", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Compensatory") 
  } 
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("E", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Pathogenic") 
  }
  
  else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Compensatory") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Pathogenic") 
  } else if (grepl("UP", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Pathogenic") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("OE", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("DUP", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Compensatory") 
  } else if (grepl("DOWN", intersectedModifierData[[i, 'DEG_Direction_AD_B42']], fixed = TRUE) && grepl("RNAi", intersectedModifierData[[i, 'Modifier Construct']], fixed = TRUE) && grepl("L", intersectedModifierData[[i, 'elav>B42']], fixed = TRUE)){
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "Causal Change - Pathogenic") 
  }
  
  else{
    Gene_Change_AD_B42 = c(Gene_Change_AD_B42, "EXCEPTION!!!!!!!!") 
  }
  i = i+1
}
intersectedModifierData$`Gene_Change_AD_B42` <- Gene_Change_AD_B42
intersectedModifierData = intersectedModifierData %>%relocate(`Gene_Change_AD_B42`, .after = `elav>B42`)

#Export Results
intersectedModifierData <- apply(intersectedModifierData,2,as.character)
write.csv(intersectedModifierData, file="/Users/ven/Documents/PhD/Rotation Botas Lab/Data/Analyzed_Intersected_Modifier_Dataset.csv", row.names = FALSE)

#Plot Interesting Genes

FBgn0002940 <- normalized_counts[rownames(normalized_counts)=="FBgn0002940", ]
PlotGroups (FBgn0002940, edesign = edesign)

FBgn0003137 <- normalized_counts[rownames(normalized_counts)=="FBgn0003137", ]
PlotGroups (FBgn0003137, edesign = edesign)

FBgn0033494 <- normalized_counts[rownames(normalized_counts)=="FBgn0033494", ]
PlotGroups (FBgn0033494, edesign = edesign)

FBgn0033763 <- normalized_counts[rownames(normalized_counts)=="FBgn0033763", ]
PlotGroups (FBgn0033763, edesign = edesign)

FBgn0039645 <- normalized_counts[rownames(normalized_counts)=="FBgn0039645", ]
PlotGroups (FBgn0039645, edesign = edesign)

FBgn0040524 <- normalized_counts[rownames(normalized_counts)=="FBgn0040524", ]
PlotGroups (FBgn0040524, edesign = edesign)

FBgn0051345 <- normalized_counts[rownames(normalized_counts)=="FBgn0051345", ]
PlotGroups (FBgn0051345, edesign = edesign)

FBgn0086347 <- normalized_counts[rownames(normalized_counts)=="FBgn0086347", ]
PlotGroups (FBgn0086347, edesign = edesign)

FBgn0260632 <- normalized_counts[rownames(normalized_counts)=="FBgn0260632", ]
PlotGroups (FBgn0260632, edesign = edesign)

FBgn0261612<- normalized_counts[rownames(normalized_counts)=="FBgn0261612", ]
PlotGroups (FBgn0261612, edesign = edesign)
