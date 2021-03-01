#Load required packages
library(DESeq2)
library(readxl)

#Initialize annotation matrix and format w/ correct Sample ID naming + column data (metafile)
status = read_excel('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/results_Sep2120_cmp_same_age/Ismael-QC-SHU_14219_B01_EXS_RNA_batchIdentifiers.xlsx')
status[status=="Negative Control (driver)"] = 'Neg_Driver'
status[status=="Negative Control (no driver)"] = 'Neg_No_Driver'
status[status=="AD-model – Secreted β42" ] = "AD_B42"
status[status=="AD-model - Tauopathy"] = "AD_Tau"
status[status=="Parkinson’s model - Synuclein"] = 'PD'
status[status=="Huntington’s Model – NT-HTT-128Q"] = "HD_128Q"
status[status=="Huntington’s Model – FL-HTT-200Q"] = "HD_200Q"
status$"SampleID" <- paste(status$`Comments (Genotype)`,status$`Age (days)`)
metafile <- na.omit(status) #remove NAs - there are duplicate IDs from failed sequencing

#Initialize and format count matrix 
counts = read.delim('/Users/ven/Documents/PhD/Rotation Botas Lab/Data/results_Sep2120_cmp_same_age/featureCounts_count_matrix.txt')
rownames(counts) = counts$Gene
counts$Gene = NULL

#Create DESeq Data Set from count matrix
data <- DESeqDataSetFromMatrix(countData = counts,colData = metafile, design = ~SampleID)

#Normalize data
data <- estimateSizeFactors(data)
normalized_counts <- counts(data, normalized=TRUE)

#Export normalized data
write.table(normalized_counts, file="/Users/ven/Documents/PhD/Rotation Botas Lab/Data/normalized_counts.txt", sep="\t", quote=F, col.names=NA)