library(readxl)
library(tidyverse)
library(maSigPro)

#GO Analysis= https://maayanlab.cloud/Enrichr/enrich#

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
intersectedModifierData$`gmr>Tau`= NULL
intersectedModifierData$`repo>FL-HTT`= NULL
intersectedModifierData$`repo>NT-HTT`= NULL
intersectedModifierData$`Any Human Homolog`= NULL
intersectedModifierData$Gene.Symbol= NULL
intersectedModifierData$Stock.Number= NULL
intersectedModifierData$Allele.Class= NULL
#intersectedModifierData = intersectedModifierData[,c(1,2,9,10,3,4,5,6,7,8)]

#remember to make results object independently available to this script!
cluster_membership <- vector(mode = "list")
for (i in intersectedModifierData[["Fly.Base.GeneID"]]){
  cluster_membership <- rbind (cluster_membership, results[["cut"]][[i]])
}
intersectedModifierData$`Cluster Membership` <- as.numeric(cluster_membership)


FBgn0000044 <- normalized_counts[rownames(normalized_counts)=="FBgn0000044",]
PlotGroups(FBgn0000044, edesign = edesign, summary.mode = "representative")

#parse out Log2FC values per condition per day and add to intersectedModifierData

intersected_modifer_df <- data.frame("FB Gene ID" = intersectedModifierData[["Fly.Base.GeneID"]], 
                                     "Drosophila Gene" = intersectedModifierData[["Drosophila.Gene"]],
                                     "Cluster Membership" = cluster_membership,
                                     "Log2FC HD200Q Day 2",
                                     "Log2FC HD200Q Day 5",
                                     "Log2FC HD200Q Day 7",
                                     "Log2FC HD200Q Day 10",
                                     "Log2FC HD200Q Day 14",
                                     "Log2FC HD128Q Day 2",
                                     "Log2FC HD128Q Day 5",
                                     "Log2FC HD128Q Day 7",
                                     "Log2FC HD128Q Day 10",
                                     "Log2FC HD128Q Day 14",
                                     "Log2FC PD Day 2",
                                     "Log2FC PD Day 5",
                                     "Log2FC PD Day 7",
                                     "Log2FC PD Day 10",
                                     "Log2FC PD Day 14",
                                     "Log2FC AD-Tau Day 2",
                                     "Log2FC AD-Tau Day 5",
                                     "Log2FC AD-Tau Day 7",
                                     "Log2FC AD-Tau Day 10",
                                     "Log2FC AD-Tau Day 14",
                                     "Log2FC AD-B42 Day 2",
                                     "Log2FC AD-B42 Day 5",
                                     "Log2FC AD-B42 Day 7",
                                     "Log2FC AD-B42 Day 10",
                                     "Log2FC AD-B42 Day 14",
                                     "Modifier Construct+Stock Number",
                                     "Modifier Direction - HD",
                                     "Modifier Direction - PD",
                                     "Modifier Direction - Alz")



