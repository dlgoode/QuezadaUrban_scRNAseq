module load pySCENIC
library(SCopeLoomR)
library(Seurat)
library(tidyverse)
#### This code chunk used for processing the SCENIC output files..
loom<-open_loom(file.path="integratedALL_aucell_mtx_step3_withAR.loom",mode="r+")
regulonsMat<-get_regulons(loom,column.attr.name="Regulons")
regulons <- SCENIC::regulonsToGeneLists(regulonsMat)
regulons_AUC <- get_regulons_AUC(loom,column.attr.name="RegulonsAUC")
integratedALL_subset_withAR_regulons<-regulonsMat[which(rowSums(regulonsMat)>20),]
step2_regulons<-read.csv("integratedALL_CTX_regulons_Step2_withAR.csv",sep="")
integratedALL_subset_withAR_regulons<-regulonsMat[which(rowSums(regulonsMat)>20),]
integratedALL_subset_withAR_regulons_set<-integratedALL_subset_withAR_regulons
row.names(integratedALL_subset_withAR_regulons_set)<-str_replace(rownames(head(integratedALL_subset_withAR_regulons)),"\\(\\+\\)", "")
names(regulons_AUC)<-str_replace(names(regulons_AUC),"\\(\\+\\)", "")
#################################################
cellInfo<-readRDS("../integratedALLcellInfo.rds")
#For generating heatmaps directly use the below code along with this data file.
regulons_AUC_integrated_ALL_withAR<-readRDS("regulons_AUC_integrated_ALL_withAR.rds")
regulons_filtered<-read.delim("integrated_ALL_regulons_filtered_list.txt",header=FALSE)
regulons_AUC_integrated_ALL_withAR_subset<-regulons_AUC_integrated_ALL_withAR[regulons_filtered$V1,]
regulonActivity_byclusters <- sapply(split(rownames(cellInfo), cellInfo$clusters),
                                     function(cells) rowMeans(getAUC(regulons_AUC_integrated_ALL_withAR_subset)[,cells]))
regulonActivity_byclusters_Scaled_subset <- t(scale(t(regulonActivity_byclusters), center = T, scale=T))
ComplexHeatmap::Heatmap(regulonActivity_byclusters_Scaled_subset, name="Regulon activity",row_names_gp = gpar(fontsize =5),column_names_gp = gpar(fontsize = 8),width=1500,height=1500)
dev.off()