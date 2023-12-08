#use this code to infer CNV subclones in PDX samples using single-cell RNASeq data
#this code is modified (for original refer:https://github.com/gabrielakinker/CCLE_heterogeneity) to fit our PDX samples
#module load R
#R version 4.0.2 (2020-06-22) -- "Taking Off Again"
library(Seurat)
library(mclust)
library(tidyverse)
pdxsample_obj<-readRDS("PDXseuratObjectfile.rds")
#####using entire genes
pdxsample_obj_umicount<-as.matrix(pdxsample_obj[['SCT']]@counts)
# read gene locus annotation
gene_locus <- readRDS("CCLE_scRNAseq_github/CCLE_heterogeneity_Rfiles/gene_locus.RDS")# Refer original sub clonal anlaysis code from Kinker et al.,
#### select common genes between datasets
common_genes<-Reduce(intersect,c(lapply(list(pdxsample_obj_umicount), rownames), list(gene_locus$HGNC.symbol)))
gene_locus <- gene_locus[match(common_genes, gene_locus$HGNC.symbol),]

pdxsample_obj_umicount_filtered <- pdxsample_obj_umicount[common_genes,]
# select top expressed genes
ave_expr<- apply(pdxsample_obj_umicount_filtered,1,mean)#this code also works and gives same output                
pdxsample_obj_umicount_filtered<-pdxsample_obj_umicount_filtered[order(ave_expr,decreasing=TRUE),]
gene_locus <- gene_locus[match(row.names(pdxsample_obj_umicount_filtered), gene_locus$HGNC.symbol),]                   

# process data   
#expr_ccle <- lapply(expr_ccle, function(x) log2((x/10) + 1))
pdxsample_obj_umicount_filtered_log<-log2((pdxsample_obj_umicount_filtered/10) + 1)
#expr_ccle <- lapply(expr_ccle, function(x) t(t(x)-colMeans(x)))
pdxsample_obj_umicount_filtered_log<-t(t(pdxsample_obj_umicount_filtered_log)-colMeans(pdxsample_obj_umicount_filtered_log))
#ave_expr_log <- rowMeans(sapply(expr_ccle, rowMeans))
ave_expr_log<-rowMeans(pdxsample_obj_umicount_filtered_log)
#expr_ccle <- lapply(expr_ccle, function(x) x - ave_expr_log)
pdxsample_obj_umicount_filtered_log<-pdxsample_obj_umicount_filtered_log-ave_expr_log
# **************************************************************************
# Infer large scale copy number aberrations in chromosome arms

# truncate expression data
for(i in names(expr_ccle)) {
  expr_ccle[[i]][expr_ccle[[i]] > 3] <- 3
  expr_ccle[[i]][expr_ccle[[i]] < -3] <- -3
}
hist(colSums(pdxsample_obj_umicount_filtered_log))                    
dev.off()
#******************************************** Didn't use this step()
# genes by chromosome location                     
gene_order <- gene_locus$HGNC.symbol[order(gene_locus$Chromosome.scaffold.name, gene_locus$Karyotype.band)]      
# calculate running average expression in windows of 100 genes 
a <- pdxsample_obj_umicount_filtered_log[gene_order,]
b <- data.frame(matrix(ncol = nrow(a) - 99, nrow = ncol(a)), row.names = colnames(a))
for(j in 1:ncol(b)) {
    b[,j] <- colMeans(a[j:(j+99),])
  }
cna_infer<-b
# determine the limits of each chromosome arm in the inferred cna matrix
chr_arms_size <- table(gene_locus$chr.arm) # get the number of genes in each chr arm 
chr_arms_size <- chr_arms_size[order(as.numeric(gsub('.{1}$', '', names(chr_arms_size))), names(chr_arms_size))] # order by chr location

window_vs_arm <- c(rep(names(chr_arms_size)[1],chr_arms_size[1]-50)) # annotate the chormosome arm of each window in the inferred cna matrix
for(i in 2:length(chr_arms_size)) {
  window_vs_arm <- c(window_vs_arm, rep(names(chr_arms_size)[i], (chr_arms_size)[i]))
}
window_vs_arm <- window_vs_arm[1:<add specific number based on the PDX sample here>]
# average the inferred cna matrix by chromossome arm                                        
#cna_infer_arm <- lapply(cna_infer, function(x) data.frame(aggregate(t(x), by=list(window_vs_arm), mean), row.names = 1))
cna_infer_arm<-data.frame(aggregate(t(cna_infer), by = list(window_vs_arm), FUN = "mean"), row.names = 1)
# **************************************************************************
# Detection of arm-level CNA subclones                 
                        
# fit a bimodal Gaussian mixture (via EM algorithm) for each arm in each cell line
#gmm_test <- lapply(cna_infer_arm, function(x) apply(x, 1, function(y) Mclust(y, G=1:2)))
gmm_test<-apply(cna_infer_arm, 1, function(y) Mclust(y, G=1:2))
# confindently define cell lines that have subclones (i.e. those with > 20 cells classified into a second mode with > 99% confidence)                                                   
####cna_subclones <- lapply(gmm_test, function(x) lapply(x, function(y) data.frame("class"=y$classification, "uncer" = y$uncertainty)))
cna_subclones<-lapply(gmm_test, function(y) data.frame("class"=y$classification, "uncer" = y$uncertainty))
###for(i in names(cna_subclones)) {
for(j in names(cna_subclones)) {
    a <- cna_subclones[[j]]
    a <- a[a$uncer<0.01,]
    if(length(which(table(a$class) > 20)) <2) cna_subclones[[j]] <- NULL
  }
#for some of the PDX samples this step will fail if there are no subclones.
cna_subclones_final<-cna_subclones[sapply(cna_subclones, function(x) length(x)!=0)]
# confindently assign cells to clones (only consider cells classified with >90% confidence)
###cna_subclones<-readRDS("pdxsample_obj_all_genes_cna_subclones_final_confidence_results.RDS")
####(pdxsample_obj_umicount_filtered_log) is nothing but expr_ccle
clone_assignment <- 
  if(length(cna_subclones) > 1) { # more than one clone - we considered all combinations of modes with at least 5 cells.
    a <- colnames(pdxsample_obj_umicount_filtered_log)
    b <- do.call(paste, lapply(cna_subclones, function(x) x[["class"]]))
    b <- data.frame("cell"=a, "class"=b, stringsAsFactors = F)
    c <- sapply(cna_subclones, function(x) x[["uncer"]])
    b <- b[apply(c,1,function(x) length(which(x<=0.1))==length(x)),]
    b <- b[is.element(b$class, names(which(table(b$class)>=5))),]
    d <- data.frame(names(sort(table(b$class), decreasing = T)), 1:length(unique(b$class)), stringsAsFactors = F)
    b$final_class <- d[,2][match(b$class, d[,1])]
    b$final_class <- paste("PDXID_clone", "_", b$final_class, sep="")
    clone_assignment <- b
  }
####mapping clone assigned class as meta data to seurat objects
df_clone<-clone_assignment %>% as_tibble() %>% column_to_rownames(var = "cell")
####using only cells/clones only i.e., 1507 cells/clones
pdxsample_obj[,row.names(df_clone)]
subset_pdxsample_obj<-pdxsample_obj[,row.names(df_clone)]
subset_pdxsample_obj_withClones<-AddMetaData(object=subset_pdxsample_obj,metadata=df_clone$final_class, col.name="cna_clone_class")
DimPlot(subset_pdxsample_obj_withClones,reduction="umap",group.by="orig.ident",pt.size=0.2,label=FALSE,label.size=2)
DimPlot(subset_pdxsample_obj_withClones,reduction="umap",split.by="orig.ident",pt.size=0.2,label=TRUE,label.size=2)
DimPlot(subset_pdxsample_obj_withClones,reduction="umap",group.by="seurat_clusters",pt.size=0.2,label=TRUE,label.size=2)
DimPlot(subset_pdxsample_obj_withClones,reduction="umap",group.by="cna_clone_class",pt.size=0.2,label=TRUE,label.size=2)
DimPlot(subset_pdxsample_obj_withClones,reduction="umap",split.by="cna_clone_class",pt.size=0.2,label=FALSE,label.size=2)
dev.off()
DimPlot(subset_pdxsample_obj_withClones,reduction="umap",group.by="Phase",pt.size=0.2,label=TRUE,label.size=2)
DimPlot(subset_pdxsample_obj_withClones,reduction="umap",split.by="Phase",pt.size=0.2,label=TRUE,label.size=2)
dev.off()
subset_cna_infer_arm_df<-t(cna_infer_arm)[row.names(df_clone),]
new_int_obj_clones<-subset_pdxsample_obj_withClones
new_new_int_obj_clones_with_cna_infer_arm<-AddMetaData(object=new_int_obj_clones,metadata=subset_cna_infer_arm_df[,c("7q","21q","8q","9p")],col.name=c("cna_arm_7q","cna_arm_21q","cna_arm_8q","cna_arm_9p"))
VlnPlot(new_new_int_obj_clones_with_cna_infer_arm,features=c("cna_arm_7q","cna_arm_21q","cna_arm_8q","cna_arm_9p"),ncol=4)
VlnPlot(new_new_int_obj_clones_with_cna_infer_arm,features=c("cna_arm_7q","cna_arm_21q"),ncol=2,pt.size=0)
VlnPlot(new_new_int_obj_clones_with_cna_infer_arm,features=c("cna_arm_8q","cna_arm_9p"),ncol=2,pt.size=0)
VlnPlot(new_new_int_obj_clones_with_cna_infer_arm,features=c("cna_arm_7q","cna_arm_21q"),ncol=2,pt.size=0,sort="increasing")
VlnPlot(new_new_int_obj_clones_with_cna_infer_arm,features=c("cna_arm_8q","cna_arm_9p"),ncol=2,pt.size=0,sort="increasing")
#proportion graph
pro.table<-table(Idents(subset_pdxsample_obj_withClones),subset_pdxsample_obj_withClones$cna_clone_class)
pt<-as.data.frame(pro.table)
library(ggplot2)
library(RColorBrewer)
###x-axis cluster with custom green gradient color
ggplot(pt, aes(x = Var1, y = Freq, fill = Var2)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Clusters") +
  ylab("Proportion") +
  scale_fill_manual(values = c("#808000","darkgreen", "#00FF00", "#90EE90",  "#00FA9A", "#20B2AA","green yellow")) +
  theme(legend.title = element_blank())
dev.off()
##x-axis copy number clone class
ggplot(pt, aes(x = Var2, y = Freq, fill = Var1)) +
  theme_bw(base_size = 15) +
  geom_col(position = "fill", width = 0.5) +
  xlab("Clone class") +
  ylab("Proportion") +
  scale_fill_manual(values = brewer.pal(12, "Paired")) +
  theme(legend.title = "Clusters")
dev.off()
###for plotting all chromosomes (for 470B) use this code for all the PDX samples....
cna_infer_arm_pdxID<-cna_infer_arm %>% rownames_to_column(var="chr_arm")
library(reshape2)
cna_infer_arm_pdxID_long=melt(cna_infer_arm_pdxID,id="chr_arm")
str_sort(cna_infer_arm_pdxID$chr_arm,numeric=TRUE)
 [1] "1p"  "1q"  "2p"  "2q"  "3p"  "3q"  "4p"  "4q"  "5p"  "5q"  "6p"  "6q"
[13] "7p"  "7q"  "8p"  "8q"  "9p"  "9q"  "10p" "10q" "11p" "11q" "12p" "12q"
[25] "13q" "14q" "15q" "16p" "16q" "17p" "17q" "18p" "18q" "19p" "19q" "20p"
[37] "20q" "21q" "22q" "23p" "23q"
ggplot(cna_infer_arm_pdxID_long,aes(x=value,y=factor(chr_arm,level=str_sort(chrm_arm,numeric=TRUE))))+geom_violin()+ylab('chrm_arm')
dev.off()