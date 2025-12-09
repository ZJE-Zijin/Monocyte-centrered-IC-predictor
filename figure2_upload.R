
# conda activate seurat5-r43
library(tidyr)
library(Seurat) # install.packages('https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_5.0.2.tar.gz',repos=NULL)
library(dplyr)
library(patchwork)
library(ggsignif)
library(ggplot2)
setwd('../figure2_output/cellchat')
monocyte_all_scanvi <- LoadH5Seurat("../figure2_output/monocytes_all.h5seurat")

# Read umap coordinates
umap_coords <- read.csv("anndata_scaled_mde2.csv", row.names = 1)
class(monocyte_all_scanvi@reductions)
names(monocyte_all_scanvi@reductions)
head(monocyte_all_scanvi)
class(monocyte_all_scanvi@tools)

monocyte_all_scanvi[["umap"]] <- CreateDimReducObject(
  embeddings = as.matrix(umap_coords),
  key = "UMAP_",  # Key must match column name prefix (e.g., "UMAP_1" → "UMAP_")
  assay = "RNA"    # Specify the corresponding Assay (usually "RNA")
)


# Read Leiden clustering results (ensure cell IDs match Seurat)
leiden_df <- read.csv("leiden_clusters2.csv", row.names = "cell_id")

# Check if cell IDs are consistent
if (!all(rownames(leiden_df) %in% colnames(monocyte_all_scanvi))) {
  stop("Cell IDs do not match! Please check data source.")
}

# Add Leiden clusters to Seurat's meta.data
monocyte_all_scanvi@meta.data$leiden_clusters2 <- leiden_df[colnames(monocyte_all_scanvi), "leiden"]

# Convert to factor variable (ensure clusters are treated as categorical)
monocyte_all_scanvi@meta.data$leiden_clusters2 <- factor(monocyte_all_scanvi@meta.data$leiden_clusters2)

# Verify results
head(monocyte_all_scanvi@meta.data)

# Read Leiden clustering results (ensure cell IDs match Seurat)
leiden_df <- read.csv("leiden_clusters15.csv", row.names = "cell_id")

# Check if cell IDs are consistent
if (!all(rownames(leiden_df) %in% colnames(monocyte_all_scanvi))) {
  stop("Cell IDs do not match! Please check data source.")
}

# Add Leiden clusters to Seurat's meta.data
monocyte_all_scanvi@meta.data$leiden_clusters15 <- leiden_df[colnames(monocyte_all_scanvi), "leiden"]

# Convert to factor variable (ensure clusters are treated as categorical)
monocyte_all_scanvi@meta.data$leiden_clusters15 <- factor(monocyte_all_scanvi@meta.data$leiden_clusters15)


Idents(monocyte_all_scanvi)<-'leiden_clusters15'
monocyte_all_scanvi<-RenameIdents(monocyte_all_scanvi,'0'='non-classical',
'14'='intermediated','9'='intermediated','8'='intermediated','11'='intermediated','23'='intermediated'
)
monocyte_all_scanvi[['leiden_annotation']]<-Idents(monocyte_all_scanvi)
tmp_id<-colnames(subset(monocyte_all_scanvi,idents=c('non-classical','intermediated'),invert=TRUE))
monocyte_all_scanvi<-SetIdent(monocyte_all_scanvi,cells=tmp_id,value='classical')
monocyte_all_scanvi[['leiden_annotation']]<-Idents(monocyte_all_scanvi)

Idents(monocyte_all_scanvi)<-'treat_result'
monocyte_all_scanvi<-RenameIdents(monocyte_all_scanvi,'CPR'='pCR/MPR','MPR'='pCR/MPR',
'CPR-after'='pCR/MPR-after','MPR-after'='pCR/MPR-after','non-MPR'='NMPR',
'non-MPR-after'='NMPR-after')
monocyte_all_scanvi[['treat_result_coarse']]<-Idents(monocyte_all_scanvi)


for(i in unique(monocyte_all_scanvi$treat_result_coarse)){
    for(j in unique(monocyte_all_scanvi$sample)){
        cell_id<-colnames(subset(monocyte_all_scanvi,idents=i,sample==j))
        monocyte_all_scanvi<-SetIdent(monocyte_all_scanvi,cells=cell_id,value=paste0(i,'_',j))
    }
}
monocyte_all_scanvi[['treat_result_coarse_sample']]<-Idents(monocyte_all_scanvi)



tmp<-monocyte_all_scanvi
Idents(tmp)<-'leiden_clusters15'
tmp<-subset(tmp,idents=c('15','18'),invert=TRUE)
Idents(tmp)<-'treat_result_coarse_sample'
tmp$treat_result_coarse_sample<-factor(tmp$treat_result_coarse_sample,levels=c(
    "pCR/MPR_blood","NMPR_blood","pCR/MPR-after_blood","NMPR-after_blood",
    "pCR/MPR_tumor","NMPR_tumor","pCR/MPR-after_tumor","NMPR-after_tumor"
))
pdf('monocyte_all_scanvi_dimplot_leiden_annotation_nolabel_figure2A.pdf',12,6)
tmp_plot<-DimPlot(tmp,pt.size=1,
        reduction = paste0('umap'),
        group.by = c('leiden_annotation'),ncol=4,split.by = 'treat_result_coarse_sample',
        combine = FALSE, label.size = 5,label=FALSE,repel=TRUE,raster=FALSE
        )
print(tmp_plot)
dev.off()
saveRDS(tmp,'monocyte_all_scanvi.rds')




monocyte_all_scanvi_formal<-readRDS('../figure2_output/monocyte_all_scanvi_formal.rds')

Idents(monocyte_all_scanvi_formal)<-'leiden_clusters'
monocyte_all_scanvi_formal<-RenameIdents(monocyte_all_scanvi_formal,'0'='classical',
'1'='classical','2'='classical','3'='non-classical','4'='classical',
'5'='intermediated','6'='intermediated',
'7'='intermediated','8'='classical','9'='classical','10'='classical',
'11'='non-classical','14'='classical','15'='classical'
)
monocyte_all_scanvi_formal[['leiden_annotation']]<-Idents(monocyte_all_scanvi_formal)

Idents(monocyte_all_scanvi_formal)<-'treat_result'
monocyte_all_scanvi_formal<-RenameIdents(monocyte_all_scanvi_formal,'CPR'='pCR/MPR-native','MPR'='pCR/MPR-native',
'CPR-after'='pCR/MPR-treated','MPR-after'='pCR/MPR-treated','non-MPR'='NMPR-native',
'non-MPR-after'='NMPR-treated')
monocyte_all_scanvi_formal[['treat_result_coarse']]<-Idents(monocyte_all_scanvi_formal)


for(i in unique(monocyte_all_scanvi_formal$treat_result_coarse)){
    for(j in unique(monocyte_all_scanvi_formal$sample)){
        cell_id<-colnames(subset(monocyte_all_scanvi_formal,idents=i,sample==j))
        monocyte_all_scanvi_formal<-SetIdent(monocyte_all_scanvi_formal,cells=cell_id,value=paste0(i,'_',j))
    }
}
monocyte_all_scanvi_formal[['treat_result_coarse_sample']]<-Idents(monocyte_all_scanvi_formal)

monocyte_all_scanvi_formal$treat_result_coarse_sample<-factor(monocyte_all_scanvi_formal$treat_result_coarse_sample,
levels=c(paste0(c('pCR/MPR-native','NMPR-native','pCR/MPR-treated','NMPR-treated'),'_','blood'),
paste0(c('pCR/MPR-native','NMPR-native','pCR/MPR-treated','NMPR-treated'),'_','tumor')))

pdf('monocyte_all_scanvi_formal_featureplot_intermediated_coexpression_figure2B.pdf',17,6)
p1<-FeaturePlot(subset(monocyte_all_scanvi_formal,condition=='native'),
cols=c("black", "#ff0000", "#00ff00"),
 features = c('CD14','FCGR3A'),
 blend=TRUE,
# split.by='treat_result_coarse_sample',
# pt.size=0.1,
reduction="umap",label=F,raster=FALSE,ncol=4,max.cutoff=3)
print(p1)
dev.off()


ident1<-'leiden_annotation'
ident2<-'orig.ident'
ident3<-'treat_result_coarse_sample'
ident1_value<-unique(monocyte_all_scanvi_formal@meta.data[[ident1]])
Idents(monocyte_all_scanvi_formal)<-ident1
ratio03 <- data.frame(table(monocyte_all_scanvi_formal@meta.data[[ident1]],monocyte_all_scanvi_formal@meta.data[[ident2]],monocyte_all_scanvi_formal@meta.data[[ident3]]),row.names = NULL)
sum<-data.frame(table(monocyte_all_scanvi_formal@meta.data[[ident2]]))
ratio03<-ratio03[ratio03$Freq!=0,]
Freq_nor <- c()
for (i in (1:nrow(ratio03))){
  sumi<-sum[sum$Var1==ratio03[i,'Var2'],'Freq']
  Freq_nor <- c(Freq_nor, ratio03[i, 'Freq']/sumi)
}
ratio03 <- cbind(ratio03, Freq_nor)



tmp_function<-function(ident1,ident2){
        freq1<-ratio_tmp[ratio_tmp$Var3==ident1,'Freq_nor']
        freq2<-ratio_tmp[ratio_tmp$Var3==ident2,'Freq_nor']
        p_value<-wilcox.test(freq1,freq2,)$p.value
}
p_df<-data.frame(matrix(ncol=4,nrow=length(ident1_value)))
count=0
for(i in ident1_value){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    p_df[count,1]<-tmp_function("pCR/MPR-native_blood","NMPR-native_blood")
    p_df[count,2]<-tmp_function("pCR/MPR-treated_blood","NMPR-treated_blood")
    p_df[count,3]<-tmp_function("pCR/MPR-native_blood","pCR/MPR-treated_blood")
    p_df[count,4]<-tmp_function("NMPR-native_blood","NMPR-treated_blood")
}

colnames(p_df)<-c('Blood native pCR/MPR vs. NMPR','Blood treated pCR/MPR vs. NMPR',
'Blood pCR/MPR native vs. treated','Blood NMPR native vs. treated')
rownames(p_df)<-ident1_value

write.csv(p_df,paste0(ident1,'_blood_pvalue_table_figure1C.csv'),quote=FALSE)

p_df<-data.frame(matrix(ncol=4,nrow=length(ident1_value)))
count=0
for(i in ident1_value){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    p_df[count,1]<-tmp_function("pCR/MPR-native_tumor","NMPR-native_tumor")
    p_df[count,2]<-tmp_function("pCR/MPR-treated_tumor","NMPR-treated_tumor")
    p_df[count,3]<-tmp_function("pCR/MPR-native_tumor","pCR/MPR-treated_tumor")
    p_df[count,4]<-tmp_function("NMPR-native_tumor","NMPR-treated_tumor")
}

colnames(p_df)<-c('Tumor native pCR/MPR vs. NMPR','Tumor treated pCR/MPR vs. NMPR',
'Tumor pCR/MPR native vs. treated','Tumor NMPR native vs. treated')
rownames(p_df)<-ident1_value

write.csv(p_df,paste0(ident1,'_tumor_pvalue_table_figure1C.csv'),quote=FALSE)







library(tidyr)
library(Seurat) # install.packages('https://cran.r-project.org/src/contrib/Archive/Seurat/Seurat_5.0.2.tar.gz',repos=NULL)
library(dplyr)
library(patchwork)
library("batchelor")
library(reticulate)
library(SeuratWrappers)
library(ggplot2) #Plot_Density_Joint_Only functionality is currently restricted to ggplot v3.4.4 or lower.
# devtools::install_version("ggplot2", version = "3.5.0",
#  repos = "http://cran.us.r-project.org")
# library(ggpubr)
library(SeuratDisk)
library(rhdf5)

count<-0
for(i in c('_blood','-after_blood','_tumor','-after_tumor')){
    count=count+1
    blood_native_cellchat_cmpr<-readRDS(paste0('CMPR',i,'_bt_all_annotation1_mono_cellchat.rds'))
    blood_native_cellchat_nmpr<-readRDS(paste0('NMPR',i,'_bt_all_annotation1_mono_cellchat.rds'))
    object.list <- list(CMPR=blood_native_cellchat_cmpr, NMPR=blood_native_cellchat_nmpr)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    assign(paste0('p',as.character(count)),netVisual_heatmap(cellchat,comparison = c(2, 1), measure = "weight"))
}

pdf('netVisual_heatmap_cmpr_nmpr_mono_blood_figureS6DE.pdf',7,4)
p1+p2
dev.off()

pdf('netVisual_heatmap_cmpr_nmpr_mono_tumor_figureS6DE.pdf',7,4)
p3+p4
dev.off()
