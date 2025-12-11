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
# library(Matrix)
library(scCustomize)
# devtools::install_version("scCustomize", version = "2.0.0",
#  repos = "http://cran.us.r-project.org")
library(viridis)
# library(RColorBrewer)
library(ggsignif)
library(SingleR)
hpca.se <- celldex::HumanPrimaryCellAtlasData()
options(future.globals.maxSize = 3e+09)
set.seed(1)
setwd('../figure1_output')
bt_all<-readRDS('../figure1_ouput/bt_all_formal_sample_id.rds')
bt_all<-NormalizeData(bt_all)%>%
    FindVariableFeatures(nfeatures = 2000)%>%
    ScaleData()
    
bt_all<-RunPCA(bt_all,
  npcs = 50, 
  verbose = FALSE
)
bt_all<-RunUMAP(bt_all,
  reduction = "pca",    
  dims = 1:50              # Select the first 50 Harmony dimensions
)%>%FindNeighbors(
  reduction = "pca",
  dims = 1:50
)%>%FindClusters(
  resolution = 0.1,         # Adjust resolution (higher value = more clusters)
  algorithm = 1             # 1 = Louvain, 2 = Leiden
)%>%RunTSNE(
    reduction = "pca", 
    dims = 1:50)

Idents(bt_all)<-'annotation1'
tmp_seurat<-subset(bt_all,idents=c('Neutrophil','Unkown1','Unkown2','Unkown3'),invert=TRUE)

pdf('figure1B_annotation1_dimplot_ds_fine9_nolabel_noneutro_without_rmbatch_labeled_figureS1A.pdf',width=14,height=7)
p1<-DimPlot(tmp_seurat, 
group.by = 'annotation1',split.by='treat_result_coarse_sample',ncol=4,
reduction="umap",label=T,raster=NULL,repel=T)+
guides(color = guide_legend(ncol = 1, byrow = TRUE, override.aes = list(size=2)))+
theme(legend.title= element_blank(),title = element_blank(),text = element_text(size = 20))+ggtitle("")
print(p1)
dev.off()



bt_all_noneutrophil<-subset(bt_all,annotation1=='Neutrophil',invert=TRUE)


# Generate files required for infercnv
set.seed(1)
library('data.table')
setwd('../figure1_output/infercnv')
for(i in c('tumor','blood')){
  print(i)
tmp<-subset(bt_all_noneutrophil,sample==i)
cellids<-colnames(tmp)
cellids_downsample<-cellids[sample(1:length(cellids),40000,replace=FALSE)]
tmp<-subset(tmp,cells=cellids_downsample)
table(tmp$annotation1)
cell_gene_matrix <- tmp[['RNA']]@counts
cell_gene_matrix.df <- as.data.frame(cell_gene_matrix)
fwrite(x = cell_gene_matrix.df, row.names=TRUE,file = paste0(i,"_cell_gene_matrix.matrix"),sep='\t')
cell_id_annotation <- data.frame(rownames(tmp[[]]),tmp[['annotation1']])
head(cell_id_annotation)
fwrite(x=cell_id_annotation,row.names=F,col.names=F,file=paste0(i,'_cell_id_annotation.txt'),sep='\t')
}


tmp<-bt_all_noneutrophil
tmp$annotation1 <- as.character(tmp$annotation1)
tmp$treat_result_coarse_sample<-as.character(tmp$treat_result_coarse_sample)
tmp$formal_sample_id<-as.character(tmp$formal_sample_id)

SaveH5Seurat(
  tmp,
  filename = "bt_all_noneutrophil_annotation1.h5seurat",
  overwrite = TRUE,
  assay = "RNA",        # Specify RNA Assay to save
  layers = "counts"     # Explicitly save counts layer
)

# SaveH5Seurat(bt_all, filename = "bt_all.h5seurat", overwrite = TRUE)
Convert("bt_all_noneutrophil_annotation1.h5seurat", dest = "h5ad", overwrite = TRUE)

markers1 <- c('AGER','CLDN18',#AT1
'SFTPC','SFTPB','SFTPA1',#AT2
'PIFO','FOXJ1','HYDIN','CFAP299',#Cliated
'SCGB3A1','SCGB3A2',#Club
'CLEC9A','XCR1','CLNK',#cDC1
'CD1C','CLEC10A',#cDC2
'CCR7','CCL22',#DC matue
'IL3RA','CLEC4C',#pDC
'VWF','CDH5','SELE',#Endothelial
'CCL21',#Endothelial lymphatic
'PDGFRA','FAP','COL1A1',#Fibroblast
'MFAP5','SCARA5',#Fibroblast adventitial
'ITGA8','SCN7A',#Fibroblast avlveolar
'APOE','CD5L','MARCO','C1QB','TREM2',#Macrophage
'FABP4',#Macrophage alveolar
'TPSB2',#Mast cell
'MSLN','CALB2',#Mesothlial
'CD14','VCAN','FCN1',#Monocyte
'S100A12',#Monocyte conventional
'LILRB1','LILRB2',#non-Monocyte conventional
'FCGR3B','CSF3R','CXCR2',#Neutrophils
'COX4I2','PDGFRB',#Pericyte
'CD19','CD79A','MS4A1',#B cell
'SDC1','MZB1',#Plasma cell
'TAGLN','MYH11',#Smooth mscule cell
'CD3E',#T
'CD4',#CD4 T
'CD8A',#CD8 T
'FOXP3','IL2RA','CTLA4',#T regulatory
'KLRD1','GNLY',#NK 
'MKI67','CDK1',#dividing
'TP63','KRT5','DSG3','KRT6A',# LUSC
'TTF1'#LUAD
)
Idents(bt_all_noneutrophil)<-'annotation1'
tmp<-subset(bt_all_noneutrophil,idents=c('Unkown1','Unkown2','Unkown3'),invert=TRUE)


pdf('bt_all_noneutrophil_markers1_figureS1B.pdf', width = 10, height = 20) # Appropriately increase size to prevent oversized points
p1<-DotPlot(tmp, features = markers1, group.by = 'annotation1') + coord_flip()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,size=15),text = element_text(size=15))
print(p1)
dev.off()


Idents(bt_all_noneutrophil)<-'annotation1'
tmp<-subset(bt_all_noneutrophil,idents=c('Tumor','Unkown1','Unkown2','Unkown3','Ciliated','Normal epithelial'),sample=='blood',invert=TRUE)
tmp<-subset(tmp,idents=c('Unkown1','Unkown2','Unkown3'),invert=TRUE)

write.csv(table(tmp$annotation1,tmp$treat_result_coarse_sample),'bt_all_noneutrophil_annotation1_table_figure1D.csv',quote=FALSE)





tmp_seurat<-subset(bt_all_noneutrophil,sample=='blood')
ident1<-'annotation1'
ident2<-'orig.ident'
ident3<-'treat_result_coarse_sample'
Idents(tmp_seurat)<-ident1
tmp_seurat<-subset(tmp_seurat,idents=c('Tumor','Unkown1','Unkown2','Unkown3','Ciliated','Normal epithelial'),invert=TRUE)
ratio03 <- data.frame(table(tmp_seurat@meta.data[[ident1]],tmp_seurat@meta.data[[ident2]],tmp_seurat@meta.data[[ident3]]),row.names = NULL)
sum<-data.frame(table(tmp_seurat@meta.data[[ident2]]))
ratio03<-ratio03[ratio03$Freq!=0,]
Freq_nor <- c()
for (i in (1:nrow(ratio03))){
  sumi<-sum[sum$Var1==ratio03[i,'Var2'],'Freq']
  Freq_nor <- c(Freq_nor, ratio03[i, 'Freq']/sumi)
}
ratio03 <- cbind(ratio03, Freq_nor)

# Retain monocyte for clinical cutoff line
ratio_monocyte<-ratio03[ratio03$Var1=='Monocyte',]
saveRDS(ratio_monocyte,'ratio_monocyte.rds')
# 0.36, take high point of NMPR, low point of pCR

tmp_function<-function(ident1,ident2){
        freq1<-ratio_tmp[ratio_tmp$Var3==ident1,'Freq_nor']
        freq2<-ratio_tmp[ratio_tmp$Var3==ident2,'Freq_nor']
        p_value<-wilcox.test(freq1,freq2,)$p.value
}
p_df<-data.frame(matrix(ncol=4,nrow=length(unique(tmp_seurat$annotation1))))
count=0
for(i in unique(tmp_seurat$annotation1)){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    p_df[count,1]<-tmp_function("pCR/MPR-native_blood","NMPR-native_blood")
    p_df[count,2]<-tmp_function("pCR/MPR-treated_blood","NMPR-treated_blood")
    p_df[count,3]<-tmp_function("pCR/MPR-native_blood","pCR/MPR-treated_blood")
    p_df[count,4]<-tmp_function("NMPR-native_blood","NMPR-treated_blood")
}

colnames(p_df)<-c('Blood-native pCR/MPR vs. NMPR','Blood-treated pCR/MPR vs. NMPR',
'Blood pCR/MPR native vs. treated','Blood NMPR native vs. treated')
rownames(p_df)<-unique(tmp_seurat$annotation1)

write.csv(p_df,'annotation1_blood_pvalue_table_figure1E.csv',quote=FALSE)




# Adjust color palette
mycolor<-ggsci::pal_npg("nrc", alpha = 1)(10) #Extract 8 colors, transparency 80%
mycolor
c("#E64B35FF", "#4DBBD5FF", "#00A087FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF", "#DC0000FF","#7E6148FF", "#B09C85FF")


draw_list<-c()
count=0
for(i in unique(tmp_seurat$annotation1)){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    tmp<-ggplot(ratio_tmp,aes(x=Var3,y=Freq_nor,color=Var3))+
    geom_boxplot(outlier.shape = NA)+
    scale_color_manual(values=c("#3C5488FF","#DC0000FF","#00A087FF","#F39B7FFF"))+
    theme_bw()+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "false")+
    ylab("")+xlab("")+ggtitle(i)+
    coord_cartesian(ylim=c(0, 1.0))+

    geom_signif(comparisons = list(c("pCR/MPR-native_blood","NMPR-native_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.55),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-treated_blood","NMPR-treated_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.65),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-native_blood","pCR/MPR-treated_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.75),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("NMPR-native_blood","NMPR-treated_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.85),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)
    assign(paste0('p',count),tmp)
}
print(count)




pdf('Figure1C_blood_annotation1_without_tumor_without_neutrophil_by_treat_result_coarse_split_formal.pdf',width = 9.5,height = 6)
print(tmp)
dev.off()





tmp_seurat<-subset(bt_all_noneutrophil,sample=='tumor')
ident1<-'annotation1'
ident2<-'orig.ident'
ident3<-'treat_result_coarse_sample'
Idents(tmp_seurat)<-ident1
tmp_seurat<-subset(tmp_seurat,idents=c('Unkown1','Unkown2','Unkown3'),invert=TRUE)
cpr_after_id<-unique(subset(tmp_seurat,idents=c('Tumor'),treat_result=='CPR-after')$orig.ident)

ratio03 <- data.frame(table(tmp_seurat@meta.data[[ident1]],tmp_seurat@meta.data[[ident2]],tmp_seurat@meta.data[[ident3]]),row.names = NULL)
sum<-data.frame(table(tmp_seurat@meta.data[[ident2]]))
ratio03<-ratio03[ratio03$Freq!=0,]
ratio03[ratio03$Var2%in%cpr_after_id&ratio03$Var1=='Tumor','Freq']<-0
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
p_df<-data.frame(matrix(ncol=4,nrow=length(unique(tmp_seurat$annotation1))))
count=0
for(i in unique(tmp_seurat$annotation1)){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    p_df[count,1]<-tmp_function("pCR/MPR-native_tumor","NMPR-native_tumor")
    p_df[count,2]<-tmp_function("pCR/MPR-treated_tumor","NMPR-treated_tumor")
    p_df[count,3]<-tmp_function("pCR/MPR-native_tumor","pCR/MPR-treated_tumor")
    p_df[count,4]<-tmp_function("NMPR-native_tumor","NMPR-treated_tumor")
}

colnames(p_df)<-c('Tumor native pCR/MPR vs. NMPR','Tumor treated pCR/MPR vs. NMPR',
'Tumor pCR/MPR native vs. treated','Tumor NMPR native vs. treated')
rownames(p_df)<-unique(tmp_seurat$annotation1)

write.csv(p_df,'annotation1_tumor_pvalue_table_figure1E.csv',quote=FALSE)
for(i in c(1:nrow(p_df))){
  for(j in c(1:ncol(p_df))){
    if(p_df[i,j]<0.05){
      p_df[i,j]<-'TRUE'
    }
  }
}



draw_list<-c()
count=0
for(i in unique(tmp_seurat$annotation1)){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    tmp<-ggplot(ratio_tmp,aes(x=Var3,y=Freq_nor,color=Var3))+
    geom_boxplot(outlier.shape = NA)+
    scale_color_manual(values=c("#3C5488FF","#DC0000FF","#00A087FF","#F39B7FFF"))+
    theme_bw()+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "false")+
    ylab("")+xlab("")+ggtitle(i)+
    coord_cartesian(ylim=c(0, 1.0))+
    geom_signif(comparisons = list(c("pCR/MPR-native_tumor","NMPR-native_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.5),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-treated_tumor","NMPR-treated_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.6),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-native_tumor","pCR/MPR-treated_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.7),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("NMPR-native_tumor","NMPR-treated_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.8),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)
    assign(paste0('p',count),tmp)
}


pdf('Figure1C_tumor_annotation1_with_tumor_without_neutrophil_by_treat_result_coarse_split_formal.pdf',width = 9.50,height = 9.00)
tmp<-wrap_plots(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13),ncol=5)
print(tmp)
dev.off()






bt_ratio<-readRDS('../figure1_output/bt_all_ds_fine9.rds')
bt_all_noneutrophil<-subset(bt_ratio,annotation1=='Neutrophil',invert=TRUE)
tmp_seurat<-subset(bt_all_noneutrophil,sample=='blood')
tmp_seurat<-subset(bt_all_noneutrophil,sample=='blood')
ident1<-'annotation1'
ident2<-'orig.ident'
ident3<-'treat_result_coarse_sample'
Idents(tmp_seurat)<-ident1
tmp_seurat<-subset(tmp_seurat,idents=c('Tumor','Unkown1','Unkown2','Unkown3','Ciliated','Normal epithelial'),invert=TRUE)
ratio03 <- data.frame(table(tmp_seurat@meta.data[[ident1]],tmp_seurat@meta.data[[ident2]],tmp_seurat@meta.data[[ident3]]),row.names = NULL)
sum<-data.frame(table(tmp_seurat@meta.data[[ident2]]))
ratio03<-ratio03[ratio03$Freq!=0,]
Freq_nor <- c()
for (i in (1:nrow(ratio03))){
  sumi<-sum[sum$Var1==ratio03[i,'Var2'],'Freq']
  Freq_nor <- c(Freq_nor, ratio03[i, 'Freq']/sumi)

}
ratio03 <- cbind(ratio03, Freq_nor)


draw_list<-c()
count=0
for(i in unique(tmp_seurat$annotation1)){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    tmp<-ggplot(ratio_tmp,aes(x=Var3,y=Freq_nor,color=Var3))+
    geom_boxplot(outlier.shape = NA)+
    # scale_color_manual(values=c("black", "red"))+
    theme_bw()+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "false")+
    ylab("")+xlab("")+ggtitle(i)+
    coord_cartesian(ylim=c(0, 1.0))+
    geom_signif(comparisons = list(c("pCR/MPR_blood","pCR/MPR-after_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.55),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR_blood","NMPR_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.65),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("NMPR_blood","NMPR-after_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.75),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-after_blood","NMPR-after_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.85),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)
    assign(paste0('p',count),tmp)
}
print(count)


pdf('blood_annotation1_without_tumor_without_neutrophil_by_treat_result_coarse_split_formal_figure1C_figureS2A.pdf',width = 9.5,height = 6)
tmp<-wrap_plots(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10),ncol=5)
print(tmp)
dev.off()


tmp_seurat<-subset(bt_all_noneutrophil,sample=='tumor')
ident1<-'annotation1'
ident2<-'orig.ident'
ident3<-'treat_result_coarse_sample'
Idents(tmp_seurat)<-ident1
tmp_seurat<-subset(tmp_seurat,idents=c('Unkown1','Unkown2','Unkown3'),invert=TRUE)
ratio03 <- data.frame(table(tmp_seurat@meta.data[[ident1]],tmp_seurat@meta.data[[ident2]],tmp_seurat@meta.data[[ident3]]),row.names = NULL)
sum<-data.frame(table(tmp_seurat@meta.data[[ident2]]))
ratio03<-ratio03[ratio03$Freq!=0,]
Freq_nor <- c()
for (i in (1:nrow(ratio03))){
  sumi<-sum[sum$Var1==ratio03[i,'Var2'],'Freq']
  Freq_nor <- c(Freq_nor, ratio03[i, 'Freq']/sumi)

}
ratio03 <- cbind(ratio03, Freq_nor)


draw_list<-c()
count=0
for(i in unique(tmp_seurat$annotation1)){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    tmp<-ggplot(ratio_tmp,aes(x=Var3,y=Freq_nor,color=Var3))+
    geom_boxplot(outlier.shape = NA)+
    # scale_color_manual(values=c("black", "red"))+
    theme_bw()+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "false")+
    ylab("")+xlab("")+ggtitle(i)+
    coord_cartesian(ylim=c(0, 1.0))+
    geom_signif(comparisons = list(c("pCR/MPR_tumor","pCR/MPR-after_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.5),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR_tumor","NMPR_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.6),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("NMPR_tumor","NMPR-after_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.7),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-after_tumor","NMPR-after_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.8),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)
    assign(paste0('p',count),tmp)
}


pdf('tumor_annotation1_with_tumor_without_neutrophil_by_treat_result_coarse_split_formal_figureS2A.pdf',width = 9.50,height = 9.00)
tmp<-wrap_plots(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13),ncol=5)
print(tmp)
dev.off()



# TCR diversity analysis
# conda activate seurat4-r40
# suppressMessages(library(scRepertoire))
library(scRepertoire)
library(ggplot2)
library(ggpubr)
# suppressMessages(library(Seurat))
library(ggsci)
library(dplyr)
library(stringr)
library(ggsignif)
library(Seurat) 
library(patchwork)


seurat_all_tcr<-readRDS('../figure1_output/seurat_all_tcr.rds')
Idents(seurat_all_tcr)<-'treat_result_coarse2_sample'
seurat_all_tcr<-RenameIdents(seurat_all_tcr,'CMPR_blood'='pCR/MPR-native_blood',
'CMPR-after_blood'='pCR/MPR-treated_blood','CMPR_tumor'='pCR/MPR-native_tumor','CMPR-after_tumor'='pCR/MPR-treated_tumor',
"NMPR_blood"='NMPR-native_blood',"NMPR-after_blood"='NMPR-treated_blood',"NMPR_tumor"='NMPR-native_tumor',"NMPR-after_tumor"='NMPR-treated_tumor')
seurat_all_tcr[['treat_result_coarse2_sample']]<-Idents(seurat_all_tcr)

seurat_all_tcr$treat_result_coarse2_sample<-factor(seurat_all_tcr$treat_result_coarse2_sample,
levels = c(paste0(c('pCR/MPR-native','NMPR-native','pCR/MPR-treated','NMPR-treated'),'_','blood'),
paste0(c('pCR/MPR-native','NMPR-native','pCR/MPR-treated','NMPR-treated'),'_','tumor')))

meta<-seurat_all_tcr@meta.data
meta_tmp<-meta
meta_tmp<-meta_tmp[,c('orig.ident','treat_result_coarse2','sample','treat_result_coarse2_sample')]
meta_tmp<-unique(meta_tmp)
contig_list <- list()
sample_names<-c()
for(var_name in meta_tmp$orig.ident){
    tmp<-get(var_name)
    # colnames(tmp)[ncol(tmp)]<-zty_colnames[length(zty_colnames)]
    contig_list[[length(contig_list)+1]]<-tmp
    # sample_names<-c(sample_names,var_name)
}

head(contig_list[[1]])

combined <- combineTCR(contig_list, 
                samples = meta_tmp$orig.ident)
saveRDS(combined,'combined.rds')
clonalDiversity_table<-clonalDiversity(combined, 
                cloneCall = "gene",  
                group.by = "sample",exportTable=TRUE)

colnames(clonalDiversity_table)[7]<-'orig.ident'

clonesize_table<-data.frame(table(seurat_all_tcr$cloneSize,seurat_all_tcr$orig.ident))
sum<-data.frame(table(seurat_all_tcr@meta.data[['orig.ident']]))
Freq_nor <- c()
for (i in (1:nrow(clonesize_table))){
  sumi<-sum[sum$Var1==clonesize_table[i,'Var2'],'Freq']
  Freq_nor <- c(Freq_nor, clonesize_table[i, 'Freq']/sumi)

}
clonesize_table <- cbind(clonesize_table, Freq_nor)
clonesize_table<-clonesize_table[clonesize_table$Var1=='Single (0 < X <= 1)',]
clonesize_table<-clonesize_table[,c('Var2','Freq_nor')]
colnames(clonesize_table)<-c('orig.ident','unique_clone')
clonalDiversity_table<-left_join(clonalDiversity_table,clonesize_table,by='orig.ident')
clonalDiversity_table<-right_join(clonalDiversity_table,meta_tmp,by='orig.ident')



clonalDiversity_table_tmp<-clonalDiversity_table[clonalDiversity_table$sample=='tumor',]
diversity_list<-c('unique_clone','chao1','shannon','inv.simpson','norm.entropy','gini.simpson')
draw_list<-c()
count=0
for(i in unique(diversity_list)){
    count=count+1
    ratio_tmp<-clonalDiversity_table_tmp[,c('orig.ident','treat_result_coarse2_sample',i)]
    colnames(ratio_tmp)[3]<-'Freq_nor'
    ymax<-max(ratio_tmp$Freq_nor)
    tmp<-ggplot(ratio_tmp,aes(x=treat_result_coarse2_sample,y=Freq_nor,color=treat_result_coarse2_sample))+
    geom_boxplot(outlier.shape = NA)+
    # scale_color_manual(values=c("black", "red"))+
    theme_bw()+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "false")+
    ylab("")+xlab("")+ggtitle(i)+
    coord_cartesian(ylim=c(0, 1.5*ymax))+
    geom_signif(comparisons = list(c("pCR/MPR-native_tumor","NMPR-native_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(1.1*ymax),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-treated_tumor","NMPR-treated_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(1.2*ymax),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-native_tumor","pCR/MPR-treated_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(1.3*ymax),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("NMPR-native_tumor","NMPR-treated_tumor")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(1.4*ymax),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)
    assign(paste0('p',count),tmp)
}

pdf('FigureS3C_tumor_clonalDiversity_table_formal.pdf',width = 1.90*2,height = 6)
tmp<-wrap_plots(list(p1,p2,p3,p6),ncol=2)
print(tmp)
dev.off()





clonalDiversity_table_tmp<-clonalDiversity_table[clonalDiversity_table$sample=='blood',]
diversity_list<-c('unique_clone','chao1','shannon','inv.simpson','norm.entropy','gini.simpson')
draw_list<-c()
count=0
for(i in unique(diversity_list)){
    count=count+1
    ratio_tmp<-clonalDiversity_table_tmp[,c('orig.ident','treat_result_coarse2_sample',i)]
    colnames(ratio_tmp)[3]<-'Freq_nor'
    ymax<-max(ratio_tmp$Freq_nor)
    tmp<-ggplot(ratio_tmp,aes(x=treat_result_coarse2_sample,y=Freq_nor,color=treat_result_coarse2_sample))+
    geom_boxplot(outlier.shape = NA)+
    # scale_color_manual(values=c("black", "red"))+
    theme_bw()+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "false")+
    ylab("")+xlab("")+ggtitle(i)+
    coord_cartesian(ylim=c(0, 1.5*ymax))+
    geom_signif(comparisons = list(c("pCR/MPR-native_blood","NMPR-native_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(1.1*ymax),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-treated_blood","NMPR-treated_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(1.2*ymax),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("pCR/MPR-native_blood","pCR/MPR-treated_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(1.3*ymax),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c("NMPR-native_blood","NMPR-treated_blood")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(1.4*ymax),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = TRUE,
                textsize = 3)
    assign(paste0('p',count),tmp)
}


pdf('FigureS3C_blood_clonalDiversity_table_formal.pdf',width = 1.90*2,height = 6)
tmp<-wrap_plots(list(p1,p2,p3,p6),ncol=2)
print(tmp)
dev.off()







# CNV result visualization
library(ggpubr)
setwd('../figure1_output/infercnv')
for(i in c('tumor','blood')){
meta<-readRDS(paste0(i,'_cnv_plotdata.rds'))
meta$annotation1<-as.character(meta$annotation1)
tmp_idents<-as.character(unique(meta$annotation1))
tmp_levels<-c(
  'Tumor','Ciliated','T','NK','Endothelial',tmp_idents[!tmp_idents%in%c('Tumor','Ciliated','T','Endothelial','NK')]
)
if(i=='tumor'){
  tmp_levels<-tmp_levels[!tmp_levels%in%c('Unkown1','Unkown2','Unkown3')]
  meta<-meta[meta$annotation1%in%tmp_levels,]
  # meta$annotation1<-factor(meta$annotation1,levels=tmp_levels)
}
if(i=='blood'){
  meta[meta$annotation1=='Tumor','annotation1']<-'Tumor?'
  tmp_idents<-as.character(unique(meta$annotation1))
  tmp_levels<-c(
  'Tumor?','Ciliated','T','NK','Endothelial',tmp_idents[!tmp_idents%in%c('Tumor?','Ciliated','T','NK','Endothelial')]
  )
  tmp_levels<-tmp_levels[!tmp_levels%in%c('Unkown1','Unkown2','Unkown3','Ciliated','Fibroblast','Endothelial','Mast')]
  meta<-meta[meta$annotation1%in%tmp_levels,]
  # meta$annotation1<-factor(meta$annotation1,levels=tmp_levels)
}
if(i=='tumor'){
p <- ggboxplot(meta, "annotation1", "cnv_score",
 fill = "annotation1",order = tmp_levels,outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, 2000)) +  
  xlab("Cluster") +  
  ylab("CNV Score") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
  axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")+
    geom_signif(comparisons = list(c('Tumor','Ciliated')),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = 1300,#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.5,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c('Tumor','T')),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = 1400,#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.5,color="black",map_signif_level = TRUE,
                textsize = 3)+    
    geom_signif(comparisons = list(c('Tumor','NK')),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = 1500,#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.5,color="black",map_signif_level = TRUE,
                textsize = 3)+
    geom_signif(comparisons = list(c('Tumor?','T')),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = 1600,#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.5,color="black",map_signif_level = TRUE,
                textsize = 3)
}

if(i=='blood'){
p <- ggboxplot(meta, "annotation1", "cnv_score",
 fill = "annotation1",order = tmp_levels,outlier.shape = NA) + 
  scale_y_continuous(limits = c(0, 2000)) +  
  xlab("Cluster") +  
  ylab("CNV Score") +
  theme(panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
  axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")+
    geom_signif(comparisons = list(c('Tumor?','Plasma')),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = 1300,#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.5,color="black",map_signif_level = TRUE,
                textsize = 3)
}


  pdf(paste0(i,'_cnv_score_figureS1C.pdf'), width = 7, height = 4)
  print(p)
  dev.off()
}







published_lusc<-readRDS('../figure1_output/lusc.rds')
head(published_lusc)
unique(published_lusc$origin)
unique(published_lusc$cell_type_major)
Idents(published_lusc)<-'cell_type_major'
published_lusc<-RenameIdents(published_lusc,'Macrophage alveolar'='Macrophage','cDC2'='DC',
'T cell regulatory'='T cell','Alveolar cell type 1'='Normal epithelial','T cell CD4'='T cell',
'cDC2'='DC','T cell CD8'='T cell','Ciliated'='Normal epithelial','transitional club/AT2'='Normal epithelial',
'Club'='Normal epithelial',"Alveolar cell type 2"='Normal epithelial','cDC1'='DC','pDC'='DC','DC mature'='DC',
'Tumor cells'='Tumor cell')
published_lusc[['cell_type_major_compare']]<-Idents(published_lusc)
unique(published_lusc$cell_type_major_compare)

Idents(published_lusc)<-'published'
published_lusc[['data_source']]<-Idents(published_lusc)
tmp<-subset(published_lusc,origin=='normal',invert=TRUE)
Idents(tmp)<-'donor_id'
tmp[['orig.ident']]<-Idents(tmp)
t_seurat<-merge(subset(annotation05,condition=='native'),tmp)

t_seurat_notumor<-subset(t_seurat,cell_type_major_compare=='Tumor cell',invert=TRUE)

t_seurat_notumor_noneutro<-subset(t_seurat_notumor,cell_type_major_compare=='Neutrophils',invert=TRUE)
t_seurat_notumor_noneutro$cell_type_major_compare<-droplevels(t_seurat_notumor_noneutro$cell_type_major_compare)


draw_list<-c()
count=0
for(i in unique(t_seurat_notumor_noneutro$cell_type_major_compare)){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    png(paste0('cell_type_major_',i,'_compare_without_tumor_without_neutro.png'),width = 300,height=300,res=100)
    tmp<-ggplot(ratio_tmp,aes(x=Var3,y=Freq_nor,color=Var3))+
    geom_boxplot(outlier.shape = NA)+
    # scale_color_manual(values=c("black", "red"))+
    theme_bw()+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "right")+
    ylab("")+xlab("")+
    coord_cartesian(ylim=c(0, 1))+
    geom_signif(comparisons = list(c("published","this_paper")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.8),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.5,color="black",
                textsize = 3)
    print(tmp)
    dev.off()
    tmp<-ggplot(ratio_tmp,aes(x=Var3,y=Freq_nor,color=Var3))+
    geom_boxplot(outlier.shape = NA)+
    # scale_color_manual(values=c("black", "red"))+
    theme_bw()+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")+
    ylab("")+xlab("")+ggtitle(i)+
    coord_cartesian(ylim=c(0, 1))+
    geom_signif(comparisons = list(c("published","this_paper")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.8),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.5,color="black",
                textsize = 3)
    assign(paste0('p',count),tmp)
}




pdf('cell_type_major_compare_without_tumor_without_neutrophil_figureS1F.pdf',width = 10,height = 5)
tmp<-wrap_plots(list(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12),ncol=6)
print(tmp)
dev.off()



# peripheral blood boxplot
pbmc_data_bybatch_present_merge<-readRDS('../figure1_output/pbmc_data_bybatch_present2.rds')

pbmc_data_bybatch_present_merge<-pbmc_data_bybatch_present
pbmc_data_bybatch_present_merge$response<-as.character(pbmc_data_bybatch_present_merge$response)
pbmc_data_bybatch_present_merge[pbmc_data_bybatch_present_merge$response%in%c('pCR','MPR'),'response']<-'pCR/MPR'
pbmc_data_bybatch_present_merge$response<-
factor(pbmc_data_bybatch_present_merge$response,levels=c('pCR/MPR','NMPR'))

pdf('cbc_boxplot_figure1G.pdf',width = 4,height=4)
tmp_plot<-ggplot(pbmc_data_bybatch_present_merge,aes(x=response,y=monocyte_ratio_without_granulocyte,color=response))+
geom_boxplot(outlier.shape = NA)+
geom_jitter(width = 0.2, alpha = 0.6)+
scale_color_manual(values=c('#00db0b','#d80404'))+
theme_bw()+
RotatedAxis()+ylim(0,0.6)+
theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")+
    geom_signif(comparisons = list(c("pCR/MPR","NMPR")),#Set groups to compare
                test = wilcox.test, ##Calculation method
                y_position = c(0.5),#Horizontal line position setting in the graph
                tip_length = c(c(0.03,0.03)),#Vertical line setting below the horizontal line
                size=0.15,color="black",map_signif_level = FALSE,
                textsize = 3)

print(tmp_plot)
dev.off()


# conda activate seurat4-r42
# CellChat heatmap for each group
library(CellChat)
# annotation1
setwd('../figure1_output/cellchat')
count<-0
for(i in c('_blood','-after_blood','_tumor','-after_tumor')){
    count=count+1
    blood_native_cellchat_cmpr<-readRDS(paste0('CMPR',i,'_bt_all_annotation1_cellchat.rds'))
    blood_native_cellchat_nmpr<-readRDS(paste0('NMPR',i,'_bt_all_annotation1_cellchat.rds'))
    object.list <- list(CMPR=blood_native_cellchat_cmpr, NMPR=blood_native_cellchat_nmpr)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    assign(paste0('p',as.character(count)),netVisual_heatmap(cellchat,comparison = c(2, 1), measure = "weight"))
}

pdf('netVisual_heatmap_cmpr_nmpr_blood_figureS3B.pdf',7,4)
p1+p2
dev.off()

pdf('netVisual_heatmap_cmpr_nmpr_tumor_figureS3B.pdf',7,4)
p3+p4
dev.off()
