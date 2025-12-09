# conda activate seurat5-r43
# rm(list=ls())
gc()
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
# library(clusterProfiler)
# library("org.Mm.eg.db")
# library(org.Hs.eg.db)
# library(data.table)
# library(biomaRt)
library(rhdf5)
# library(Matrix)
library(scCustomize)
# devtools::install_version("scCustomize", version = "2.0.0",
#  repos = "http://cran.us.r-project.org")
library(viridis)
library(RColorBrewer)
# library(dittoSeq)
# library(plot1cell)
# library(CellChat)
library(SingleR)
library(ggpubr)
library(SCP)
# hpca.se <- celldex::HumanPrimaryCellAtlasData()
options(future.globals.maxSize = 3e+09)
set.seed(1)# Set seed for reproducibility
setwd('../figure3_output')

monocyte_all_scanvi_formal<-readRDS('../figure2_output/monocyte_all_scanvi_formal.rds')
unique(monocyte_all_scanvi_formal$treat_result_coarse_sample)



# Because the y-axis limits for blood-after are too narrow, causing missing points, plot the y-axis separately.
# Some data have abnormal p-values, suspected to be due to incomplete red blood cell lysis.
tmp<-subset(monocyte_all_scanvi_formal,idents=c('pCR/MPR-after_blood','NMPR-after_blood'))
# # Rename NMPR so that pCR/MPR appears in the second column, adjusting y-axis range.
# Idents(tmp)<-'treat_result_coarse_sample'
# tmp<-RenameIdents(tmp,'NMPR-after_blood'='aNMPR-after_blood')
# tmp[['treat_result_coarse_sample']]<-Idents(tmp)
tmp_function2<-function(t_seurat_blood,file_name){
t_seurat_blood<-subset(t_seurat_blood,orig.ident=='B_15949890_WNG',invert=TRUE)
t_seurat_blood <- RunDEtest(srt = t_seurat_blood, group_by = "treat_result_coarse_sample", fc.threshold = 1, only.pos = FALSE,
BPPARAM = BiocParallel::MulticoreParam(workers = 4))
DEGs_b <- t_seurat_blood@tools$DEtest_treat_result_coarse$AllMarkers_wilcox
DEGs_b <- DEGs_b[with(DEGs_b, avg_log2FC > 0 & p_val_adj < 0.05), ]
tmp<-t_seurat_blood@tools[[1]][[1]]
t_seurat_blood@tools[[1]][[1]]<-tmp[tmp$p_val_adj>10e-500,]# Remove genes with abnormal p-values in tools
png(paste0(file_name,'.png'),width=1200,height=700)
tmp_plot<-VolcanoPlot(srt = t_seurat_blood, 
#  DE_threshold = "avg_log2FC > 0 & p_val_adj < 0.05 & p_val_adj > 10e-300",
group_by = "treat_result_coarse_sample",nlabel = 5)# To modify y-axis display range, add +ylim()
print(tmp_plot)
dev.off()
pdf(paste0(file_name,'.pdf'),width=12,height=7)
tmp_plot<-VolcanoPlot(srt = t_seurat_blood, group_by = "treat_result_coarse_sample",nlabel = 5)
print(tmp_plot)
dev.off()
}
tmp_function2(subset(monocyte_all_scanvi_formal,idents=c('pCR/MPR_blood','NMPR_blood')),'blood_native_monocyte_volcanplot_figure3A')
tmp_function2(subset(monocyte_all_scanvi_formal,idents=c('pCR/MPR-after_blood','NMPR-after_blood')),'blood_after_monocyte_volcanplot_figure3A')
tmp_function2(subset(monocyte_all_scanvi_formal,idents=c('pCR/MPR_tumor','NMPR_tumor')),'tumor_native_monocyte_volcanplot_figure3A')
tmp_function2(subset(monocyte_all_scanvi_formal,idents=c('pCR/MPR-after_tumor','NMPR-after_tumor')),'tumor_after_monocyte_volcanplot_figure3A')


setwd('../figure3_output/cellchat')
count<-0
for(i in c('_blood','-after_blood','_tumor','-after_tumor')){
    count=count+1
    blood_native_cellchat_cmpr<-readRDS(paste0('CMPR',i,'_bt_all_annotation1_apobec3a_cellchat.rds'))
    blood_native_cellchat_nmpr<-readRDS(paste0('NMPR',i,'_bt_all_annotation1_apobec3a_cellchat.rds'))
    object.list <- list(CMPR=blood_native_cellchat_cmpr, NMPR=blood_native_cellchat_nmpr)
    cellchat <- mergeCellChat(object.list, add.names = names(object.list))
    assign(paste0('p',as.character(count)),netVisual_heatmap(cellchat,comparison = c(2, 1), measure = "weight"))
}

pdf('netVisual_heatmap_cmpr_nmpr_a3a_blood_figureS8B.pdf',7,4)
p1+p2
dev.off()

pdf('netVisual_heatmap_cmpr_nmpr_a3a_tumor_figureS8B.pdf',7,4)
p3+p4
dev.off()




# Explanation for why WNB-B-n is removed
setwd('../figure3_output')
monocyte_all_scanvi_formal<-readRDS('../figure2_output/monocyte_all_scanvi_formal.rds')
orig_sample_id_df<-read.csv('../figure1_output/filename_samplename_df.csv',sep=',',header = T)

Idents(monocyte_all_scanvi_formal)<-'orig.ident'
for(i in orig_sample_id_df$file_name){
    cell_id<-colnames(subset(monocyte_all_scanvi_formal,orig.ident==i))
    monocyte_all_scanvi_formal<-SetIdent(monocyte_all_scanvi_formal,cells=cell_id,value=orig_sample_id_df[orig_sample_id_df$file_name==i,'sample_id'])
}
monocyte_all_scanvi_formal[['formal_sample_id']]<-Idents(monocyte_all_scanvi_formal)

monocyte_all_scanvi_formal_n<-subset(monocyte_all_scanvi_formal,condition=='native')
all_genes<-rownames(monocyte_all_scanvi_formal_n)

for(s in c('tumor','blood')){
out_deg_df<-data.frame(gene_sym=all_genes)
monocyte_all_scanvi_formal_b_n<-subset(monocyte_all_scanvi_formal_n,sample==s)
for(i in unique(monocyte_all_scanvi_formal_b_n$formal_sample_id)){
tmp<-subset(monocyte_all_scanvi_formal_b_n,formal_sample_id==i,invert=TRUE)
tmp_deg<-deg_sort2(tmp,ident.i='pCR/MPR',ident.j='NMPR',idents='treat_result_coarse')
tmp_deg<-tmp_deg[,c(1,2)]
tmp_deg<-tmp_deg[abs(tmp_deg$avg_log2FC)>0.5,]
colnames(tmp_deg)[2]<-paste0(i,'_FC')
out_deg_df<-left_join(out_deg_df,tmp_deg,by='gene_sym')
}

out_deg_df[is.na(out_deg_df)]<-0
out_deg_df_filter<-filter(out_deg_df,!rowSums(out_deg_df[,-1])==0)
combinded_df.t <- t(out_deg_df_filter[,2:ncol(out_deg_df_filter)])


#PCA calculation
comdf.prc <- prcomp(combinded_df.t, center = FALSE,scale. = FALSE)# scale. = TRUE normalizes data before analysis
summary(comdf.prc)
#Extract PCA scores
prc_score <- comdf.prc$x
#Merge sample classification into data frame
prc_score <- data.frame(prc_score, sample_id = rownames(combinded_df.t))
# annotation_col_df <- readRDS('/home/zijin/ouyang_OA/code/output/annotation_col_df.rds')# Load sample_id classification info
# annotation_col_df['sample_id'] <- row.names(annotation_col_df)
# prc_df <- left_join(prc_score,annotation_col_df)
# prc_df <- prc_df[,!names(prc_df) %in% c('sample_id','OA_sample')]

#Extract variance contribution rate of principal components to generate axis titles
summ<-summary(comdf.prc)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
library(ggrepel)
p2<-ggplot(data = prc_score,aes(x=PC1,y=PC2,color=sample_id))+
stat_ellipse(aes(fill=sample_id),
type = "norm", geom ="polygon",alpha=0.2,color=NA)+
geom_point(label=TRUE)+labs(x=xlab,y=ylab,color="")+
guides(fill="none")+ggrepel::geom_text_repel(
  aes(label=sample_id)
  )+theme_bw()

setwd('/public/home/JianLiu/zijin/lusc_patients/code2/formal_code/output/figure3_new')

pdf(paste0(s,'_CDD_PCA_figureS7A.pdf'))
print(p2)
# +scale_fill_manual(values = c("purple","orange","blue","red","#00ff15","#03c7f8","brown"))+
# scale_colour_manual(values = c("purple","orange","blue","red","#00ff15","#03c7f8","brown"))
dev.off()

}

head(out_deg_df)
out_deg_df[out_deg_df$gene_sym=='APOBEC3A',]


# The above results indicate potential issues with WNG-B-native and DHL-T-native, requiring further validation.

tmp_function<-function(t_seurat_blood,invert_ident=NULL,file_name){
    if(!is.null(invert_ident)){
    t_seurat_blood<-subset(t_seurat_blood,orig.ident==invert_ident,invert=TRUE)
    }
t_seurat_blood <- RunDEtest(srt = t_seurat_blood, group_by = "treat_result_coarse_sample", fc.threshold = 1, only.pos = FALSE,
BPPARAM = BiocParallel::MulticoreParam(workers = 1))
DEGs_b <- t_seurat_blood@tools$DEtest_treat_result_coarse$AllMarkers_wilcox
DEGs_b <- DEGs_b[with(DEGs_b, p_val_adj < 0.05), ]

pdf(paste0(file_name,'.pdf'),width=12,height=7)
tmp_plot<-VolcanoPlot(srt = t_seurat_blood, group_by = "treat_result_coarse_sample",nlabel = 5)
print(tmp_plot)
dev.off()
    changes<-c()
    for(i in c(1:nrow(DEGs_b))){
      if(DEGs_b[i,'avg_log2FC']>0){
        changes<-c(changes,paste0(DEGs_b[i,'gene'],'+'))
      }else{
        changes<-c(changes,paste0(DEGs_b[i,'gene'],'-'))
      }
    }
    DEGs_b['changes']<-changes
return(DEGs_b[DEGs_b$group1==unique(DEGs_b$group1)[1],])
}

monocyte_all_scanvi_formal_b_n<-subset(monocyte_all_scanvi_formal_n,sample=='blood')
control_b_n_deg<-tmp_function(monocyte_all_scanvi_formal_b_n,file_name='blood_invert_deg_control')
test_b_n_deg<-tmp_function(monocyte_all_scanvi_formal_b_n,'B_15949890_WNG',file_name='blood_invert_deg_test')
intersect_top<-function(control_b_n_deg,test_b_n_deg,n_top){
    control_b_n_deg<-control_b_n_deg[order(control_b_n_deg$avg_log2FC,decreasing=TRUE),]
    test_b_n_deg<-test_b_n_deg[order(test_b_n_deg$avg_log2FC,decreasing=TRUE),]
    top_deg1<-control_b_n_deg[c(1:n_top),'changes']
    bot_deg1<-control_b_n_deg[c((nrow(control_b_n_deg)+1-n_top):nrow(control_b_n_deg)),'changes']
    top_deg2<-test_b_n_deg[c(1:n_top),'changes']
    bot_deg2<-test_b_n_deg[c((nrow(test_b_n_deg)+1-n_top):nrow(test_b_n_deg)),'changes']
    return(c(length(intersect(top_deg1,top_deg2)),length(intersect(bot_deg1,bot_deg2))))
}

intersect_top(control_b_n_deg,test_b_n_deg,100)


monocyte_all_scanvi_formal_b_n<-subset(monocyte_all_scanvi_formal_n,sample=='tumor')
control_t_n_deg<-tmp_function(monocyte_all_scanvi_formal_b_n,file_name='tumor_invert_deg_control')
test_t_n_deg<-tmp_function(monocyte_all_scanvi_formal_b_n,file_name='T_15464852_DHL','tumor_invert_deg_test')

intersect_top(control_t_n_deg,test_t_n_deg,100)


# The above results indicate that WNG-B-native is problematic, while DHL-T-native is not. Therefore, peripheral blood needs to exclude WNG-B-native and run CDD again to confirm no issues.
monocyte_all_scanvi_formal_n_filter<-subset(monocyte_all_scanvi_formal_n,orig.ident=='B_15949890_WNG',invert=TRUE)
for(s in c('blood')){
out_deg_df<-data.frame(gene_sym=all_genes)
monocyte_all_scanvi_formal_b_n<-subset(monocyte_all_scanvi_formal_n_filter,sample==s)
for(i in unique(monocyte_all_scanvi_formal_b_n$formal_sample_id)){
tmp<-subset(monocyte_all_scanvi_formal_b_n,formal_sample_id==i,invert=TRUE)
tmp_deg<-deg_sort2(tmp,ident.i='pCR/MPR',ident.j='NMPR',idents='treat_result_coarse')
tmp_deg<-tmp_deg[,c(1,2)]
tmp_deg<-tmp_deg[abs(tmp_deg$avg_log2FC)>0.5,]
colnames(tmp_deg)[2]<-paste0(i,'_FC')
out_deg_df<-left_join(out_deg_df,tmp_deg,by='gene_sym')
}

out_deg_df[is.na(out_deg_df)]<-0
out_deg_df_filter<-filter(out_deg_df,!rowSums(out_deg_df[,-1])==0)
combinded_df.t <- t(out_deg_df_filter[,2:ncol(out_deg_df_filter)])


#PCA calculation
comdf.prc <- prcomp(combinded_df.t, center = FALSE,scale. = FALSE)# scale. = TRUE normalizes data before analysis
summary(comdf.prc)
#Extract PCA scores
prc_score <- comdf.prc$x
#Merge sample classification into data frame
prc_score <- data.frame(prc_score, sample_id = rownames(combinded_df.t))
# annotation_col_df <- readRDS('/home/zijin/ouyang_OA/code/output/annotation_col_df.rds')# Load sample_id classification info
# annotation_col_df['sample_id'] <- row.names(annotation_col_df)
# prc_df <- left_join(prc_score,annotation_col_df)
# prc_df <- prc_df[,!names(prc_df) %in% c('sample_id','OA_sample')]

#Extract variance contribution rate of principal components to generate axis titles
summ<-summary(comdf.prc)
xlab<-paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab<-paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
library(ggrepel)
p2<-ggplot(data = prc_score,aes(x=PC1,y=PC2,color=sample_id))+
stat_ellipse(aes(fill=sample_id),
type = "norm", geom ="polygon",alpha=0.2,color=NA)+
geom_point(label=TRUE)+labs(x=xlab,y=ylab,color="")+
guides(fill="none")+ggrepel::geom_text_repel(
  aes(label=sample_id)
  )+theme_bw()


pdf(paste0(s,'_CDD_PCA_rmP12_figureS7A.pdf'))
print(p2)
# +scale_fill_manual(values = c("purple","orange","blue","red","#00ff15","#03c7f8","brown"))+
# scale_colour_manual(values = c("purple","orange","blue","red","#00ff15","#03c7f8","brown"))
dev.off()

}





tmp_function<-function(t_seurat,file_name,comparison,colors){
ident1<-'annotation1_apobec3a'
ident2<-'orig.ident'
ident3<-'treat_result_coarse_sample'
ratio03 <- data.frame(table(t_seurat@meta.data[[ident1]],t_seurat@meta.data[[ident2]],t_seurat@meta.data[[ident3]]),row.names = NULL)
sum<-data.frame(table(t_seurat@meta.data[[ident2]]))
ratio03<-ratio03[ratio03$Freq!=0,]
Freq_nor <- c()
for (i in (1:nrow(ratio03))){
  sumi<-sum[sum$Var1==ratio03[i,'Var2'],'Freq']
  Freq_nor <- c(Freq_nor, ratio03[i, 'Freq']/sumi)

}
ratio03 <- cbind(ratio03, Freq_nor)


draw_list<-c()
count=0
for(i in c('APOBEC3A+','APOBEC3A-')){
    count=count+1
    ratio_tmp<-ratio03[ratio03$Var1==i,]
    pi<-ggplot(ratio_tmp,aes(x=Var3,y=Freq_nor,color=Var3))+
    geom_boxplot(outlier.shape = NA)+
    scale_color_manual(values=colors)+
    theme_bw()+
    RotatedAxis()+
    theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.position = "none")+
    ylab("")+xlab("")+ggtitle(i)+
    coord_cartesian(ylim=c(0, 1))+
    geom_signif(comparisons = list(comparison),#Set groups for comparison
                test = wilcox.test, ##Calculation method
                y_position = c(0.93),#Set horizontal line position in plot
                tip_length = c(c(0.03,0.03)),#Set vertical lines below horizontal line
                size=0.5,color="black",
                textsize = 3)
    assign(paste0('p',count),pi)
}



png(paste0(file_name,'.png'),width = 400,height = 300,res=100)
tmp<-wrap_plots(list(p1,p2),ncol=2)
print(tmp)
dev.off()

pdf(paste0(file_name,'.pdf'),width = 4,height = 3)
tmp<-wrap_plots(list(p1,p2),ncol=2)
print(tmp)
dev.off()

}

tmp<-subset(bt_all,orig.ident=='B_15949890_WNG',invert=TRUE)

colors1=c('blue','red')
Idents(bt_all)<-'treat_result_coarse_sample'
tmp_function(subset(tmp,idents=c('pCR/MPR-native_blood','NMPR-native_blood')),'blood_native_a3amonocyte_boxplot_figure3D',c('pCR/MPR-native_blood','NMPR-native_blood'),colors1)
tmp_function(subset(tmp,idents=c('pCR/MPR-treated_blood','NMPR-treated_blood')),'blood_after_a3amonocyte_boxplot_figure3D',c('pCR/MPR-treated_blood','NMPR-treated_blood'),colors1)
tmp_function(subset(tmp,idents=c('pCR/MPR-native_tumor','NMPR-native_tumor')),'tumor_native_a3amonocyte_boxplot_figure3D',c('pCR/MPR-native_tumor','NMPR-native_tumor'),colors1)
tmp_function(subset(tmp,idents=c('pCR/MPR-treated_tumor','NMPR-treated_tumor')),'tumor_after_a3amonocyte_boxplot_figure3D',c('pCR/MPR-treated_tumor','NMPR-treated_tumor'),colors1)