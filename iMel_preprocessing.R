library(Seurat)
library(Signac)
library(cellranger)
library(magrittr)
library(data.table)
library(patchwork)
library(dplyr)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(GGally)
library(gplots)
library(plotly)
library(cowplot)
library(scales)

# Read and merge timecourse data 
#iMel_iPS quality control

#iMel_iPS_D23
D23.data <- Read10X(data.dir = '')
colnames(D23.data) <- paste(colnames(D23.data),"D23",sep = "_")
D23 <- CreateSeuratObject(counts = D23.data, project = "iMel_D23", assay = 'RNA',min.cells = 3)
D23[["percent.mt"]] <- PercentageFeatureSet(D23, pattern = "^MT-")
rm(D23.data)

D23 <- subset(
  x = D23, 
  subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & 
	   nCount_RNA >500 & nCount_RNA <5e4 & 
           percent.mt < 10)
 
VlnPlot(D23, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0)


#iMel_iPS_D32
D32.data <- Read10X(data.dir = '')
colnames(D32.data) <- paste0("D32_",colnames(D32.data))
D32 <- CreateSeuratObject(counts = D32.data, project = "iMel_D32", assay = 'RNA',min.cells = 3)
D32[["percent.mt"]] <- PercentageFeatureSet(D32, pattern = "^MT-")
rm(D32.data)

D32 <- subset(
  x = D32, 
  subset = nFeature_RNA > 3000 &nFeature_RNA < 9000 & 
           nCount_RNA >500 & nCount_RNA <1e5 & 
           percent.mt < 10)

VlnPlot(D32, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0)


##iMel_iPS_D32
D44.data <- Read10X(data.dir = '')
colnames(D44.data) <- paste0("D44_",colnames(D44.data))
D44 <- CreateSeuratObject(counts = D44.data, project = "iMel_D44", assay = 'RNA',min.cells = 3)
D44[["percent.mt"]] <- PercentageFeatureSet(D44, pattern = "^MT-")
rm(D44.data)

D44 <- subset(
  x = D44, 
  subset = nFeature_RNA > 3000 & nFeature_RNA < 9000 &
           nCount_RNA > 500 & nCount_RNA <1e5 & 
           percent.mt < 10)

VlnPlot(D44, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), ncol = 3, pt.size = 0)


# merge timecourse data
iMel.merged <- merge(D23, y = c(D32, D44), project = "iMel_",merge.data = T)
common_gene = rownames(D23)[rownames(D23) %in% rownames(D32) ==T];
common_gene = common_gene[common_gene %in% rownames(D44) == T ]
iMel.merged.comgene <- subset(iMel.merged,feature = common_gene)


# get gene name lists for cell cycle scoring and regression (s and g2m phases)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
iMel.merged.comgene <- CellCycleScoring(iMel.merged.comgene, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

iMel.merged.comgene <- iMel.merged.comgene %>%
  NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData() %>%
  SCTransform(vars.to.regress = c("S.Score","G2M.Score"))

iMel.merged.comgene <- RunPCA(iMel.merged.comgene, assay = "SCT", npcs = 50)


#harmony
library(harmony)
iMel_seu <- RunHarmony(iMel.merged.comgene, 
                       group.by.vars = c("orig.ident", "batch"), 
                       reduction = "pca", assay.use = "SCT", reduction.save = "harmony")
iMel_seu <- RunUMAP(iMel_seu, reduction = "harmony", assay = "SCT", dims = 1:50)
iMel_seu <- FindNeighbors(iMel_seu, reduction = "harmony")
iMel_seu <- FindClusters(iMel_seu, resolution = 0.4)


#UMAP plot
#color setting
color_iMel_time <- c("#B0CB9A","#D0A364","#A90000")
color_iMel_group_1 <- c("#516B97","#9DBD67","#DF99B8")
color_iMel_res0.4 <- c("#132D94","#50619F","#54B9F5","#81C1C4","#583BF1","#4B9BC9","#68E2DA","#DF99B8","#9DBD67","#E5828E")

FigS1C <- DimPlot(iMel_seu, group.by = "orig.ident",reduction = "umap",label.size = 5, cols = color_iMel_time)

DimPlot(iMel_seu, group_by = "SCT_snn_res.0.4", reduction = "umap",label = TRUE,label.size = 5, cols = cluster.col1)

#add metadata
iMel_seu $group_1 = NA
iMel_seu $group_1[iMel_seu $SCT_snn_res.0.4 %in% c(0,1,2,3,4,5) ==T]  = "C1"
iMel_seu $group_1[iMel_seu $SCT_snn_res.0.4 %in% c(7) ==T]  = "C2"
iMel_seu $group_1[iMel_seu $SCT_snn_res.0.4 %in% c(6,8) ==T]  = "C3"

FigS1D <-  DimPlot(iMel_seu, group.by = "group_1",reduction = "umap",label.size = 5,cols = color_iMel_group_1)

Fig2A <- DimPlot(iMel_seu, group.by = "group_1",split.by = "orig.ident", reduction = "umap",label.size = 5,cols = color_iMel_group_1)

##cell cluster percentage plot
data <- as.data.frame(table(iMel_seu $orig.ident, iMel_seu $group_1))

Fig2B <- ggplot(data,aes(x =Var1, y= Freq, fill = Var2)) + 
  geom_bar(stat = "identity",position = "fill",width = 0.7,size = 0.5,colour = '#222222')+ 
  scale_y_continuous(labels = scales::percent) + 
  theme_bw(base_line_size = 0) +
  labs(x='',y = 'Percentage')+
  scale_fill_manual(name='clusters', 
                    labels=c("C1","C2","C3"),
                    values= color_iMel_group_1)+
  theme(legend.title = element_blank(),
        axis.text = element_text(hjust = 0.5,vjust = 1,angle = 0,size = 12),
        axis.title = element_text(size = 12)) + 
  NoLegend()

##qc violin plot
FigS1B <- VlnPlot(iMel_seu,group.by = "orig.ident",features = 'nCount_RNA', #'percent.mt', 'nFeature_RNA'
                    pt.size = 0, cols=sample_id_cols) +
  geom_boxplot(width=0.1, fill="gray", color='black') + 
  NoLegend() +
  ggtitle("nCount_RNA") +
  theme(axis.title = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5,size = 12),
        title = element_text(size = 12)) 

#C1 sub-cluster UMAP plot
#color setting
color_iMel_group_2 = c("#132D94","#50619F","#54B9F5","#81C1C4","#583BF1","#4B9BC9","#9DBD67","#DF99B8")

#add metadata
iMel_seu $group_2 =NA
iMel_seu $group_2[iMel_seu $SCT_snn_res.0.4=="0"] = "C1a"
iMel_seu $group_2[iMel_seu $SCT_snn_res.0.4=="1"] = "C1b"
iMel_seu $group_2[iMel_seu $SCT_snn_res.0.4=="2"] = "C1c"
iMel_seu $group_2[iMel_seu $SCT_snn_res.0.4=="3"] = "C1d"
iMel_seu $group_2[iMel_seu $SCT_snn_res.0.4=="4"] = "C1e"
iMel_seu $group_2[iMel_seu $SCT_snn_res.0.4=="5"] = "C1f"
iMel_seu $group_2[iMel_seu $SCT_snn_res.0.4 %in% c(7) ==T]  = "C2"
iMel_seu $group_2[iMel_seu $SCT_snn_res.0.4 %in% c(6,8) ==T]  = "C3"

iMel_seu $group_2 <- factor(iMel_seu $group_2, levels = c("C1a","C1b","C1c","C1d","C1e","C1f","C2","C3")) 

C1_sub <- subset(iMel_seu,idents = "C1")

Fig3C <- DimPlot(C1_sub, group.by = "group_2", reduction = "umap", label.size = 5, cols = color_iMel_group_2, label = T, label.box = T, repel = T, label.color = "white")
