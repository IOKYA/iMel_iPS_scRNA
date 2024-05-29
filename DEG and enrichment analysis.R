library(clusterProfiler)
library(openxlsx)
library(ggplot2)

# find DEGs in iMel 

allmarkers_group_1 <- FindAllMarkers(iMel_seu, logfc.threshold = 0.25, min.pct = 0.25, only.pos = T, test.use = "wilcox")

#Dot plot of top20 DEGs 

top20_by_group_1 <- allmarkers_group_1 %>% group_by(cluster) %>% top_n(n = 20,  wt = avg_log2FC)
gene = rev(unique(top20_by_group_1$gene))
DoHeatmap(subset(iMel_seu, downsample = 100), features = unique(top20_by_group_1$gene), size = 3,label = F,group.colors = color_group_1) +
    theme(axis.text.y = element_text(face = "italic")) +
  scale_fill_gradient2( low = rev(c('#d1e5f0','#67a9cf','#2166ac')), mid = "white", high = rev(c('#b2182b','#ef8a62','#fddbc7')), midpoint = 0, guide = "colourbar", aesthetics = "fill")


#Feature plot of melanocyte genes

custom_red <- scale_color_gradient(low = alpha('grey', 0.1), high = alpha('red', 1), limits = c(0, 4))  
mel.gene <- c("MITF", "DCT", "PMEL", "TYR", "TYRP1", "PAX3", "SOX10", "KIT", "MLANA")
for (i in 1:length(mel.gene)) {
  p <- FeaturePlot(iMel_seu, mel.gene[i], pt.size = 1) + 
    custom_red + 
    theme_bw(base_line_size = 0) +
    theme(title = element_text(face = 'italic'),
          axis.text = element_blank(),  
          axis.title = element_blank()) +
    NoLegend() 
  ggsave(p, filename = paste0(mel.gene[i], "_featureplot.jpg"), width = 2.8, height = 3)
}


# GO-bp of DEGs in group_1

df_sig <- allmarkers_group_1[allmarkers_group_1$p_val_adj <0.05,]

group <- data.frame(gene = df_sig$gene, group = df_sig$cluster)
Gene_ID <- bitr(df_sig$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
data <- merge(Gene_ID,group,by.x = "SYMBOL",by.y = "gene")

go_analysis <- compareCluster(ENTREZID~group,
                          data = data,
                          fun = "enrichGO",
                          OrgDb =  "org.Hs.eg.db",
                          ont = "BP",
                          pAdjustMethod = "bonferroni", 
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = T)

go_analysis_sim <- simplify(go_analysis,
                        cutoff =0.7,
                        by = "p.adjust",
                        select_fun = min)

dotplot(go_analysis_sim,showCategory = 10,font.size =10) +  
  theme(panel.grid = element_blank(), 
    axis.title.x = element_text(size = 12))


# GSEA plot (C1 vs rest)

geneset <- read.gmt("enrichment_gmt/msigdb_v2022.1.Hs_GMTs/c2.cp.kegg.v2022.1.Hs.entrez.gmt")  

markers <- FindMarkers(iMel_seu, ident.1 = "C1", min.pct = 0.25, logfc.threshold = 0.25)
gs <-bitr(rownames(markers), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
markers1<-cbind(markers[gs[,1],],gs)
geneList = markers1$avg_log2FC
names(geneList) = markers1$ENTREZID
geneList = sort(geneList,decreasing = T)

library(enrichplot)

gseaplot2(egmt, geneSetID = i, title = egmt@result$ID[4]) 


#GO-bp of DEGs in C1 sub clusters 

C1_sub <- subset(iMel_seu,idents = "C1")
Idents(C1_sub) <- "group_2”;levels(C1_sub)

allmarkers_C1_sub <- FindAllMarkers(C1_sub,l ogfc.threshold = 0.25, min.pct = 0.25, only.pos = T)

group <- data.frame(gene = allmarkers_C1_sub$gene, group = allmarkers_C1_sub$cluster)
Gene_ID <- bitr(group$gene ,fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = 'org.Hs.eg.db')
data <- merge(Gene_ID,group, by.x = "SYMBOL", by.y = "gene")

data_GO <- compareCluster(ENTREZID~group,
                          data = data,
                          fun = "enrichGO",
                          OrgDb =  "org.Hs.eg.db",
                          ont = "BP",
                          pAdjustMethod = "bonferroni",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = T)
data_GO_result <- data_GO@compareClusterResult

dotplot(data_GO_sim, showCategory = 3, font.size =10) +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12))


#GO-bp of DEGs in time point

allmarkers_time <- FindAllMarkers(harmonized_seurat,logfc.threshold = 0.25,min.pct = 0.25,only.pos = T)
group <- data.frame(gene = allmarkers_time$gene, group = allmarkers_time$cluster)
Gene_ID <- bitr(group$gene,fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = 'org.Hs.eg.db')
data <- merge(Gene_ID,group,by.x = "SYMBOL",by.y = "gene")

data_GO <- compareCluster(ENTREZID~group,
                          data = data,
                          fun = "enrichGO",
                          OrgDb =  "org.Hs.eg.db",
                          ont = "BP",
                          pAdjustMethod = "bonferroni",
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.05,
                          readable = T)
data_GO_sim <- clusterProfiler::simplify(data_GO, cutoff =0.7,by = "p.adjust",select_fun = min)

dotplot(data_GO_sim,showCategory = 8,font.size =10) +
  theme(panel.grid = element_blank(), 
        axis.title.x = element_text(size = 12))


