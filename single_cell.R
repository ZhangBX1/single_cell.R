library(Seurat)
library(patchwork) 
library(dplyr)
early=readRDS("D://scavrrpm1//GSE226826_AvrRpm1_6h_peak.rds")
late=readRDS("D://scavrrpm1//GSE226826_AvrRpm1_9h_peak.rds")

early@active.assay <- "RNA"
late@active.assay <- "RNA"

early[["ATAC"]] <- NULL
late[["ATAC"]] <- NULL


early_new <- CreateSeuratObject(counts = GetAssayData(early, slot = "counts", assay = "RNA"))
late_new <- CreateSeuratObject(counts = GetAssayData(late, slot = "counts", assay = "RNA"))
early_new$orig.ident <- "6h"
late_new$orig.ident <- "9h"
pbmc <- merge(early_new, y = late_new, add.cell.ids = c("6h", "9h"))


pbmc@active.assay <- "RNA"

pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Scale the data
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

pbmc <- JoinLayers(pbmc)


# 查看最高变的10个基因
top10 <- head(VariableFeatures(pbmc), 10)

# 画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) #1个PC 500个细胞

DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) #15个PC

pbmc <- FindNeighbors(pbmc, dims = 1:50)
pbmc <- FindClusters(pbmc, resolution = 0.5)

pbmc <- RunUMAP(pbmc, dims = 1:50)
DimPlot(pbmc, 
        reduction = "umap", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) 


DimPlot(pbmc, reduction = "umap")
# 显示在聚类标签

# 使用TSNE聚类
pbmc <- RunTSNE(pbmc, dims = 1:50)
DimPlot(pbmc, reduction = "tsne")
# 显示在聚类标签
DimPlot(pbmc, reduction = "tsne", label = TRUE)
pbmc <- JoinLayers(pbmc)


# cluster 1的标记基因
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
head(cluster1.markers, n = 5)
#找出区分cluster 5与cluster 0和cluster 3的所有标记
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
head(cluster5.markers, n = 5)
# 找出每个cluster的标记与所有剩余的细胞相比较，只报告阳性细胞
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

VlnPlot(pbmc, features = c("MS4A1", "CD79A"))

VlnPlot(pbmc, features = c("RIN4"))

FeaturePlot(pbmc, features = c("RIN4"))
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

VlnPlot(pbmc, features = c("RRTF1"),split.by = "orig.ident")

new_cluster_ids <- c("Mesophyll", "Epidermis", "Mesophyll", "Stomatal pore","Epidermis","Mesophyll","Bundle sheath","Other cells","Mesophyll","Epidermis")
names(new_cluster_ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new_cluster_ids)







