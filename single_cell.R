setwd("/home/bingxu/scrna_paper/aa/")

library(Seurat)
library(patchwork) 
library(dplyr)
library(future)
library(ggplot2)
library(harmony)
library(DoubletFinder)
inte@meta.data$orig.ident[inte@meta.data$orig.ident == "Col0"] <- "Col-0"

# 确认修改是否成功
unique(inte@meta.data$orig.ident)



plan(multisession, workers = 20)
options(future.globals.maxSize = 200 * 1024^3)
plan(sequential)

wt=Read10X("/home/bingxu/pap1_scrna/all/result/wt/output/filter_matrix/",gene.column=1)
mutant=Read10X("/home/bingxu/pap1_scrna/all/result/mutant/output/filter_matrix/",gene.column=1)
chl=c("accD",
"ArthCp001",
"ArthCp002",
"ArthCp003",
"ArthCp004",
"ArthCp005",
"ArthCp006",
"ArthCp007",
"ArthCp008",
"ArthCp009",
"ArthCp010",
"ArthCp011",
"ArthCp012",
"ArthCp013",
"ArthCp014",
"ArthCp015",
"ArthCp016",
"ArthCp017",
"ArthCp018",
"ArthCp019",
"ArthCp020",
"ArthCp021",
"ArthCp022",
"ArthCp023",
"ArthCp024",
"ArthCp025",
"ArthCp026",
"ArthCp027",
"ArthCp028",
"ArthCp029",
"ArthCp030",
"ArthCp031",
"ArthCp032",
"ArthCp033",
"ArthCp034",
"ArthCp035",
"ArthCp036",
"ArthCp037",
"ArthCp038",
"ArthCp039",
"ArthCp040",
"ArthCp041",
"ArthCp042",
"ArthCp043",
"ArthCp044",
"ArthCp045",
"ArthCp047",
"ArthCp048",
"ArthCp049",
"ArthCp050",
"ArthCp051",
"ArthCp052",
"ArthCp053",
"ArthCp054",
"ArthCp055",
"ArthCp056",
"ArthCp057",
"ArthCp058",
"ArthCp059",
"ArthCp060",
"ArthCp061",
"ArthCp062",
"ArthCp063",
"ArthCp064",
"ArthCp065",
"ArthCp066",
"ArthCp068",
"ArthCp069",
"ArthCp070",
"ArthCp071",
"ArthCp072",
"ArthCp073",
"ArthCp074",
"ArthCp075",
"ArthCp076",
"ArthCp077",
"ArthCp078",
"ArthCp079",
"ArthCp080",
"ArthCp081",
"ArthCp083",
"ArthCp084",
"ArthCp085",
"ArthCp086",
"Arthcp087",
"ArthCp088",
"ArthCr088",
"ArthCr089",
"ArthCr090",
"ArthCr091",
"ArthCt100",
"ArthCt112",
"atpA",
"atpB",
"atpE",
"atpF",
"atpH",
"atpI",
"ccsA",
"cemA",
"clpP",
"matK",
"ndhA",
"ndhB",
"ndhC",
"ndhD",
"ndhE",
"ndhF",
"ndhG",
"ndhH",
"ndhI",
"ndhJ",
"ndhK",
"petA",
"petB",
"petD",
"petG",
"petL",
"petN",
"psaA",
"psaB",
"psaC",
"psaI",
"psaJ",
"psbA",
"psbB",
"psbC",
"psbD",
"psbE",
"psbF",
"psbH",
"psbI",
"psbJ",
"psbK",
"psbL",
"psbM",
"psbN",
"psbT",
"psbZ",
"rbcL",
"rpl14",
"rpl16",
"rpl2",
"rpl20",
"rpl22",
"rpl23",
"rpl32",
"rpl33",
"rpl36",
"rpoA",
"rpoB",
"rpoC1",
"rpoC2",
"rps11",
"rps12",
"rps14",
"rps15",
"rps16",
"rps18",
"rps19",
"rps2",
"rps3",
"rps4",
"rps7",
"rps8",
"rrn16S",
"rrn23S",
"rrn4.5S",
"trnA",
"trnC",
"trnD",
"trnE",
"trnF",
"trnfM",
"trnG",
"trnH",
"trnI",
"trnK",
"trnL",
"trnM",
"trnN",
"trnP",
"trnQ",
"trnR",
"trnS",
"trnT",
"trnV",
"trnW",
"trnY",
"ycf1",
"ycf2",
"ycf3",
"ycf4")
chl=unique(chl)

wt <- wt[!rownames(wt) %in% chl, ]
mutant <- mutant[!rownames(mutant) %in% chl, ]

rownames(wt) <- gsub("_", "-", rownames(wt))  
rownames(mutant) <- gsub("_", "-", rownames(mutant))  

wt_object <- CreateSeuratObject(counts =wt, project = "Col-0", min.cells = 5, min.features = 1000)
mutant_object <- CreateSeuratObject(counts = mutant, project = "PAP1-D", min.cells =2, min.features = 1000)

wt_object <- subset(wt_object , features = setdiff(rownames(wt_object), chl))
mutant_object <- subset(mutant_object, features = setdiff(rownames(mutant_object), chl))
object.list <- list(wt_object, mutant_object)

# 获取每个细胞的orig.ident和cluster信息
data <- data.frame(orig.ident = inte$orig.ident, cluster = Idents(inte))

# 统计每个orig.ident和cluster的细胞数
cell_counts <- data %>%
  group_by(orig.ident, cluster) %>%
  summarise(cell_count = n())

# 画合并后的柱状图并标上数字
p <- ggplot(cell_counts, aes(x = cluster, y = cell_count, fill = orig.ident)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = cell_count), vjust = -0.3, position = position_dodge(0.9), size = 3.5) +
  labs(title = "Cell Counts per Cluster for Each orig.ident", x = "Cluster", y = "Cell Count") +
  theme_minimal()

ggsave(filename = "cells_cluster_number.pdf", plot = p, width = 12, height = 8, dpi = 300)
ggsave(filename = "cells_cluster_number.tiff", plot = p, width = 12, height = 8, dpi = 300)





# 2. 对每个数据集标准化和识别高变基因
object.list <- lapply(X = object.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

features <- SelectIntegrationFeatures(object.list = object.list)

# 3. 选择整合锚点(anchors)
anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)

# 4. 整合数据
integrated_object <- IntegrateData(anchorset = anchors)

inte=integrated_object
# 5. 切换到整合后的assay并进行标准分析
DefaultAssay(inte) <- "integrated"

inte <- FindVariableFeatures(inte, selection.method = "vst", nfeatures = 1000)
top10 <- head(VariableFeatures(inte), 10)


plot1 <- VariableFeaturePlot(inte)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

inte <- ScaleData(inte, verbose = FALSE)
inte <- RunPCA(inte, npcs = 40)
inte <- FindNeighbors(inte, reduction = "pca", dims = 1:40)

inte$source <- inte$orig.ident

# 运行harmony去除批次效应
inte <- RunHarmony(inte, "source")

inte2 <- RunPCA(inte, npcs = 50)
inte2 <- FindNeighbors(inte2, reduction = "pca", dims = 1:50)
a=ElbowPlot(inte2, ndims = 50)
ggsave(filename = "PCs.pdf", plot = a, width = 8, height = 8, dpi = 300)
ggsave(filename = "PCs.tiff", plot = a, width = 8, height = 8, dpi = 300)

inte2 <- RunPCA(inte, npcs = 40)
inte2 <- FindNeighbors(inte2, reduction = "pca", dims = 1:40)
resolutions <- seq(0.1, 1.2, by = 0.1)
for(res in resolutions) {
    inte <- FindClusters(inte2, resolution = res, 
                        algorithm = 1, 
                        random.seed = 42,
                        verbose = FALSE)
}
library(clustree)
reso=clustree(inte@meta.data, prefix = "integrated_snn_res.")

ggsave(filename = "resolution.pdf", plot = reso, width = 10, height = 10, dpi = 300)
ggsave(filename = "resolution.tiff", plot = reso, width = 10, height = 10, dpi = 300)






inte <- RunUMAP(inte, reduction = "pca", dims = 1:40)
inte <- RunTSNE(inte, reduction = "pca", dims = 1:40)
inte <- FindClusters(inte, resolution = 0.7)
DimPlot(inte, 
        reduction = "umap", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("UMAP") +theme_minimal()
umap=DimPlot(inte, 
        reduction = "umap", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("UMAP")
ggsave(filename = "umap.pdf", plot = umap, width = 8, height = 8, dpi = 300)
ggsave(filename = "umap.tiff", plot = umap, width = 8, height = 8, dpi = 300)
umap_fen=DimPlot(inte, 
        reduction = "umap", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("UMAP")
ggsave(filename = "umap_fen.pdf", plot = umap_fen, width = 14, height = 8, dpi = 300)
ggsave(filename = "umap_fen.tiff", plot = umap_fen, width = 14, height = 8, dpi = 300)

DimPlot(inte, 
        reduction = "tsne", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("TSNE") +theme_minimal()
tsne=DimPlot(inte, 
        reduction = "tsne", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("TSNE")
ggsave(filename = "tsne.pdf", plot = tsne, width = 8, height = 8, dpi = 300)
ggsave(filename = "tsne.tiff", plot = tsne, width = 8, height = 8, dpi = 300)
tsne_fen=DimPlot(inte, 
        reduction = "tsne", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) +
	ggtitle("TSNE")
ggsave(filename = "tsne_fen.pdf", plot = tsne_fen, width = 14, height = 8, dpi = 300)
ggsave(filename = "tsne_fen.tiff", plot = tsne_fen, width = 14, height = 8, dpi = 300)


library(Seurat)
library(ggplot2)
library(dplyr)
library(reshape2)

# 创建一个数据框，用于存放UMI counts、genes counts和reads counts
df <- data.frame(
  ident = inte$orig.ident,
  nCount_RNA = inte$nCount_RNA,
  nFeature_RNA = inte$nFeature_RNA
  # 如果有reads counts，可以添加相应的列
)

# 将数据转为长格式，便于绘图
df_long <- melt(df, id.vars = "ident")
df_long$type <- ifelse(df_long$variable == "nCount_RNA", "UMI counts per cell", "Genes per cell")
df_long$variable <- NULL

# 生成UMI和Genes的小提琴图
p1 <- ggplot(subset(df_long, type == "UMI counts per cell"), aes(x = ident, y = value, fill = ident)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  scale_y_continuous(trans = 'log10') +
  theme_classic() +
  labs(title = "Distribution of UMI counts per cell",
       y = "UMI Counts (log10 scale)",
       x = "") +
  theme(legend.position = "none")+
theme(
    plot.title = element_text(face = "bold", size = 16)  # 加粗并调大标题字体
  )

p2 <- ggplot(subset(df_long, type == "Genes per cell"), aes(x = ident, y = value, fill = ident)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, outlier.shape = NA) +
  theme_classic() +
  labs(title = "Distribution of Genes per cell",
       y = "Genes Counts",
       x = "") +
  theme(legend.position = "none")+
theme(
    plot.title = element_text(face = "bold", size = 16)  # 加粗并调大标题字体
  )

# 生成细胞数目柱形图
cell_count_df <- df %>%
  count(ident) %>%
  rename(cell_count = n)

p3 <- ggplot(cell_count_df, aes(x = ident, y = cell_count, fill = ident)) +
  geom_bar(stat = 'identity', position = position_dodge()) +
  theme_classic() +
  labs(title = "Cell Counts per Group",
       y = "Cell Counts",
       x = "") +
  theme(legend.position = "none")+
theme(
    plot.title = element_text(face = "bold", size = 16)  # 加粗并调大标题字体
  )

ggsave(filename = "UMI.pdf", plot = p1, width = 8, height = 8, dpi = 300)
ggsave(filename = "UMI.tiff", plot = p1, width = 8, height = 8, dpi = 300)
ggsave(filename = "UMI.tiff", plot = p1, width = 8, height = 8, dpi = 300)
ggsave(filename = "genes.pdf", plot = p2, width = 8, height = 8, dpi = 300)
ggsave(filename = "genes.tiff", plot = p2, width = 8, height = 8, dpi = 300)
ggsave(filename = "genes.tiff", plot = p2, width = 8, height = 8, dpi = 300)
ggsave(filename = "cells.pdf", plot = p3, width = 8, height = 8, dpi = 300)
ggsave(filename = "cells.tiff", plot = p3, width = 8, height = 8, dpi = 300)
ggsave(filename = "cells.tiff", plot = p3, width = 8, height = 8, dpi = 300)


cell_metadata <- inte@meta.data

# 统计每个 cluster 中各个 orig.ident 的细胞数目
cell_counts <- cell_metadata %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarize(cell_count = n(), .groups = "drop")

# 创建条形图并在每个柱子上添加数值标签
p <- ggplot(cell_counts, aes(x = factor(seurat_clusters), y = cell_count, fill = orig.ident)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = cell_count), 
            position = position_dodge(width = 0.9), 
            vjust = -0.5, 
            size = 3) +
  labs(title = "Number of Cells per Cluster by orig.ident",
       x = "Cluster",
       y = "Number of Cells") +
  theme_minimal()

ggsave(filename = "cluster_cell_nu.pdf", plot = p, width = 14, height = 8, dpi = 300)
ggsave(filename = "cluster_cell_nu.tiff", plot = p, width = 14, height = 8, dpi = 300)

metadata <- inte@meta.data

# 计算每细胞平均UMI数
average_umi_per_cell <- aggregate(metadata$nCount_RNA, by=list(metadata$orig.ident), FUN=mean)
colnames(average_umi_per_cell) <- c("orig.ident", "average_UMI_per_cell")

# 计算每个orig.ident的细胞总数
cell_counts <- table(metadata$orig.ident)
cell_counts <- as.data.frame(cell_counts)
colnames(cell_counts) <- c("orig.ident", "cell_count")

# 计算每细胞平均基因数
average_genes_per_cell <- aggregate(metadata$nFeature_RNA, by=list(metadata$orig.ident), FUN=mean)
colnames(average_genes_per_cell) <- c("orig.ident", "average_genes_per_cell")

# 输出结果
print("每细胞平均UMI数:")
print(average_umi_per_cell)

print("\n每个orig.ident的细胞总数:")
print(cell_counts)

print("\n每细胞平均基因数:")
print(average_genes_per_cell)































all_clusters_markers <- FindAllMarkers(inte)

# 将结果导出为CSV文件
write.csv(all_clusters_markers, file = "markers.csv")

DefaultAssay(inte) <- "RNA"

inte=JoinLayers(inte)

library(dplyr)
# 提取所有的 cluster
all_clusters <- unique(Idents(inte))
diff_genes_list <- list()

# 获取所有unique的cluster
clusters <- unique(Idents(inte))

# 对每个cluster进行循环分析
for(cluster in clusters) {
    # 子集化特定cluster的细胞
    inte_subset <- subset(inte, idents = cluster)
    
    # 在这个cluster中进行PAP1-D vs Col-0的差异分析
    diff_genes <- FindMarkers(inte_subset, 
                            ident.1 = "PAP1-D", 
                            ident.2 = "Col-0",
                            group.by = "orig.ident",
                            min.pct = 0.25,
                            logfc.threshold = 0.25)
    
    # 添加cluster信息列
    diff_genes$cluster <- cluster
    # 添加gene名称列
    diff_genes$gene <- rownames(diff_genes)
    
    # 只保留显著上调的基因（p_val_adj < 0.05 且 avg_log2FC > 0）
    up_genes <- diff_genes[diff_genes$p_val_adj < 0.05 & diff_genes$avg_log2FC > 0, ]
    
    # 将结果存储在列表中
    diff_genes_list[[paste0("cluster_", cluster)]] <- up_genes
}

# 合并所有cluster的结果
all_diff_genes <- do.call(rbind, diff_genes_list)

# 重新排序列，使gene和cluster在前面
all_diff_genes <- all_diff_genes[, c("gene", "cluster", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]

# 将合并后的结果保存为CSV文件
write.csv(all_diff_genes, 
          file = "all_clusters_PAP1D_vs_Col0_up.csv", 
          row.names = FALSE)

# 输出每个cluster中上调基因的数量统计
cluster_stats <- table(all_diff_genes$cluster)
print("Each cluster up-regulated genes count:")
print(cluster_stats)
















VlnPlot(inte, features = c("rna_AAP2"), split.by = "orig.ident", split.plot = TRUE,log=TRUE)

FeaturePlot(inte, features = c("rna_CHS"),split.by = "orig.ident", cols = c("grey", "red"))

VlnPlot(inte, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


VlnPlot(subset(inte, subset = orig.ident == "PAP1-D"), 
        features = c("rna_CHI"), 
        log = TRUE)

genes_of_interest <- c("rna_AT2G37040",
"rna_AT3G21240",
"rna_AT1G20440",
"rna_AT1G51680",
"rna_AT1G05680",
"rna_AT2G36790",
"rna_AT4G01070",
"rna_AT2G19450",
"rna_AT5G05010",
"rna_AT1G02930",
"rna_AT1G17170",
"rna_AT1G78380",
"rna_AT4G02520",
"rna_AT1G30400",
"rna_AT2G34660",
"rna_AT2G36910",
"rna_AT1G56650",
"rna_AT1G59870",
"rna_AT4G22880",
"rna_AT5G17220",
"rna_AT5G24520")



celltype <- c(
"0" = "Mesophyll", 
"1" = "Mesophyll", 
"2" = "Mesophyll", 
"3" = "Vasculature", 
"4" = "Mesophyll", 
"5" = "Mesophyll", 
"6" = "Mesophyll", 
"7" = "Epidermis", 
"8" = "Mesophyll",
"9" = "Mesophyll",
"10" = "Mesophyll",
"11" = "Epidermis",
"12" = "Companion")
inte_ident <- RenameIdents(inte, celltype)
inte_ident@meta.data$celltype <- Idents(inte_ident)


umap=DimPlot(inte_ident, 
        reduction = "umap", 
        repel = FALSE,
        label = TRUE) 
ggsave(filename = "ident.pdf", plot = umap, width = 8, height = 8, dpi = 300)
ggsave(filename = "ident.tiff", plot = umap, width = 8, height = 8, dpi = 300)
umap_fen=DimPlot(inte_ident, 
        reduction = "umap", 
        split.by = "orig.ident", 
        repel = FALSE,
        label = TRUE) 
ggsave(filename = "ident_fen.pdf", plot = umap_fen, width = 14, height = 8, dpi = 300)
ggsave(filename = "ident_fen.tiff", plot = umap_fen, width = 14, height = 8, dpi = 300)


celltype_detailed <- c(
    "0" = "Others-1",
    "1" = "Mesophyll-1",
    "2" = "Others-2",
    "3" = "Mesophyll-2",
    "4" = "Mesophyll-3",
    "5" = "Others-3",
    "6" = "Mesophyll-4",
    "7" = "Mesophyll-5",
    "8" = "Epidermis-1",
    "9" = "Vasculature-1",
    "10" = "Mesophyll-7",
    "11" = "Mesophyll-6",
    "12" = "Vasculature-2",
    "13" = "Vasculature-3",
    "14" = "Mesophyll-8",
    "15" = "Epidermis-2",
    "16" = "Companion cells",
    "17" = "Guard cells"
)

# 设置细胞类型顺序
celltype_levels <- c(
   "Others-1",
"Others-2",
"Others-3",
    "Mesophyll-1",
    
    "Mesophyll-2",
    "Mesophyll-3",
    
     "Mesophyll-4",
     "Mesophyll-5",
"Mesophyll-6",
"Mesophyll-7",
"Mesophyll-8",
     "Epidermis-1",
"Epidermis-2",
    
    
     "Vasculature-1",
     "Vasculature-2",
     "Vasculature-3",
     
     
     "Companion cells",
     "Guard cells"
)

# 重命名细胞类型
inte_ident <- RenameIdents(inte, celltype_detailed)
Idents(inte_ident) <- factor(Idents(inte_ident), levels = celltype_levels)
inte_ident@meta.data$celltype <- Idents(inte_ident)

# 设置配色方案
color_scheme3 <- c(
    # 间充质细胞 - 渐变的红色系
    "Mesophyll-1" = "#FF7F7F",
    "Mesophyll-2" = "#FF6666",
    "Mesophyll-3" = "#FF4D4D",
    "Mesophyll-4" = "#FF3333",
    "Mesophyll-5" = "#FF1A1A",
    "Mesophyll-6" = "#FF0000",
    "Mesophyll-7" = "#CC0000",
    "Mesophyll-8" = "#B20000",
    # 表皮细胞 - 绿色系
    "Epidermis-1" = "#90EE90",
    "Epidermis-2" = "#228B22",

    "Others-1" = "grey",
"Others-2" = "grey",
"Others-3" = "grey",
    
    
    # 维管组织 - 蓝色系
    "Vasculature-1" = "#87CEEB",
    "Vasculature-2" = "#4169E1",
# 韧皮部 - 紫色
    "Vasculature-3" = "#000080",
    
    # 伴细胞 - 棕色
    "Companion cells" = "#8B4513",
    
    # 保卫细胞 - 青色
    "Guard cells" = "#20B2AA"
)

# 绘制整体UMAP图
umap <- DimPlot(inte_ident, 
        reduction = "umap", 
        repel = TRUE,  # 使用repel
        label = TRUE,
        label.size = 4,
  # 添加文字背景框
        label.color = "black") + # 文字颜色
        scale_color_manual(values = color_scheme3) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right"
        )

# 绘制分组UMAP图
umap_fen <- DimPlot(inte_ident, 
        reduction = "umap", 
        split.by = "orig.ident", 
        repel = TRUE,
        label = TRUE,
        label.size = 4,

        label.color = "black") +
        scale_color_manual(values = color_scheme3) +
        theme_bw() +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "right"
        )

# 保存图片
ggsave(filename = "ident.pdf", plot = umap, width = 10, height = 8, dpi = 300)
ggsave(filename = "ident.tiff", plot = umap, width = 10, height = 8, dpi = 300)
ggsave(filename = "ident_fen.pdf", plot = umap_fen, width = 14, height = 8, dpi = 300)
ggsave(filename = "ident_fen.tiff", plot = umap_fen, width = 14, height = 8, dpi = 300)







# 绘制UMI在UMAP图上的密度图

library(ggExtra)

# 假设你的Seurat对象是inte
# 计算UMI数量并取log10
inte$logUMI <- log10(inte$nCount_RNA)

# 生成UMAP坐标
umap_coords <- Embeddings(inte, "umap")

# 创建数据框，包含UMAP坐标、logUMI数量和orig.ident
umap_df <- data.frame(UMAP_1 = umap_coords[,1], UMAP_2 = umap_coords[,2], logUMI = inte$logUMI, orig.ident = inte$orig.ident)

# 分别提取Col-0和PAP1-D的数据
umap_df_col0 <- umap_df[umap_df$orig.ident == "Col-0", ]
umap_df_pap1d <- umap_df[umap_df$orig.ident == "PAP1-D", ]

# 找到Col-0和PAP1-D的UMI最小值的点
min_umi_col0 <- umap_df_col0[which.min(umap_df_col0$logUMI), ]
min_umi_pap1d <- umap_df_pap1d[which.min(umap_df_pap1d$logUMI), ]

# 绘制Col-0的UMAP图
umap_plot_col0 <- ggplot(umap_df_col0, aes(x = UMAP_1, y = UMAP_2, color = logUMI)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "UMAP with log10(UMI) Density (Col-0)", x = "UMAP 1", y = "UMAP 2", color = "log10(UMI)") +
  geom_point(data = min_umi_col0, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 2) +  # 标记Col-0最小值的点
  geom_text(data = min_umi_col0, aes(x = UMAP_1, y = UMAP_2, label = paste("Col-0 Min UMI:", round(logUMI, 2))), vjust = -1, color = "red", size = 7)+
 theme(
    plot.title = element_text(face = "bold", size = 16)  # 加粗并调大标题字体
  )

# 绘制PAP1-D的UMAP图
umap_plot_pap1d <- ggplot(umap_df_pap1d, aes(x = UMAP_1, y = UMAP_2, color = logUMI)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "UMAP with log10(UMI) Density (PAP1-D)", x = "UMAP 1", y = "UMAP 2", color = "log10(UMI)") +
  geom_point(data = min_umi_pap1d, aes(x = UMAP_1, y = UMAP_2), color = "red", size = 2) +  # 标记PAP1-D最小值的点
  geom_text(data = min_umi_pap1d, aes(x = UMAP_1, y = UMAP_2, label = paste("PAP1-D Min UMI:", round(logUMI, 2))), vjust = -1, color = "red", size = 7)+
 theme(
    plot.title = element_text(face = "bold", size = 16)  # 加粗并调大标题字体
  )



ggsave(filename = "Col0_UMI.tiff", plot = umap_plot_col0 , width = 10, height = 8, dpi = 300)
ggsave(filename = "PAP1-D_UMI.tiff", plot = umap_plot_pap1d , width = 10, height = 8, dpi = 300)

















a <- "rna_UGT76E11"
gly=VlnPlot(inte, features = c(a), split.by = "orig.ident", split.plot = TRUE,log=TRUE)

# 保存为PDF
b <- paste(a, "vln.pdf", sep = "")
ggsave(b, gly, width = 8, height = 5, dpi = 300)
# 保存为PNG
c <- paste(a, "vln.tiff", sep = "")
ggsave(c, gly, width = 8, height = 5, dpi = 300)

gly <- FeaturePlot(inte,
                   features = c(a),
                   split.by = "orig.ident",
                   cols = c("grey", "red")) 
# 保存为PDF
b <- paste(a, ".pdf", sep = "")
ggsave(b, gly, width =10, height = 5, dpi = 300)
# 保存为PNG
c <- paste(a, ".tiff", sep = "")
ggsave(c, gly, width = 10, height = 5, dpi = 300)


a <- "rna_AT3G51240"
 FeaturePlot(inte,
                   features = c(a),
                   split.by = "orig.ident",
                   cols = c("grey", "red")) 





genes_df <- read.table(text = "cluster	gene
1	PAP17
1	NRT1.1
1	G3Pp1
3	RNS1
3	RBL14
3	ATSDI1
4	CPN60A
4	AT4G23493
4	AT4G10450
6	XTH24
6	CNI1
6	AT5G62150
7	FSD1
7	RBCS1A
7	COR15A
8	KCS1
8	GPAT2
8	AT3G47800
9	THA2
9	AT4G29100
9	AT1G61660
10	AT3G51730
10	AT3G47540
11	AT3G17110
11	AT1G72060
11	AT2G34810
11	PSAD-1
12	SIAR1
12	BZIP9
13	AT5G43620
13	RGLG1
13	AT1G20823
14	PETC
14	LHCB4.1
15	AT1G78830
15	AT1G05135
15	AT1G02360
16	SUC2
16	MT3
16	HA3
17	TGG2
17	AT5G44567
17	AT4G18280
", header = TRUE)

# 创建细胞类型数据框
cell_types <- data.frame(
  cluster = c(0,1,2,3, 4, 5, 6, 7, 8, 9,10, 11, 12,13,14,15,16,17),
  cell_type = c("others","mesophyll","others", "mesophyll", "mesophyll", "others", "mesophyll","mesophyll", "epidermis","vasculature", "mesophyll", "mesophyll", "vasculature","vasculature", "mesophyll","epidermis","companion cell","guard cell")
)

# 合并数据
genes_df <- merge(genes_df, cell_types, by = "cluster")

# 设置细胞类型顺序
desired_order <- c("mesophyll", "epidermis", "phloem","vasculature","companion cell", "guard cell","others")
genes_df$cell_type <- factor(genes_df$cell_type, levels = desired_order)
genes_df <- genes_df[order(genes_df$cell_type, genes_df$cluster),]
genes_to_plot <- genes_df$gene

cluster_order <- factor(c(1,3,4,6,7,11,10,14,8,15,9,12,13,16,17,0,2,5))

p <- DotPlot(inte,
    features = genes_to_plot,
    cols = c("lightgrey", "purple"),
    dot.scale = 8,
    cluster.idents = TRUE) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
    labs(x = "Genes", y = "Clusters") +
    scale_y_discrete(limits = rev(cluster_order)) +  # 添加这一行来指定y轴顺序
    coord_cartesian(ylim = c(0.5, length(unique(Idents(inte))) + 1))

# 后面的代码保持不变
cell_type_changes <- which(diff(as.integer(genes_df$cell_type)) != 0)
for(pos in cell_type_changes) {
    p <- p + geom_vline(xintercept = pos + 0.5, linetype = "dashed", color = "grey70", alpha = 0.5)
}

cell_type_labels <- tapply(1:nrow(genes_df), genes_df$cell_type, mean)
p <- p + annotate("text",
    x = unname(cell_type_labels),
    y = length(unique(Idents(inte))) + 1,
    label = names(cell_type_labels),
    size = 4)

p <- p + annotate("text",
    x = 1:nrow(genes_df),
    y = length(unique(Idents(inte))) + 0.5,
    label = paste("C", genes_df$cluster, sep=""),
    size = 4)
p
ggsave("markers_points.pdf", p, width = 18, height = 8, dpi = 300)
ggsave("markers_points.tiff", p, width = 18, height = 8, dpi = 300)



library(ggplot2)
library(cowplot)
library(ggrepel)  # 添加这个包
library(dplyr)    # 用于数据处理

# 创建一个空的列表来存储所有火山图
volcano_plots <- list()

# 获取所有cluster编号
clusters <- unique(inte@meta.data$seurat_clusters)

# 对每个cluster进行循环
for(cluster in clusters) {
    # 提取当前cluster的细胞
    cells_in_cluster <- WhichCells(inte, expression = seurat_clusters == cluster)
    cluster_obj <- subset(inte, cells = cells_in_cluster)
    
    # 计算差异基因
    markers <- FindMarkers(cluster_obj,
                          ident.1 = "PAP1-D",
                          ident.2 = "Col-0",
                          group.by = "orig.ident",
                          min.pct = 0.1,
                          logfc.threshold = 0)
    
    # 添加基因名和显著性列
    markers$gene <- rownames(markers)
    markers$significant <- ifelse(markers$p_val_adj < 0.05,
                                ifelse(markers$avg_log2FC > 0.25, "Up",
                                       ifelse(markers$avg_log2FC < -0.25, "Down", "NS")),
                                "NS")
    
    # 选择要标记的顶部基因
    top_genes <- markers %>% 
        filter(p_val_adj < 0.00001 & abs(avg_log2FC) > 1) %>%
        top_n(10, wt = abs(avg_log2FC))
    
    # 创建火山图
    p <- ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj), color = significant)) +
        geom_point(size = 1, alpha = 0.7) +
        scale_color_manual(values = c("Down" = "blue", "Up" = "red", "NS" = "grey")) +
        geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey") +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
        geom_text_repel(data = top_genes, 
                       aes(label = gene),
                       size = 3,
                       box.padding = 0.5,
                       max.overlaps = 20) +
        labs(title = paste("Cluster", cluster),
             x = "log2(Fold Change)",
             y = "-log10(adjusted p-value)") +
        theme_classic() +
        theme(plot.title = element_text(hjust = 0.5, face = "bold"),
              legend.position = "right")
    
    volcano_plots[[paste0("cluster_", cluster)]] <- p
}

# 计算需要的行数和列数
n_plots <- length(clusters)
n_cols <- 3  # 每行3个图
n_rows <- ceiling(n_plots/n_cols)

# 将所有火山图保存到一个PDF文件
pdf("all_clusters_volcano_plots.pdf", width = 15, height = 5*n_rows)
plot_grid(plotlist = volcano_plots, ncol = n_cols)
dev.off()
