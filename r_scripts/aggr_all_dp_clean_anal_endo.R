library('Seurat')
library('plyr')
library('dplyr')


read_aggr_all_matrix <- readMM('matrix_Endo_norep1.mtx')
aggr_all_matrix <- as.matrix(read_aggr_all_matrix)
Genes <- read.delim('genes_Endo_norep1.tsv', sep = '\t', header = FALSE)
Barcodes <- read.delim('barcodes_Endo_norep1.tsv', sep = '\t', header = FALSE)
colnames(aggr_all_matrix) <- (Barcodes$V1)
rownames(aggr_all_matrix) <- make.names(Genes$V1, unique = TRUE)

aggr_all_seurat <- CreateSeuratObject(raw.data = aggr_all_matrix, min.cells = 50, project = 'Endothelial cells')

mito.genes <- grep(pattern = "^mt.", x = rownames(x = aggr_all_seurat@data), value = TRUE)
percent.mito <- Matrix::colSums(aggr_all_seurat@raw.data[mito.genes, ])/Matrix::colSums(aggr_all_seurat@raw.data)

aggr_all_seurat <- AddMetaData(object = aggr_all_seurat, metadata = percent.mito, col.name = "percent.mito")

Library_ID <- read.delim('Endo_norep1.csv', sep = ",", header = TRUE)
CellsMeta = aggr_all_seurat@meta.data
Library_ID$nGene <- NULL
Library_ID$nUMI <- NULL
rownames(Library_ID) <- Library_ID$Barcode
Library_ID$Barcode <- NULL

#head(CellsMeta)
#CellsMeta["LibraryID"] <- Library_ID$LibraryID
#head(CellsMeta)
#CellsMeta["CellType"] <- Library_ID$CellType
#head(CellsMeta)
#CellsMetaTrim <- subset(CellsMeta, select = c("LibraryID", "CellType"))
#head(CellsMetaTrim)

CellsMeta2 <- merge(CellsMeta, Library_ID, by=0, sort = FALSE)
rownames(CellsMeta2) <- CellsMeta2$Row.names
CellsMeta2$Row.names <- NULL
aggr_all_seurat <- AddMetaData(aggr_all_seurat, CellsMeta2)
head(aggr_all_seurat@meta.data)

pdf(file = "endo_nGene_nUMI_percent_mito.pdf")
VlnPlot(object = aggr_all_seurat, features.plot = c('nGene', 'nUMI', 'percent.mito'), nCol = 3)
dev.off()

pdf(file = "endo_filter1.pdf")
GenePlot(object = aggr_all_seurat, gene1 = 'nUMI', gene2 = 'nGene')
dev.off()

pdf(file = "endo_filter2.pdf")
GenePlot(object = aggr_all_seurat, gene1 = 'nUMI', gene2 = 'percent.mito')
dev.off()

aggr_all_seurat <- FilterCells(object = aggr_all_seurat, subset.names = c("nGene", "percent.mito"), 
                               low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

aggr_all_seurat <- NormalizeData(object = aggr_all_seurat, normalization.method = "LogNormalize", 
                                 scale.factor = 10000)

pdf(file = "endo_variable_genes.pdf")
aggr_all_seurat <- FindVariableGenes(object = aggr_all_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, y.high.cutoff = 20)
dev.off()
length(x = aggr_all_seurat@var.genes)

aggr_all_seurat <- ScaleData(object = aggr_all_seurat, vars.to.regress = c("nUMI", 'percent.mito'))

aggr_all_seurat <- RunPCA(object = aggr_all_seurat, pcs.compute = 40, pc.genes = aggr_all_seurat@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)

pdf(file = "endo_pca_genes.pdf")
VizPCA(object = aggr_all_seurat, pcs.use = 1:2)
dev.off()

pdf(file = "endo_pca.pdf")
PCAPlot(object = aggr_all_seurat, dim.1 = 1, dim.2 = 2, no.legend = TRUE)
dev.off()

pdf(file = "endo_pca_heat.pdf")
PCHeatmap(object = aggr_all_seurat, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
dev.off()

pdf(file = "endo_pca_elbow.pdf")
PCElbowPlot(object = aggr_all_seurat, num.pc = 40)
dev.off()

aggr_all_seurat <- FindClusters(object = aggr_all_seurat, reduction.type = "pca", dims.use = 1:30, 
                                resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

aggr_all_seurat <- RunTSNE(object = aggr_all_seurat, dims.use = 1:30, do.fast = TRUE)
pdf(file = "endo_tsne.pdf")
TSNEPlot(object = aggr_all_seurat)
dev.off()

pdf(file = "endo_marker_pop.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = c('Ptprb', 'Stab2', 'Vwf'), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

all_markers <- FindAllMarkers(object = aggr_all_seurat, only.pos = TRUE, min.pct = 0.60, logfc.threshold = 2)

top10 <- all_markers %>% group_by(cluster) %>% top_n(15, avg_logFC)

heatmap_markers <- read.delim('heatmap_markers.csv', sep = ",", header = FALSE)

pdf(file = "marker_heatmap.pdf")
DoHeatmap(object = aggr_all_seurat, genes.use = heatmap_markers$V1, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top', remove.key = TRUE)
dev.off()

pdf(file = "marker_heatmap_color.pdf")
DoHeatmap(object = aggr_all_seurat, genes.use = heatmap_markers$V1, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top')
dev.off()

png(file = "marker_heatmap.png", width = 13000, height = 10000, res = 1200)
DoHeatmap(object = aggr_all_seurat, genes.use = heatmap_markers$V1, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top', remove.key = TRUE)
dev.off()

png(file = "marker_heatmap_color.png", width = 13000, height = 10000, res = 1200)
DoHeatmap(object = aggr_all_seurat, genes.use = heatmap_markers$V1, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top')
dev.off()

tiff(file = "marker_heatmap.tiff", width = 4000, height = 3000, res = 400)
DoHeatmap(object = aggr_all_seurat, genes.use = heatmap_markers$V1, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top', remove.key = TRUE, cex.row = 15)
dev.off()

tiff(file = "marker_heatmap_color.tiff", width = 4000, height = 3000, res = 400)
DoHeatmap(object = aggr_all_seurat, genes.use = heatmap_markers$V1, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top')
dev.off()

cluster_unk <- FindMarkers(object = aggr_all_seurat, ident.1 = 6, ident.2 = c(0, 1, 2, 3, 4, 5, 7, 8), min.pct = 0.25)
unk_print <- row.names(cluster_unk)[1:8]

pdf(file = "Endo_unk_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = unk_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_lrat <- FindMarkers(object = aggr_all_seurat, ident.1 = c(2, 3), ident.2 = c(0, 1, 4, 5, 8), min.pct = 0.25)
lrat_print <- row.names(cluster_lrat)[1:12]

pdf(file = "Endo_lrat_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = lrat_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_endo <- FindMarkers(object = aggr_all_seurat, ident.1 = c(1, 5), ident.2 = c(0, 2, 3, 4, 8), min.pct = 0.25)
endo_print <- row.names(cluster_endo)[1:12]

pdf(file = "Endo_endo_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = endo_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_endo <- FindMarkers(object = aggr_all_seurat, ident.1 = c(0, 4), ident.2 = c(1, 2, 3, 5, 8), min.pct = 0.25)
endo_print <- row.names(cluster_endo)[1:12]

pdf(file = "Endo_endo_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = endo_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

aggr_all_seurat <- RunUMAP(object = aggr_all_seurat, reduction.use = "pca", dims.use = 1:10, n_neighbors = 10)

pdf(file = "Endo_umap.pdf")
DimPlot(object = aggr_all_seurat, reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "Endo_marker_pop_umap.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = c('Lrat', 'Ptprb', 'Adgre1', 'Alb'), cols.use = c("grey85", "red1"), reduction.use = "umap", no.legend = FALSE)
dev.off()

pdf(file = "Endo_umap_library.pdf")
DimPlot(object = aggr_all_seurat, group.by = 'LibraryID', reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "Endo_umap_cluster.pdf")
DimPlot(object = aggr_all_seurat, group.by = 'ident', reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

save.image('aggr_Endo_clean.RData')

#Subset for slingshot analysis
aggr_all_subset <- SubsetData(aggr_all_seurat, subset.name = 'Stab2', accept.low = 2)

pdf(file = "subset_endo_variable_genes.pdf")
aggr_all_subset <- FindVariableGenes(object = aggr_all_subset, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, y.high.cutoff = 20)
dev.off()
length(x = aggr_all_subset@var.genes)

aggr_all_subset <- RunPCA(object = aggr_all_subset, pcs.compute = 40, pc.genes = aggr_all_subset@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)

pdf(file = "subset_endo_pca_elbow.pdf")
PCElbowPlot(object = aggr_all_subset, num.pc = 40)
dev.off()

aggr_all_subset <- FindClusters(object = aggr_all_subset, reduction.type = "pca", dims.use = 1:22, 
                                resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

aggr_all_subset <- RunTSNE(object = aggr_all_subset, dims.use = 1:22, do.fast = TRUE)
pdf(file = "subset_endo_tsne.pdf")
TSNEPlot(object = aggr_all_subset)
dev.off()

aggr_all_subset <- RunUMAP(object = aggr_all_subset, reduction.use = "pca", dims.use = 1:22, n_neighbors = 10)

pdf(file = "subset_endo_umap_library.pdf")
DimPlot(object = aggr_all_subset, group.by = 'LibraryID', reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Convert AnnData to Seurat
library(reticulate)
ad <- import("anndata", convert = FALSE)
endo_seurat <- ad$read_h5ad('scanpy_anal_DE.h5ad')
endo_seurat_anal <- Convert(endo_seurat, to = "seurat")
endo_seurat_anal <- SetAllIdent(endo_seurat_anal, id = 'louvain')

#Convert to SCE
seurat_sce <- Convert(from = aggr_all_seurat, to = "sce")

#Run Diffusionmap
aggr_all_seurat <- RunDiffusion(aggr_all_seurat, dims.use = 1:25, max.dim = 4)

pdf(file = 'seurat_difmap_cluster.pdf')
DMPlot(aggr_all_seurat, group.by = 'ident')
dev.off()

pdf(file = 'seurat_difmap_libraryid.pdf')
DMPlot(aggr_all_seurat, group.by = 'LibraryID')
dev.off()


aggr_diff_dm1 <- GetCellEmbeddings(object = aggr_all_seurat, reduction.type = "dm", dims.use = 1, cells.use = NULL)
aggr_diff_dm2 <- GetCellEmbeddings(object = aggr_all_seurat, reduction.type = "dm", dims.use = 2, cells.use = NULL)

pdf(file = 'seurat_difmap_libraryid.pdf')
plot(x=aggr_diff_dm1, y=aggr_diff_dm2, pch=20, col=aggr_all_seurat@meta.data$LibraryID)
dev.off()

#Slingshot
library(destiny)
dm <- DiffusionMap(seurat_sce)

pdf(file = 'destiny_diffmap.pdf')
plot(dm, 1:2, col_by ='LibraryID', legend_main ='Treatment', pal = c('#c02b2b', '#c02b2b', '#c0b82b', '#c0b82b', '#c0b82b', '#5ec02b', '#5ec02b', '#5ec02b'))
dev.off()

pdf(file = 'destiny_diffmap_3d.pdf')
plot(dm, pch = 20, col_by ='LibraryID', legend_main ='Treatment', pal = c('#c02b2b', '#c02b2b', '#c0b82b', '#c0b82b', '#c0b82b', '#5ec02b', '#5ec02b', '#5ec02b'))
dev.off()

pdf(file = 'destiny_diffmap_ggplot.pdf')
qplot(DC1, DC2, data = dm, colour = LibraryID) + scale_color_cube_helix() + 
  scale_color_manual(breaks=c('Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3', 'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3', 'CCl4 4w - 2', 'CCl4 4w - 3'), 
                     values=c('#c02b2b', '#c02b2b', '#c02b2b', '#c0b82b', '#c0b82b', '#5ec02b', '#5ec02b', '#5ec02b'))

dev.off()

library(slingshot)
library(RColorBrewer)
library(scater)
library(scales)

seurat_sce <- runDiffusionMap(seurat_sce, ncomponents = 4, ntop = 1000)
pdf(file = 'endo_scat_difmap.pdf')
plotDiffusionMap(seurat_sce, colour_by = "Col1a1")
dev.off()


louvain <- read.delim('endo_louvain.csv', sep = ",", header = TRUE)
rownames(louvain) <- louvain$barcode
louvain$barcode = NULL

CellsMeta = colData(seurat_sce)
head(CellsMeta)

CellsMeta2 <- merge(CellsMeta, louvain, by=0, sort = FALSE)
louvain <- as.data.frame(CellsMeta2$louvain)
names(louvain)[1] <- 'louvain'
seurat_sce2 <- mutate(seurat_sce, as.factor(louvain$louvain))

seurat_sub <- subset(seurat_sce2, as.factor.louvain.louvain.==c('0', '1', '2', '3'))

slingshot_sce <- slingshot(seurat_sce2, clusterLabels = 'as.factor.louvain.louvain.', reducedDim = 'DiffusionMap')
summary(slingshot_sce$slingPseudotime_1)

colors <- colorRampPalette(palette(c('#E3E3E3', '#CE7373', '#C02424', '#8F1A1A', '#681212', '#4B0C0C'))[-6])(12)
pdf(file = 'endo_slingshot_scat_color.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(colors[cut(slingshot_sce$slingPseudotime_1, breaks=12)], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
#lines(slingshot_sce@int_metadata$slingshot@curves$curve1, lwd=2, col = '#8B1A1A')
dev.off()

color <- palette(c('#A7414A', '#A37C27', '#A37C27', '#A37C27', '#A7414A', '#A7414A', '#6A8A82', '#6A8A82'))
pdf(file = 'endo_slingshot_scat_library.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(color[slingshot_sce$LibraryID], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

#color1 <- palette(c('#33689E', '#B22222', '#B47746', '#50A4A8', '#4CBBC7', '#6E508C', '#A171B0'))
pdf(file = 'endo_slingshot_scat_ident_louvain.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(color1[slingshot_sce$as.factor.louvain.louvain.], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

pdf(file = 'slingshot_seurat_cluster_scat.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = brewer.pal(9,'Set1')[slingshot_sce$ident], pch=16, asp = 1)
lines(SlingshotDataSet(slingshot_sce), lwd=2, type = 'lineages')
dev.off()

pdf(file = 'slingshot_pairs_scat.pdf')
pairs(SlingshotDataSet(slingshot_sce), type="curves", col = slingshot_sce$ident)
dev.off()

#Subset genes
out <- as.data.frame(row.names(aggr_all_seurat@meta.data))
names(out)[1] <- c("Barcode")
out['nGene'] <- aggr_all_seurat@meta.data$nGene
out['nUMI'] <- aggr_all_seurat@meta.data$nUMI
out['LibraryID'] <- aggr_all_seurat@meta.data$LibraryID
out2 <- as.data.frame(aggr_all_seurat@ident)
out['CellType'] <- out2$'aggr_all_seurat@ident'
write.table(out, file = 'aggr_all_3_pops.txt', sep='\t', col.names=NA)
write.table(out$Barcode, file = 'barcodes_3_pops.tsv', sep='\t', col.names=FALSE, row.names=FALSE)
write.table(rownames(aggr_all_matrix), file = 'genes_3_pops.tsv', sep='\t', col.names=FALSE, row.names=FALSE)


out <- as.data.frame(rownames(aggr_all_seurat@data))
names(out)[1] <- c("Genes")
aggr_all_matrix2 <- as.data.frame(aggr_all_matrix)
new <- aggr_all_matrix %>% filter(one_of(dput(as.character(out$Genes))))
aggr_all_matrix[!rownames(aggr_all_matrix) %in% out$Genes, , drop = FALSE]
new <- aggr_all_matrix2 %>% filter(rownames(aggr_all_matrix2) %in% as.character(out$Genes))
library(data.table)
new <- setDT(aggr_all_matrix)[rownames(aggr_all_matrix) %chin% rownames(out)]
new <- as.matrix(new)
rownames(new) <- c()
colnames(new) <- NULL
new1 <- as(new, "sparseMatrix")
writeMM(new1, 'matrix_3_pops.mtx')

outEndo <- out[out$CellType == 'Hepatic stellate cells',]
outEndo <- out[out$CellType == 'Endothelial cells',]
outendo <- out[out$CellType == 'endophages',]

Library_ID$LibraryID_agg <- replace(as.character(Library_ID$LibraryID_agg), Library_ID$LibraryID_agg == "Ctrl - 1", "Ctrl")
Library_ID$LibraryID_agg <- replace(as.character(Library_ID$LibraryID_agg), Library_ID$LibraryID_agg == "Ctrl - 2", "Ctrl")
Library_ID$LibraryID_agg <- replace(as.character(Library_ID$LibraryID_agg), Library_ID$LibraryID_agg == "Ctrl - 3", "Ctrl")
Library_ID$LibraryID_agg <- replace(as.character(Library_ID$LibraryID_agg), Library_ID$LibraryID_agg == "CCl4 2w - 1", "CCl4 2w")
Library_ID$LibraryID_agg <- replace(as.character(Library_ID$LibraryID_agg), Library_ID$LibraryID_agg == "CCl4 2w - 2", "CCl4 2w")
Library_ID$LibraryID_agg <- replace(as.character(Library_ID$LibraryID_agg), Library_ID$LibraryID_agg == "CCl4 2w - 3", "CCl4 2w")
Library_ID$LibraryID_agg <- replace(as.character(Library_ID$LibraryID_agg), Library_ID$LibraryID_agg == "CCl4 4w - 1", "CCl4 4w")
Library_ID$LibraryID_agg <- replace(as.character(Library_ID$LibraryID_agg), Library_ID$LibraryID_agg == "CCl4 4w - 2", "CCl4 4w")
Library_ID$LibraryID_agg <- replace(as.character(Library_ID$LibraryID_agg), Library_ID$LibraryID_agg == "CCl4 4w - 3", "CCl4 4w")

aggr_all_matrix <- as.data.frame(aggr_all_matrix)
newEndo <- aggr_all_matrix %>% select(one_of(dput(as.character(outEndo$Barcode))))
newEndo <- as.matrix(newEndo)
rownames(newEndo) <- c()
colnames(newEndo) <- NULL
newEndo1 <- as(newEndo, "sparseMatrix")
writeMM(newEndo1, 'matrix_Endo.mtx')

write.table(outEndo, file = 'aggr_all_Endo.txt', sep='\t', col.names=NA)
write.table(outEndo$Barcode, file = 'barcodes_Endo.tsv', sep='\t', col.names=FALSE, row.names=FALSE)
write.table(rownames(aggr_all_matrix), file = 'genes_Endo.tsv', sep='\t', col.names=FALSE, row.names=FALSE)

#Generate sorted (pc1/pc2) genelist with description
library("annotables")

Gene1 <- as.data.frame(cluster_unk)
Gene1[6] <- rownames(cluster_unk)
names(Gene1)[6] <- c("symbol")
Gene2 <- merge(Gene1, grcm38, by="symbol")
View(Gene2)
Gene3 <- Gene2[order(-Gene2$pct.1),]
write.table(Gene3, file = "genelist_unk.txt", sep = "\t")

Gene3[15] <- Gene3$pct.1/Gene3$pct.2
names(Gene3)[15] <- c("pct1/pct2")
Gene4 <- Gene3[order(-Gene3$`pct1/pct2`),]
write.table(Gene4, file = "genelist2.txt", sep = "\t")

Endo_markers_list2 <- as.data.frame(Endo_markers2)
Endo_markers_list2[6] <- rownames(Endo_markers_list2)
names(Endo_markers_list2)[6] <- c("symbol")
Endo_markers_list2 <- merge(Endo_markers_list2, grcm38, by="symbol")
View(Endo_markers_list2)
Endo_markers_list2 <- Endo_markers_list2[order(-Endo_markers_list2$pct.1),]
write.table(Endo_markers_list2, file = "Endo_markers_list2_pct1.txt", sep = "\t")

Endo_markers_list2[15] <- Endo_markers_list2$pct.1/Endo_markers_list2$pct.2
names(Endo_markers_list2)[15] <- c("pct1/pct2")
Endo_markers_list2 <- Endo_markers_list2[order(-Endo_markers_list2$`pct1/pct2`),]
write.table(Endo_markers_list2, file = "Endo_markers_list2.txt", sep = "\t")

Endo_markers_list2_genes <- scan("genelist2.txt", what="", sep="\n")

Endo_markers_list2_genes <- split(Endo_markers_list2_genes, ceiling(seq_along(Endo_markers_list2_genes)/15))

pdf(file = "Endo_genes.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = Endo_markers_list2_genes$`1`, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE, nCol = 7)
dev.off()

FeaturePlot(object = Lrat_tSNE, features.plot = c("Colec10", "Gas7", "Fgfr2", "Cdh2", "Ifitm1", "Tln2", "Ank3", "Tmem56", "Wt1", "Ctsk", "Htra1", "Dpt", "Aebp1", "Ltbp2"), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE, nCol = 7)

