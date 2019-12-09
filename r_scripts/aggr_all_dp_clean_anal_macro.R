library('Seurat')
library('plyr')
library('dplyr')


read_aggr_all_matrix <- readMM('matrix_Macro_norep1.mtx')
aggr_all_matrix <- as.matrix(read_aggr_all_matrix)
Genes <- read.delim('genes_Macro_norep1.tsv', sep = '\t', header = FALSE)
Barcodes <- read.delim('barcodes_Macro_norep1.tsv', sep = '\t', header = FALSE)
colnames(aggr_all_matrix) <- (Barcodes$V1)
rownames(aggr_all_matrix) <- make.names(Genes$V1, unique = TRUE)

aggr_all_seurat <- CreateSeuratObject(raw.data = aggr_all_matrix, min.cells = 50, project = 'Macrophages')

mito.genes <- grep(pattern = "^mt.", x = rownames(x = aggr_all_seurat@data), value = TRUE)
percent.mito <- Matrix::colSums(aggr_all_seurat@raw.data[mito.genes, ])/Matrix::colSums(aggr_all_seurat@raw.data)

aggr_all_seurat <- AddMetaData(object = aggr_all_seurat, metadata = percent.mito, col.name = "percent.mito")

Library_ID <- read.delim('Macro_norep1.csv', sep = ",", header = TRUE)
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

pdf(file = "macro_nGene_nUMI_percent_mito.pdf")
VlnPlot(object = aggr_all_seurat, features.plot = c('nGene', 'nUMI', 'percent.mito'), nCol = 3)
dev.off()

pdf(file = "macro_filter1.pdf")
GenePlot(object = aggr_all_seurat, gene1 = 'nUMI', gene2 = 'nGene')
dev.off()

pdf(file = "macro_filter2.pdf")
GenePlot(object = aggr_all_seurat, gene1 = 'nUMI', gene2 = 'percent.mito')
dev.off()

aggr_all_seurat <- FilterCells(object = aggr_all_seurat, subset.names = c("nGene", "percent.mito"), 
                               low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

aggr_all_seurat <- NormalizeData(object = aggr_all_seurat, normalization.method = "LogNormalize", 
                                 scale.factor = 10000)

pdf(file = "macro_variable_genes.pdf")
aggr_all_seurat <- FindVariableGenes(object = aggr_all_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, y.high.cutoff = 20)
dev.off()
length(x = aggr_all_seurat@var.genes)

aggr_all_seurat <- ScaleData(object = aggr_all_seurat, vars.to.regress = c("nUMI", 'percent.mito'))

aggr_all_seurat <- RunPCA(object = aggr_all_seurat, pcs.compute = 40, pc.genes = aggr_all_seurat@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)

pdf(file = "macro_pca_genes.pdf")
VizPCA(object = aggr_all_seurat, pcs.use = 1:2)
dev.off()

pdf(file = "macro_pca.pdf")
PCAPlot(object = aggr_all_seurat, dim.1 = 1, dim.2 = 2, no.legend = TRUE)
dev.off()

pdf(file = "macro_pca_heat.pdf")
PCHeatmap(object = aggr_all_seurat, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
dev.off()

pdf(file = "macro_pca_elbow.pdf")
PCElbowPlot(object = aggr_all_seurat, num.pc = 40)
dev.off()

aggr_all_seurat <- FindClusters(object = aggr_all_seurat, reduction.type = "pca", dims.use = 1:22, 
                                resolution = 0.02, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

aggr_all_seurat <- RunTSNE(object = aggr_all_seurat, dims.use = 1:22, do.fast = TRUE)
pdf(file = "macro_tsne.pdf")
TSNEPlot(object = aggr_all_seurat)
dev.off()

pdf(file = "macro_marker_pop.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = c('Adgre1', 'Lyz2', 'Clec4f'), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
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

pdf(file = "Macro_unk_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = unk_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_lrat <- FindMarkers(object = aggr_all_seurat, ident.1 = c(2, 3), ident.2 = c(0, 1, 4, 5, 8), min.pct = 0.25)
lrat_print <- row.names(cluster_lrat)[1:12]

pdf(file = "Macro_lrat_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = lrat_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_macro <- FindMarkers(object = aggr_all_seurat, ident.1 = c(1, 5), ident.2 = c(0, 2, 3, 4, 8), min.pct = 0.25)
macro_print <- row.names(cluster_macro)[1:12]

pdf(file = "Macro_macro_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = macro_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_endo <- FindMarkers(object = aggr_all_seurat, ident.1 = c(0, 4), ident.2 = c(1, 2, 3, 5, 8), min.pct = 0.25)
endo_print <- row.names(cluster_endo)[1:12]

pdf(file = "Macro_endo_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = endo_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

aggr_all_seurat <- RunUMAP(object = aggr_all_seurat, reduction.use = "pca", dims.use = 1:22, n_neighbors = 10)

pdf(file = "Macro_umap.pdf")
DimPlot(object = aggr_all_seurat, reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "Macro_marker_pop_umap.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = c('Adgre1', 'Lyz2', 'Clec4f'), cols.use = c("grey85", "red1"), reduction.use = "umap", no.legend = FALSE)
dev.off()

pdf(file = "Macro_umap_library.pdf")
DimPlot(object = aggr_all_seurat, group.by = 'LibraryID', reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

save.image('aggr_Macro_clean.RData')

#Heatmap of goi 
Louvain$Louvain[Louvain$Louvain==4] <- 5
Louvain$Louvain[Louvain$Louvain==3] <- 4
Louvain$Louvain[Louvain$Louvain==2] <- 3
Louvain$Louvain[Louvain$Louvain==1] <- 2
Louvain$Louvain[Louvain$Louvain==0] <- 1

aggr_all_seurat <- AddMetaData(aggr_all_seurat, Louvain)

pdf(file = "macro_markers_cluster_4.pdf")
DoHeatmap(object = aggr_all_seurat, genes.use = c('Clec4f', 'Vsig4', 'Ccr2', 'Cx3cr1', 'Hmgb2', 'Stmn1', 'Birc5', 'Ube2c', 'Hmgn2', 'Cenpa', 'Ccnb2', 'Cks2', 'Mki67', 'Racgap1', 'Ccna2'), 
          slim.col.label = TRUE, col.high = 'firebrick', col.mid = '#EBEBD2', col.low = 'steelblue', group.label.loc = 'top', group.by = 'Louvain', draw.line = TRUE, remove.key = FALSE)
dev.off()

library(gplots)
aggr_all_sub <- subset(aggr_all_seurat@scale.data, rownames(aggr_all_seurat@scale.data) %in% c('Clec4f', 'Vsig4', 'Lyz2', 'Adgre1', 'Fcgr3', 'Cd74', 'Ccr2', 'Cx3cr1', 'Hmgb2', 'Stmn1', 'Birc5', 'Ube2c', 'Hmgn2', 'Cenpa', 'Ccnb2', 'Cks2', 'Mki67', 'Racgap1', 'Ccna2'))
aggr_all_sub <- aggr_all_sub[match(c('Clec4f', 'Vsig4', 'Lyz2', 'Adgre1', 'Fcgr3', 'Cd74', 'Ccr2', 'Cx3cr1', 'Hmgb2', 'Stmn1', 'Birc5', 'Ube2c', 'Hmgn2', 'Cenpa', 'Ccnb2', 'Cks2', 'Mki67', 'Racgap1', 'Ccna2'), rownames(aggr_all_sub)),]

Louvain$cellid <- rownames(Louvain)
cluster_one_val <- Louvain[Louvain$Louvain == '1',]
cluster_two_val <- Louvain[Louvain$Louvain == '2',]
cluster_three_val <- Louvain[Louvain$Louvain == '3',]
cluster_four_val <- Louvain[Louvain$Louvain == '4',]
cluster_five_val <- Louvain[Louvain$Louvain == '5',]
heat_mean <- as.data.frame(rowMeans(subset(aggr_all_sub, select = cluster_one_val$cellid)))
names(heat_mean)[1] <- '1'
heat_mean$'2' <- rowMeans(subset(aggr_all_sub, select = cluster_two_val$cellid))
heat_mean$'3' <- rowMeans(subset(aggr_all_sub, select = cluster_three_val$cellid))
heat_mean$'4' <- rowMeans(subset(aggr_all_sub, select = cluster_four_val$cellid))
heat_mean$'5' <- rowMeans(subset(aggr_all_sub, select = cluster_five_val$cellid))

heat_mean_mat <- as.matrix(heat_mean)

col <- colorRampPalette(palette(c('steelblue', '#EBEBD2', 'firebrick')))(1000)
pdf(file = "macro_markers_cluster_more_markers.pdf")
heatmap.2(heat_mean_mat, dendrogram = 'none', Colv = FALSE, Rowv = FALSE, col = col, srtCol = 0, offsetCol = -47, adjCol = c(NA, 1.5), scale = 'row', trace = 'none', density.info = 'none')
dev.off()

#Subset for slingshot analysis
aggr_all_seurat@meta.data$Louvain <- mapvalues(aggr_all_seurat@meta.data$Louvain, from = '4', to = '5')

pdf(file = "subset_macro_variable_genes.pdf")
aggr_all_subset <- FindVariableGenes(object = aggr_all_subset, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, y.high.cutoff = 20)
dev.off()
length(x = aggr_all_subset@var.genes)

aggr_all_subset <- RunPCA(object = aggr_all_subset, pcs.compute = 40, pc.genes = aggr_all_subset@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)

pdf(file = "subset_macro_pca_elbow.pdf")
PCElbowPlot(object = aggr_all_subset, num.pc = 40)
dev.off()

aggr_all_subset <- FindClusters(object = aggr_all_subset, reduction.type = "pca", dims.use = 1:25, 
                                resolution = 0.4, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

aggr_all_subset <- RunTSNE(object = aggr_all_subset, dims.use = 1:25, do.fast = TRUE)
pdf(file = "subset_macro_tsne.pdf")
TSNEPlot(object = aggr_all_subset)
dev.off()

aggr_all_subset <- RunUMAP(object = aggr_all_subset, reduction.use = "pca", dims.use = 1:25, n_neighbors = 10)

pdf(file = "subset_macro_umap_library.pdf")
DimPlot(object = aggr_all_subset, group.by = 'LibraryID', reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

#Convert AnnData to Seurat
library(reticulate)
ad <- import("anndata", convert = FALSE)
macro_seurat <- ad$read_h5ad('scanpy_anal_DE.h5ad')
macro_seurat_anal <- Convert(macro_seurat, to = "seurat")
macro_seurat_anal <- SetAllIdent(macro_seurat_anal, id = 'louvain')

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
pdf(file = 'macro_scat_difmap.pdf')
plotDiffusionMap(seurat_sce, colour_by = "Col1a1")
dev.off()


louvain <- read.delim('macro_louvain.csv', sep = ",", header = TRUE)
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

colors <- colorRampPalette(palette(c('#E3E3E3', '#A6D6A6', '#7CCD7C', '#508D50', '#2D5A2D', '#133413'))[-6])(12)
pdf(file = 'macro_slingshot_scat_color.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(colors[cut(slingshot_sce$slingPseudotime_1, breaks=12)], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
#lines(slingshot_sce@int_metadata$slingshot@curves$curve1, lwd=2, col = '#8B1A1A')
dev.off()

color <- palette(c('#A7414A', '#A37C27', '#A37C27', '#A37C27', '#A7414A', '#A7414A', '#6A8A82', '#6A8A82'))
pdf(file = 'macro_slingshot_scat_library.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(color[slingshot_sce$LibraryID], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

color1 <- palette(c('#33689E', '#B22222', '#B47746', '#50A4A8', '#4CBBC7'))
pdf(file = 'macro_slingshot_scat_ident_louvain.pdf')
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

#Subset for specific pop
Lrat_Macro <- SubsetData(aggr_all_seurat, ident.use = c('2', '3'))

pdf(file = "Macro_variable_genes_Macro.pdf")
Lrat_Macro <- FindVariableGenes(object = Lrat_Macro, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()
length(x = Lrat_Macro@var.genes)
write.table(as.character(Lrat_Macro@var.genes), file = 'HVG_Macros.txt', sep='\n', row.names=FALSE, col.names=FALSE)

#Subset for 3 major pops
CCl4_pops <- SubsetData(aggr_all_seurat, ident.use = c('0', '1'))

#Cell type annotation and output
CCl4_pops <- FindClusters(object = CCl4_pops, reduction.type = "pca", dims.use = 1:12, 
                                resolution = 0.02, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

CCl4_pops <- RunTSNE(object = CCl4_pops, dims.use = 1:12, do.fast = TRUE)
pdf(file = "aggr_dp_clean_tsne.pdf", width = 4000, height = 3000, res = 400)
TSNEPlot(object = CCl4_pops)
dev.off()

pdf(file = "aggr_dp_clean_markers.pdf", width = 8000, height = 5000, res = 400)
FeaturePlot(object = CCl4_pops, features.plot = c('Lrat', 'Ptprb', 'Adgre1'), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

current.cluster.ids <- c(0, 1, 2)
new.cluster.ids <- c('Endothelial cells', 'Hepatic stellate cells', 'Macrophages')
CCl4_pops@ident <- plyr::mapvalues(x = CCl4_pops@ident, from = current.cluster.ids, to = new.cluster.ids)

out <- as.data.frame(row.names(CCl4_pops@meta.data))
names(out)[1] <- c("Barcode")
out['nGene'] <- CCl4_pops@meta.data$nGene
out['nUMI'] <- CCl4_pops@meta.data$nUMI
out['LibraryID'] <- CCl4_pops@meta.data$LibraryID
out2 <- as.data.frame(CCl4_pops@ident)
out['CellType'] <- out2$'CCl4_pops@ident'
write.table(out, file = 'aggr_all_3_pops.txt', sep='\t', col.names=NA)
write.table(out$Barcode, file = 'barcodes_3_pops.tsv', sep='\t', col.names=FALSE, row.names=FALSE)
write.table(rownames(aggr_all_matrix), file = 'genes_3_pops.tsv', sep='\t', col.names=FALSE, row.names=FALSE)

aggr_all_matrix <- as.data.frame(aggr_all_matrix)
new <- aggr_all_matrix %>% select(one_of(dput(as.character(out$Barcode))))
new <- as.matrix(new)
rownames(new) <- c()
colnames(new) <- NULL
new1 <- as(new, "sparseMatrix")
writeMM(new1, 'matrix_3_pops.mtx')

outMacro <- out[out$CellType == 'Hepatic stellate cells',]
outEndo <- out[out$CellType == 'Endothelial cells',]
outMacro <- out[out$CellType == 'Macrophages',]

aggr_all_matrix <- as.data.frame(aggr_all_matrix)
newMacro <- aggr_all_matrix %>% select(one_of(dput(as.character(out$Barcode))))
newMacro <- as.matrix(newMacro)
rownames(newMacro) <- c()
colnames(newMacro) <- NULL
newMacro1 <- as(newMacro, "sparseMatrix")
writeMM(newMacro1, 'matrix_Macro_edit.mtx')

write.table(out, file = 'aggr_all_Macro_edit.txt', sep='\t', col.names=NA)
write.table(out$Barcode, file = 'barcodes_Macro_edit.tsv', sep='\t', col.names=FALSE, row.names=FALSE)
write.table(rownames(aggr_all_matrix), file = 'genes_Macro_edit.tsv', sep='\t', col.names=FALSE, row.names=FALSE)

#DoubletFinder
aggr_all_seurat <- doubletFinder(aggr_all_seurat, expected.doublets = 5200)

pdf(file = "Macro_doublets_umap.pdf", width = 4000, height = 3000, res = 400)
DimPlot(object = aggr_all_seurat, group.by = 'pANNPredictions', reduction.use = "umap", no.legend = FALSE, do.return = TRUE, 
        cols.use = c("tomato", "grey85"))
dev.off()

pdf(file = "Macro_doublets_tsne.pdf", width = 4000, height = 3000, res = 400)
DimPlot(object = aggr_all_seurat, group.by = 'pANNPredictions', reduction.use = "tsne", no.legend = FALSE, do.return = TRUE, 
        cols.use = c("tomato", "grey85"))
dev.off()

write.table(aggr_all_seurat@meta.data, file = 'aggr_all_doublet_filtering.txt', sep='\t', col.names=NA)

#Filter data to remove doublets
out <- as.data.frame(row.names(aggr_all_seurat@meta.data))
names(out)[1] <- c("Barcode")
out['DP'] <- aggr_all_seurat@meta.data$pANNPredictions
out['LibraryID'] <- aggr_all_seurat@meta.data$LibraryID
out['nGene'] <- aggr_all_seurat@meta.data$nGene
out['nUMI'] <- aggr_all_seurat@meta.data$nUMI
out2 <- out[out$DP == 'Singlet',]
write.table(out2, file = 'aggr_all_dp_removal.txt', sep='\t', header=TRUE, col.names=NA)
read.delim('aggr_all_dp_removal.txt', sep = '\t', header = TRUE)
new <- aggr_all_matrix %>% select(one_of(dput(as.character(filter$X))))
new1 <- as(new, "sparseMatrix")
writeMM(new1, 'matrix_dp.mtx')

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

Macro_markers_list2 <- as.data.frame(Macro_markers2)
Macro_markers_list2[6] <- rownames(Macro_markers_list2)
names(Macro_markers_list2)[6] <- c("symbol")
Macro_markers_list2 <- merge(Macro_markers_list2, grcm38, by="symbol")
View(Macro_markers_list2)
Macro_markers_list2 <- Macro_markers_list2[order(-Macro_markers_list2$pct.1),]
write.table(Macro_markers_list2, file = "Macro_markers_list2_pct1.txt", sep = "\t")

Macro_markers_list2[15] <- Macro_markers_list2$pct.1/Macro_markers_list2$pct.2
names(Macro_markers_list2)[15] <- c("pct1/pct2")
Macro_markers_list2 <- Macro_markers_list2[order(-Macro_markers_list2$`pct1/pct2`),]
write.table(Macro_markers_list2, file = "Macro_markers_list2.txt", sep = "\t")

Macro_markers_list2_genes <- scan("genelist2.txt", what="", sep="\n")

Macro_markers_list2_genes <- split(Macro_markers_list2_genes, ceiling(seq_along(Macro_markers_list2_genes)/15))

pdf(file = "Macro_genes.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = Macro_markers_list2_genes$`1`, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE, nCol = 7)
dev.off()

FeaturePlot(object = Lrat_tSNE, features.plot = c("Colec10", "Gas7", "Fgfr2", "Cdh2", "Ifitm1", "Tln2", "Ank3", "Tmem56", "Wt1", "Ctsk", "Htra1", "Dpt", "Aebp1", "Ltbp2"), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE, nCol = 7)

