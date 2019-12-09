library('Seurat')
library('plyr')
library('dplyr')

matrix <- as.matrix(read.delim('matrix_HSC.csv', sep = ','))
Genes <- read.delim('genes_HSC.tsv', sep = '\t', header = FALSE)
Barcodes <- read.delim('barcodes_order_HSC.tsv', sep = '\t', header = FALSE)

write.table(colnames(aggr_all_matrix2), file = 'barcodes_order_HSC.tsv', sep = '\t', col.names = FALSE, row.names = FALSE)
newHSC <- as.matrix(aggr_all_matrix2)
rownames(newHSC) <- c()
colnames(newHSC) <- NULL
newHSC1 <- as(newHSC, "sparseMatrix")
writeMM(newHSC1, 'matrix_order_HSC.mtx')

aggr_all_matrix2 <- as.data.frame(aggr_all_matrix)
aggr_all_matrix2 <- aggr_all_matrix2 %>% select(ends_with('-1'), ends_with('-5'), ends_with('-6'), ends_with('-2'), ends_with('-3'), ends_with('-4'), ends_with('-7'), ends_with('-8'), ends_with('-9'))



read_aggr_all_matrix <- readMM('matrix_order_HSC.mtx')
aggr_all_matrix <- as.matrix(read_aggr_all_matrix)
Genes <- read.delim('genes_HSC_norep1.tsv', sep = '\t', header = FALSE)
Barcodes <- read.delim('barcodes_HSC_norep1.tsv', sep = '\t', header = FALSE)
colnames(aggr_all_matrix) <- (Barcodes$V1)
rownames(aggr_all_matrix) <- make.names(Genes$V1, unique = TRUE)

aggr_all_seurat <- CreateSeuratObject(raw.data = aggr_all_matrix, min.cells = 50, project = 'Aggr all clean')

mito.genes <- grep(pattern = "^mt.", x = rownames(x = aggr_all_seurat@data), value = TRUE)
percent.mito <- Matrix::colSums(aggr_all_seurat@raw.data[mito.genes, ])/Matrix::colSums(aggr_all_seurat@raw.data)

aggr_all_seurat <- AddMetaData(object = aggr_all_seurat, metadata = percent.mito, col.name = "percent.mito")

Library_ID <- read.delim('HSC_norep1.csv', sep = ",", header = TRUE)
CellsMeta = aggr_all_seurat@meta.data
head(CellsMeta)
CellsMeta["LibraryID"] <- Library_ID$LibraryID
head(CellsMeta)
CellsMetaTrim <- subset(CellsMeta, select = c("LibraryID"))
head(CellsMetaTrim)
aggr_all_seurat <- AddMetaData(aggr_all_seurat, CellsMetaTrim)
head(aggr_all_seurat@meta.data)

Library_ID$nGene <- NULL
Library_ID$nUMI <- NULL
rownames(Library_ID) <- Library_ID$Barcode
Library_ID$Barcode <- NULL
Library_ID$CellType <- NULL
write.table(Library_ID, file = 'Library_ID.csv', sep = ',', col.names = NA)

CellsMeta2 <- merge(CellsMeta, Library_ID, by=0, sort = FALSE)
rownames(CellsMeta2) <- CellsMeta2$Row.names
CellsMeta2$Row.names <- NULL
aggr_all_seurat <- AddMetaData(aggr_all_seurat, CellsMeta2)
head(aggr_all_seurat@meta.data)


pdf(file = "HSC_nGene_nUMI_percent_mito.pdf")
VlnPlot(object = aggr_all_seurat, features.plot = c('nGene', 'nUMI', 'percent.mito'), nCol = 3)
dev.off()

pdf(file = "HSC_filter1.pdf")
GenePlot(object = aggr_all_seurat, gene1 = 'nUMI', gene2 = 'nGene')
dev.off()

pdf(file = "HSC_filter2.pdf")
GenePlot(object = aggr_all_seurat, gene1 = 'nUMI', gene2 = 'percent.mito')
dev.off()

aggr_all_seurat <- FilterCells(object = aggr_all_seurat, subset.names = c("nGene", "percent.mito"), 
                               low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

aggr_all_seurat <- NormalizeData(object = aggr_all_seurat, normalization.method = "LogNormalize", 
                                 scale.factor = 10000)

pdf(file = "HSC_variable_genes.pdf")
aggr_all_seurat <- FindVariableGenes(object = aggr_all_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, y.high.cutoff = 20)
dev.off()
length(x = aggr_all_seurat@var.genes)

aggr_all_seurat <- ScaleData(object = aggr_all_seurat, vars.to.regress = c("nUMI", 'percent.mito'))

aggr_all_seurat <- RunPCA(object = aggr_all_seurat, pcs.compute = 40, pc.genes = aggr_all_seurat@var.genes, do.print = TRUE, pcs.print = 1:5, 
                          genes.print = 5)

pdf(file = "HSC_pca_genes.pdf")
VizPCA(object = aggr_all_seurat, pcs.use = 1:2)
dev.off()

pdf(file = "HSC_pca.pdf")
PCAPlot(object = aggr_all_seurat, dim.1 = 1, dim.2 = 2, no.legend = TRUE)
dev.off()

pdf(file = "HSC_pca_heat.pdf")
PCHeatmap(object = aggr_all_seurat, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
dev.off()

pdf(file = "HSC_pca_elbow.pdf")
PCElbowPlot(object = aggr_all_seurat, num.pc = 40)
dev.off()

aggr_all_seurat <- FindClusters(object = aggr_all_seurat, reduction.type = "pca", dims.use = 1:25, 
                                resolution = 0.3, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

aggr_all_seurat <- RunTSNE(object = aggr_all_seurat, dims.use = 1:25, do.fast = TRUE)
pdf(file = "HSC_tsne.pdf")
TSNEPlot(object = aggr_all_seurat)
dev.off()

pdf(file = "HSC_marker_pop.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = c('Lrat', 'Acta2', 'Col1a1'), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

all_markers <- FindAllMarkers(object = aggr_all_seurat, only.pos = TRUE, min.pct = 0.50, logfc.threshold = 0.8)

top10 <- all_markers %>% group_by(cluster) %>% top_n(15, avg_logFC)

pdf(file = "marker_hsc_heatmap.pdf")
DoHeatmap(object = aggr_all_seurat, genes.use = top10$gene, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top', remove.key = TRUE)
dev.off()

pdf(file = "marker_hsc_heatmap_color.pdf")
DoHeatmap(object = aggr_all_seurat, genes.use = top10$gene, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top')
dev.off()

tiff(file = "marker_hsc_heatmap.tiff", width = 4000, height = 3000, res = 400)
DoHeatmap(object = aggr_all_seurat, genes.use = top10$gene, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top', remove.key = TRUE, cex.row = 15)
dev.off()

tiff(file = "marker_hsc_heatmap_color.tiff", width = 4000, height = 3000, res = 400)
DoHeatmap(object = aggr_all_seurat, genes.use = top10$gene, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top')
dev.off()

cluster_unk <- FindMarkers(object = aggr_all_seurat, ident.1 = 6, ident.2 = c(0, 1, 2, 3, 4, 5, 7, 8), min.pct = 0.25)
unk_print <- row.names(cluster_unk)[1:8]

pdf(file = "hsc_unk_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = unk_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_lrat <- FindMarkers(object = aggr_all_seurat, ident.1 = c(2), ident.2 = c(1), min.pct = 0.70, logfc.threshold = 1.5)
lrat_print <- row.names(cluster_lrat)[1:12]

pdf(file = "hsc_lrat_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = lrat_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_macro <- FindMarkers(object = aggr_all_seurat, ident.1 = c(1, 5), ident.2 = c(0, 2, 3, 4, 8), min.pct = 0.25)
macro_print <- row.names(cluster_macro)[1:12]

pdf(file = "hsc_macro_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = macro_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_endo <- FindMarkers(object = aggr_all_seurat, ident.1 = c(0, 4), ident.2 = c(1, 2, 3, 5, 8), min.pct = 0.25)
endo_print <- row.names(cluster_endo)[1:12]

pdf(file = "hsc_endo_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = endo_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

aggr_all_seurat <- RunUMAP(object = aggr_all_seurat, reduction.use = "pca", dims.use = 1:25, n_neighbors = 10)

pdf(file = "hsc_umap.pdf")
DimPlot(object = aggr_all_seurat, reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "hsc_marker_pop_umap.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = c('Lrat', 'Acta2', 'Col1a1'), cols.use = c("grey85", "red1"), reduction.use = "umap", no.legend = FALSE)
dev.off()

pdf(file = "hsc_umap_library.pdf")
DimPlot(object = aggr_all_seurat, group.by = 'LibraryID', reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

save.image('HSC_norep1.RData')

#Convert AnnData to Seurat
library(reticulate)
ad <- import("anndata", convert = FALSE)
HSC_seurat <- ad$read_h5ad('scanpy_anal_DE.h5ad')
HSC_seurat_anal <- Convert(HSC_seurat, to = "seurat")
HSC_seurat_anal <- SetAllIdent(HSC_seurat_anal, id = 'louvain')

#Convert to SCE
seurat_sce <- Convert(from = aggr_all_seurat, to = "sce")

library(reticulate)
ad <- import("anndata", convert = FALSE)
HSC_seurat <- ad$read_h5ad('scanpy_anal_15n.h5ad')
HSC_seurat_anal <- Convert(HSC_seurat, to = "seurat")
HSC_seurat_anal <- SetAllIdent(HSC_seurat_anal, id = 'louvain')

#Convert to SCE
seurat_sce <- Convert(from = HSC_seurat_anal, to = "sce")

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
plot(dm, 1:2, col_by ='LibraryID', legend_main ='Treatment', pal = c('#5ec02b', '#c0b82b', '#c0b82b', '#c0b82b', '#5ec02b', '#5ec02b', '#c02b2b', '#c02b2b'))
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
pdf(file = 'HSC_scat_difmap.pdf')
plotDiffusionMap(seurat_sce, colour_by = "Col1a1")
dev.off()


louvain <- read.delim('HSC_louvain.csv', sep = ",", header = TRUE)
rownames(louvain) <- louvain$index
louvain$index = NULL

CellsMeta = colData(seurat_sce)
head(CellsMeta)

CellsMeta2 <- merge(CellsMeta, louvain, by=0, sort = FALSE)
louvain <- as.data.frame(CellsMeta2$louvain)
names(louvain)[1] <- 'louvain'
seurat_sce2 <- mutate(seurat_sce, as.factor(louvain$louvain))


slingshot_sce <- slingshot(seurat_sce2, clusterLabels = 'as.factor.louvain.louvain.', reducedDim = 'DiffusionMap')
summary(slingshot_sce$slingPseudotime_1)

colors <- colorRampPalette(palette(c('#E3E3E3', '#8CB4D6', '#4F94CD', '#407AAA', '#275070', '#142E42'))[-6])(12)
colors <- viridis(12)
pdf(file = 'HSC_slingshot_scat_color.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(colors[cut(slingshot_sce$slingPseudotime_1, breaks=12)], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

#Legend out
legend_image <- as.raster(rev(color_gene), ncol=1)
pdf(file = 'legend_blue_1.pdf')
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

color_gene <- colorRampPalette(palette(c('#E3E3E3', '#99BDD9', '#69A1D1', '#4F94CC', '#5075A1', '#386396', '#20518B', '#183D68')))(12)
pdf(file = 'HSC_slingshot_scat_col1a1_new_scale.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(color_gene[cut(exma['Col1a1',], breaks=12)], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

pdf(file = 'HSC_slingshot_scat_fn1_new_scale.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(color_gene[cut(exma['Fn1',], breaks=12)], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

pdf(file = 'HSC_slingshot_scat_angptl6_new_scale.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(color_gene[cut(exma['Angptl6',], breaks=12)], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

pdf(file = 'HSC_slingshot_scat_mfap4_new_scale.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(color_gene[cut(exma['Mfap4',], breaks=12)], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

#Plotting in ggplot
slingshot_sce4 <- slingshot_sce3
slingshot_sce4$DC1 <- reducedDim(slingshot_sce4, 'DiffusionMap')[,1]
slingshot_sce4$DC2 <- reducedDim(slingshot_sce4, 'DiffusionMap')[,2]

pdf(file = 'HSC_slingshot_scat_ggplot_col1a1.pdf')
ggplot() + geom_point(data=as.data.frame(colData(slingshot_sce4)), mapping=aes(x = DC1, y = DC2, color = exma['Col1a1',]), size=.5) + scale_color_gradientn(colours = color_gene)
dev.off()

pdf(file = 'HSC_slingshot_scat_ggplot_fn1.pdf')
ggplot() + geom_point(data=as.data.frame(colData(slingshot_sce4)), mapping=aes(x = DC1, y = DC2, color = exma['Fn1',]), size=.5) + scale_color_gradientn(colours = color_gene)
dev.off()

pdf(file = 'HSC_slingshot_scat_ggplot_angptl6.pdf')
ggplot() + geom_point(data=as.data.frame(colData(slingshot_sce4)), mapping=aes(x = DC1, y = DC2, color = exma['Angptl6',]), size=.5) + scale_color_gradientn(colours = color_gene)
dev.off()

pdf(file = 'HSC_slingshot_scat_ggplot_mfap4.pdf')
ggplot() + geom_point(data=as.data.frame(colData(slingshot_sce4)), mapping=aes(x = DC1, y = DC2, color = exma['Mfap4',]), size=.5) + scale_color_gradientn(colours = color_gene)
dev.off()


color <- c('#BCE1FF', '#4F94CD', '#274763')
pdf(file = 'HSC_slingshot_scat_library.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(color[slingshot_sce3$short_lib.V1], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

pdf(file = 'HSC_slingshot_scat_col1a1.pdf')
plot(reducedDims(slingshot_sce)$DiffusionMap, col = alpha(colors[exma['Col1a1',]], 0.8), pch=16, asp = 1, cex=.5, xaxt = 'n', yaxt = 'n', ann=FALSE)
lines(SlingshotDataSet(slingshot_sce), lwd=2, col = '#2D2D2D')
dev.off()

color1 <- palette(c('#B22222', '#B47746', '#50A85E', '#33689E', '#50A4A8', '#674CC7', '#C64CC7', '#C74C83'))
pdf(file = 'HSC_slingshot_scat_ident_louvain_new_color.pdf')
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


#Plot genes in pseudotime with smoothened line
c('Plvap','Col1a1','Acta2','Mfap4','Abcc9','Rgs5')

exma_sub <- subset(exma_raw, rownames(exma_raw) %in% c('Plvap', 'Col1a1', 'Acta2', 'Mfap4', 'Abcc9', 'Rgs5'))
exma_sub <- as.data.frame(t(exma_sub))
exma_scale_sub <- subset(exma_scale, rownames(exma_scale) %in% c('Plvap', 'Col1a1', 'Acta2', 'Mfap4', 'Abcc9', 'Rgs5'))
exma_scale_sub <- as.data.frame(t(exma_scale_sub))
exma_scale_sub$DC1 <- as.data.frame(colData(slingshot_sce4))$DC1

pdf(file = 'HSC_slingshot_ggplot_col1a1_smooth.pdf')
ggplot(as.data.frame(colData(slingshot_sce4)), aes(x = DC1, y = exma_sub)) + geom_point(size=.5) + geom_smooth()
dev.off()

exma_sub <- as.data.frame(t(exma_sub))
exma_sub2 <- exma_sub
exma_sub2$DC1 = NULL
exma_sub2 <- apply(exma_sub2, MARGIN = 2, FUN = function(X) (X - min(X))/diff(range(X)))
exma_sub2 <- as.data.frame(exma_sub2)
exma_sub2$DC1 <- exma_sub$DC1

curve_agg <- ggplot(exma_sub2, aes(x = DC1, group = 1)) + geom_smooth(aes(y = Plvap, colour = 'Plvap'), se = FALSE) + geom_smooth(aes(y = Col1a1, colour = 'Col1a1'), se = FALSE) +
              geom_smooth(aes(y = Acta2, colour = 'Acta2'), se = FALSE) + geom_smooth(aes(y = Mfap4, colour = 'Mfap4'), se = FALSE) + geom_smooth(aes(y = Abcc9, colour = 'Abcc9'), se = FALSE) + 
              geom_smooth(aes(y = Rgs5, colour = 'Rgs5'), se = FALSE) + scale_colour_manual('', breaks = c('Plvap', 'Col1a1', 'Acta2', 'Mfap4', 'Abcc9', 'Rgs5'), values = c('#951818', '#8F9518', '#18951A', '#18958C', '#183A95', '#95185A'))
ggsave('HSC_slingshot_ggplot_smooth.pdf', curve_agg, width = 10, height = 6)

#Density plot from pseudotime
dens <- as.data.frame(reducedDims(slingshot_sce)$DiffusionMap)
dens$DC2 <- NULL
dens$DC3 <- NULL
dens$DC4 <- NULL
lib_dens <- as.data.frame(slingshot_sce$LibraryID)
names(lib_dens)[1] <- 'treatment'
dens$treatment <- lib_dens$treatment

pseu <- as.data.frame(slingshot_sce$slingPseudotime_1)
names(pseu)[1] <- 'pseudo'
dens$pseudo <- pseu$pseudo

pdf(file = 'lib_dens.pdf')
ggplot(dens, aes(x=DC1, fill = treatment_2)) + geom_density(alpha = 0.8) + scale_shape_manual(values = c(16, 16, 16)) + scale_fill_manual(breaks=c('Ctrl', 'CCl4 2w', 'CCl4 4w'), values=c('#4F94CD', '#274763', '#A6D6FF')) + labs(fill = 'Treatment')
dev.off()

pdf(file = 'lib_dens_pseu.pdf')
ggplot(dens, aes(x=pseudo, fill = treatment_2)) + geom_density(alpha = 0.8) + scale_shape_manual(values = c(16, 16, 16)) + scale_fill_manual(breaks=c('Ctrl', 'CCl4 2w', 'CCl4 4w'), values=c('#4F94CD', '#274763', '#A6D6FF')) + labs(fill = 'Treatment')
dev.off()

dens$treatment_3 <- factor(dens$treatment_2, levels = c('Ctrl', 'CCl4 2w', 'CCl4 4w'))
g <- ggplot(dens, aes(x = pseudo, y = treatment_3, fill = treatment_3)) + geom_density_ridges(scale = 4, alpha=0.7) + scale_fill_manual(breaks=c('Ctrl', 'CCl4 2w', 'CCl4 4w'), values=c('#A6D6FF', '#4F94CD', '#274763')) + labs(fill = 'Treatment') + theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0))
ggsave('lib_dens_pseu_ridge_ctrl_up.pdf', g, width = 9, height = 7)

dens$treatment_3 <- factor(dens$treatment_2, levels = c('CCl4 4w', 'CCl4 2w', 'Ctrl'))
h <- ggplot(dens, aes(x = pseudo, y = treatment_3, fill = treatment_3)) + geom_density_ridges(scale = 4, alpha=0.7) + scale_fill_manual(breaks=c('Ctrl', 'CCl4 2w', 'CCl4 4w'), values=c('#274763', '#4F94CD', '#A6D6FF')) + labs(fill = 'Treatment') + theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0))
ggsave('lib_dens_pseu_ridge_ctrl_down.pdf', h, width = 9, height = 7)

h <- ggplot(dens_stage, aes(x = pseudo, y = Stage, fill = treatment_3)) + geom_density_ridges(scale = 4, alpha=0.7) + scale_fill_manual(breaks=c('Ctrl', 'CCl4 2w', 'CCl4 4w'), values=c('#274763', '#4F94CD', '#A6D6FF')) + labs(fill = 'Treatment') + theme_ridges(grid = FALSE, center_axis_labels = TRUE) +
  scale_x_continuous(expand = c(0.01, 0)) +
  scale_y_discrete(expand = c(0.01, 0))
ggsave('lib_dens_pseu_ridge_stage.pdf', h, width = 9, height = 7)


#Subset data based on density peaks
dens_ctrl <- dens[dens$treatment_2 == 'Ctrl',]
dens_2w <- dens[dens$treatment_2 == 'CCl4 2w',]
dens_4w <- dens[dens$treatment_2 == 'CCl4 4w',]

which.max(density(dens_ctrl$pseudo)$y)
density(dens_ctrl$pseudo)$x[194]
#max for ctrl = 0.03526632

second_peak <- max(density(dens_ctrl$pseudo)$y[density(dens_ctrl$pseudo)$x > 0.05])
which(density(dens_ctrl$pseudo)$y == second_peak)
density(dens_ctrl$pseudo)$x[346]
#second max for ctrl = 0.06841782

which.max(density(dens_2w$pseudo)$y)
density(dens_2w$pseudo)$x[211]
#max for 2w = 0.05301848

second_peak <- max(density(dens_2w$pseudo)$y[density(dens_2w$pseudo)$x > 0.06])
which(density(dens_2w$pseudo)$y == second_peak)
density(dens_2w$pseudo)$x[284]
#second max for 2w = 0.06895075

which.max(density(dens_4w$pseudo)$y)
density(dens_4w$pseudo)$x[318]
#max for 4w = 0.09145961

second_peak <- max(density(dens_4w$pseudo)$y[density(dens_4w$pseudo)$x < 0.08])
which(density(dens_4w$pseudo)$y == second_peak)
density(dens_4w$pseudo)$x[218]
#second max for 4w = 0.06794788

#Subset list of barcodes for analysis app 500 cells
ini <- dens[dens$pseudo > 0.03276632 & dens$pseudo < 0.03776632,]
#(0.06841782 + 0.06895075 + 0.06794788)/3 = 0.06843882
mid_2w <- dens[dens$pseudo > 0.05250848 & dens$pseudo < 0.05352848,]
mid <- dens[dens$pseudo > 0.06895075 & dens$pseudo < 0.06895075,]
late <- dens[dens$pseudo > 0.09045961 & dens$pseudo < 0.09245961,]

#Subset list for stage assignment
ini <- dens[dens$pseudo > 0.03026632 & dens$pseudo < 0.04026632,]
ini$Stage_1 <- 1
ini$DC1 <- ini$treatment <- ini$treatment_2 <- ini$pseudo <- ini$treatment_3 <- NULL
#(0.06841782 + 0.06895075 + 0.06794788)/3 = 0.06843882
mid_2w <- dens[dens$pseudo > 0.04801848 & dens$pseudo < 0.05801848,]
mid_2w$Stage_2 <- 2
mid_2w$DC1 <- mid_2w$treatment <- mid_2w$treatment_2 <- mid_2w$pseudo <- mid_2w$treatment_3 <- NULL
mid <- dens[dens$pseudo > 0.06395075 & dens$pseudo < 0.07395075,]
mid$Stage_3 <- 3
mid$DC1 <- mid$treatment <- mid$treatment_2 <- mid$pseudo <- mid$treatment_3 <- NULL
late <- dens[dens$pseudo > 0.08645961 & dens$pseudo < 0.09645961,]
late$Stage_4 <- 4
late$DC1 <- late$treatment <- late$treatment_2 <- late$pseudo <- late$treatment_3 <- NULL

dens2 <- dens
dens2 <- merge(dens2, ini, by=0, all.x = TRUE, sort = FALSE)
rownames(dens2) <- dens2$Row.names
dens2$Row.names <- NULL
dens2 <- merge(dens2, mid_2w, by=0, all.x = TRUE, sort = FALSE)
rownames(dens2) <- dens2$Row.names
dens2$Row.names <- NULL
dens2 <- merge(dens2, mid, by=0, all.x = TRUE, sort = FALSE)
rownames(dens2) <- dens2$Row.names
dens2$Row.names <- NULL
dens2 <- merge(dens2, late, by=0, all.x = TRUE, sort = FALSE)
rownames(dens2) <- dens2$Row.names
dens2$Row.names <- NULL

dens3 <- unite(dens2, 'Stage', c('Stage_1', 'Stage_2', 'Stage_3', 'Stage_4'))

dens2 <- dens2[with(dens2, order(pseudo)), ]

library(gam)
library(clusterExperiment)
t <- slingshot_sce$slingPseudotime_1

# for time, only look at the 1,000 most variable genes
Y <- assays(seurat_sce)$logcounts
var1K <- names(sort(apply(Y,1,var),decreasing = TRUE))[1:1000]
Y <- Y[var1K,]

# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

topgenes <- names(sort(gam.pval, decreasing = FALSE))[1:200]
heatdata <- assays(seurat_sce)$norm[rownames(assays(seurat_sce)$norm) %in% topgenes, 
                             order(t, na.last = NA)]
heatdata_scale <- assays(seurat_sce)$scale[rownames(assays(seurat_sce)$scale) %in% topgenes, 
                                           order(t, na.last = NA)]

heatclus <- slingshot_sce$as.factor.louvain.louvain.[order(t, na.last = NA)]
heatlib <- slingshot_sce$LibraryID[order(t, na.last = NA)]
heatlib <- mapvalues(heatlib, c('Ctrl - 1', 'Ctrl - 2', 'Ctrl - 3',
                                'CCl4 2w - 1', 'CCl4 2w - 2', 'CCl4 2w - 3',
                                'CCl4 4w - 2', 'CCl4 4w - 3'), 
                              c('Ctrl', 'Ctrl', 'Ctrl',
                                'CCl4 2w', 'CCl4 2w', 'CCl4 2w',
                                'CCl4 4w', 'CCl4 4w'))

ce <- ClusterExperiment(heatdata, heatclus, transformation = log1p)
ce_scale <- ClusterExperiment(heatdata_scale, heatclus)


library(colorRamps)
hmcol <- rev(brewer.pal(11,"RdBu"))
hmcol2 <- matlab.like(20)

new_names <- c('0', '1', '2', '3')
color_clus <- c('#B22222', '#B47746', '#50A85E', '#33689E')
color_lib <- c('#4F94CD', '#274763', '#A6D6FF')
names(color_lib)<-names(new_names)<-clusterLegend(ce)[[1]][,'name']
ce <- recolorClusters(ce, whichCluster=1, value=color_clus)
ce_scale <- recolorClusters(ce_scale, whichCluster=1, value=color_clus)
ce <- recolorClusters(ce, whichCluster=1, value=color_lib)
ce_scale <- recolorClusters(ce_scale, whichCluster=1, value=color_lib)

ce_scale <- mutate(ce_scale, heatlib)

col_clust = c('#50A85E', '#B22222', '#33689E', '#B47746')
col_lib = c('#4F94CD', '#274763', '#A6D6FF')
ann_colors = list(heatlib = c('clusterIds' = 1, 2, 3, 4, 'color' = col_lib, 'name' = c('Ctrl', 'CCl4 2w', 'CCl4 4w')))

ann_colors = matrix(c('1', '2', '3', '#A6D6FF', '#4F94CD', '#274763', 'Ctrl', 'CCl4 2w', 'CCl4 4w'), nrow=3, ncol=3) 
colnames(ann_colors) <- c('clusterIds', 'color', 'name')
ann_colors = list(heatlib = ann_colors)

pdf(file = 'pseudo_heatmap_200.pdf')
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed', colorScale = hmcol2, breaks = 0.99, capBreaksLegend = TRUE)
dev.off()

png(file = "pseudo_heatmap_200.png", width = 8000, height = 10000, res = 1200)
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue", visualizeData = 'transformed', colorScale = hmcol2, breaks = 0.99, capBreaksLegend = TRUE)
dev.off()

pdf(file = 'pseudo_heatmap_scale_200_heatlib.pdf')
plotHeatmap(ce_scale, clusterSamplesData = "orderSamplesValue", colData = 'heatlib', clusterLegend = ann_colors, colorScale = hmcol, breaks = 0.99, capBreaksLegend = TRUE, symmetricBreaks = TRUE)
dev.off()

png(file = "pseudo_heatmap_scale_200_heatlib.png", width = 8000, height = 10000, res = 1200)
plotHeatmap(ce_scale, clusterSamplesData = "orderSamplesValue", colData = 'heatlib', clusterLegend = ann_colors, colorScale = hmcol, breaks = 0.99, capBreaksLegend = TRUE, symmetricBreaks = TRUE)
dev.off()


#Subset for specific pop
Lrat_HSC <- SubsetData(aggr_all_seurat, ident.use = c('2', '3'))

pdf(file = "hsc_variable_genes_HSC.pdf")
Lrat_HSC <- FindVariableGenes(object = Lrat_HSC, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
dev.off()
length(x = Lrat_HSC@var.genes)
write.table(as.character(Lrat_HSC@var.genes), file = 'HVG_HSCs.txt', sep='\n', row.names=FALSE, col.names=FALSE)

#Subset for 3 major pops
CCl4_pops <- SubsetData(aggr_all_seurat, ident.use = c('0', '1', '2', '3', '4', '5'))

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

outHSC <- out[out$CellType == 'Hepatic stellate cells',]
outEndo <- out[out$CellType == 'Endothelial cells',]
outMacro <- out[out$CellType == 'Macrophages',]

aggr_all_matrix <- as.data.frame(aggr_all_matrix)
newHSC <- aggr_all_matrix %>% select(one_of(dput(as.character(outHSC$Barcode))))
newHSC <- as.matrix(newHSC)
rownames(newHSC) <- c()
colnames(newHSC) <- NULL
newHSC1 <- as(newHSC, "sparseMatrix")
writeMM(newHSC1, 'matrix_HSC.mtx')

write.table(outHSC, file = 'aggr_all_HSC.txt', sep='\t', col.names=NA)
write.table(outHSC$Barcode, file = 'barcodes_HSC.tsv', sep='\t', col.names=FALSE, row.names=FALSE)
write.table(rownames(aggr_all_matrix), file = 'genes_HSC.tsv', sep='\t', col.names=FALSE, row.names=FALSE)

#DoubletFinder
aggr_all_seurat <- doubletFinder(aggr_all_seurat, expected.doublets = 5200)

pdf(file = "hsc_doublets_umap.pdf", width = 4000, height = 3000, res = 400)
DimPlot(object = aggr_all_seurat, group.by = 'pANNPredictions', reduction.use = "umap", no.legend = FALSE, do.return = TRUE, 
        cols.use = c("tomato", "grey85"))
dev.off()

pdf(file = "hsc_doublets_tsne.pdf", width = 4000, height = 3000, res = 400)
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

HSC_markers_list2 <- as.data.frame(HSC_markers2)
HSC_markers_list2[6] <- rownames(HSC_markers_list2)
names(HSC_markers_list2)[6] <- c("symbol")
HSC_markers_list2 <- merge(HSC_markers_list2, grcm38, by="symbol")
View(HSC_markers_list2)
HSC_markers_list2 <- HSC_markers_list2[order(-HSC_markers_list2$pct.1),]
write.table(HSC_markers_list2, file = "HSC_markers_list2_pct1.txt", sep = "\t")

HSC_markers_list2[15] <- HSC_markers_list2$pct.1/HSC_markers_list2$pct.2
names(HSC_markers_list2)[15] <- c("pct1/pct2")
HSC_markers_list2 <- HSC_markers_list2[order(-HSC_markers_list2$`pct1/pct2`),]
write.table(HSC_markers_list2, file = "HSC_markers_list2.txt", sep = "\t")

HSC_markers_list2_genes <- scan("genelist2.txt", what="", sep="\n")

HSC_markers_list2_genes <- split(HSC_markers_list2_genes, ceiling(seq_along(HSC_markers_list2_genes)/15))

pdf(file = "hsc_genes.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = HSC_markers_list2_genes$`1`, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE, nCol = 7)
dev.off()

FeaturePlot(object = Lrat_tSNE, features.plot = c("Colec10", "Gas7", "Fgfr2", "Cdh2", "Ifitm1", "Tln2", "Ank3", "Tmem56", "Wt1", "Ctsk", "Htra1", "Dpt", "Aebp1", "Ltbp2"), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE, nCol = 7)

#Convert pilot to seurat
seurat_sce <- Convert(from = aggr_all_seurat, to = "sce")

library(Seurat)
library(reticulate)
ad <- import("anndata", convert = FALSE)
pilot_seurat <- ad$read_h5ad('C:/Users/mikekrogh/Desktop/single_cell_ccl4/processed_data/pilot_exp/scanpy_anal.h5ad')
pilot_seurat <- Convert(pilot_seurat, to = "seurat")
pilot_seurat <- SetAllIdent(pilot_seurat, id = 'leiden')


