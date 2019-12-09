library('Seurat')
library('plyr')
library('dplyr')


read_aggr_all_matrix <- readMM('matrix_norep1.mtx')
aggr_all_matrix <- as.matrix(read_aggr_all_matrix)
Genes <- read.delim('genes_norep1.tsv', sep = '\t', header = FALSE)
Barcodes <- read.delim('barcodes_norep1.tsv', sep = '\t', header = FALSE)
colnames(aggr_all_matrix) <- (Barcodes$V1)
rownames(aggr_all_matrix) <- make.names(Genes$V1, unique = TRUE)

aggr_all_matrix2 <- as.data.frame(aggr_all_matrix)
aggr_all_matrix2 <- aggr_all_matrix2 %>% select(ends_with('-1'), ends_with('-5'), ends_with('-6'), ends_with('-2'), ends_with('-3'), ends_with('-4'), ends_with('-7'), ends_with('-8'), ends_with('-9'))

aggr_all_seurat <- CreateSeuratObject(raw.data = aggr_all_matrix2, min.cells = 50, project = 'Aggr all clean')

mito.genes <- grep(pattern = "^mt.", x = rownames(x = aggr_all_seurat@data), value = TRUE)
percent.mito <- Matrix::colSums(aggr_all_seurat@raw.data[mito.genes, ])/Matrix::colSums(aggr_all_seurat@raw.data)

aggr_all_seurat <- AddMetaData(object = aggr_all_seurat, metadata = percent.mito, col.name = "percent.mito")

Library_ID <- read.delim('norep1.csv', sep = ",", header = TRUE)
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

pdf(file = "aggr_all_nGene_nUMI_percent_mito.pdf")
VlnPlot(object = aggr_all_seurat, features.plot = c('nGene', 'nUMI', 'percent.mito'), nCol = 3)
dev.off()

pdf(file = "aggr_all_filter1.pdf")
GenePlot(object = aggr_all_seurat, gene1 = 'nUMI', gene2 = 'nGene')
dev.off()

pdf(file = "aggr_all_filter2.pdf")
GenePlot(object = aggr_all_seurat, gene1 = 'nUMI', gene2 = 'percent.mito')
dev.off()

aggr_all_seurat <- FilterCells(object = aggr_all_seurat, subset.names = c("nGene", "percent.mito"), 
                               low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

aggr_all_seurat <- NormalizeData(object = aggr_all_seurat, normalization.method = "LogNormalize", 
                                 scale.factor = 10000)

pdf(file = "aggr_all_variable_genes.pdf")
aggr_all_seurat <- FindVariableGenes(object = aggr_all_seurat, mean.function = ExpMean, dispersion.function = LogVMR, 
                                     x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, y.high.cutoff = 20)
dev.off()
length(x = aggr_all_seurat@var.genes)

aggr_all_seurat <- ScaleData(object = aggr_all_seurat, vars.to.regress = c("nUMI", 'percent.mito'))

aggr_all_seurat <- RunPCA(object = aggr_all_seurat, pcs.compute = 40, pc.genes = aggr_all_seurat@var.genes, do.print = TRUE, pcs.print = 1:5, 
                        genes.print = 5)

pdf(file = "aggr_all_pca_genes.pdf")
VizPCA(object = aggr_all_seurat, pcs.use = 1:2)
dev.off()

pdf(file = "aggr_all_pca.pdf")
PCAPlot(object = aggr_all_seurat, dim.1 = 1, dim.2 = 2, no.legend = TRUE)
dev.off()

pdf(file = "aggr_all_pca_heat.pdf")
PCHeatmap(object = aggr_all_seurat, pc.use = 1:20, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)
dev.off()

pdf(file = "aggr_all_pca_elbow.pdf")
PCElbowPlot(object = aggr_all_seurat, num.pc = 40)
dev.off()

aggr_all_seurat <- FindClusters(object = aggr_all_seurat, reduction.type = "pca", dims.use = 1:35, 
                              resolution = 0.02, print.output = 0, save.SNN = TRUE, force.recalc = TRUE)

aggr_all_seurat <- RunTSNE(object = aggr_all_seurat, dims.use = 1:35, do.fast = TRUE)
pdf(file = "aggr_all_tsne.pdf")
TSNEPlot(object = aggr_all_seurat)
dev.off()

pdf(file = "aggr_all_marker_pop.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = c('Lrat', 'Ptprb', 'Adgre1'), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

all_markers <- FindAllMarkers(object = aggr_all_seurat, only.pos = TRUE, min.pct = 0.60, logfc.threshold = 2)

top10 <- all_markers %>% group_by(cluster) %>% top_n(15, avg_logFC)

heatmap_markers <- read.delim('heatmap_markers.csv', sep = ",", header = FALSE)

pdf(file = "marker_heatmap_ord.pdf")
DoHeatmap(object = aggr_all_seurat2, cells.use = heat_nam2, genes.use = heatmap_markers$V1, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top', remove.key = TRUE)
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


all_markers_low <- FindAllMarkers(object = aggr_all_seurat, only.pos = TRUE, min.pct = 0.60, logfc.threshold = 1.6)

top20 <- all_markers_low %>% group_by(cluster) %>% top_n(30, avg_logFC)

pdf(file = "marker_heatmap_low.pdf")
DoHeatmap(object = aggr_all_seurat, genes.use = top20$gene, slim.col.label = TRUE, col.high = '#ff0000', col.mid = '#ffffff', col.low = '#0066cc', group.label.loc = 'top', remove.key = TRUE, cex.row = 6)
dev.off()


cluster_unk <- FindMarkers(object = aggr_all_seurat, ident.1 = 6, ident.2 = c(0, 1, 2, 3, 4, 5, 7, 8), min.pct = 0.25)
unk_print <- row.names(cluster_unk)[1:8]

pdf(file = "aggr_all_unk_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = unk_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_lrat <- FindMarkers(object = aggr_all_seurat, ident.1 = c(2, 3), ident.2 = c(0, 1, 4, 5, 8), min.pct = 0.25)
lrat_print <- row.names(cluster_lrat)[1:12]

pdf(file = "aggr_all_lrat_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = lrat_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_macro <- FindMarkers(object = aggr_all_seurat, ident.1 = c(1, 5), ident.2 = c(0, 2, 3, 4, 8), min.pct = 0.25)
macro_print <- row.names(cluster_macro)[1:12]

pdf(file = "aggr_all_macro_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = macro_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

cluster_endo <- FindMarkers(object = aggr_all_seurat, ident.1 = c(0, 4), ident.2 = c(1, 2, 3, 5, 8), min.pct = 0.25)
endo_print <- row.names(cluster_endo)[1:12]

pdf(file = "aggr_all_endo_print.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = endo_print, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

aggr_all_seurat <- RunUMAP(object = aggr_all_seurat, reduction.use = "pca", dims.use = 1:35, n_neighbors = 15)

pdf(file = "aggr_all_umap.pdf")
DimPlot(object = aggr_all_seurat, reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

pdf(file = "aggr_all_marker_pop_umap.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = c('Lrat', 'Ptprb', 'Adgre1', 'Alb'), cols.use = c("grey85", "red1"), reduction.use = "umap", no.legend = FALSE)
dev.off()

pdf(file = "aggr_all_umap_library.pdf")
DimPlot(object = aggr_all_seurat, group.by = 'LibraryID', reduction.use = "umap", no.legend = FALSE, do.return = TRUE) + ggtitle("UMAP") + theme(plot.title = element_text(hjust = 0.5))
dev.off()

save.image('aggr_all_new_dp.RData')

#Heatmap of clustered cells
load('aggr_norep1_seurat.RData')
heat_mat <- aggr_all_seurat@scale.data
heat_order <- read.delim('heat_barcodes_300.csv', sep = ',', header = TRUE)
rownames(heat_order) <- heat_order$X
heat_order$X = NULL
heat_order <- as.vector(colnames(heat_order))
heat_order <- gsub('[.]', '-', heat_order)
heat_mat <- as.data.frame(heat_mat)
heat_mat2 <- heat_mat %>% select(one_of(heat_order))
heat_mat2 <- tibble::rownames_to_column(heat_mat2)
heat_mat_clean <- filter(heat_mat2, heat_mat2$rowname %in% heatmap_markers$V1)
heat_mat_clean <- heat_mat_clean %>% slice(match(heatmap_markers$V1, heat_mat_clean$rowname))
rownames(heat_mat_clean) <- heat_mat_clean$rowname
heat_mat_clean$rowname = NULL
write.table(heat_mat_clean, file = 'heat_mat_300.csv', sep = ',', col.names = NA)

library(gplots)
library(RColorBrewer)
hmcol <- brewer.pal(11,"RdBu")
plot_heat <- as.matrix(heat_mat_clean)
png(file = "new_heat.png", width = 13000, height = 10000, res = 1200)
heatmap.2(plot_heat, col=hmcol, dendrogram='none', Rowv=FALSE, Colv=FALSE, trace='none')
dev.off()

#Subset for specific pop
Lrat_HSC <- SubsetData(aggr_all_seurat, ident.use = c('2', '3'))

pdf(file = "aggr_all_variable_genes_HSC.pdf")
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

pdf(file = "aggr_dp_clean_markers.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = c('Lrat', 'Ptprb', 'Adgre1'), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE)
dev.off()

current.cluster.ids <- c(0, 1, 2)
new.cluster.ids <- c('Endothelial cells', 'Hepatic stellate cells', 'Macrophages')
aggr_all_seurat@ident <- plyr::mapvalues(x = aggr_all_seurat@ident, from = current.cluster.ids, to = new.cluster.ids)

out <- as.data.frame(row.names(aggr_all_seurat@meta.data))
names(out)[1] <- c("Barcode")
out['nGene'] <- aggr_all_seurat@meta.data$nGene
out['nUMI'] <- aggr_all_seurat@meta.data$nUMI
out['LibraryID'] <- aggr_all_seurat@meta.data$LibraryID
out2 <- as.data.frame(aggr_all_seurat@ident)
out['CellType'] <- out2$'aggr_all_seurat@ident'
write.table(out, file = 'C.csv', sep=',', col.names=NA)
write.table(out$Barcode, file = 'barcodes_4w_reps.tsv', sep='\t', col.names=FALSE, row.names=FALSE)
write.table(rownames(aggr_all_matrix), file = 'genes_4w_reps.tsv', sep='\t', col.names=FALSE, row.names=FALSE)

out <- out[out$LibraryID == c('CCl4 4w - 1', 'CCl4 4w - 2', 'CCl4 4w - 3'),]
out <- out[out$Cluster %in% c('0', '1', '2', '3', '4', '5'),]

aggr_all_matrix <- as.data.frame(aggr_all_matrix)
new <- aggr_all_matrix %>% select(one_of(dput(as.character(out$Barcode))))
new <- as.matrix(new)
rownames(new) <- c()
colnames(new) <- NULL
new1 <- as(new, "sparseMatrix")
writeMM(new1, 'CCl4_4w_reps.mtx')

outHSC <- out[out$CellType == 'Hepatic stellate cells',]
outEndo <- out[out$CellType == 'Endothelial cells',]
outMacro <- out[out$CellType == 'Macrophages',]

aggr_all_matrix <- as.data.frame(aggr_all_matrix)
newall <- aggr_all_matrix %>% select(one_of(dput(as.character(out$Barcode))))
newall <- as.matrix(newall)
rownames(newall) <- c()
colnames(newall) <- NULL
newall1 <- as(newall, "sparseMatrix")
writeMM(newall1, 'matrix_all_clean_clean.mtx')

write.table(out, file = 'all_clean_clean.txt', sep='\t', col.names=NA)
write.table(out$Barcode, file = 'barcodes_all_clean_clean.tsv', sep='\t', col.names=FALSE, row.names=FALSE)
write.table(rownames(aggr_all_matrix), file = 'genes_all_clean_clean.tsv', sep='\t', col.names=FALSE, row.names=FALSE)

#CellPhoneDB
matrix <- as.matrix(read.delim('X.csv', sep = ',', header = FALSE))
genes <- read.delim('var.csv', sep = ',')
barcodes <- read.delim('obs.csv', sep = ',')
genes <- as.data.frame(gsub('X(.*Rik)', '\\1', genes$X0))
names(genes)[1] <- 'X0'
genes <- as.data.frame(gsub('[.]', '-', genes$X0))
names(genes)[1] <- 'X0'
genes <- read.delim('genes_var.csv', sep = ',', header = FALSE)


matrix2 <- t(matrix)
colnames(matrix2) <- barcodes$X0
rownames(matrix2) <- genes$V1
write.table(matrix2, file = 'cellphone_matrix.csv', sep=',')

gtf <- rtracklayer::import('genes.gtf')
gtf_df=as.data.frame(gtf)
gtf_df2 <- as.data.frame(gtf_df$gene_id)
names(gtf_df2)[1] <- 'gene_id'
gtf_df2['gene_name'] <- gtf_df$gene_name
gtf_df3 <- distinct(gtf_df2)
gtf_df4 <- gtf_df3[!duplicated(gtf_df3$gene_name), ]

row_matrix <- as.data.frame(rownames(matrix3))
names(row_matrix)[1] <- 'gene_name'
row_matrix <- as.data.frame(gsub('X(.*Rik)', '\\1', row_matrix$gene_name))
names(row_matrix)[1] <- 'gene_name'
row_matrix <- as.data.frame(gsub('[.]', '-', row_matrix$gene_name))
names(row_matrix)[1] <- 'gene_name'
row_matrix2 <- merge(row_matrix, gtf_df4, by.all="gene_name", sort=FALSE, all.x=T)
#row_matrix2 <- inner_join(row_matrix, gtf_df3, by='gene_name')

row_names_remove <- c("CR974586-5", "4930556M19Rik-1", "Gt-ROSA-26Sor", "AC168977-2", "Sept15", "Pcdhga8-1", "Smim20-1", "Gbp6-1", "Gm16701-1", "AC168977-1", "AC125149-3", "CAAA01147332-1", "Atn1-1", "Dancr-1", "Schip1-1")
matrix3 <- matrix2[!(row.names(matrix2) %in% row_names_remove), ]

rownames(matrix3) <- row_matrix2$gene_id
matrix3[1:10,1:10]

require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
genes_con <- row_matrix2$gene_id
genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = genes_con, mart = mouse, attributesL = c("ensembl_gene_id"), martL = human)
names(genesV2)[1] <- 'gene_id'
names(genesV2)[2] <- 'gene_id_human'

row_matrix3 <- merge(row_matrix2, genesV2, by.all="gene_id", sort=FALSE, all.x=T)
row_matrix3 <- left_join(row_matrix2, genesV3, by='gene_id')

genesV3 <- genesV2[!duplicated(genesV2$gene_id), ]
row_matrix3 <- left_join(row_matrix2, genesV3, by='gene_id')
rownames(matrix3) <- row_matrix3$gene_id_human
matrix4 <- matrix3[!is.na(rownames(matrix3)), ]

#DoubletFinder
aggr_all_seurat <- doubletFinder(aggr_all_seurat, expected.doublets = 5200)

pdf(file = "aggr_all_doublets_umap.pdf", width = 4000, height = 3000, res = 400)
DimPlot(object = aggr_all_seurat, group.by = 'pANNPredictions', reduction.use = "umap", no.legend = FALSE, do.return = TRUE, 
        cols.use = c("tomato", "grey85"))
dev.off()

pdf(file = "aggr_all_doublets_tsne.pdf", width = 4000, height = 3000, res = 400)
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

pdf(file = "aggr_all_genes.pdf")
FeaturePlot(object = aggr_all_seurat, features.plot = HSC_markers_list2_genes$`1`, cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE, nCol = 7)
dev.off()

FeaturePlot(object = Lrat_tSNE, features.plot = c("Colec10", "Gas7", "Fgfr2", "Cdh2", "Ifitm1", "Tln2", "Ank3", "Tmem56", "Wt1", "Ctsk", "Htra1", "Dpt", "Aebp1", "Ltbp2"), cols.use = c("grey85", "red1"), reduction.use = "tsne", no.legend = FALSE, nCol = 7)

#Subset count matrix based on expression of Hbb-bs
aggr_all_matrix2 <- aggr_all_matrix
row.names(aggr_all_matrix) <- gsub('-', '_', row.names(aggr_all_matrix))
aggr_hbbbs_subset <- aggr_all_matrix[,aggr_all_matrix['Hbb-bs',]<1,drop=FALSE]

aggr_all_seurat <- SubsetData(aggr_all_seurat, subset.name = 'Hbb.bs', accept.high = 1)

