#Workflow adapted from https://bioconductor.org/packages/devel/workflows/vignettes/maEndToEnd/inst/doc/MA-Workflow.html#2_workflow_package_installation

#Unix command to extract meta data from GEO soft file
#zgrep -P '^(\^SAMPLE = GSM|\!Sample_characteristics_ch1 = other_id:)' GSE48452_family.soft.gz | paste - - | sed -r 's/(\^SAMPLE = |\!Sample_characteristics_ch1 = other_id: )//g'
#zgrep -P '^(\^SAMPLE = GSM|\!Sample_title)' GSE48452_family.soft.gz | paste - - | sed -r 's/(\^SAMPLE = |\!Sample_title)//g'

#zgrep -P '^(\^SAMPLE = GSM|\!Sample_characteristics_ch1 = tissue:)' GSE49541_family.soft.gz | paste - - | sed -r 's/(\^SAMPLE = |\!Sample_characteristics_ch1 = tissue: )//g'
#zgrep -P '^(\^SAMPLE = GSM|\!Sample_title)' GSE49541_family.soft.gz | paste - - | sed -r 's/(\^SAMPLE = |\!Sample_title)//g'


library(oligo)
library(ggplot2)
library(arrayQualityMetrics)
library(pheatmap)
library(RColorBrewer)
library(stringr)
library(plyr)
library(dplyr)
library(limma)
library(pd.hg.u133.plus.2)
library(hgu133plus2.db)

celpath <- 'ccl4_array_anal/Murphy/cel_files/'
cellist = list.files(celpath, full.names = TRUE)
meta_data <- read.delim('ccl4_array_anal/Murphy/human_array_oligo_meta_new.csv', header = TRUE, sep = ',', row.names = 1)
meta_data <- AnnotatedDataFrame(meta_data)

data = read.celfiles(cellist, verbose = FALSE, phenoData = meta_data)

#Log2 transform and perform PCA
exp_raw <- log2(exprs(data))
PCA_raw <- prcomp(t(exp_raw), scale. = FALSE)

percentVar <- round(100*PCA_raw$sdev^2/sum(PCA_raw$sdev^2),1)
sd_ratio <- sqrt(percentVar[2] / percentVar[1])

dataGG <- data.frame(PC1 = PCA_raw$x[,1], PC2 = PCA_raw$x[,2],
                     tissue = pData(data)$Tissue,
                     stage = pData(data)$Stage,
                     Individual = rownames(pData(data)))

PCA_plot <- ggplot(dataGG, aes(PC1, PC2)) +
            geom_point(aes(shape = tissue, colour = stage)) + ggtitle('PCA of log-transformed expression data') +
            xlab(paste0('PC1, VarExp: ', percentVar[1], '%')) + ylab(paste0('PC2, VarExp: ', percentVar[2], '%')) +
            theme(plot.title = element_text(hjust = 0.5))+ coord_fixed(ratio = sd_ratio) + scale_shape_manual(values = c(19)) + 
            scale_color_manual(values = c('firebrick', 'goldenrod3', 'darkorange2', 'chartreuse3', 'darkolivegreen', 'dodgerblue4', 'cyan4', 'darkorchid4')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.key=element_blank())
ggsave('human_array_pca.pdf', PCA_plot, width = 10, height = 6)

#Boxplot to visualize intensity distribution (whether or not to normalize)
pdf('human_array_boxplot.pdf', width = 7, height = 4)
boxplot(data, target = 'core', main = 'Boxplot of log2-intensitites for the raw data', xaxt='n')
dev.off()

#Report of quality metrics
arrayQualityMetrics(expressionset = data, outdir = 'ccl4_array_anal/Murphy', force = TRUE, do.logtransform = TRUE, 
                    intgroup = c('Tissue', 'Stage'))

#Relative Log Expression
rle_data <- oligo::rma(data, normalize = FALSE)

row_medians_assayData <- Biobase::rowMedians(as.matrix(Biobase::exprs(rle_data)))

RLE_data <- sweep(Biobase::exprs(rle_data), 1, row_medians_assayData)

RLE_data <- as.data.frame(RLE_data)
RLE_data_gathered <- tidyr::gather(RLE_data, patient_array, log2_expression_deviation)

rle_boxplot <- ggplot(RLE_data_gathered, aes(patient_array, log2_expression_deviation)) + geom_boxplot(outlier.shape = NA) + 
                ylim(c(-2, 2)) + theme(axis.text.x = element_text(colour = 'aquamarine4', angle = 60, size = 6.5, hjust = 1 ,
                face = 'bold'))
ggsave('human_array_rle_boxplot.pdf', rle_boxplot, width = 12, height = 5)

#RMA calibration of data
data_norm <- oligo::rma(data)

#PCA analysis
exp_data <- Biobase::exprs(data_norm)
PCA <- prcomp(t(exp_data), scale = FALSE)

percentVar_new <- round(100*PCA$sdev^2/sum(PCA$sdev^2),1)
sd_ratio_new <- sqrt(percentVar_new[2] / percentVar_new[1])

dataGG <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     tissue = Biobase::pData(data_norm)$Tissue,
                     stage = Biobase::pData(data_norm)$Stage)


PCA_plot2 <- ggplot(dataGG, aes(PC1, PC2)) + geom_point(aes(shape = tissue, colour = stage)) + ggtitle('PCA plot of the calibrated, summarized data') +
                    xlab(paste0('PC1, VarExp: ', percentVar_new[1], '%')) + ylab(paste0('PC2, VarExp: ', percentVar_new[2], '%')) +
                    theme(plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio = sd_ratio_new) + scale_shape_manual(values = c(19)) + 
                    scale_color_manual(values = c('firebrick', 'goldenrod3', 'darkorange2', 'chartreuse3', 'darkolivegreen', 'dodgerblue4', 'cyan4', 'darkorchid4')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.key=element_blank())
ggsave('human_array_pca_calibrated_summarized.pdf', PCA_plot2, width = 10, height = 6)


dataGG2 <- data.frame(PC1 = PCA$x[,1], PC2 = PCA$x[,2],
                     fibrosis = Biobase::pData(data_norm)$fibrosis,
                     inflammation = Biobase::pData(data_norm)$inflammation)


PCA_plot3 <- ggplot(dataGG2, aes(PC1, PC2)) + geom_point(aes(shape = as.factor(inflammation), colour = as.factor(fibrosis))) + ggtitle('PCA plot of the calibrated, summarized data') +
                    xlab(paste0('PC1, VarExp: ', percentVar_new[1], '%')) + ylab(paste0('PC2, VarExp: ', percentVar_new[2], '%')) +
                    theme(plot.title = element_text(hjust = 0.5)) + coord_fixed(ratio = sd_ratio_new) + scale_shape_manual(values = c(19, 4, 15, 2)) + 
                    scale_color_manual(values = c('firebrick', 'goldenrod3', 'darkorange2', 'chartreuse3', 'darkolivegreen', 'dodgerblue4', 'cyan4')) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = 'black'), legend.key=element_blank())
ggsave('human_array_pca_calibrated_summarized_other.pdf', PCA_plot3, width = 10, height = 6)

#Heatmap clustering analysis
group_names <- ifelse(str_detect(pData(data_norm)$group, 'Control'), 'Control', 'Nash')
fibrosis_names <- ifelse(str_detect(pData(data_norm)$fibrosis, '0'), '0', '4')
nas_names <- ifelse(str_detect(pData(data_norm)$nas, '0'), '0', '7')
inflammation_names <- ifelse(str_detect(pData(data_norm)$inflammation, '0'), '0', '3')

annotation_for_heatmap <- data.frame(Group = group_names, Fibrosis_score = fibrosis_names, NAS = nas_names, Inflammation = inflammation_names)

row.names(annotation_for_heatmap) <- row.names(pData(data_norm))

dists <- as.matrix(dist(t(exp_data), method = 'manhattan'))

rownames(dists) <- row.names(pData(data_norm))
hmcol <- rev(colorRampPalette(RColorBrewer::brewer.pal(9, 'RdBu'))(255))
colnames(dists) <- NULL
diag(dists) <- NA

#ann_colors <- list(Group = c('' = 'chartreuse4', 'NASH' = 'burlywood3'), 
#                   Fibrosis_score = c('female' = 'blue4', 'male' = 'cadetblue2'),
#                   NAS = c('BS' = 'darkorange2', 'Diet' = 'dodgerblue4'),
#                   Inflammation = c('baseline' = 'firebrick', 'follow_up' = 'gold4'))

pdf('Clustering_heatmap_calibrated_sampels.pdf', width = 10, height = 15)
pheatmap(dists, col = (hmcol), annotation_row = annotation_for_heatmap,
         legend = TRUE, treeheight_row = 0, legend_breaks = c(min(dists, na.rm = TRUE), max(dists, na.rm = TRUE)), 
         legend_labels = (c('small distance', 'large distance')), main = 'Clustering heatmap for the calibrated samples')
dev.off()

#Filter lowly expressed genes based on intensity
palmieri_medians <- rowMedians(Biobase::exprs(data_norm))

pdf('intensity_hist.pdf')
hist_res <- hist(palmieri_medians, 100, col = "cornsilk1", freq = FALSE, main = "Histogram of the median intensities", 
                 border = "antiquewhite4", xlab = "Median intensities")
dev.off()

man_threshold <- 3

pdf('intensity_hist_thre.pdf')
hist_res <- hist(palmieri_medians, 100, col = "cornsilk", freq = FALSE, main = "Histogram of the median intensities",
                 border = "antiquewhite4", xlab = "Median intensities")
abline(v = man_threshold, col = "coral4", lwd = 2)
dev.off()

#Number of samples in exp groups
no_of_samples <- table(paste0(pData(data_norm)$Stage))
no_of_samples 

to_remove <- c(3, 4)
no_of_samples <- no_of_samples[!no_of_samples %in% to_remove]
samples_cutoff <- min(no_of_samples)

idx_man_threshold <- apply(Biobase::exprs(data_norm), 1,
                           function(x){
                             sum(x > man_threshold) >= samples_cutoff})
                             table(idx_man_threshold)

data_norm_filt <- subset(data_norm, idx_man_threshold)

#Annotating the transcript clusters
anno_data <- AnnotationDbi::select(hgu133plus2.db,
                                       keys = (featureNames(data_norm_filt)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")

anno_data <- subset(anno_data, !is.na(SYMBOL))

#Remove multiple mappings
anno_grouped <- group_by(anno_data, PROBEID)
anno_summarized <- dplyr::summarize(anno_grouped, no_of_matches = n_distinct(SYMBOL))
head(anno_summarized)
anno_filtered <- filter(anno_summarized, no_of_matches > 1)
head(anno_filtered)
probe_stats <- anno_filtered 
nrow(probe_stats)


ids_to_exlude <- (featureNames(data_norm_filt) %in% probe_stats$PROBEID)
table(ids_to_exlude)
data_norm_final <- subset(data_norm_filt, !ids_to_exlude)
validObject(data_norm_final)
head(anno_data)

fData(data_norm_final)$PROBEID <- rownames(fData(data_norm_final))
fData(data_norm_final) <- left_join(fData(data_norm_final), anno_data)

#Restore rownames after left_join
rownames(fData(data_norm_final)) <- fData(data_norm_final)$PROBEID 
validObject(data_norm_final)

#Fitting linear model to the data
individual <- as.character(Biobase::pData(data_norm_final)$index)

tissue <- str_replace_all(Biobase::pData(data_norm_final)$Liver_status, ' ', '_')
tissue <- ifelse(str_detect(Biobase::pData(data_norm_final)$Liver_status, 'no_NASH'), 'no_NASH', 'NASH')
tissue <- replace(tissue, tissue=='NASH', 'Liver')
tissue <- replace(tissue, tissue=='no_NASH', 'Liver')

disease <- str_replace_all(Biobase::pData(data_norm_final)$Liver_status, ' ', '_')
disease <- ifelse(str_detect(Biobase::pData(data_norm_final)$Liver_status, 'no_NASH'), 'no_NASH', 'NASH')

NASH_status <- individual[disease == 'no_NASH']
design_palmieri_CD <- model.matrix(~ 0 + tissue[disease == 'no_NASH'] + NASH_status)
colnames(design_palmieri_CD)[1:2] <- c('Liver')
rownames(design_palmieri_CD) <- NASH_status 


#Alternative (working)
fibrosis_score <- factor(pData(data_norm_final)$Stage, levels=c('Mild', 'Advanced'))
design <- model.matrix(~ 0 + fibrosis_score)
colnames(design) <- c('Mild', 'Advanced')
palmieri_fit <- lmFit(exprs(data_norm_final), design = design)

contrast_matrix <- makeContrasts(Advanced-Mild, levels = design)
palmieri_fit <- eBayes(contrasts.fit(lmFit(data_norm_final, design = design), contrast_matrix))
ext <- topTable(palmieri_fit, number = Inf)
head(ext)

ext2 <- as.data.frame(ext)
ext2 <- ext2[order(ext2$logFC),]

ext3 <- ext3[seq(dim(ext3)[1],1),]
write.table(ext, file = 'human_array_exprs.csv', sep = ',', col.names = NA)

#Exp for further anaylsis
data_for_exp <- exprs(data_norm_final)
rownames(data_for_exp) <- fData(data_norm_final)$SYMBOL
data_for_exp2 <- subset(data_for_exp, rownames(data_for_exp) %in% c('COL1A1' ,'GPX3', 'COL1A2', 'DPEP1', 'DPT', 'GAS6', 'MFAP4'))
data_for_exp2 <- as.data.frame(t(data_for_exp2))
data_for_exp2$Stage <- pData(data_norm_final)$Stage
write.table(data_for_exp2, file = 'human_array_corr_Murphy.csv', sep = ',', col.names = NA)

#HUM analysis
library(HUM)
human_array_corr <- read.delim('human_array_corr_Murphy.csv', sep = ',', row.names = 1)
human_array_corr <- na.omit(human_array_corr)
indexF=c(1:16);
indexClass=17;
label=unique(human_array_corr[,indexClass])
indexLabel=label
out=CalculateHUM_seq(human_array_corr, indexF, indexClass, indexLabel)

#Coordinates for plot
indexF=names(human_array_corr[,c(4),drop = FALSE])
indexClass=9
label=unique(human_array_corr[,indexClass])
indexLabel=label
out=CalculateHUM_seq(human_array_corr, indexF, indexClass, indexLabel)
HUM<-out$HUM
seq<-out$seq
out=CalculateHUM_ROC(human_array_corr, indexF, indexClass, indexLabel, seq)

out2 <- as.data.frame(out)
library(ggplot2)
ggplot(data = out2, mapping = aes(x = COL1A1, y = COL1A1.1)) + geom_point()

library(pROC)
df <- read.delim('C:/Users/mikekrogh/Desktop/single_cell_ccl4/human_micro_array/human_array_Murphy/human_array_corr_Murphy_roc - Copy - Copy.csv', sep = ',')

pdf('C:/Users/mikekrogh/Desktop/AUC_Murphy.pdf')
par(mfrow=c(3,3))
pROC_obj <- roc(df, 'Class', 'Marker1',
                smoothed = TRUE,
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)
dev.off()

pROC_M1 <- roc(df, 'Class', 'Marker1')
pROC_M2 <- roc(df, 'Class', 'Marker2')
pROC_M3 <- roc(df, 'Class', 'Marker3')
pROC_M4 <- roc(df, 'Class', 'Marker4')
pROC_M5 <- roc(df, 'Class', 'Marker5')
pROC_M6 <- roc(df, 'Class', 'Marker6')
pROC_M7 <- roc(df, 'Class', 'Marker7')

ROC <- data.frame(CI_low = c(paste(ci(pROC_M1))[1], paste(ci(pROC_M2))[1], paste(ci(pROC_M3))[1], paste(ci(pROC_M4))[1], paste(ci(pROC_M5))[1], paste(ci(pROC_M6))[1], paste(ci(pROC_M7))[1]), 
                  CI_high = c(paste(ci(pROC_M1))[3], paste(ci(pROC_M2))[3], paste(ci(pROC_M3))[3], paste(ci(pROC_M4))[3], paste(ci(pROC_M5))[3], paste(ci(pROC_M6))[3], paste(ci(pROC_M7))[3]),
                  SE = c(var(pROC_M1), var(pROC_M2), var(pROC_M3), var(pROC_M4), var(pROC_M5), var(pROC_M6), var(pROC_M7)))
rownames(ROC) <- c('DPEP1', 'MFAP4', 'GPX3', 'GAS6', 'COL1A1', 'COL1A2', 'DPT')


library(pROC)
df <- read.delim('C:/Users/mikekrogh/Desktop/single_cell_ccl4/human_micro_array/human_array_Murphy/human_array_corr_Murphy_roc - Copy - Copy.csv', sep = ',')

comb_roc <- function(df, combs, meta, markers){
require(pROC)
require(dplyr)
datalist = list()
for (i in 1:combs){
  df_temp <- data.frame(combn(colnames(df)[3:(markers+2)], i))
  names(df_temp)[1:ncol(df_temp)] <- paste0(rep('rep', each=ncol(df_temp)), '_', i, '_', 1:ncol(df_temp))
  datalist[[i]] <- df_temp
}

out_dat <- setNames(data.frame(matrix(ncol = 2, nrow = 0)), c('comb_auc', 'comb'))
for(w in 1:combs){
  for(m in 1:ncol(datalist[[w]])){
    m_comb <- paste0(datalist[[w]][[paste0('rep_', w, '_', m)]], collapse = ' + ')
    form <- as.formula(paste(meta, m_comb, sep = ' ~ '))
    mod1 <- glm(form, data = df, family = 'binomial')
    form1 <- as.formula(paste(meta, paste('predict', '(', 'mod1', ',', 'type', '=', "'",'response', "'", ')', sep = ''), sep = ' ~ '))
    suppressMessages(out_roc <- roc(form1, data = df))
    comb_roc <- round(as.numeric(paste0(out_roc$auc)), 4)
    comb_roc <- append(comb_roc, m_comb)
    comb_roc <- as.data.frame(t(comb_roc))
    names(comb_roc) <- names(out_dat) 
    suppressWarnings(out_dat <- bind_rows(out_dat, comb_roc))
    }}
return(out_dat)
}


mod1 <- glm(Class ~ Marker1 + Marker2 + Marker4 + Marker6 + Marker7, data = df, family = 'binomial')
out_roc <- roc(Class ~ predict(mod1, type = 'response'), data = df)
coords(pROC_M6, 'best', ret = 'all', transpose = FALSE, best.method='youden')

pdf('marker_comb_4_roc.pdf', width = 2.5, height = 3)
plot.roc(out_roc, asp = FALSE, legacy.axes = TRUE, col = '#F44D4D', grid.lty = 'solid', grid = c(0.20, 0.25), identity.lty = 'dashed', print.thres='best', print.thres.best.method='youden')
dev.off()

pdf('C:/Users/mikekrogh/Desktop/single_cell_ccl4/human_micro_array/human_array_Murphy/Dpt_Col1a2_Gas6_Mfap4_Dpep1_auc_0.990.pdf', width = 2.5, height = 3)
plot.roc(out_roc, asp = FALSE, legacy.axes = TRUE, col = '#F44D4D', grid.lty = 'solid', grid = c(0.20, 0.25), identity.lty = 'dashed')
dev.off()

#Calc se and ci for auc
se_auc <- function(auc, n1, n2, ci=FALSE, ci_per=95) {
  se <- sqrt(((auc*(1-auc))+(n1-1)*((auc/(2-auc))-auc^(2))+(n2-1)*(((2*auc^(2))/(1+auc))-auc^(2)))/(n1*n2))
  newData2 <- data.frame('Standard_error' = c(round(se, digits = 3)))
  if (ci==FALSE) return(print(newData2, row.names = FALSE))
  percen <- (100-ci_per)/100
  ci_low <- auc-qnorm(1-percen/2)*se
  ci_high <- auc+qnorm(1-percen/2)*se
  newData <- data.frame('Standard_error' = c(round(se, digits = 3)), 'Confidence_interval_lower' = c(round(ci_low, digits = 3)), 
                        'Confidence_interval_upper' = c(round(ci_high, digits = 3)))
  if (ci==TRUE) return(print(newData, row.names = FALSE))
}

se_auc(0.526, 40, 32, ci = TRUE, ci_per = 95)

t.test(df$Marker1~df$Class)
t.test(df$Marker2~df$Class)
t.test(df$Marker3~df$Class)
t.test(df$Marker4~df$Class)
t.test(df$Marker5~df$Class)
t.test(df$Marker6~df$Class)
t.test(df$Marker7~df$Class)

p_val <- data.frame(Marker = c('DPEP1', 'MFAP4', 'GPX3', 'GAS6', 'COL1A1', 'COL1A2', 'DPT'), 
                    pval = c(paste(t.test(df$Marker1~df$Class))[3], paste(t.test(df$Marker2~df$Class))[3], paste(t.test(df$Marker3~df$Class))[3],
                             paste(t.test(df$Marker4~df$Class))[3], paste(t.test(df$Marker5~df$Class))[3], paste(t.test(df$Marker6~df$Class))[3],
                             paste(t.test(df$Marker7~df$Class))[3]))

p_val$BH <- p.adjust(pval, method = 'BH')

#Plot distribution
df <- read.delim('C:/Users/mikekrogh/Desktop/single_cell_ccl4/human_micro_array/human_array_Murphy/human_array_corr_Murphy_roc_plot.csv', sep = ',')

library(ggplot2)
library(ggpubr)
# p <- ggplot(df, aes(x=Class, y=DPT)) + geom_dotplot(binaxis='y', stackdir='center', dotsize = 0.4) + stat_summary(fun.data=mean_sdl, fun.args = list(mult=1), geom="pointrange", color='grey', shape = 18, size = 0.5) +
#       scale_fill_grey() + theme_classic() + scale_x_discrete(limits=c('Mild', 'Advanced'))


p <- ggdotplot(df, x = 'Class', y = 'DPEP1', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040') + 
  stat_compare_means(aes(group = Class), method = 't.test', label.y.npc = 1) + stat_compare_means(aes(group = Class), method = 't.test', label = 'p.signif', label.y.npc = 0.9)
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('DPEP1_dist_plot.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'MFAP4', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040') + 
  stat_compare_means(aes(group = Class), method = 't.test', label.y.npc = 1) + stat_compare_means(aes(group = Class), method = 't.test', label = 'p.signif', label.y.npc = 0.9)
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('MFAP4_dist_plot.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'GPX3', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040') + 
  stat_compare_means(aes(group = Class), method = 't.test', label.y.npc = 1) + stat_compare_means(aes(group = Class), method = 't.test', label = 'p.signif', label.y.npc = 0.9)
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('GPX3_dist_plot.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'GAS6', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040') + 
  stat_compare_means(aes(group = Class), method = 't.test', label.y.npc = 1) + stat_compare_means(aes(group = Class), method = 't.test', label = 'p.signif', label.y.npc = 0.9)
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('GAS6_dist_plot.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'COL1A1', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040') + 
  stat_compare_means(aes(group = Class), method = 't.test', label.y.npc = 1) + stat_compare_means(aes(group = Class), method = 't.test', label = 'p.signif', label.y.npc = 0.9)
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('COL1A1_dist_plot.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'COL1A2', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040') + 
  stat_compare_means(aes(group = Class), method = 't.test', label.y.npc = 1) + stat_compare_means(aes(group = Class), method = 't.test', label = 'p.signif', label.y.npc = 0.9)
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('COL1A2_dist_plot.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'DPT', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040') + 
  stat_compare_means(aes(group = Class), method = 't.test', label.y.npc = 1) + stat_compare_means(aes(group = Class), method = 't.test', label = 'p.signif', label.y.npc = 0.9)
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('DPT_dist_plot.pdf', p, width = 3, height = 4)


p <- ggdotplot(df, x = 'Class', y = 'DPEP1', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040')
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('DPEP1_dist_plot_WO.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'MFAP4', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040')
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('MFAP4_dist_plot_WO.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'GPX3', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040')
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('GPX3_dist_plot_WO.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'GAS6', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040')
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('GAS6_dist_plot_WO.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'COL1A1', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040')
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('COL1A1_dist_plot_WO.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'COL1A2', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040')
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('COL1A2_dist_plot_WO.pdf', p, width = 3, height = 4)
p <- ggdotplot(df, x = 'Class', y = 'DPT', size = 0.9, color = '#F44D4D', fill = '#F44D4D', order = c('Mild', 'Advanced'), font.x = '#404040', font.y = '#404040', font.tickslab = '#404040')
p <- p + stat_summary(fun.y = mean, geom = "errorbar", aes(ymax = ..y.., ymin = ..y..), width = .60, linetype = 'solid') + stat_summary(fun.data = mean_se, geom = "errorbar", fun.args = list(mult=1), width = .40) + theme(axis.line = element_line(color = '#404040'), axis.ticks = element_line(colour='#404040'))
ggsave('DPT_dist_plot_WO.pdf', p, width = 3, height = 4)



