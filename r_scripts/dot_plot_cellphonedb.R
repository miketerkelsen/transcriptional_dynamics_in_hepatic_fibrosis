#Subsetting
setwd('C:/Users/mikekrogh/Google Drive/SDU/Master_PhD/CellPhoneDB/out_louvain_2nd')
means <- read.delim('means_new_new.txt', header = TRUE, sep = '\t')

means2 <- means[rowSums(means[12:75] > 1) > 0, ]
write.table(means2$interacting_pair, file = 'row_output_thr_3.txt', sep = '\t', col.names = FALSE, row.names = FALSE)


#Dotplot
#height 24 interactions = 7
#height 74 interactions = 21

library(ggplot2)
dot_plot = function(selected_rows = NULL,
                    selected_columns = NULL,
                    filename = 'dot_plot_HSC_thr_3.pdf',
                    width = 15,
                    height = 7,
                    means_path = './means.txt',
                    pvalues_path = './pvalues.txt',
                    means_separator = '\t',
                    pvalues_separator = '\t',
                    output_extension = '.pdf'
){
  
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  rows_use = read.delim('./row_output_thr_3.txt', sep = '\t', header = FALSE)
  rows_use = as.vector(rows_use$V1)
  cols_use = read.delim('./columns.txt', sep = '\t', header = FALSE)
  cols_use = as.vector(cols_use$V1)
  
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = rows_use
  }
  
  if(is.null(selected_columns)){
    selected_columns = cols_use
  }
  
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  pr[pr==0] = 1
  plot.data = cbind(plot.data,log2(pr))
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  
  my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
  
  ggplot(plot.data,aes(x=clusters,y=pair)) +
    geom_point(aes(size=-log10(pvalue),color=mean)) +
    scale_color_gradientn('Log2 mean (Molecule 1, Molecule 2)', colors=my_palette) +
    theme_bw() +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          axis.text=element_text(size=14, colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, family = 'Arial'),
          axis.text.y = element_text(size=12, colour = "black", family = 'Arial'),
          axis.title=element_blank(),
          text = element_text('Arial'),
          panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
  
  if (output_extension == '.pdf') {
    ggsave(filename, width = width, height = height, device = cairo_pdf, limitsize=F)
  }
  else {
    ggsave(filename, width = width, height = height, limitsize=F)
  }
}
