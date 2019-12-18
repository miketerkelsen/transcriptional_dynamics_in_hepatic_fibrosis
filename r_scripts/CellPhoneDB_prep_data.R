#Replace genenames with Ensembl ID
#Import 10x gtf file and remove duplicate genes
gtf <- rtracklayer::import('genes.gtf')
gtf_df=as.data.frame(gtf)
gtf_df2 <- as.data.frame(gtf_df$gene_id)
names(gtf_df2)[1] <- 'gene_id'
gtf_df2['gene_name'] <- gtf_df$gene_name
gtf_df3 <- distinct(gtf_df2)
gtf_df4 <- gtf_df3[!duplicated(gtf_df3$gene_name), ]

#Export own gene list from expression matrix, convert some names which have been changed and 
#merge to get gene/Ensembl ID dataframe
row_matrix <- as.data.frame(rownames(matrix3))
names(row_matrix)[1] <- 'gene_name'
row_matrix <- as.data.frame(gsub('X(.*Rik)', '\\1', row_matrix$gene_name))
names(row_matrix)[1] <- 'gene_name'
row_matrix <- as.data.frame(gsub('[.]', '-', row_matrix$gene_name))
names(row_matrix)[1] <- 'gene_name'
row_matrix2 <- merge(row_matrix, gtf_df4, by.all="gene_name", sort=FALSE, all.x=T)
#row_matrix2 <- inner_join(row_matrix, gtf_df3, by='gene_name')

#Remove conflicting rows not having Ensembl ID
row_names_remove <- c("CR974586-5", "4930556M19Rik-1", "Gt-ROSA-26Sor", "AC168977-2", "Sept15", "Pcdhga8-1", "Smim20-1", "Gbp6-1", "Gm16701-1", "AC168977-1", "AC125149-3", "CAAA01147332-1", "Atn1-1", "Dancr-1", "Schip1-1")
matrix3 <- matrix2[!(row.names(matrix2) %in% row_names_remove), ]

#Replace rownames in expression matrix to check compliance
rownames(matrix3) <- row_matrix2$gene_id
matrix3[1:10,1:10]

#Import human/mouse Ensembl gene set
require("biomaRt")
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

#Do biomaRt comparison
genes_con <- row_matrix2$gene_id
genesV2 = getLDS(attributes = c("ensembl_gene_id"), filters = "ensembl_gene_id", values = genes_con, mart = mouse, attributesL = c("ensembl_gene_id"), martL = human)
names(genesV2)[1] <- 'gene_id'
names(genesV2)[2] <- 'gene_id_human'

#Merge new Ensembl comparison with preivous export from expression matrix to get the right order(Use only one below one will work better)
row_matrix3 <- merge(row_matrix2, genesV2, by.all="gene_id", sort=FALSE, all.x=T)
row_matrix3 <- left_join(row_matrix2, genesV3, by='gene_id')

#If dupicated gene names present remove duplicates(To remove the above 2 lines)
genesV3 <- genesV2[!duplicated(genesV2$gene_id), ]
row_matrix3 <- left_join(row_matrix2, genesV3, by='gene_id')

#Import human Ensembl back into expression matrix as rownames and remove rows with not imported Ensembl
rownames(matrix3) <- row_matrix3$gene_id_human
matrix4 <- matrix3[!is.na(rownames(matrix3)), ]
