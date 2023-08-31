library(ribor)
library(data.table)
library(corrplot)
# Import data 

original.ribo <- Ribo('/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/celegans_duplicate_allstages.ribo')

new.ribo  <- Ribo('/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/20230720_additional_replicates/celegans.ribo') 


ribo_cell_or<- get_region_counts(original.ribo,
                                 range.lower = 25,
                                 range.upper = 35,
                                 experiment = c('1cell_A', '1cell_B','2cell_A', '2cell_B','4cell_A', '4cell_B','8cell_A', '8cell_B'),
                                 length      = TRUE,
                                 transcript  = FALSE,
                                 tidy = F,
                                 region      = c("CDS") )
ribo_cell_new <-get_region_counts(new.ribo,
                                                 range.lower = 25,
                                                 range.upper = 35,
                                                 experiment = c('1cell_C', '2cell_C', '4cell_C','4cell_D', '8cell_C', '8cell_D'),
                                                 length      = TRUE,
                                                 transcript  = FALSE,
                                                 tidy = F,
                                                 region      = c("CDS") )


# Create a correlation heatmap of ribosome occupancy data 
# 1st column should be all transcript names and the next columns are data from each cell replicate 

ribo_cell_or<- as.data.table(ribo_cell_or)

ribo_cell_or [,CDS := (CDS+1)]
ribo_cell_new <-as.data.table(ribo_cell_new)
ribo_cell_new [,CDS := (CDS+1)]

rcw_diff_or = dcast(ribo_cell_or, transcript ~ experiment)  

rcw_diff_new = dcast(ribo_cell_new , transcript ~ experiment)  

rcw_diff_ribo = merge.data.table(rcw_diff_or,rcw_diff_new,by='transcript')
rcw_diff_ribo_corr<- rcw_diff_ribo [,-1]
data <- cor(rcw_diff_ribo_corr, method = "spearman") 


corrplot(data, method = 'color',type = 'upper',col = COL1('Oranges', 500))# Create a correlation heatmap of RNA-seq data 

data_outlier_removed <- rcw_diff_ribo_corr [,-c('8cell_B','1cell_C','8cell_C')]
data_out <- cor(data_outlier_removed, method = "spearman") 
corrplot(data_out, method = 'color',type = 'upper',col = COL1('Oranges', 500), col.lim = c(0.7,1), is.corr = FALSE)# Create a correlation heatmap of RNA-seq data 

# RNA-seq to RNA-seq correlation plot  

rnaseq_diff <- get_rnaseq(ribo.object = original.ribo,
                          tidy        = F,
                          compact = F,
                          experiment = c('1cell_A', '1cell_B','2cell_A', '2cell_B','4cell_A', '4cell_B','8cell_A', '8cell_B'),
                          region = "CDS")
rnaseq_diff <- as.data.table(rnaseq_diff)

rnaseq_diff [,CDS := (CDS+1)]
rnaseq_w_diff = dcast(rnaseq_diff, transcript ~ experiment)                          

rnaseq_w_diff_corr<- rnaseq_w_diff [,-1]
data_rnaseq <- cor(rnaseq_w_diff_corr, method = "spearman") 
corrplot(data_rnaseq, method = 'color',type = 'upper',col = COL1('Oranges',200), col.lim = c(0.7,1), is.corr = FALSE)#


#RNA-seq vs Ribo correlation plot  

all_counts_diff = merge(rcw_diff_ribo,rnaseq_w_diff, by= "transcript")
all_counts_diff_corr <- all_counts_diff[,-1]
data_all <-  cor(all_counts_diff_corr,method = "spearman")
data_rna_ribo <- data_all[(1:14),-(1:14)]
data_rna_ribo_outlier_remove <- data_rna_ribo [c('1cell_A.x','1cell_B.x','2cell_A.x','2cell_B.x','4cell_A.x','4cell_B.x','8cell_A.x','2cell_C','4cell_C','4cell_D','8cell_D'),]
corrplot(data_rna_ribo_outlier_remove , method = 'color',col = COL1('Oranges',200), col.lim = c(0.5,1), is.corr = FALSE)#
