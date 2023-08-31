# first import the file - get both RNAseq and Ribodata 

library(ribor)
library(tidyverse)
library(edgeR) 
library(data.table)

color.palette0 = colorRampPalette(c("#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"), space="Lab")
original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo")
original.ribo

#get the required region counts for 1 and 2-cell data 

ribo_one_two_cell<- get_region_counts(original.ribo,
                        range.lower = 25,
                        range.upper = 35,
                        experiment = c('1cell_A','1cell_B','2cell_B'),
                        length      = TRUE,
                        transcript  = FALSE,
                        tidy = F,
                        region      = c("CDS") )
ribo_two_cell_A <-  get_region_counts(original.ribo,
                                              range.lower = 28,
                                              range.upper = 38,
                                              experiment = c('2cell_A'),
                                              length      = TRUE,
                                              transcript  = FALSE,
                                              tidy = F,
                                              region      = c("CDS") )
ribo_one_two_cell<- as.data.table(ribo_one_two_cell)
ribo_two_cell_A <-as.data.table(ribo_two_cell_A)

rcw_diff_one_two = dcast(ribo_one_two_cell, transcript ~ experiment)  

rcw_diff_two_A = dcast(ribo_two_cell_A, transcript ~ experiment)  

rcw_diff_ribo = inner_join(rcw_diff_one_two,rcw_diff_two_A, by = "transcript")

#processing RNA-seq 
rnaseq_diff <- get_rnaseq(ribo.object = original.ribo,
                          tidy        = F,
                          compact = F,
                          experiment = c('1cell_A','1cell_B','2cell_A','2cell_B'),
                          region = "CDS")
rnaseq_diff <- as.data.table(rnaseq_diff)
rnaseq_w_diff = dcast(rnaseq_diff, transcript ~ experiment)                          

all_counts_diff = merge(rnaseq_w_diff, rcw_diff_ribo, by= "transcript")
all_counts_diff


exptype <- factor(c( "twocell.Ribo","twocell.Ribo","fourcell.Ribo",
                    "fourcell.Ribo","twocell.Ribo", "fourcell.Ribo","fourcell.Ribo",
                    "twocell.RNA","twocell.RNA","fourcell.RNA","fourcell.RNA")  )
y <- DGEList(counts=all_counts_diff[,-1],
             group=exptype, genes = all_counts_diff[,1])

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~0+exptype)
colnames(design) <- levels(exptype)
y <- estimateDisp(y,design)
plotBCV(y)
plotMDS(y)

fit <- glmQLFit(y,design)

my.contrasts <- makeContrasts(
  TE_onecellvstwocell = (twocell.Ribo - twocell.RNA) - (onecell.Ribo - onecell.RNA) ,
  Ribo_onevstwo =  twocell.Ribo - onecell.Ribo  ,
  RNA_onevstwo = twocell.RNA - onecell.RNA ,
  levels = design
)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"RNA_onevstwo"])    
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"Ribo_onevstwo"])    
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"TE_onecellvstwocell"])
qlf

summary(decideTests(qlf, p.value = 0.05, adjust.method = "fdr"))
plotMD(qlf, p.value = 0.05, hl.cex = 0.75, main = "",
       hl.col = color.palette0(2), legend = F)
head (qlf$genes)
 

# to get top hits 

top_value <-topTags(qlf, n =1000 , p.value = 0.01, adjust.method = "fdr")
tables <-top_value$table

#to point out genes that are predominately present in the p-cell 

p_cell_list <- fread("/Users/yashshukla/Documents/Computational/Celegans_riboITP/Projects/Pairwise_tranlsation_analysis/data/exported_genes.csv")
logFC <- y$logFC
p_cell_highlight <- p_cell_list[,1]

colors <- rep("black", nrow(qlf$table))
colors[y$genes %in% p_cell_highlight] <- "red"


plotMD(qlf, col = colors, cex = 0.8)



# Create a vector to store the colors for each gene
colors <- rep("black", nrow(qlf$table))
colors[rownames(qlf$table) %in% p_cell_highlight] <- "red"  # Set the color for genes of interest

# Create the plot with highlighted genes of interest
plotMD(qlf, col = colors, p.value = 0.05, hl.cex = 100000, main = "",
       hl.col = colors, legend = FALSE)


