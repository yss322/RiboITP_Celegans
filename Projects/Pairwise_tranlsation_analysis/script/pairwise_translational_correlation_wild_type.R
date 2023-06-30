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


exptype <- factor(c("onecell.RNA","onecell.RNA",
                    "twocell.RNA","twocell.RNA",
                    "onecell.Ribo","onecell.Ribo","twocell.Ribo",
                    "twocell.Ribo" )  )
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

top_value <-topTags(qlf, n =1000 , p.value = 0.01, adjust.method = "fdr")
tables <-top_value$table
summary(decideTests(qlf, p.value = 0.05, adjust.method = "fdr"))
plotMD(qlf, p.value = 0.005, hl.cex = 0.75, main = "",
       hl.col = color.palette0(2), legend = F)
head (qlf$genes)

