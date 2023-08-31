# first import the file - get both RNAseq and Ribodata 

library(ribor)
library(tidyverse)
library(edgeR) 
library(data.table)

color.palette0 = colorRampPalette(c("#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"), space="Lab")
new.ribo <- Ribo("/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/20230720_additional_replicates/celegans.ribo")
original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo")


#get the required region counts for 1 and 2-cell data 

ribo_cell_or<- get_region_counts(original.ribo,
                                      range.lower = 25,
                                      range.upper = 35,
                                      experiment = c('4cell_A','4cell_B','8cell_A'),
                                      length      = TRUE,
                                      transcript  = FALSE,
                                      tidy = F,
                                      region      = c("CDS") )
ribo_cell_new <-  get_region_counts(new.ribo,
                                      range.lower = 25,
                                      range.upper = 35,
                                      experiment = c('4cell_C','4cell_D','8cell_D'),
                                      length      = TRUE,
                                      transcript  = FALSE,
                                      tidy = F,
                                      region      = c("CDS") )
ribo_cell_or<- as.data.table(ribo_cell_or)

ribo_cell_or [,CDS := (CDS+1)]
ribo_cell_new <-as.data.table(ribo_cell_new)
ribo_cell_new [,CDS := (CDS+1)]

rcw_diff_or = dcast(ribo_cell_or, transcript ~ experiment)  

rcw_diff_new = dcast(ribo_cell_new , transcript ~ experiment)  

rcw_diff_ribo = merge.data.table(rcw_diff_or,rcw_diff_new,by='transcript')

setcolorder(rcw_diff_ribo, c('transcript','4cell_A','4cell_B','4cell_C','4cell_D','8cell_A','8cell_D'))
#processing RNA-seq 
rnaseq_diff <- get_rnaseq(ribo.object = original.ribo,
                          tidy        = F,
                          compact = F,
                          experiment = c('4cell_A','4cell_B','8cell_A','8cell_B'),
                          region = "CDS")
rnaseq_diff <- as.data.table(rnaseq_diff)
rnaseq_w_diff = dcast(rnaseq_diff, transcript ~ experiment)   

colnames(rnaseq_w_diff) <- c('transcript','4cell_A_RNA','4cell_B_RNA','8cell_A_RNA','8cell_B_RNA')

all_counts_diff = merge(rcw_diff_ribo,rnaseq_w_diff, by= "transcript")
all_counts_diff

exptype <- factor(c( "fourcell.Ribo",
                     "fourcell.Ribo",
                     "fourcell.Ribo",
                     "fourcell.Ribo",
                     "eightcell.Ribo",
                     "eightcell.Ribo",
                     "fourcell.RNA",
                     "fourcell.RNA",
                     "eightcell.RNA",
                     "eightcell.RNA")  )
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
  TE_eightcellvsfourcell = (eightcell.Ribo - eightcell.RNA) - (fourcell.Ribo - fourcell.RNA) ,
  Ribo_eightvsfour =  eightcell.Ribo - fourcell.Ribo  ,
  RNA_eightvsfour = eightcell.RNA - fourcell.RNA ,
  levels = design
)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"RNA_eightvsfour"])    
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"Ribo_eightvsfour"])    
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"TE_eightcellvsfourcell"])


summary(decideTests(qlf, p.value = 0.99, adjust.method = "fdr"))
plotMD(qlf, 
       ylim =c(-4,4), 
       p.value = 0.05,
       hl.cex = 0.75, 
       main = "",
       hl.col = color.palette0(2), 
       legend = F)








