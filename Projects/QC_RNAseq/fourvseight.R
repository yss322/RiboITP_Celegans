
library(ribor)
library(tidyverse)
library(edgeR) 
library(data.table)

color.palette0 = colorRampPalette(c("#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"), space="Lab")
original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo")
original.ribo 


#comparing 2-cell to 4-cell stage 


rnaseq_diff <- get_rnaseq(ribo.object = original.ribo,
                          tidy        = F,
                          compact = F,
                          experiment = c('4cell_A','4cell_B','8cell_A','8cell_B'),
                          region = "CDS")
rnaseq_diff <- as.data.table(rnaseq_diff)
rnaseq_w_diff = dcast(rnaseq_diff, transcript ~ experiment)     
rnaseq_w_diff 
exptype <- factor(c("fourcell.RNA","fourcell.RNA",
                    "eightcell.RNA","eightcell.RNA"))
y <- DGEList(counts=rnaseq_w_diff[,-1],
             group=exptype, genes = rnaseq_w_diff[,1])
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
  RNA_eightvsfour = eightcell.RNA - fourcell.RNA ,
  levels = design)
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"RNA_eightvsfour"])   
summary(decideTests(qlf, p.value = 0.05, adjust.method = "fdr"))
plotMD(qlf, p.value = 0.05, hl.cex = 0.75, main = "",
       hl.col = color.palette0(2), legend = F)

