# first import the file - get both RNAseq and Ribodata 

library(ribor)
library(tidyverse)
library(edgeR) 
library(data.table)

color.palette0 = colorRampPalette(c("#98abc5", "#8a89a6", "#7b6888", "#6b486b", "#a05d56", "#d0743c", "#ff8c00"), space="Lab")
new.ribo <- Ribo("/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/20230720_additional_replicates/celegans.ribo")
original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo")


#get the required region counts for 2 and 4-cell data 

ribo_cell_or<- get_region_counts(original.ribo,
                                      range.lower = 25,
                                      range.upper = 35,
                                      experiment = c('4cell_B', '2cell_B'),
                                      length      = TRUE,
                                      transcript  = FALSE,
                                      tidy = F,
                                      region      = c("CDS") )
ribo_cell_new <-  get_region_counts(new.ribo,
                                      range.lower = 25,
                                      range.upper = 35,
                                      experiment = c('2cell_C','4cell_C','4cell_D'),
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

colnames(rcw_diff_ribo)
setcolorder(rcw_diff_ribo, c("transcript",    "2cell_B",  "2cell_C"   ,"4cell_B","4cell_C"  ,  "4cell_D"   ))
colnames(rcw_diff_ribo)<- c("transcript",    "2cell_B_ribo",  "2cell_C_ribo"   ,"4cell_B_ribo","4cell_C_ribo"  ,  "4cell_D_ribo")
#processing RNA-seq 
rnaseq_diff <- get_rnaseq(ribo.object = original.ribo,
                          tidy        = F,
                          compact = F,
                          experiment = c('2cell_A','2cell_B','4cell_A','4cell_B'),
                          region = "CDS")
rnaseq_diff <- as.data.table(rnaseq_diff)

rnaseq_diff [,CDS := (CDS+1)]
rnaseq_w_diff = dcast(rnaseq_diff, transcript ~ experiment)                          

all_counts_diff = merge(rcw_diff_ribo,rnaseq_w_diff, by= "transcript")


exptype <- factor(c("twocell.Ribo",
                     "twocell.Ribo",
                     "fourcell.Ribo",
                     "fourcell.Ribo",
                     "fourcell.Ribo",
                     "twocell.RNA",
                     "twocell.RNA",
                     "fourcell.RNA",
                     "fourcell.RNA")  )
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
  TE_fourcellvstwocell = (fourcell.Ribo - fourcell.RNA) - (twocell.Ribo - twocell.RNA) ,
  Ribofourvstwo =  fourcell.Ribo - twocell.Ribo  ,
  RNA_fourvstwo = fourcell.RNA - twocell.RNA ,
  levels = design
)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"RNA_fourvstwo"])    
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"Ribofourvstwo"])    
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"TE_fourcellvstwocell"])


summary(decideTests(qlf, p.value = 0.05, adjust.method = "fdr"))
plotMD(qlf, 
       ylim =c(-8,8), 
       p.value = 0.05,
       hl.cex = 0.75, 
       main = "",
       hl.col = color.palette0(2), 
       legend = F)
top <-topTags(qlf)
head (qlf$genes)

view(top)


x_coord <- qlf$table$logCPM  # Replace with your actual x-coordinate
y_coord <- qlf$table$logFC  # Replace with your actual y-coordinate
# Replace with the label you want to display

# Run the identify function to interactively label the point
identify(x_coord, y_coord )

qlf$table[('4241'),]
