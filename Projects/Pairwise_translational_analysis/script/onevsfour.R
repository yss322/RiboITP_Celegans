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
                                 experiment = c('1cell_A', '1cell_B','4cell_B'),
                                 length      = TRUE,
                                 transcript  = FALSE,
                                 tidy = F,
                                 region      = c("CDS") )
ribo_cell_new <-  get_region_counts(new.ribo,
                                    range.lower = 25,
                                    range.upper = 35,
                                    experiment = c('4cell_C','4cell_D'),
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
setcolorder(rcw_diff_ribo, c("transcript", "1cell_A",    "1cell_B"  ,"4cell_B","4cell_C"  ,  "4cell_D"   ))
colnames(rcw_diff_ribo)<- c("transcript", "1cell_A_ribo",    "1cell_B_ribo" ,"4cell_B_ribo","4cell_C_ribo"  ,  "4cell_D_ribo")
#processing RNA-seq 
rnaseq_diff <- get_rnaseq(ribo.object = original.ribo,
                          tidy        = F,
                          compact = F,
                          experiment = c('1cell_A','1cell_B','4cell_A','4cell_B'),
                          region = "CDS")
rnaseq_diff <- as.data.table(rnaseq_diff)

rnaseq_diff [,CDS := (CDS+1)]
rnaseq_w_diff = dcast(rnaseq_diff, transcript ~ experiment)                          

all_counts_diff = merge(rcw_diff_ribo,rnaseq_w_diff, by= "transcript")
all_counts_diff <- all_counts_diff[, transcript := as.character(transcript)]
all_counts_diff[, gene_name := sapply(strsplit(transcript, "gene_symbol:"), function(x) tail(x, n = 1))]

all_counts_diff [,"transcript":= NULL ]

all_counts_diff <- setcolorder(all_counts_diff, "gene_name") 


exptype <- factor(c( "onecell.Ribo",
                     "onecell.Ribo",
                    
                     "fourcell.Ribo",
                     "fourcell.Ribo",
                     "fourcell.Ribo",
                     "onecell.RNA",
                     "onecell.RNA",
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
  TE_fourcellvsonecell = (fourcell.Ribo - fourcell.RNA) - (onecell.Ribo - onecell.RNA) ,
  Ribofourvsone =  fourcell.Ribo - onecell.Ribo  ,
  RNA_fourvsone = fourcell.RNA - onecell.RNA ,
  levels = design
)

qlf <- glmQLFTest(fit, contrast=my.contrasts[,"RNA_fourvsone"])    
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"Ribofourvsone"])    
qlf <- glmQLFTest(fit, contrast=my.contrasts[,"TE_fourcellvsonecell"])


summary(decideTests(qlf, p.value = 0.001, adjust.method = "fdr"))
plotMD(qlf, p.value = 0.001, hl.cex = 0.75, main = "",
       hl.col = color.palette0(2), legend = F)
top <-topTags(qlf,n = Inf)
head (qlf$genes)

view(top)




#first convert my list to Entrez gene IDS

ensembl <- useEnsembl(biomart = "genes", dataset = "celegans_gene_ensembl")
listFilters(ensembl) 

listofgenes<- qlf$genes

entrez_IDs <- getBM(attributes = c('entrezgene_id','wormbase_gene','external_gene_name'),
                    filters = 'external_gene_name',
                    values = listofgenes ,
                    mart = ensembl) 

colnames(entrez_IDs) <- c('entrezgene_id','wormbase_gene','gene_name')
entrez_IDs


qlf_entrez <- merge(qlf,entrez_IDs, by = 'gene_name')

head(qlf_entrez)

toplist <- topTags(qlf,p.value = 0.015, n = Inf)

qlf_entrez <-as.data.table( merge(toplist,entrez_IDs, by = 'gene_name'))
qlf_entrez_up <- qlf_entrez[logFC>0]
qlf_entrez_down <- qlf_entrez[logFC<0]
# now onto go_annotation 
genelist <- qlf_entrez$entrezgene_id
genelist <- qlf_entrez_up$entrezgene_id
genelist <- qlf_entrez_down$entrezgene_id


GO <-goana(genelist, species = "Ce")
view(topGO(GO, n = Inf, p.value = 0.05))




