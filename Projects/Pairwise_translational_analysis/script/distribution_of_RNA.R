
# first import the file - get both RNAseq and Ribodata 

library(ribor)
library(tidyverse)
library(edgeR) 
library(data.table)
library(ggpubr)
library(stack)
library(EnhancedVolcano)
library(biomaRt)


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

#adding pseudocounts of +1 

ribo_one_two_cell [,CDS := (CDS+1)]
ribo_two_cell_A [,CDS := (CDS+1)]

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


# Extract the gene name using strsplit and data.table
all_counts_diff <- all_counts_diff[, transcript := as.character(transcript)]
all_counts_diff[, gene_name := sapply(strsplit(transcript, "gene_symbol:"), function(x) tail(x, n = 1))]

all_counts_diff [,"transcript":= NULL ]

all_counts_diff <- setcolorder(all_counts_diff, "gene_name") 


# Print the updated data.table


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
plotMD(qlf, 
       ylim =c(-10,10), 
       p.value = 0.05,
       hl.cex = 0.5, 
       main = "",
       hl.col = color.palette0(2), 
       legend = F)

#assessing p-value distribution 
out <- topTags(qlf, n = "Inf")$table
ggplot(out, aes(x = PValue)) + 
  geom_histogram(bins = 70)

check <- topTags(qlf)
summary(decideTests(qlf, p.value = 0.05, adjust.method = "fdr"))
plotMD(qlf, p.value = 0.05, hl.cex = 0.75, main = "",
       hl.col = color.palette0(2), legend = F)
ggplot(qlf$table, aes(x = logFC))+
  geom_histogram(bins = 55, ) + 
  scale_x_continuous(breaks = -20:10)

#using ggplot 

AB_cell_list <- fread("/Users/yashshukla/Documents/Computational/Celegans_riboITP/Projects/Pairwise_tranlsation_analysis/data/exported_genes.csv")

AB_cell_highlight <- AB_cell_list[,1]



overlap_AB_cell <- qlf$genes$gene_name %in% AB_cell_highlight$V1 

highlight_AB_cell <- rownames(qlf$genes)[overlap_AB_cell]

P_cell_list <- fread("/Users/yashshukla/Documents/Computational/Celegans_riboITP/Projects/Pairwise_tranlsation_analysis/data/exported_genes_P-cell.csv")
p_cell_highlight <- P_cell_list[,1]
overlap_P_cell <-qlf$genes$gene_name %in% p_cell_highlight$V1
highlight_p_cell <- rownames(qlf$genes)[overlap_P_cell]
ggplot(qlf$table, 
       aes(
         x = logFC , 
         y = -log(PValue))) + 
  geom_jitter(data = subset(qlf$table,rownames(qlf$table) %in% highlight_AB_cell), color = "red",size = 1.5) +
  geom_jitter(data = subset(qlf$table,rownames(qlf$table) %in% highlight_p_cell), color = "blue",size = 1.5) +
  geom_jitter(color ="black",size = 0.8, alpha = 0.3)

Oma_pull_down <- read_excel("/Users/yashshukla/Documents/Computational/Celegans_riboITP/Projects/Pairwise_tranlsation_analysis/data/FPKM_OMA-1_Pull_down_wbgene.xlsx")
oma_highlight <- Oma_pull_down[,2]
overlap_oma_pull <-qlf$genes$gene_name %in% oma_highlight$gene_id
highlight_oma_pull <- rownames(qlf$genes)[overlap_oma_pull]
ggplot(qlf$table, aes(x = logFC , y = -log(PValue))) + 
  geom_jitter(data = subset(qlf$table,rownames(qlf$table) %in% highlight_oma_pull), color = "red",size = 1.5) +
  geom_jitter(color ="black",size = 0.8, alpha = 0.3)

qlf$table



keyvals = ifelse(out$logFC < 0 & out$FDR < 0.05, '#6e005f',
                 ifelse(out$logFC > 0 & out$FDR < 0.05, '#045275', 'grey60'))
names(keyvals)[keyvals == '#6e005f'] <- 'Down'
names(keyvals)[keyvals == '#045275'] <- 'Up'
names(keyvals)[keyvals == 'grey60'] <- 'NS'



EnhancedVolcano(out,
                xlim = c(-10, 7),
                ylim = c(0, 3.5),
                lab = out$gene_name,
                x = 'logFC',
                y = 'FDR',
                title = 'Translational efficiency',
                subtitle = NULL ,
                axisLabSize = 14,
                titleLabSize = 16,
                subtitleLabSize = 14,
                legendIconSize = 2.0,
                legendLabSize = 12,
                pCutoff = 0.05,
                FCcutoff = 0, 
                cutoffLineType = 'blank',
                vline = c(0),
                #                vlineCol = c('grey50', 'grey0','grey50'),
                vlineType = 'blank',
                vlineWidth = 0.5,
                drawConnectors = TRUE,
                pointSize = c(ifelse(out$FDR< 0.05, 1, 0.2)),
                caption = "",
                legendLabels = c("NS", "NS", "NS", "5% FDR"),
                widthConnectors = 0.5,
                colConnectors = 'black',
                colCustom = keyvals, 
                #                col=c('grey60', 'grey60', 'grey60', color.palette0(2)[2]),
                colAlpha = 0.9,
                gridlines.minor = FALSE,
                gridlines.major = T,
                legendPosition = 'right')

# highlight AB and P-cell data 
keyvals = ifelse(out$gene_name %in% AB_cell_highlight$V1   , 'red',ifelse(out$gene_name %in% p_cell_highlight$V1,'blue','grey60'))
                 
names(keyvals)[keyvals == 'red'] <- 'AB'
names(keyvals)[keyvals == 'blue'] <- 'P'




EnhancedVolcano(out,
                xlim = c(-10, 7),
                ylim = c(0, 3.5),
                lab = out$gene_name,
                x = 'logFC',
                y = 'FDR',
                title = 'Translational efficiency',
                subtitle = NULL ,
                axisLabSize = 14,
                titleLabSize = 16,
                subtitleLabSize = 14,
                legendIconSize = 2.0,
                legendLabSize = 12,
                pCutoff = 0.05,
                FCcutoff = 0, 
                cutoffLineType = 'blank',
                vline = c(0),
                #                vlineCol = c('grey50', 'grey0','grey50'),
                labSize = 0,
                vlineType = 'blank',
                vlineWidth = 0.5,
                drawConnectors = FALSE,
                pointSize = c(ifelse(out$FDR< 0.05, 1, 0.2)),
                caption = "",
                legendLabels = c("NS", "NS", "NS", "5% FDR"),
                widthConnectors = 0.5,
                colConnectors = 'black',
                colCustom = keyvals, 
                #                col=c('grey60', 'grey60', 'grey60', color.palette0(2)[2]),
                colAlpha = 0.9,
                gridlines.minor = FALSE,
                gridlines.major = T,
                legendPosition = 'right')
  
              
keyvals = ifelse(out$gene_name %in% AB_cell_highlight$V1   , 'red',ifelse(out$gene_name %in% p_cell_highlight$V1,'blue','grey60'))

names(keyvals)[keyvals == 'red'] <- 'AB'
names(keyvals)[keyvals == 'blue'] <- 'P'



#Oma-1 and AB-cell 
somatic_oma_1 <- oma_highlight [oma_highlight$gene_id %in% AB_cell_highlight$V1,]
somatic_oma_1$gene_id <- c("F14H3.6", "C48B4.7" ,  "C44B9.3", "dhhc-3"  ,"ent-5","pigv-1","dtmk-1","C05C8.5" , "F40G12.11"  )
only_oma_1 <- oma_highlight [-(oma_highlight$gene_id %in% somatic_oma_1$gene_id),]
keyvals = ifelse(out$gene_name %in% somatic_oma_1$gene_id  , '#ffba39',ifelse(out$gene_name %in% only_oma_1$gene_id,'green','grey60'))

names(keyvals)[keyvals == 'green'] <- 'OMA-1 Bound RNA'
names(keyvals)[keyvals == '#ffba39' ] <- 'AB-cell transcripts'

selected_labels <- c("F14H3.6", "C48B4.7" ,  "C44B9.3", "dhhc-3"  ,"ent-5","pigv-1","dtmk-1","C05C8.5" , "F40G12.11"  )
EnhancedVolcano(out,
                xlim = c(-10, 7),
                ylim = c(0, 3.5),
                lab = out$gene_name,
                x = 'logFC',
                y = 'FDR',
                title = 'Translational efficiency',
                subtitle = NULL ,
                axisLabSize = 14,
                titleLabSize = 16,
                subtitleLabSize = 14,
                legendIconSize = 2.0,
                legendLabSize = 12,
                pCutoff = 0.05,
                FCcutoff = 0, 
                cutoffLineType = 'blank',
                vline = c(0),
                #                vlineCol = c('grey50', 'grey0','grey50'),
                labSize = 3.5,
                vlineType = 'blank',
                vlineWidth = 0.5,
                drawConnectors = TRUE,
                pointSize = c(ifelse(out$FDR< 0.05, 1, 0.2)),
                caption = "",
                legendLabels = c("NS", "NS", "NS", "5% FDR"),
                widthConnectors = 0.5,
                colConnectors = 'black',
                colCustom = keyvals, 
                #                col=c('grey60', 'grey60', 'grey60', color.palette0(2)[2]),
                selectLab = selected_labels,
                colAlpha = 0.9,
                gridlines.minor = FALSE,
                gridlines.major = T,
                legendPosition = 'right')


# barchart combining all the results 

# let us make a table with all FDR values: 
out <- as.data.table(out)


ab_cell <- out[out$gene_name %in% AB_cell_highlight$V1,c("logFC","FDR")] 
ab_cell [,regulation := ifelse(
  ab_cell$FDR > 0.05, 
                "NS", 
  ifelse(ab_cell$logFC > 0,
         "Up",
         "Down")
)]


p_cell <- out [out$gene_name %in% p_cell_highlight$V1,c("logFC","FDR")] 

ggplot() + 
  geom_boxplot(data = ab_cell,
               aes(
                 y = "AB-cell",
                 x = logFC,
                 fill = "#FF787D"
               )
  ) + 
  geom_boxplot(data = p_cell,
               aes(
                 y = "P-cell",
                 x = logFC,
                 fill = "#92ABFF"
               )
  )+
  scale_fill_manual(values = c( "#92ABFF","#FF787D"))+
theme_pubr()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20))

oma_1_all <-  out [out$gene_name %in% oma_highlight$gene_id, c("logFC", "FDR")]

oma_1_ab_cell <- out [out$gene_name %in% somatic_oma_1$gene_id, c("logFC", "FDR")]


ggplot() + 
  geom_boxplot(data = oma_1_all ,
               aes(
                 y = "All transcripts",
                 x = logFC,
                 fill = '#67FA5D'
               )
  ) + 
  geom_boxplot(data = oma_1_ab_cell,
               aes(
                 y = "AB-cell transcripts",
                 x = logFC,
                 fill = '#695736'
               )
  )+
  scale_fill_manual(values = c("green", "#695736"))+
  theme_pubr()+
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 20)) 



# doing GO analysis 

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

     