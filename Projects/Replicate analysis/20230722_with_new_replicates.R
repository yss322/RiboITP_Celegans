
library(ribor) 
library(tidyverse) 
library(ggpubr)
library(edgeR)
library(DGEobj)
library(data.table)

BURNT_ORANGE = "#bf5700"
UT_BLUE      = "#005f86"

MOUSE_MIN_LENGTH = 29
MOUSE_MAX_LENGTH = 35

FONT_LABEL_SIZE = 16
FONT_TITLE_SIZE = 20

PDF_resolution = 600
FIGURE_FONT    = "helvetica"


ribo_orange = rgb(228,88,10 , maxColorValue = 255)
rna_blue    = rgb(55,135,192, maxColorValue = 255)

# function for pairwise correlation 

plot_pairwise_relationships = function (counts_w, 
                                        RNA, RIBO, 
                                        xlab    = "RNA-seq (Counts)", 
                                        ylab    = "Ribosome-Profiling (Counts)",
                                        main = "",
                                        num_bin = 30, 
                                        xrange  = 100000, 
                                        yrange  = 100000  ) { 
  
  sp = ggscatter(counts_w, x = RNA, y = RIBO ,title = main,
                 #                add = "reg.line", conf.int = FALSE,     
                 #                add.params = list(color = "blue", size = 0.5),
                 font.family = "Helvetica", 
                 size        = 0.2,
                 color       = "gray", 
                 alpha       = 0.3, 
                 ggtheme     = theme_bw()) 
  
  formatted =   sp +   
    scale_x_log10(labels = scales::label_number_si(), limits = c(0.3, xrange)) +   
    scale_y_log10(labels = scales::label_number_si(), limits = c(0.3, yrange)) + 
    labs (x=xlab, y = ylab) +
    stat_cor(method        = "spearman", 
             aes(label     = ..r.label..), 
             cor.coef.name = "rho", 
             digits        = 2,
             size = 10,
    )  + 
    geom_hex(bins= num_bin, aes(alpha=log10(..count..) ), fill="#bf5700" ,
             show.legend = FALSE) +
    theme( axis.text.x      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           axis.title.x     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           axis.text.y      = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           axis.title.y     = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE),
           plot.title       = element_text(family = FIGURE_FONT, face = "plain", size = FONT_LABEL_SIZE)
    )
  return (formatted)  
}


#extracting data 

original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo") 
new.ribo <- Ribo("/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/20230720_additional_replicates/celegans.ribo")

#riboseq 
ribo_rc_or <- get_region_counts(
  original.ribo,
  range.lower = 25,
  range.upper = 35,
  length      = TRUE,
  transcript  = FALSE,
  tidy = F ,
  region      = c("CDS"), 
  compact = F)

ribo_rc_new <- get_region_counts(
  new.ribo,
  range.lower = 25,
  range.upper = 35,
  length      = TRUE,
  transcript  = FALSE,
  tidy = F ,
  region      = c("CDS"), 
  compact = F)
#adding pseudocount
ribo_rc_or <- as.data.table(ribo_rc_or)
ribo_rc_or [,CDS := (CDS+1)]

ribo_rc_new <- as.data.table(ribo_rc_new)
ribo_rc_new [,CDS := (CDS+1)]

#transform
rcw_or = dcast(ribo_rc_or, transcript ~ experiment) 
colnames(rcw_or)[2:9] <-c("one_cell_A","one_cell_B","two_cell_A","two_cell_B","four_cell_A","four_cell_B","eight_cell_A","eight_cell_B")
rcw_new = dcast(ribo_rc_new, transcript ~ experiment) 
colnames(rcw_new)[2:7] <-c("one_cell_C","two_cell_C","four_cell_C","four_cell_D","eight_cell_C","eight_cell_D")
rcw_ribo_combine <- merge(rcw_or,rcw_new, by = "transcript")
#RNAseq 
rna_rc <- get_rnaseq(original.ribo,
                     tidy = F ,
                     region      = c("CDS"), 
                     compact = F)
#adding pseudocount
rna_rc <- as.data.table(rna_rc)
rna_rc [,CDS := (CDS+1)]
rna_rcw= dcast(rna_rc, transcript ~ experiment) 
colnames(rna_rcw)[2:9] <-c("one_cell_A_RNA","one_cell_B_RNA","two_cell_A_RNA","two_cell_B_RNA","four_cell_A_RNA","four_cell_B_RNA","eight_cell_A_RNA","eight_cell_B_RNA")
# doing by raw counts 

rcw_combined <- as.data.table(merge(rcw_ribo_combine,rna_rcw, by = "transcript"))
exptype <- factor(c( "one_cell.ribo",
                     "one_cell.ribo",
                     "two_cell.ribo",
                     "two_cell.ribo",
                     "four_cell.ribo",
                     "four_cell.ribo",
                     "eight_cell.ribo",
                     "eight_cell.ribo",
                     "one_cell.ribo",
                     "two_cell.ribo",
                     "four_cell.ribo",
                     "four_cell.ribo",
                     "eight_cell.ribo",
                     "eight_cell.ribo",
                     "one_cell.RNA",
                     "one_cell.RNA",
                     "two_cell.RNA",
                     "two_cell.RNA",
                     "four_cell.RNA",
                     "four_cell.RNA",
                     "eight_cell.RNA",
                     "eight_cell.RNA")  )
y <- DGEList(counts=rcw_combined[,-1],
             genes = rcw_combined[,1],group = exptype)

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method = "TMM")
#Ribo to Ribocorrelation graphs 
y$counts[,1]

plot_pairwise_relationships(y$counts  ,
                            RNA =   y$counts[,1],
                            RIBO = y$counts[,1],
                            main =  "1-cell Ribo-RNA correlation")

plot_pairwise_relationships(rcw_combined  ,
                            id1 = "two_cell_A" ,
                            id2 = "two_cell_B",
                            main =  "2-cell vs 2-cell Ribo-RNA")
plot_pairwise_relationships(rcw_combined  ,
                            id1 = "four_cell_A" ,
                            id2 = "four_cell_B",
                            main =  "4-cell Ribo-RNA correlation")
plot_pairwise_relationships(rcw_combined  ,
                            id1 = "four_cell_C" ,
                            id2 = "four_cell_D",
                            main =  "4-cell Ribo-RNA correlation")
plot_pairwise_relationships(rcw_combined  ,
                            id1 = "eight_cell_A" ,
                            id2 = "eight_cell_D",
                            main =  "8-cell Ribo-RNA correlation")


#Ribo to RNA correlation graphs 


plot_pairwise_relationships(rcw_combined  ,
                            RNA =    "one_cell_A_RNA",
                            RIBO =  "one_cell_C",
                            main =  "1-cell Ribo-RNA correlation")

plot_pairwise_relationships(rcw_combined  ,
                            RNA =   "eight_cell_B_RNA" ,
                            RIBO = "eight_cell_D",
                            main =  "8-cell Ribo-RNA replicate D")
plot_pairwise_relationships(rcw_combined  ,
                            RNA =   "four_cell_B_RNA" ,
                            RIBO = "four_cell_B",
                            main =  "4-cell Ribo-RNA correlation")

plot_pairwise_relationships(rcw_combined  ,
                            RNA = "eight_cell_A_RNA" ,
                            RIBO  = "eight_cell_A",
                            main =  "8-cell Ribo-RNA correlation")








#it seems like for riboITP 1_cell_A and B is good
# 2-cell A,B,C is good 
# 4-cells A,B,C,D is good 
#8-cell A and D is good 
# all RNA-seq duplicate data is good 
