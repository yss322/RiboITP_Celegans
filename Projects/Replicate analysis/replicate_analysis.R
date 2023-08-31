
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
                                        id1, id2, 
                                        xlab    = "Replicate 2 (Counts)", 
                                        ylab    = "Replicate 1 (Counts)",
                                        main = "",
                                        num_bin = 52, 
                                        xrange  = 100000, 
                                        yrange  = 100000  ) { 
  
  sp = ggscatter(counts_w, x = id1, y = id2, title = main,
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

#import file 

original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo") 



ribo_rc <- get_region_counts(
                            original.ribo,
                             range.lower = 25,
                             range.upper = 35,
                             length      = TRUE,
                             transcript  = FALSE,
                             tidy = F ,
                             region      = c("CDS"), 
                             compact = F)
ribo_rc <-as.data.table(ribo_rc)

# adding pseudocount 
ribo_rc [,CDS := (CDS+1)]

rcw = dcast(ribo_rc, transcript ~ experiment) 
colnames(rcw)[2:9] <-c("one_cell_A","one_cell_B","two_cell_A","two_cell_B","four_cell_A","four_cell_B","eight_cell_A","eight_cell_B")
rcw <- as.data.table(rcw)
# filter transcripts less than 2 



  
plot(log2(one_cell_A), log2(rcw_filtered$one_cell_B), xlab = "Replicate 2 (counts)", ylab = "Replicate 1 (counts)", main = "1-cell replicate RNA correlation", pch = 19, cex = 0.5) + 
  stat_cor(method        = "spearman", 
           aes(label     = ..r.label..), 
           cor.coef.name = "rho", 
           digits        = 2)


plot_pairwise_relationships(rcw, id1 = "eight_cell_A", id2 = "eight_cell_", main =  "1-cell replicate RiboITP correlation")
plot_pairwise_relationships(rcw, id1 = "two_cell_A", id2 = "two_cell_B", main =  "2-cell replicate RiboITP correlation")
plot_pairwise_relationships(rcw_filtered, id1 = "four_cell_A", id2 = "four_cell_B", main =  "4-cell replicate RiboITP correlation")
plot_pairwise_relationships(rcw_filtered, id1 = "eight_cell_A", id2 = "eight_cell_B", main =  "8-cell replicate RiboITP correlation") 

# RNA correlation 

rna_rc <- get_rnaseq(original.ribo,
             tidy = F ,
             region      = c("CDS"),
             compact = F)
rna_rc <- as.data.table(rna_rc)
rna_rcw= dcast(rna_rc, transcript ~ experiment) 
colnames(rna_rcw)[2:9] <-c("one_cell_A","one_cell_B","two_cell_A","two_cell_B","four_cell_A","four_cell_B","eight_cell_A","eight_cell_B")
rna_rcw <- as.data.table(rna_rcw)
filtered_transcripts_rna <-rna_rcw[one_cell_A > 2 & one_cell_B >2]


plot_pairwise_relationships(filtered_transcripts_rna, id1 = "one_cell_A", id2 = "one_cell_B", main =  "1-cell replicate RNA-seq correlation")
plot_pairwise_relationships(filtered_transcripts_rna, id1 = "two_cell_A", id2 = "two_cell_B", main =  "2-cell replicate RNA-seq correlation")
plot_pairwise_relationships(filtered_transcripts_rna, id1 = "four_cell_A", id2 = "four_cell_B", main =  "4-cell replicate RNA-seq correlation")

plot_pairwise_relationships(filtered_transcripts_rna, id1 = "eight_cell_A", id2 = "eight_cell_B", main =  "8-cell replicate RNA-seq correlation")





