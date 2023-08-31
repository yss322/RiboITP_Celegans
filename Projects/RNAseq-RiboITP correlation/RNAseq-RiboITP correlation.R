
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
                                        xlab    = "RNA-seq (Counts)", 
                                        ylab    = "Ribosome-Profiling (Counts)",
                                        main = "",
                                        num_bin = 30, 
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


#extracting data 

original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo") 


#riboseq 
ribo_rc <- get_region_counts(
  original.ribo,
  range.lower = 25,
  range.upper = 37,
  length      = TRUE,
  transcript  = FALSE,
  tidy = F ,
  region      = c("CDS"), 
  compact = F)

#adding pseudocount
ribo_rc <- as.data.table(ribo_rc)
ribo_rc [,CDS := (CDS+1)]

rcw = dcast(ribo_rc, transcript ~ experiment) 
colnames(rcw)[2:9] <-c("one_cell_A","one_cell_B","two_cell_A","two_cell_B","four_cell_A","four_cell_B","eight_cell_A","eight_cell_B")

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

rcw_combined <- as.data.table(merge(rcw,rna_rcw, by = "transcript"))
rcw_filtered    <- rcw_combined[one_cell_A_RNA > 2 & one_cell_B_RNA > 2]

      


#Ribo to RNA correlation graphs 


plot_pairwise_relationships(rcw_combined  ,
                             id1 =  "one_cell_A_RNA",
                             id2 = "one_cell_A",
                            main =  "1-cell Ribo-RNA correlation")

plot_pairwise_relationships(rcw_combined  ,
                             id1 = "two_cell_A_RNA" ,
                             id2 = "two_cell_B",
                             main =  "2-cell vs 2-cell Ribo-RNA")
plot_pairwise_relationships(rcw_combined  ,
                            id1 = "four_cell_A_RNA" ,
                            id2 = "four_cell_A",
                            main =  "4-cell Ribo-RNA correlation")

plot_pairwise_relationships(rcw_combined  ,
                            id1 = "eight_cell_A_RNA" ,
                            id2 = "eight_cell_A",
                            main =  "8-cell Ribo-RNA correlation")
plot_pairwise_relationships(rcw_combined  ,
                            id1 = "four_cell_A_RNA" ,
                            id2 = "two_cell_A",
                            main =  "4-cell Ribo- 2-cell RNA correlation")
plot_pairwise_relationships(rcw_combined  ,
                            id1 = "two_cell_A" ,
                            id2 = "two_cell_A_RNA",
                            main =  "4-cell Ribo- 2-cell RNA correlation")











#combine information - if you want to use CPM


rcw_CPM <- cbind(
  rcw[,1],RawCPM
)

rcw_CPM[rcw_CPM==0] <- NA
rcw_CPM_filtered<-rcw_CPM[complete.cases(rcw_CPM),] 


RawCPM <- cpm(rcw_CPM_filtered[2:17],
              unit      = "CPM",
              log       = FALSE,
              normalize = "none") 


rcw_CPM_filtered <- rcw_CPM_filtered[, (2:17) := lapply(.SD, as.numeric), .SDcols = 2:17]


