#Set the working directory to 20230301_RiboITP_RNAseq_analysis file 
library(ribor)
library(tidyverse)
library(data.table)
library(readxl)
original.ribo <- Ribo("/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/20230720_additional_replicates/celegans.ribo")

# Replotting length distrubution after using the two 3-embryo data as replicates. (Sanity Check)
plot_length_distribution(x           = original.ribo,
                         region      = "CDS",
                         range.lower = 20,
                         range.upper = 38,
                         fraction    = FALSE)
                         
get_info(original.ribo)$attributes$metagene_radius



plot_metagene(original.ribo,
              site        = "start",
              
              range.lower = 28,
              range.upper = 32
             )

plot_metagene(original.ribo,
              site        = "stop",
              normalize   = TRUE,
              title       = "Stop Site Coverage",
              range.lower = 28,
              range.upper = 32)


plot_region_counts(x           = original.ribo,
                   range.lower = 28,
                   range.upper = 32)

rcw <- get_region_counts(original.ribo, 
                   range.lower = 28,
                   range.upper = 32,
                   region = "CDS",
                   transcript = FALSE,
                   tidy = TRUE)
rcw <- as.data.table(rcw)

rcw_counts <- dcast(rcw,transcript ~ experiment )
colnames(rcw_counts) <- c("transcript","one_cell_c","two_cell_c","four_cell_c","four_cell_d","eight_cell_c","eight_cell_d")

rcw_counts[,sum(eight_cell_d > 0)]
