#Set the working directory to 20230301_RiboITP_RNAseq_analysis file 
library(ribor)
library(tidyverse)
library(data.table)
library(readxl)
original.ribo <- Ribo("/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/celegans_duplicate_allstages.ribo")

# Replotting length distrubution after using the two 3-embryo data as replicates. (Sanity Check)
plot_length_distribution(x           = original.ribo,
                         region      = "CDS",
                         experiment = c("1cell_A","1cell_B","2cell_B","4cell_A","4cell_B","8cell_A"),
                         range.lower = 20,
                         range.upper = 40,
                         fraction    = TRUE)

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
                   experiment = c("1cell_A","1cell_B","2cell_B","4cell_A","4cell_B","8cell_A"),
                   range.lower = 28,
                   range.upper = 32)

rcw <- get_region_counts(original.ribo, 
                         range.lower = 28,
                         range.upper = 32,
                         region = "CDS",
                         experiment = c("1cell_A","1cell_B","2cell_B","4cell_A","4cell_B","8cell_A"),
                         transcript = FALSE,
                         tidy = TRUE)
rcw <- as.data.table(rcw)

rcw_counts <- dcast(rcw,transcript ~ experiment )
colnames(rcw_counts) <- c("transcript","one_cell_c","two_cell_c","four_cell_c","four_cell_d","eight_cell_c","eight_cell_d")

rcw_counts[,sum(rcw_counts$'1cell_A' > 0)]
rcw_counts[,sum(rcw_counts$'2cell_A' > 0)]
rcw_counts[,sum(rcw_counts$'4cell_A' > 0)]
rcw_counts[,sum(rcw_counts$'8cell_A' > 0)]
rcw_counts[,sum(rcw_counts$'1cell_B' > 0)]
rcw_counts[,sum(rcw_counts$'2cell_B' > 0)]
rcw_counts[,sum(rcw_counts$'4cell_B' > 0)]
rcw_counts[,sum(rcw_counts$'8cell_B' > 0)]
