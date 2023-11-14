library(ribor)
library(tidyverse)
library(data.table)

original.ribo <- Ribo("/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/20230929_WT_ITP_reseq1/output_umi/ribo/all.ribo")


plot_length_distribution(x           = original.ribo,
                         region      = "CDS",
                         range.lower = 22,
                         range.upper = 40,
                         fraction    = TRUE)





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

# getting region counts and distrubution 


ribo_cell_or<- get_region_counts(original.ribo,
                                 range.lower = 25,
                                 range.upper = 35,
                                 transcript  = FALSE,
                                 tidy = F,
                                 region      = c("CDS") )



ribo_cell_or <- as.data.table(ribo_cell_or)

ggplot(ribo_cell_or, aes(x = log2(CDS), color = experiment)) + 
  geom_density(linewidth = 0.75) 
