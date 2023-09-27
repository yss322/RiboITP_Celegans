library(ribor)
library(data.table)
library(ggplot2) 

# We are importing data over here
original.ribo <- Ribo('/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/celegans_duplicate_allstages.ribo')

new.ribo  <- Ribo('/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/20230720_additional_replicates/celegans.ribo') 


ribo_cell_or<- get_region_counts(original.ribo,
                                 range.lower = 25,
                                 range.upper = 35,
                                 experiment = c('1cell_A', '1cell_B','2cell_A', '2cell_B','4cell_A', '4cell_B','8cell_A', '8cell_B'),
                                 length      = TRUE,
                                 transcript  = FALSE,
                                 tidy = F,
                                 region      = c("CDS") )
ribo_cell_new <-get_region_counts(new.ribo,
                                  range.lower = 25,
                                  range.upper = 35,
                                  experiment = c('1cell_C', '2cell_C', '4cell_C','4cell_D', '8cell_C', '8cell_D'),
                                  length      = TRUE,
                                  transcript  = FALSE,
                                  tidy = F,
                                  region      = c("CDS") )

# let us make into data.table because it is easier to handle 

ribo_cell_or <- as.data.table(ribo_cell_or)
ribo_cell_new <- as.data.table(ribo_cell_new)
ribo_merged <- rbind(ribo_cell_or,ribo_cell_new)
ribo_cell_display<- ribo_merged[(experiment == "8cell_A" | experiment == "8cell_B" | experiment == "8cell_C" | experiment == "8cell_D")]
ggplot(ribo_cell_display, aes(x = log2(CDS), color = experiment)) + 
  geom_density(linewidth = 0.75) 



