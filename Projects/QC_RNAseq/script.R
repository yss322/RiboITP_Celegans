# checking number of genes observed in RNAseq data 
# first import the file 

library(ribor)
library(tidyverse)

original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo")
original.ribo
# get the tidytable of the RNAseq data 
rna_seq <- get_rnaseq(ribo.object = original.ribo, 
                      tidy = TRUE,
                      region = c("CDS")
)
# filter and then summarise the number of genes with count >1 
rna_seq <- as_tibble(rna_seq)
rna_seq |> filter(count >2 ) |> 
  group_by(experiment) |> 
  summarise(n = n())

# conclusion : got enough number of genes detected. 

