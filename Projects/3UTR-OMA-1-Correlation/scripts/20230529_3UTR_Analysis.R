#Here I'll be writing a script for analysing the 3'UTR of 1_cell stage 
#see if it matches that of OMA-1 pulldown 
library(ribor)
library(tidyverse)
library(readxl)



# import the ribo file 
original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo")
original.ribo

# obtaining 3UTR data from ribo file at the 1-cells stage and fitering for top100 genes 

rc_UTR_1cell_A <- get_region_counts(ribo.object    = original.ribo,
                            range.lower = 25,
                            range.upper = 35,
                            tidy       = TRUE,
                            transcript = FALSE,
                            region     = "UTR3",
                            compact    = FALSE,
                            experiment = "1cell_A" )
rc_UTR_1cell_A_top_100 <- rc_UTR_1cell_A |> 
  arrange(desc(count)) |>  
  filter (count > 119)
# cleaning the data to extract the wbgeneIDs for each: 

wbgene_rc_UTR_1cell_A_top_100 <- rc_UTR_1cell_A_top_100 |> 
  separate(transcript, sep = "gene:",  into = c('a','b'), extra = "merge") |> 
  separate(b,sep = ".1\\|gene_biotype:", extra = "merge", into = c ('Wbgene_id','d')) |> 
  select(Wbgene_id)
  

# import the excel file of oma-1 pull down 
oma_1_pull_down <- read_excel("data/FPKM_OMA-1_Pull_down_wbgene.xlsx")

# Let us first clean up the oma_1_pull_down data 
oma_1_pull_down_top_200 <- oma_1_pull_down |>  
  filter(FPKM > 480) |>  
  rename(Wbgene_id = converted_alias ) |>  
  select(Wbgene_id)

# finding number of overlap
overlapping_genes <- wbgene_rc_UTR_1cell_A_top_100 |> 
  semi_join(oma_1_pull_down_top_200 )
number_overlapping <- nrow(overlapping_genes)

fmat = matrix (nrow =2, ncol = 2)
fmat[1,] = c(number_overlapping, nrow(wbgene_rc_UTR_1cell_A_top_100 ) - number_overlapping)
fmat[2,] = c(nrow(oma_1_pull_down_top_200) - number_overlapping, 24244- nrow(oma_1_pull_down_top_200) - fmat[1,2])
statistical_test_results = fisher.test(fmat)
statistical_test_results
print(statistical_test_results$p.value)
