#Here I'll be writing a script for analysing the 3'UTR of 1_cell stage 
#see if it matches that of OMA-1 pulldown 
library(ribor)
library(tidyverse)
library(readxl)
library("biomaRt")


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
  

g# import the excel file of oma-1 pull down 
oma_1_pull_down <- read_excel("data/FPKM_OMA-1_Pull_down.xlsx")

# Let us first clean up the oma_1_pull_down data 
oma_1_pull_down_top_200 <- oma_1_pull_down |>  
  select(tracking_id,gene_id,FPKM) |> 
  filter(FPKM > 480)



