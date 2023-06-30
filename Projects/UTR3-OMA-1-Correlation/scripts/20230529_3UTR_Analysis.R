#Here I'll be writing a script for analysing the 3'UTR of 1_cell stage 
#see if it matches that of OMA-1 pulldown 
library(ribor)
library(tidyverse)
library(readxl)
library(biomaRt)
library(conflicted)
library(tibble)

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

rc_CDS_UTR3_1cell_A <- get_region_counts(ribo.object    = original.ribo,
                                                     range.lower = 25,
                                                     range.upper = 35,
                                                     tidy       = TRUE,
                                                     transcript = FALSE,
                                                     region     = c("CDS","UTR3"),
                                                     compact    = FALSE,
                                                     experiment = "1cell_A" )


rc_UTR_1cell_A_top_100 <- rc_UTR_1cell_A |> 
  arrange(desc(count)) |>  
  filter (count > 119)
# cleaning the data to extract the wbgeneIDs for each: 

#let us make a cleanup function first 
cleanup <- function(x) {
 x |> 
    separate(transcript, sep = "gene:",  into = c('a','b'), extra = "merge") |> 
    separate(b,sep = ".1\\|gene_biotype:", extra = "merge", into = c ('Wbgene_id','d')) |> 
    dplyr::select(Wbgene_id)
}

wbgene_rc_UTR_1cell_A_top_100  <- cleanup(rc_UTR_1cell_A_top_100)

  

# import the excel file of oma-1 pull down 
oma_1_pull_down_toplist <- read_excel("data/FPKM_OMA-1_Pull_down_wbgene.xlsx")
oma_1_pull_down_complete_list <- read_excel("data/FPKM_OMA-1_Pull_down.xlsx")

# Let us first clean up the oma_1_pull_down data 
oma_1_pull_down_top_200 <- oma_1_pull_down_toplist |>  
  filter(FPKM > 480) |>  
  rename(Wbgene_id = converted_alias ) |>  
  select(Wbgene_id)

# finding number of overlap
overlapping_genes <- wbgene_rc_UTR_1cell_A_top_100 |> 
  semi_join(oma_1_pull_down_top_200 )
number_overlapping <- nrow(overlapping_genes) 
overlapping_genes
# We would like to obtain the superset of all the genes detected by riboITP and the OMA-1 pull down
 
count_greater_than_2_riboITP <- rc_CDS_UTR3_1cell_A |>  
  pivot_wider(
    names_from = region,
    values_from =count
  ) |> 
  mutate(Sum_CDS_UTR3 = CDS + UTR3) |>  
  filter(Sum_CDS_UTR3 > 1 ) 
  
  Wbgene_count_greater_than_2_riboITP <-cleanup(count_greater_than_2_riboITP)
  
  #cleaning up the oma-1 list and making sure all the gene names are in wbgene format
  
  celegans_mart <- useMart(host = "https://metazoa.ensembl.org", 
                           biomart = "metazoa_mart",
                           dataset = "celegans_eg_gene")




refseq_oma_1 <- oma_1_pull_down_complete_list |> 
  filter(FPKM > 1) |> 
  select (tracking_id)
ref_seq_oma_1_list <- pull(refseq_oma_1,tracking_id)


celegans_mart <- useMart(host = "https://metazoa.ensembl.org", 
                         biomart = "metazoa_mart",
                         dataset = "celegans_eg_gene")

 all_wbgene_oma_1 <- getBM(mart = celegans_mart,
        filters = "refseq_mrna",
        values = ref_seq_oma_1_list,
        attributes = c("refseq_mrna", "ensembl_gene_id")) |> 
   select(ensembl_gene_id) |> 
   rename(Wbgene_id = ensembl_gene_id) 
  
# finally we will get what is common between everything detected by OMA-1 and RiboITP
 overlapping_genes_all <-  all_wbgene_oma_1 |> 
   semi_join(  Wbgene_count_greater_than_2_riboITP )
 number_overlapping_all <- nrow(overlapping_genes_all) 
 number_overlapping_all
#Running a fisher's test 

fmat = matrix (nrow =2, ncol = 2)
fmat[1,] = c(number_overlapping, nrow(wbgene_rc_UTR_1cell_A_top_100 ) - number_overlapping)
fmat[2,] = c(nrow(oma_1_pull_down_top_200) - number_overlapping, number_overlapping_all - nrow(oma_1_pull_down_top_200) - fmat[1,2])
statistical_test_results = fisher.test(fmat)
statistical_test_results
print(statistical_test_results$p.value)
