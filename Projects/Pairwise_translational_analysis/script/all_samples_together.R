#combining all comparisions 
library(ribor)
library(tidyverse)
library(edgeR) 
library(data.table)
library(ggpubr)
library(stack)
library(EnhancedVolcano)
library(biomaRt)

subsampleMatrix <- function(counts, desiredDepth) {
  rns <- rownames(counts)
  n <- nrow(counts)
  proportion <- desiredDepth / colSums(counts)
  
  if (sum(proportion > 1)) {
    stop("DesiredDepth must be less than library size for all!")
  }
  
  ret <- counts
  for (experiment in 1:ncol(counts)) {
    ret[, experiment] <- rbinom(n, counts[, experiment], proportion[experiment])
  }
  
  rownames(ret) <- rns
  return(ret)
}


original.ribo <- Ribo("/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/celegans_duplicate_allstages.ribo")
reseq.ribo <- Ribo("/Users/yashshukla/Documents/Computational/Celegans_riboITP/data/20230929_WT_ITP_reseq1/output_umi/ribo/all.ribo")

# we are treating two-cell_A replicate seperately because the peaks are biased towards
#higher sizes
ribo_original <- get_region_counts(original.ribo,
                    range.lower = 25,
                    range.upper = 35,
                    experiment = c('1cell_A','1cell_B','2cell_B','4cell_A','4cell_B','8cell_A'),
                    length      = TRUE,
                    transcript  = FALSE,
                    tidy = F,
                    region      = c("CDS") )

ribo_cell_new <-  get_region_counts(reseq.ribo,
                                    range.lower = 25,
                                    range.upper = 35,
                                    experiment = c('2-cell_C','4-cell_C','4-cell_D','8-cell_D'),
                                    length      = TRUE,
                                    transcript  = FALSE,
                                    tidy = F,
                                    region      = c("CDS") )
ribo_two_cell_A <-  get_region_counts(original.ribo,
                                      range.lower = 28,
                                      range.upper = 38,
                                      experiment = c('2cell_A'),
                                      length      = TRUE,
                                      transcript  = FALSE,
                                      tidy = F,
                                      region      = c("CDS") )
ribo_original<- as.data.table(ribo_original)
ribo_cell_new <-as.data.table(ribo_cell_new)
ribo_two_cell_A <-as.data.table(ribo_two_cell_A)

ribo_original[,CDS := (CDS+1)]
ribo_cell_new[,CDS := (CDS+1)]
ribo_two_cell_A[,CDS := (CDS+1)]

rcw_diff_or = dcast(ribo_original, transcript ~ experiment) 
rcw_diff_new =  dcast(ribo_cell_new, transcript ~ experiment) 
rcw_diff_two_cell_A = dcast(ribo_two_cell_A, transcript ~ experiment) 

rcw_diff = inner_join(rcw_diff_or,rcw_diff_new, by = "transcript")
rcw_diff_ribo = inner_join(rcw_diff,rcw_diff_two_cell_A , by = "transcript")

setcolorder(rcw_diff_ribo,c( "transcript", 
                               "1cell_A"  , 
                               "1cell_B", 
                               "2cell_A",   
                              "2cell_B"  ,
                              "2-cell_C" ,
                              "4cell_A" ,
                              "4cell_B", 
                              "4-cell_C" ,
                              "4-cell_D" ,
                             "8cell_A"  ,     
                               "8-cell_D" ))
colnames(rcw_diff_ribo) <- c( "transcript", 
                              "1cell_A.RIBO"  , 
                              "1cell_B.RIBO", 
                              "2cell_A.RIBO",   
                              "2cell_B.RIBO"  ,
                              "2cell_C.RIBO" ,
                              "4cell_A.RIBO" ,
                              "4cell_B.RIBO", 
                              "4cell_C.RIBO" ,
                              "4cell_D.RIBO" ,
                              "8cell_A.RIBO"  ,     
                              "8cell_D.RIBO" )
rnaseq_diff <- get_rnaseq(ribo.object = original.ribo,
                          tidy        = F,
                          compact = F,
                          experiment = c('1cell_A','1cell_B','2cell_A','2cell_B','4cell_A','4cell_B','8cell_A','8cell_B'),
                          region = "CDS")
rnaseq_diff <- as.data.table(rnaseq_diff)
rnaseq_w_diff = dcast(rnaseq_diff, transcript ~ experiment) 


all_counts_diff = merge(rnaseq_w_diff, rcw_diff_ribo, by= "transcript")

colnames(all_counts_diff)<-c("transcript", 
                             "1cell_A.RNA"  , 
                             "1cell_B.RNA", 
                             "2cell_A.RNA",   
                             "2cell_B.RNA"  ,
                             "4cell_A.RNA" ,
                             "4cell_B.RNA", 
                             "8cell_A.RNA"  ,     
                             "8cell_B.RNA", 
                             "1cell_A.RIBO"  , 
                             "1cell_B.RIBO", 
                             "2cell_A.RIBO",   
                             "2cell_B.RIBO"  ,
                             "2cell_C.RIBO" ,
                             "4cell_A.RIBO" ,
                             "4cell_B.RIBO", 
                             "4cell_C.RIBO" ,
                             "4cell_D.RIBO" ,
                             "8cell_A.RIBO"  ,     
                             "8cell_D.RIBO")


exptype <- factor(c("1cell.RNA"  , 
                    "1cell.RNA", 
                    "2cell.RNA",   
                    "2cell.RNA"  ,
                    "4cell.RNA" ,
                    "4cell.RNA", 
                    "8cell.RNA"  ,     
                    "8cell.RNA", 
                    "1cell.RIBO"  , 
                    "1cell.RIBO", 
                    "2cell.RIBO",   
                    "2cell.RIBO"  ,
                    "2cell.RIBO" ,
                    "4cell.RIBO" ,
                    "4cell.RIBO", 
                    "4cell.RIBO" ,
                    "4cell.RIBO" ,
                    "8cell.RIBO"  ,     
                    "8cell.RIBO")  )



y <- DGEList(counts=(all_counts_diff[,-1]),
             group=exptype, genes = all_counts_diff[,1])


keep <- filterByExpr(y,min.total.count = 15 )
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~0+exptype)
colnames(design) <- levels(exptype)
y <- estimateDisp(y,design)
plotBCV(y)
plotMDS(y[,(9:19)])


######trail and error###########

all_counts_diff_matrix<- as.data.frame(all_counts_diff[,-(1:9)])
subsample_all_count_diff <- subsampleMatrix(all_counts_diff_matrix,71000)
summarized_df <- subsample_all_count_diff |> 
  summarise_all(sum)

exptype <- factor(c(  "1cell.RIBO"  , 
                    "1cell.RIBO", 
                    "2cell.RIBO",   
                    "2cell.RIBO"  ,
                    "2cell.RIBO" ,
                    "4cell.RIBO" ,
                    "4cell.RIBO", 
                    "4cell.RIBO" ,
                    "4cell.RIBO" ,
                    "8cell.RIBO"  ,     
                    "8cell.RIBO")  )
y <- DGEList(counts=(subsample_all_count_diff ),
             group=exptype, genes = all_counts_diff[,1])


keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y, method = "TMM")
design <- model.matrix(~0+exptype)
colnames(design) <- levels(exptype)
y <- estimateDisp(y,design)
plotBCV(y)
plotMDS(y)

UTR3_counts <- get_region_counts(original.ribo,
                  range.lower = 25,
                  range.upper = 35,
                  length      = TRUE,
                  transcript  = FALSE,
                  tidy = F, 
                  experiment = c('1cell_A','1cell_B','2cell_A','2cell_B','4cell_A','4cell_B','8cell_A'),
                  region      = c("UTR3"))
UTR3_counts <-as.data.table(UTR3_counts)
UTR3_counts_rcw= dcast(UTR3_counts, transcript ~ experiment) 
plot_length_distribution(original.ribo,
                         range.lower = 25,
                         range.upper = 38,
                         
                         experiment = c('1cell_A','1cell_B','2cell_A','2cell_B','4cell_A','4cell_B','8cell_A'),
                         region      = c("UTR3"))

plot_region_counts(original.ribo,
                   range.lower = 25,
                   range.upper = 35,
                   
                   experiment = c('1cell_A','1cell_B','2cell_A','2cell_B','4cell_A','4cell_B','8cell_A')
)
ribo_cell_new <-  get_region_counts(reseq.ribo,
                                    range.lower = 25,
                                    range.upper = 35,
                                    experiment = c('2-cell_C','4-cell_C','4-cell_D','8-cell_D'),
                                    length      = TRUE,
                                    transcript  = FALSE,
                                    tidy = F,
                                    region      = c("CDS") )


plot_length_distribution(reseq.ribo,
                         range.lower = 25,
                         range.upper = 35,
                         experiment = c('2-cell_C','4-cell_C','4-cell_D','8-cell_D'),
                         region      = c("CDS"))

plot_region_counts(reseq.ribo,
                   range.lower = 25,
                   range.upper = 35,
                   experiment = c('2-cell_C','4-cell_C','4-cell_D','8-cell_D')
)
