library(data.table)
library(ggplot2)
library(biomaRt)
# Define a function to extract the first part
extract_first_part <- function(gene_id) {
  split_gene_id <- unlist(strsplit(gene_id, "/"))
  return(split_gene_id[1])
}
#import and convert the poly-a-tail count file 

raw_poly_a <- as.data.table (read.csv(file = "data/GSE68002_tail-count.csv"))

raw_embryo_poly_a<- raw_poly_a[,.(Name,X1.2.cell.egg.rep1,X1.2.cell.egg.rep2)]
raw_embryo_poly_a[,Name := sapply(Name, extract_first_part)]

raw_embryo_poly_avg <- raw_embryo_poly_a[, average := ((X1.2.cell.egg.rep1 + X1.2.cell.egg.rep2)/2) ]

geneIDs <- raw_embryo_poly_avg[,1]
#converting to wormbase annotations using biomart 
listMarts()
ensembl<-useMart("ensembl")
ensembl <- useDataset("celegans_gene_ensembl",mart=ensembl)

filters<-listFilters(ensembl) 
attributes <- listAttributes(ensembl)

results<- as.data.table(getBM(attributes = c("external_gene_name","refseq_mrna"),
      filters = "refseq_mrna",
      values = geneIDs,
      mart = ensembl))

results[,Name:=refseq_mrna]
view(results)
embyro_polya_wormbaseIDs <- merge(results,raw_embryo_poly_avg)
embyro_polya_wormbaseIDs <- embyro_polya_wormbaseIDs[,!c("Name")]
embyro_polya_wormbaseIDs[, gene_name := external_gene_name]
#filter out all the low values  

embyro_polya_wormbaseIDs  <- embyro_polya_wormbaseIDs [average > 5,]


#plot historgram to understand the distrubition of the data 

ggplot(embyro_polya_wormbaseIDs, 
       aes(
         x = log2(
           average)
         )) +
  geom_histogram(bins = 75)

#Let us get the fold change and ribosome occupancy values 

one_cell_ribo <- all_counts_diff[,c(1,6,7)]#checking with one cell stage first 
one_cell_ribo_avg <- one_cell_ribo[,average_ribo := ((one_cell_ribo[,2]+one_cell_ribo[,3])/2)]
one_cell_ribo_avg<- one_cell_ribo_avg[average_ribo >5,]

#merge the tables of polyatail and ribosome occupancy 

table <- inner_join(one_cell_ribo_avg,embyro_polya_wormbaseIDs)

ggplot(table, aes( x = log2(average) , y = log2(average_ribo)))  + 
  geom_jitter()
cor(table$average_ribo,table$average,method = "spearman")



# does the log fold change correlate with average   
  
