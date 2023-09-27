#in order to sift out the gene ontologies of things changing from 1cell to 2cell 
#I will be using the following guide https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO.pdf
if (!requireNamespace("BiocManager", quietly=TRUE))
  + install.packages("BiocManager")
BiocManager::install("topGO")
library(topGO)
library(ALL)
data(ALL)
data(geneList)

View(geneList)
affyLib <- paste(annotation(ALL), "db", sep = ".")

library(package = affyLib, character.only = TRUE)

sum(topDiffGenes(geneList))
?topDiffGenes()
sampleGOdata <- new("topGOdata",
                     description = "Simple session", ontology = "BP",
                     allGenes = geneList, geneSel = topDiffGenes,
                     nodeSize = 10,
                     annot = annFUN.db)
