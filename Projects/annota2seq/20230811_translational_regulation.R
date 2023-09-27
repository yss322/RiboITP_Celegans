library(anota2seq)
library(tidyverse)
library(gridExtra)
library(ribor)
library(limma) 
library(edgeR) 
library(data.table)

original.ribo <- Ribo("data/celegans_duplicate_allstages.ribo")
ribo_one_two_cell<- get_region_counts(original.ribo,
                                      range.lower = 25,
                                      range.upper = 35,
                                      experiment = c('1cell_A','1cell_B','2cell_B'),
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
ribo_one_two_cell<- as.data.table(ribo_one_two_cell)
ribo_two_cell_A <-as.data.table(ribo_two_cell_A)

ribo_one_two_cell[,CDS := CDS+1]
ribo_two_cell_A [,CDS := CDS+1]

rcw_diff_one_two = dcast(ribo_one_two_cell, transcript ~ experiment)  

rcw_diff_two_A = dcast(ribo_two_cell_A, transcript ~ experiment)  

rcw_diff_ribo = inner_join(rcw_diff_one_two,rcw_diff_two_A, by = "transcript")

simulated_rep <- rcw_diff_ribo
colnames(simulated_rep) <- c("transcript","onecell_C","onecell_D","twocell_C","twocell_D")
rcw_diff_ribo_sim <- inner_join(rcw_diff_ribo,simulated_rep, by = "transcript")
#processing RNA-seq 
rnaseq_diff <- get_rnaseq(ribo.object = original.ribo,
                          tidy        = F,
                          compact = F,
                          experiment = c('1cell_A','1cell_B','2cell_A','2cell_B'),
                          region = "CDS")

rnaseq_diff <- as.data.table(rnaseq_diff)

rnaseq_w_diff = dcast(rnaseq_diff, transcript ~ experiment) 
sim_rep_RNA <- rnaseq_w_diff

colnames(sim_rep_RNA) <- c("transcript","onecell_C_RNA","onecell_D_RNA","twocell_C_RNA","twocell_D_RNA")
merge_sim_rep_RNA<- inner_join(rnaseq_w_diff,sim_rep_RNA,by = "transcript")
all_counts_diff = merge(merge_sim_rep_RNA, rcw_diff_ribo_sim , by= "transcript")






anot <- data.frame(sample= c("1cell_A.x", 
                             "1cell_B.x",
                             "2cell_A.x", 
                             "2cell_B.x",
                             "onecell_C_RNA",
                             "onecell_D_RNA",
                            "twocell_C_RNA",
                            "twocell_D_RNA",
                             "1cell_A.y",
                             "1cell_B.y", 
                             "2cell_A.y", 
                             "2cell_B.y",
                            "onecell_C",
                            "onecell_D",
                            "twocell_C",
                            "twocell_D"
                            ),
                  RNA = c("cyto","cyto","cyto","cyto", "cyto","cyto","cyto","cyto","poly","poly","poly","poly","poly","poly","poly","poly"), 
                  condition = c("1cell","1cell","2cell","2cell","1cell","1cell", "2cell","2cell","1cell","1cell","2cell","2cell","1cell", "1cell","2cell","2cell"),
                  replicate=c(1,2,1,2,3,4,3,4,1,2,1,2,3,4,3,4))
rownames(anot) <- anot$sample

# basic QC using PCA analysis 
#the +1 that I did for ribosome profiling will not work here. Rememnber to r
#remove zero value

tmpCounts <- as.data.frame(all_counts_diff)
rownames(tmpCounts) <- tmpCounts[,1]
tmpCounts <- tmpCounts[,-1]
tmpPCA <- as.matrix(tmpCounts)

#remove non-zero values 
tmpPCANoZero <- tmpPCA[!apply(tmpPCA,1,FUN=function(x) any(x == 0)),]

#normalize the values 
tmpPCANorm <- voom(calcNormFactors(DGEList(tmpPCANoZero)))$E

#4. Calculate standard deviations and remove mRNAs with an SD in the first three quartiles

tmpSD <- apply(tmpPCANorm,1,sd) 
sdQuantiles <- quantile(tmpSD) 
tmpPCASDFilt <- tmpPCANorm[tmpSD > sdQuantiles[4],]

# Transpose data
tmpPCATrans <- t(tmpPCASDFilt)

# Perform the PCA 
pcaOut <- prcomp(tmpPCATrans)

# Setup a plotting window
par(mfrow=c(2,2))
# generate a barplot of the proportion of variance
barplot(summary(pcaOut)$importance[2,],ylab = "Proportion of variance")
abline(h=0.05,col="red")
# generate a barplot of the cumulative proportion of variance
barplot(summary(pcaOut)$importance[3,],ylab = "Cumulative proportion ofvariance")
abline(h=0.7,col="red")

# merge the PCA output with the annotation file 
pcaPlot <- merge(pcaOut$x,anot,by="row.names")

# plot the principal components 
grid.arrange( 
  
  ggplot(data=pcaPlot,aes(x=PC1,y=PC2,shape=RNA,col=condition))+ 
    geom_point()+ 
    geom_text(aes(label=replicate,vjust=1,hjust=1)),
  
  ggplot(data=pcaPlot,aes(x=PC1,y=PC3,shape=RNA,col=condition))+ 
    geom_point()+ 
    geom_text(aes(label=replicate,vjust=1,hjust=1)),
  
  ggplot(data=pcaPlot,aes(x=PC1,y=PC4,shape=RNA,col=condition))+ 
    geom_point()+ 
    geom_text(aes(label=replicate,vjust=1,hjust=1)), 
  nrow=2,ncol=2 
) 

# extract dataP and dataT from tmpCounts 
tmpDataP <- tmpCounts[,anot[anot$RNA == "poly","sample"]] 
tmpDataT <- tmpCounts[,anot[anot$RNA == "cyto","sample"]]

# Check if you extracted the correct samples 
colnames(tmpDataP)
colnames(tmpDataT)
dim(tmpDataP)
dim(tmpDataT)


#generate phenoVec
phenoVec <- anot[anot$RNA == "poly","condition"]

#initialize anota2seq object needs 3 replicates 
ads <- anota2seqDataSetFromMatrix(dataP=tmpDataP, dataT=tmpDataT, phenoVec = phenoVec, dataType = "RNAseq", normalize = T)


#annota2seq QC 
ads <- anota2seqPerformQC(ads)

ads <- anota2seqResidOutlierTest(ads)

myContrast <- matrix(nrow=length(levels(as.factor(phenoVec))), ncol=length(levels(as.factor(phenoVec)))-1) 
rownames(myContrast) <- levels(as.factor(phenoVec)) 
myContrast[,1] <- c(-1,1)

ads <- anota2seqAnalyze(ads, contrasts=myContrast)
anota2seqPlotPvalues(ads, selContrast = 1, plotToFile = FALSE) 


# Apply significance filters 
ads <- anota2seqSelSigGenes(ads, maxPAdj = 0.01, 
                            selDeltaPT = log2(2),
                            selDeltaTP = log2(2) )

anota2seqPlotGenes(ads,
                   selContrast = 1,
                   analysis = 'translation',
                   plotToFile = F)

ads <-anota2seqRegModes(ads)
anota2seqPlotFC(ads, selContrast = 1, plotToFile = FALSE)

#during the run it seems like the zero counts are removed , which something I don't what
#one strategy is to add +1 or remove the ones which have zero value in the rnaseq. Unsure 
#what is the best way to go about analyzing the data. 
anota2seqPerformQC(Anota2seqDataSet = ads,
                   generateSingleGenePlots = TRUE)

