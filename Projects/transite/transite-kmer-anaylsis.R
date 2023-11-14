library(Biostrings)
library(transite)
library(ggplot2)
library(data.table)
# Read in the sequences into R
fore_seqs<-as.data.frame( Biostrings::readRNAStringSet("/Users/yashshukla/Documents/Computational/Celegans_riboITP/Projects/transite/output/up_signal.foreground.fa", format="fasta", use.names=TRUE))
back_seqs<-as.data.frame( Biostrings::readRNAStringSet("/Users/yashshukla/Documents/Computational/Celegans_riboITP/Projects/transite/output/up_signal.background.fa", format="fasta", use.names=TRUE))

fore_seqs_input <-c(fore_seqs$x)
back_seqs_input <- c(back_seqs$x)

test_fore_seqs <- c("CAGUCAAGACUCC", "AAUUGGUGUCUGGAUACUUCCCUGUACAU",
                     "AGAU", "CCAGUAA")
test_back_seqs <- c(test_fore_seqs, "CAACAGCCUUAAUU", "CUUUGGGGAAU",
                     "UCAUUUUAUUAAA", "AUCAAAUUA", "GACACUUAAAGAUCCU",
                     "UAGCAUUAACUUAAUG", "AUGGA", "GAAGAGUGCUCA",
                     "AUAGAC", "AGUUC")
# Run transite


results <- transite::run_kmer_tsma(foreground_sets = fore_seqs_input, 
                        background_set= back_seqs_input , k=7, n_cores=2)
results_kmer<- calculate_kmer_enrichment(
  foreground_sets =fore_seqs_input,
  background_set=back_seqs_input,
  k = 7 ,
  permutation = FALSE,
  chisq_p_value_threshold = 0.05,
  p_adjust_method = "BH",
  n_cores = 4
)

enrichment_table_kmer <- (results[[1]][["enrichment_df"]])

ggplot(enrichment_table_kmer,aes(x = p_value)) + 
  geom_histogram()
 
ggplot(enrichment_table_kmer,aes(x = enrichment)) + 
  geom_histogram() +
  xlim(0,3)

compute_kmer_enrichment(
  foreground_kmers = list(as.character(fore_seqs)) ,
  background_kmers = generate_kmers(list(as.character(back_seqs)),7),
  permutation = FALSE,
  chisq_p_value_threshold = 0.05,
  p_adjust_method = "BH"
)

generate_kmers(list(as.character(fore_seqs), 7))



################################# Trail and error space 

test_sequence <- as.data.frame(fore_seqs)

rownames(test_sequence) <- NULL
               
View(test_sequence)

test_sequence$x

foreground_sets=list(as.character(fore_seqs))

test <-  c(test_sequence$x)
?as.data.frame()
foreground_sets=as.list(as.character(fore_seqs)
                     

                        