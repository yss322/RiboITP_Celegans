

run_transite<- function(go_term, species, te_dt, foreseqs, backseqs, expressed_genes){
  
  select_genes<- te_dt[term==go_term & Species==species, unlist(strsplit(gene_list_reduced, split=","))]
  
  # Currently foreseqs contains all genes across the GO terms
  foreseqs_select<- foreseqs[names(foreseqs)%in%select_genes]
  
  # Remove the background genes with no expression
  backseqs_select<- backseqs[names(backseqs)%in%expressed_genes]
  
  transite_res<- run_kmer_tsma(foreground_sets=list(as.character(foreseqs_select)), 
                               background_set=as.character(backseqs_select), k=7, n_cores=2, produce_plot=FALSE)
  
  res<- as.data.table(transite_res[[1]]$enrichment_df)
  return(res)
}
  
# Annotate the RBPs for the heptamers
transite_kmer_motifs_file<- "motifs.rda"

load(transite_kmer_motifs_file)
# RBPMotif classes in each element

transite_kmer_heptamers_list<- lapply(motifs, get_heptamers)

get_transite_motifs<- function(x){
  rbp_names<- get_rbps(x)
  rbp_collapse<- paste(rbp_names, collapse=",")
  return(rbp_collapse)
}

# Get the RBP names associated with k-mers
names(transite_kmer_heptamers_list)<- sapply(motifs, get_transite_motifs)

is_in<- function(e_set, x){
  return(x%in%e_set)
}

transite_kmer_to_rbps<- function(kmer, transite_kmers_list){
  match_logical<- sapply(transite_kmers_list, is_in, kmer)
  rbp_names<- names(transite_kmers_list[match_logical])
  rbp_names_split<- unlist(strsplit(rbp_names, ","))
  rbp_names_collapse<- paste(sort(unique(rbp_names_split)), collapse=",")
  return(rbp_names_collapse)
}


# Also annotate with oRNAment RBPs

oRNAment_RBP_pwms_dir<- "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/PWMs.tgz"
oRNAment_RBP_ids_file<- "http://rnabiology.ircm.qc.ca/BIF/oRNAment/static/Caenorhabditis_elegans_string_to_int_ID_conversion.csv.gz"
  
oRNAment_RBP_pwms_files<- list.files(oRNAment_RBP_pwms_dir, pattern="PWM", full.names=TRUE)
ornament_pwms_list<- lapply(oRNAment_RBP_pwms_files, fread, sep="\t")
ornament_pwms_all_dt<- rbindlist(ornament_pwms_list)
ornament_pwms_all_dt[,RBP:=rep(seq(1,453), each=7)]
rm(ornament_pwms_list)

oRNAment_RBP_ids_dt<- fread(oRNAment_RBP_ids_file, sep=",")
names(oRNAment_RBP_ids_dt)<- ornament_rbp_ids_header

# Annotate with RBP name
ornament_pwms_all_merge_dt<- merge(ornament_pwms_all_dt, oRNAment_RBP_ids_dt, by="RBP")
ornament_pwms_all_merge_dt[,RBP_short:=sapply(base::strsplit(RBP_name, " ", fixed=T), "[", 1)]

# Annotate the consensus motif for each RBP PWM
get_consensus_motif<- function(pwm){
  pwm_letters<- colnames(pwm)
  max_probs_indices<- apply(pwm, 1, which.max)
  consensus_seq<- NULL
  for(i in 1:nrow(pwm)){
    consensus_seq<- c(consensus_seq, pwm_letters[max_probs_indices[i]])
  }
  consensus_seq_res<- paste0(consensus_seq, collapse="")
  return(consensus_seq_res)
}

ornament_pwms_all_merge_consensus_dt<- ornament_pwms_all_merge_dt[
  ,get_consensus_motif(.SD), by=.(RBP_name), .SDcols=c("A", "C", "G", "U")]

names(ornament_pwms_all_merge_consensus_dt)<- c("RBP", "Consensus")

ornament_pwms_all_merge_consensus_dt[,RBP_short:=sapply(strsplit(RBP, " "), "[", 1)]

ornament_pwms_all_merge_consensus_key_dt<- ornament_pwms_all_merge_consensus_dt[
  ,.(N=.N, Key=paste0(RBP, collapse=",")), by=.(Consensus)]


# Allows match of kmer to oRNAment RBP database within a Hamming distance threshold
ornament_kmer_to_rbps<- function(kmer, ornament_pwms_all_merge_consensus_key_dt, hd=1){
  
  consensus_hd<- sapply(
    ornament_pwms_all_merge_consensus_key_dt[,strsplit(Consensus, "")], 
    e1071::hamming.distance, y=unlist(strsplit(kmer, "")))
  
  match_logical<- (consensus_hd <= hd)
  consensus_matches<- ornament_pwms_all_merge_consensus_key_dt[match_logical, Key]
  res<- paste(consensus_matches, collapse=";")
  
  if(res==""){
    return(NA)
  }
  
  return(res)
}

# Better yet score the heptamer against all possible RBP PWMs and select those RBPs above a certain MSS
# MSS defined by oRNAment db
compute_mss<- function(curr_score, min_score, max_score){
  max_min_diff<- max_score-min_score
  if(max_min_diff == 0){
    return(NaN)
  } else{
    res<- (curr_score - min_score) / max_min_diff
    return(res)
  }
}

# Returns MSS for a seq
compute_mss_from_probs<- function(pwm, seq){
  
  seq_split<- unlist(strsplit(seq, ""))
  pwm_letters<- colnames(pwm)
  pwm_matrix<- as.matrix(pwm)
  
  probs<- NULL
  max_scores<- NULL
  min_scores<- NULL
  
  for(i in 1:length(seq_split)){
    seq_i<- seq_split[i]
    
    # We need the product of the probs in the pwm for the kmer in questions
    j<- which(pwm_letters==seq_i)
    pwm_i_j<- pwm_matrix[i,j]
    probs<- c(probs, pwm_i_j)
    
    # For oRNAment MSS score, we normalize to the max and min probs
    max_prob_j<- which.max(pwm_matrix[i,])
    min_prob_j<- which.min(pwm_matrix[i,])
    max_scores<- c(max_scores, pwm_matrix[i, max_prob_j])
    min_scores<- c(min_scores, pwm_matrix[i, min_prob_j])
  }
  
  curr_score<- exp(sum(log(probs)))
  max_score<- exp(sum(log(max_scores)))
  min_score<- exp(sum(log(min_scores)))
  mss<- compute_mss(curr_score, min_score, max_score)
  
  return(mss)
}

# Redefine to use PWMs
ornament_kmer_to_rbps<- function(kmer, pwm_dt=ornament_pwms_all_merge_dt, mss_thresh=0.8){
  
  mss_dt<- pwm_dt[,.(MSS=compute_mss_from_probs(.SD, kmer)), 
                  by=.(RBP_name, RBP_short), .SDcols=c("A", "C", "G", "U")]
  
  rbps<- mss_dt[MSS>=mss_thresh, RBP_short]
  res<- paste(unique(rbps), collapse=",")
  return(res)
}

transite_all_res_filt_dt$Transite_RBP<- sapply(transite_all_res_filt_dt[,kmer], transite_kmer_to_rbps, transite_kmers_list=transite_kmer_heptamers_list)

# This may take some time as all heptamers are checked against hundreds of RBP PWMs
transite_all_res_filt_dt$oRNAment_RBP<- sapply(transite_all_res_filt_dt[,kmer], ornament_kmer_to_rbps)

combine_annots<- function(x, y){
  x_split<- unlist(strsplit(x, ","))
  y_split<- unlist(strsplit(y, ","))
  res<- unique(c(x_split, y_split))
  res<- paste(res, collapse=",")
  return(res)
}

transite_all_res_filt_dt[,Transite_oRNAment_RBP:=combine_annots(Transite_RBP, oRNAment_RBP), by=seq_len(nrow(transite_all_res_filt_dt))]
transite_all_res_filt_dt[,kmer_RBP:=paste(kmer, Transite_oRNAment_RBP, sep="; ")]

