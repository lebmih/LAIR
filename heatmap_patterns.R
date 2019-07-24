### PATTERN FINDER
#############  LIBRARIES  ##################
## Libraries
x <- c('tidyr', 'dplyr', 'ggplot2', 'ggrepel', 'reshape2', 'stringr', 'Biostrings', 'tcR')
lapply(x, require, character.only = TRUE)
rm(x)
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("BSgenome.Hsapiens.NCBI.GRCh38")

library(BSgenome.Hsapiens.NCBI.GRCh38)

#############  FUNCTIONS  ##################
## Functions
# randseq() gives a string of specified length with random DNA sequence
randseq <- function(length){
  paste(sample(c("A","C","G","T"), length, replace = TRUE), collapse = '')
}
# The same but element-wise
randseq_elwise <- function(lengths){
  seqs <- c()
  for(i in lengths){
    seqs <- c(seqs, paste(sample(c("A","C","G","T"), i, replace = TRUE), collapse = ''))
  }
  return(seqs)
}

# randseqlist() gives a vector of specified length with sequences of specified length with random DNA seq
randseqlist <- function(listlength, seqlength){
  replicate(listlength, paste(sample(c("A","C","G","T"), seqlength, replace = TRUE), collapse = ''))
}

# example: how to create a list of sequences with randomly selected length
# all sequences will be of the same length!
randseqlist(4, sample(1:10, 1, replace = TRUE))

# generating toy data function
toydata <- function(data_length, spike_in = NA){
  toy_seq_list <- c()
  if(is.na(spike_in)){
    for(i in seq(data_length)){
      toy_seq_list <- c(toy_seq_list, randseq(str_length(toy_seq_list[i])))
    }
  }
  return(toy_seq_list)
}

# spiked-in toy data generator
toydata_spikein <- function(data_length, spike_in_left = "A", spike_in_right = "A",
                            distance_min = 1, distance_max = 15,
                            left_shoulder_min_length = 1, left_shoulder_max_length = 50,
                            right_shoulder_min_length = 1, right_shoulder_max_length = 50){
  toy_seq_list <- c()
  for(i in seq(data_length)){
    toy_seq_list <- c(toy_seq_list,
                      paste(randseq(
                        sample(
                          left_shoulder_min_length:left_shoulder_max_length, 1, replace = TRUE)),
                        spike_in_left,
                        randseq(sample(distance_min:distance_max)),
                        spike_in_right,
                        randseq(sample(
                          right_shoulder_min_length:right_shoulder_max_length, 1, replace = TRUE)),
                        sep = ''))
  }
  return(toy_seq_list)
}
 
# patternization makes the regex out of kmer sequence
patternization <- function(kmer, distances){
  kmer %>% str_split('') %>% unlist() -> split.mer
  k <- length(split.mer)
  pattern <- paste(split.mer[1],"[A|C|G|T]{",distances[1],"}",sep='')
  for(i in 2:(k-1)){
    pattern <- paste(pattern,split.mer[i],"[A|C|G|T]{",distances[i],"}",sep='')
  }
  pattern <- paste(pattern,split.mer[k],sep='')
  return(pattern)
}

# Colocalization matrix construction
# SLOW FUNCTION, OBVIOUSLY
colocalization <- function(seq_data, mode = 'count', kmer_length = 2){
  maxlength <- max(str_length(seq_data))
  N = 2
  K = 4^kmer_length
  kmers <- generate.kmers(.k = kmer_length)
  matrix_names<- list(c(as.character(seq(maxlength)-1)), kmers)
  coloc_matrix <- matrix(data = NA, nrow = maxlength, ncol = K, dimnames = matrix_names)
  if(mode == 'count'){
    for(i in seq(maxlength)){
      for(j in seq(K)){
        nt_left <- str_sub(kmers[j],1,1)
        nt_right<- str_sub(kmers[j],2,2)
        coloc_matrix[i,j] <- sum(str_count(seq_data,
                                           pattern = paste(nt_left,
                                                           "[A|C|G|T]{",i-1,"}",
                                                           nt_right,
                                                           sep='')),
                                 na.rm = TRUE)
      }
    }
  }
  if(mode == 'detect'){
    for(i in seq(maxlength)){
      for(j in seq(K)){
        nt_left <- str_sub(kmers[j],1,1)
        nt_right<- str_sub(kmers[j],2,2)
        coloc_matrix[i,j] <- sum(str_detect(seq_data,
                                            pattern = paste(nt_left,
                                                            "[A|C|G|T]{",i-1,"}",
                                                            nt_right,
                                                            sep='')),
                                 na.rm = TRUE)
      }
    }
  }
  return(coloc_matrix)
}

## Background subtraction
# Background synthesis - gives you the vector with the same number of sequences, each sequence is
# of the same length as in your initial data, but the DNA sequence is random
# Two modes: 'tech' gives you random DNA strings, just combination of ACGTs
#           'homo' gives you DNA strings randomly picked from Homo sapiens genome
# 'tech', of course, is muuuuuuch faster
synthesize_background <- function(seq_list, mode = 'tech'){
  if(mode == 'tech'){
    seq_list %>% str_length() %>% randseq_elwise() -> back_seq_list
  }
  if(mode == 'homo'){
    number_of_chrom <- nrow(chrlen)
    seq_list %>% str_length() -> lengths_seq_list
    back_seq_list <- c()
    for(seq_length in lengths_seq_list){
      rand.row <- sample(seq(1:number_of_chrom), 1, replace = TRUE)
      rand.chrom <- chrlen$Chromosome[rand.row]
      start.point <- sample(seq(0:(chrlen$Length.bp[rand.row]-seq_length)), 1, replace = TRUE)
      back_seq <- as.character(getSeq(Hsapiens, rand.chrom,
                                      start = start.point, end = start.point + seq_length - 1))
      back_seq_list <- c(back_seq_list, back_seq)
    }
  }
  return(back_seq_list)
}

subtract_background <- function(coloc_matrix, coloc_matrix_background){
  if(!identical(dimnames(coloc_matrix), dimnames(coloc_matrix_background))){
    stop("Dimensions of the matrices are not identical")
  }
  
  coloc_matrix_back_subtracted <- matrix(data = NA, nrow = nrow(coloc_matrix), ncol = ncol(coloc_matrix),
                                         dimnames = dimnames(coloc_matrix))
  for(i in seq(nrow(coloc_matrix))){
    for(j in seq(ncol(coloc_matrix))){
      coloc_matrix_back_subtracted[i,j] <-
        coloc_matrix[i,j] - coloc_matrix_background[i,j]
    }
  }
  
  return(coloc_matrix_back_subtracted)
}



#############  DATA       ##################
## Importing the chromosome length data for GRCh38
chrlen <- read.table(file = 'C:/Users/User/Documents/Lab/LongBCRs/MethodPaper/RunTables/GRCh38_chrlength.txt',
                     sep = '\t', header = TRUE)
chrlen <- select(chrlen, Chromosome, Length.bp)
chrlen$Length.bp <- as.numeric(str_replace_all(chrlen$Length.bp, ',', ''))
chrlen$Chromosome <- as.character(chrlen$Chromosome)
chrlen$Chromosome[chrlen$Chromosome == 'M'] <- 'MT'


#############  BODY       ##################
## COLOCALIZATIONAL PATTERNS
colocalization(toy_seq_list) -> coloc_matrix_count
colocalization(toy_seq_list, mode = 'detect') -> coloc_matrix_detect

# Background data synthesis
synthesize_background(toy_seq_list, mode = 'tech') -> toy_seq_list_back_tech
synthesize_background(toy_seq_list, mode = 'homo') -> toy_seq_list_back_homo

# Colocalization matrices for background
colocalization(toy_seq_list_back_tech) -> coloc_matrix_count_back_tech
colocalization(toy_seq_list_back_homo) -> coloc_matrix_count_back_homo
colocalization(toy_seq_list_back_tech, mode = 'detect') -> coloc_matrix_detect_back_tech
colocalization(toy_seq_list_back_homo, mode = 'detect') -> coloc_matrix_detect_back_homo

## Subtracting background
# technical background from count colocalization matrix
subtract_background(coloc_matrix_count, coloc_matrix_count_back_tech) -> col_mat_count_minus_tech
# homo background from count colocalization matrix
subtract_background(coloc_matrix_count, coloc_matrix_count_back_homo) -> col_mat_count_minus_homo
# technical background from detect colocalization matrix
subtract_background(coloc_matrix_count, coloc_matrix_detect_back_tech) -> col_mat_detect_minus_tech
# homo background from detect colocalization matrix
subtract_background(coloc_matrix_count, coloc_matrix_detect_back_homo) -> col_mat_detect_minus_homo

# technical background from homo background
subtract_background(coloc_matrix_count_back_homo, coloc_matrix_count_back_tech) -> col_mat_count_homo_minus_tech
subtract_background(coloc_matrix_detect_back_homo, coloc_matrix_detect_back_tech) -> col_mat_detect_homo_minus_tech


# Melting the data for plotting
col_mat_count_minus_tech_molten <- melt(col_mat_count_minus_tech)
col_mat_count_minus_homo_molten <- melt(col_mat_count_minus_homo)
col_mat_detect_minus_tech_molten <- melt(col_mat_detect_minus_tech)
col_mat_detect_minus_homo_molten <- melt(col_mat_detect_minus_homo)
col_mat_count_homo_minus_tech_molten <- melt(col_mat_count_homo_minus_tech)
col_mat_detect_homo_minus_tech_molten <- melt(col_mat_detect_homo_minus_tech)

# Trying to make it prettier
col_mat_count_homo_minus_tech %>% abs() %>% log() -> col_mat_count_homo_minus_tech_log
col_mat_count_homo_minus_tech_log_molten <- melt(col_mat_count_homo_minus_tech_log)

g1 <- ggplot(data = col_mat_count_minus_tech_molten, aes(y=Var1,x=Var2))+
  geom_tile(aes(fill = value))+
  scale_fill_viridis_c()+
  theme_minimal()+
  xlab("Colocalized nucleotides")+
  ylab("Distance btw colocalized nucleotides, nt")+
  labs(fill = "Count")


g2 <- ggplot(data = col_mat_count_minus_homo_molten, aes(y=Var1,x=Var2, fill = value))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()+
  xlab("Colocalized nucleotides")+
  ylab("Distance btw colocalized nucleotides")+
  labs(fill = "Count")


g3 <- ggplot(data = col_mat_detect_minus_tech_molten, aes(y=Var1,x=Var2, fill = value))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()+
  xlab("Colocalized nucleotides")+
  ylab("Distance btw colocalized nucleotides")+
  labs(fill = "Sequences")


g4 <- ggplot(data = col_mat_detect_minus_homo_molten, aes(y=Var1,x=Var2, fill = value))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()+
  xlab("Colocalized nucleotides")+
  ylab("Distance btw colocalized nucleotides")+
  labs(fill = "Sequences")

g5 <- ggplot(data = col_mat_count_homo_minus_tech_molten, aes(y=Var1,x=Var2, fill = value))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()+
  xlab("Colocalized nucleotides")+
  ylab("Distance btw colocalized nucleotides")+
  labs(fill = "Count")

g6 <- ggplot(data = col_mat_detect_homo_minus_tech_molten, aes(y=Var1,x=Var2, fill = value))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()+
  xlab("Colocalized nucleotides")+
  ylab("Distance btw colocalized nucleotides")+
  labs(fill = "Sequences")

g7 <- ggplot(data = col_mat_count_homo_minus_tech_log_molten, aes(y=Var1,x=Var2, fill = value))+
  geom_tile()+
  scale_fill_viridis_c()+
  theme_minimal()+
  xlab("Colocalized nucleotides")+
  ylab("Distance btw colocalized nucleotides")+
  labs(fill = "Count")


ggsave(paste('colocalization_count_tech',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = g1)

ggsave(paste('colocalization_count_homo',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = g2)

ggsave(paste('colocalization_detect_tech',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = g3)

ggsave(paste('colocalization_detect_homo',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = g4)

ggsave(paste('colocalization_count_homovstech',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = g5)

ggsave(paste('colocalization_detect_homovstech',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = g6)

ggsave(paste('colocalization_count_homovstech_log',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = g7)

################################

get.kmers("AAGCTGTGCACGTAAG", .k = 3)

col_mat_count_homo_minus_tech_log %>% is.infinite() %>% sum()
col_mat_count_homo_minus_tech_log[col_mat_count_homo_minus_tech_log %>% is.infinite()] <- 0
