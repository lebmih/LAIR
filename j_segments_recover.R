setwd("C:/Users/User/Documents/Lab/LongBCRs/MethodPaper/RunTables")

all_inserts <- read.table(file = 'Output/inserts_allRuns_11_04_19.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

library(seqinr)
j_segments <- read.fasta(file = 'C:/Users/User/Documents/Lab/LongBCRs/Homo_IGH_J_segments.fa',
                         seqtype = 'DNA', as.string = TRUE)

###########################
##### TIMESTAMP: 19.07.19
#### Annotate J segments in the contigs of undefined type 
detect_j_segment <- function(contig, threshold = 10, j_segments){
  edits <- c()
  for (i in j_segments){
    psa <- pairwiseAlignment(i, contig, type = 'global')
    edits <- c(edits, nedit(psa))
  }
  if (min(edits) > threshold){
    print('Trying to reverse the sequence')
    edits <- c()
    for (i in j_segments){
      psa <- pairwiseAlignment(i, reverseComplement(contig), type = 'global')
      edits <- c(edits, nedit(psa))
    }
    if (min(edits) < threshold){
      print('Alignment successfull')
    } else {
      return('No alignment')
    }
    j_segment <- getName(j_segments[edits == min(edits)][[1]])[[1]]
    
  } else {
    print('Alignment successfull')
    j_segment <- getName(j_segments[edits == min(edits)][[1]])[[1]]
  }
}

contigs_of_interest <- DNAStringSet(all_inserts$contig.seq[is.na(all_inserts$insert.type) & is.na(all_inserts$Jgene.id)])

contigs_j_segments <- c()
for (i in seq(length(contigs_of_interest))){
  contigs_j_segments <- c(contigs_j_segments, detect_j_segment(contigs_of_interest[i], j_segments = j_segments))
}


all_inserts$Jgene.id[is.na(all_inserts$insert.type) & is.na(all_inserts$Jgene.id)] <- contigs_j_segments
describe(all_inserts$Jgene.id)

#### Classify the sequences further
all_inserts$insert.type[is.na(all_inserts$insert.type) & all_inserts$Jgene.id == 'No alignment'] <- 'V-CH1'

describe(all_inserts$insert.type)

all_inserts_missing <- all_inserts[is.na(all_inserts$insert.type),]
## Manually determined type (by SnapGene)
all_inserts$insert.type[is.na(all_inserts$insert.type)] <- c(rep('V-DJ',4),'VD-J',rep('V-DJ',3))

write.table(all_inserts, file = 'Output/inserts_allRuns_22_07_19.txt', sep = '\t',
            row.names = FALSE)

#############################################################
