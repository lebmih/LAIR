setwd("C:/Users/User/Documents/Lab/LongBCRs/MethodPaper/RunTables")

all_inserts <- read.table(file = 'data_cured_Kathrins.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)

all_inserts$contig.id <- seq(1:nrow(all_inserts))

library(seqinr)
library(Biostrings)
j_segments <- read.fasta(file = 'd_segments_study/database/j_segments_imgt_simple.fa',
                         seqtype = 'DNA', as.string = TRUE)


###########################
##### TIMESTAMP: 19.07.19
#### Annotate J segments in the contigs of undefined type 
detect_segment <- function(contig, threshold = 0, segments, mode = "score", enable.reversal = FALSE, align.type = "local"){
  edits <- c()
  scores <- c()
  for (i in segments){
    psa <- pairwiseAlignment(i, contig, type = align.type, gapOpening = 4)
    edits <- c(edits, nedit(psa))
    scores <- c(scores, score(psa))
  }
  if(mode == "edit"){
    if ((min(edits) > threshold) & enable.reversal){
      print('Trying to reverse the sequence')
      edits <- c()
      for (i in segments){
        psa <- pairwiseAlignment(i, reverseComplement(contig), type = align.type)
        edits <- c(edits, nedit(psa))
      }
      if (min(edits) < threshold){
        print('Alignment successfull')
      } else {
        return('No alignment')
      }
      segment <- getName(segments[edits == min(edits)][[1]])[[1]]
      min_edit <- min(edits)
      return(paste(segment, min_edit, sep = ";"))
    } else if(min(edits) > threshold) {
      return('No alignment')
    } else {
      print('Alignment successfull')
      segment <- getName(segments[edits == min(edits)][[1]])[[1]]
      min_edit <- min(edits)
      return(paste(segment, min_edit, sep = ";"))
    }
  }
  if(mode == "score"){
    if ((max(scores) < threshold) & enable.reversal){
      print('Trying to reverse the sequence')
      scores <- c()
      for (i in segments){
        psa <- pairwiseAlignment(i, reverseComplement(contig), type = align.type)
        scores <- c(scores, score(psa))
      }
      if (max(scores) > threshold){
        print('Alignment successfull')
      } else {
        return('No alignment')
        print('No alignment')
      }
      segment <- getName(segments[scores == max(scores)][[1]])[[1]]
      top_score <- max(scores)
      return(paste(segment, top_score, sep = ";"))
    } else if(max(edits) < threshold) {
      return('No alignment')
      print('No alignment')
    } else {
      print('Alignment successfull')
      segment <- getName(segments[scores == max(scores)][[1]])[[1]]
      top_score <- max(scores)
      return(paste(segment, top_score, sep = ";"))
    }
  }
}

sum(is.na(all_inserts$insert.type))
# 186 contigs with non-identified insert type
sum(is.na(all_inserts$insert.type) & is.na(all_inserts$Vgene.id))
# 0 of them don't have V gene identified (makes sense, the pipeline pivots on the V segment to define the rest of the contig)
sum(is.na(all_inserts$insert.type) & is.na(all_inserts$Dgene.id))
# 148 contigs with non-identified D segment
sum(is.na(all_inserts$insert.type) & is.na(all_inserts$Jgene.id))
# All 186 contigs don't have j segment identified

na.contigs <- all_inserts %>% filter(is.na(insert.type))
na.contigs$contig.seq[2]
na.contigs$contig.id[2]

na.contigs %>% mutate(v.ins.gap = mini.start.insert.coord - V.contig.coord.end,
                      ins.ch1.gap = const.contig.coord.start - max.end.insert.coord) -> na.contigs

na.contigs %>% mutate(ch.length = const.contig.coord.end - const.contig.coord.start) -> na.contigs


ggplot(data = na.contigs)+
  geom_point(aes(x = v.ins.gap, y = ch.length))

na.contigs %>% filter(ins.ch1.gap < 0) -> na.contigs.wrong
na.contigs.wrong$contig.id -> na.contigs.ch1gap.negative
na.contigs %>% filter(ch.length > 50) -> na.contigs.wrong
na.contigs.wrong$contig.id -> na.contigs.ch1.toolong
intersect(na.contigs.ch1.toolong, na.contigs.ch1gap.negative) %>% length()
# All contigs with ch too long are giving the negative insert-ch1 gap
# But there are two contigs with ch not too long, but still giving the negative insert-ch1 gap
# Let's check them manually
setdiff(na.contigs.ch1gap.negative, na.contigs.ch1.toolong)
# Ah, those are 797 and 798, the ones that have both V-insert and insert-CH1 gaps negative
# I already checked them, they have a problem with being unreversed
# Then let's take a closer look at contigs that have too long CH1
rm(na.contigs.ch1gap.negative, na.contigs.ch1.toolong)

constant.no.primer <- read.fasta(file = "constants_no_primer_071019.txt", as.string = TRUE)
constant.no.primer <- unlist(constant.no.primer)

contigs_of_interest <- DNAStringSet(na.contigs.wrong$contig.seq)

contigs.const <- c()
for(i in seq(nrow(na.contigs.wrong))){
  contigs.const <- c(contigs.const, detect_segment(contigs_of_interest[i], segments = constant.no.primer))
}



na.contigs.wrong$constants.alignments <- contigs.const

# Manually checking those whose score is lower than median: 997, 999, 1119 and 1212
na.contigs.wrong$contig.seq[na.contigs.wrong$contig.id == 997]
# Checked. they are all wrong, they indeed do not contain the CH1 and probably are the result of unspecific annealing 
wrong.contigs.to.filter <- c(997, 999, 1119, 1212)

# what of all others?
# let's find the real coordinates of CH1 in them
na.contigs.wrong %>% separate(col = "constants.alignments", into = c("const.new", "const.score")) -> na.contigs.wrong

new.coord <- c()
for(i in seq(nrow(na.contigs.wrong))){
  psa <- pairwiseAlignment(unlist(constant.no.primer[na.contigs.wrong$const.new[i]]), contigs_of_interest[i],
                           type = "local")
  print(score(psa))
  new.coord <- c(new.coord, start(subject(psa)))
}

na.contigs.wrong$const.contig.coord.start <- new.coord

na.contigs$const.contig.coord.start[na.contigs$contig.id %in% na.contigs.wrong$contig.id] <- new.coord
# Now there should not be negative gaps between insert and CH1

na.contigs %>% mutate(v.ins.gap = mini.start.insert.coord - V.contig.coord.end,
                      ins.ch1.gap = const.contig.coord.start - max.end.insert.coord) -> na.contigs

na.contigs %>% mutate(ch.length = const.contig.coord.end - const.contig.coord.start) -> na.contigs


ggplot(data = na.contigs)+
  geom_point(aes(x = v.ins.gap, y = ins.ch1.gap))

# There are still some
na.contigs$contig.id[na.contigs$ins.ch1.gap < 0]
# 797 and 798 should be reversed
# 997, 999, 1119, 1212 should be deleted
# 695??? why
# I checked 695 manually, mtDNA is aligned to the part of IGHM
na.contigs$max.end.insert.coord[na.contigs$contig.id == 695] <- 157

# Filtering out 997, 999, 1119, 1212
na.contigs %>% filter(!contig.id %in% c(997, 999, 1119, 1212)) -> na.contigs

# Reversing the 797 and 798
na.contigs %>% filter(contig.id %in% c(797, 798)) -> na.contigs.797.798
na.contigs.797.798 %>% mutate(const.contig.coord.start.new = contig.length - const.contig.coord.end,
                              const.contig.coord.end.new = contig.length - const.contig.coord.start,
                              V.contig.coord.start.new = contig.length - V.contig.coord.end,
                              V.contig.coord.end.new = contig.length - V.contig.coord.start,
                              mini.start.insert.coord.new = contig.length - max.end.insert.coord,
                              max.end.insert.coord.new = contig.length - mini.start.insert.coord) %>%
  mutate(const.contig.coord.start = const.contig.coord.start.new,
         const.contig.coord.end = const.contig.coord.end.new,
         V.contig.coord.start = V.contig.coord.start.new,
         V.contig.coord.end = V.contig.coord.end.new,
         mini.start.insert.coord = mini.start.insert.coord.new,
         max.end.insert.coord = max.end.insert.coord.new) %>%
  select(-c("const.contig.coord.start.new", "const.contig.coord.end.new",
            "V.contig.coord.start.new", "V.contig.coord.end.new",
            "mini.start.insert.coord.new", "max.end.insert.coord.new")) -> na.contigs.797.798

na.contigs.797.798$contig.seq <- as.character(reverseComplement(DNAStringSet(na.contigs.797.798$contig.seq)))
na.contigs.797.798$contigWithoutInsert.seq <- as.character(reverseComplement(DNAStringSet(na.contigs.797.798$contigWithoutInsert.seq)))

# Insert sequence is not even aligning to the contig
na.contigs.797.798$contig.seq[1]
pairwiseAlignment(na.contigs.797.798$insert.seq.s.[1], na.contigs.797.798$contig.seq[1])
# From 119bp insert.seq.s. only 81 bp is aligned correctly, while left part of sequence, AGATCCCAGTGGAGAACACCAAAGCGAGCGAGGAGGAG
# is also coming from this gene, but is not present in the contig at all
# WHY IS THIS THING HERE THEN?
# The insert length is correct though
# Dig deeper into that later
na.contigs %>% filter(!contig.id %in% c(797, 798)) %>% bind_rows(na.contigs.797.798) -> na.contigs
rm(na.contigs.797.798, na.contigs.ch1.toolong, na.contigs.wrong)

# Now let's check the plot again
na.contigs %>% mutate(v.ins.gap = mini.start.insert.coord - V.contig.coord.end,
                      ins.ch1.gap = const.contig.coord.start - max.end.insert.coord) -> na.contigs

na.contigs %>% mutate(ch.length = const.contig.coord.end - const.contig.coord.start) -> na.contigs


ggplot(data = na.contigs)+
  geom_point(aes(x = v.ins.gap, y = ins.ch1.gap))

# Now there is no negative gap between Insert and CH1, but still there are dots with negative V-Insert gap
na.contigs %>% filter(v.ins.gap < 0) -> na.contigs.wrong
na.contigs.wrong$contig.id
# Again those are 797, 798 and 832
# But I already corrected 797 and 798
# The V segment borders are assigned erroneously
na.contigs$V.contig.coord.end[na.contigs$contig.id %in% c(797, 798)] <-
  na.contigs$mini.start.insert.coord[na.contigs$contig.id %in% c(797, 798)] - 1
# I also checked 823 manually
na.contigs$V.contig.coord.end[na.contigs$contig.id == 832] <- 66
na.contigs$mini.start.insert.coord[na.contigs$contig.id == 832] <- 67
na.contigs$max.end.insert.coord[na.contigs$contig.id == 832] <- 236

rm(na.contigs.wrong)

# let's check the plot again
na.contigs %>% mutate(v.ins.gap = mini.start.insert.coord - V.contig.coord.end,
                      ins.ch1.gap = const.contig.coord.start - max.end.insert.coord) -> na.contigs

na.contigs %>% mutate(ch.length = const.contig.coord.end - const.contig.coord.start) -> na.contigs


ggplot(data = na.contigs)+
  geom_point(aes(x = v.ins.gap, y = ins.ch1.gap))+
  #scale_y_continuous(limits = c(0,100))+
  #scale_x_continuous(limits = c(0,100))+
  theme_classic()+
  xlab("Gap between V and Insert, bp")+
  ylab("Gap between Insert and CH1, bp")
# Great, now it's clean
# Let's try to classify the contigs now
# But first let's make the plot for all inserts
all_inserts %>% mutate(v.ins.gap = mini.start.insert.coord - V.contig.coord.end,
                       ins.ch1.gap = const.contig.coord.start - max.end.insert.coord) -> all_inserts.1
all_inserts.1 %>% mutate(ch.length = const.contig.coord.end - const.contig.coord.start) -> all_inserts.1
ggplot(data = all_inserts.1)+
  geom_point(aes(x = v.ins.gap, y = ins.ch1.gap))+
  #scale_y_continuous(limits = c(0,100))+
  #scale_x_continuous(limits = c(0,100))+
  theme_classic()+
  xlab("Gap between V and Insert, bp")+
  ylab("Gap between Insert and CH1, bp")
# There are a lot of contigs in normal set that have negative gaps in between V and Insert and Insert and CH1
# First I will have a look at those whose Insert-CH1 gap is negative
all_inserts.1 %>% filter(v.ins.gap < 0) -> all_inserts.wrong
# There are 31 of them
# Ah, wait, I didn't merge the corrected table with the all_inserts table
all_inserts.1 %>% filter(!contig.id %in% na.contigs$contig.id) %>%
  filter(!contig.id %in% c(997, 999, 1119, 1212)) %>% 
  bind_rows(na.contigs) -> all_inserts.1
all_inserts.1 %>% mutate(v.ins.gap = mini.start.insert.coord - V.contig.coord.end,
                       ins.ch1.gap = const.contig.coord.start - max.end.insert.coord) -> all_inserts.1
all_inserts.1 %>% mutate(ch.length = const.contig.coord.end - const.contig.coord.start) -> all_inserts.1
# Let's have a look again
ggplot(data = all_inserts.1)+
  geom_point(aes(x = v.ins.gap, y = ins.ch1.gap))+
  #scale_y_continuous(limits = c(0,100))+
  #scale_x_continuous(limits = c(0,100))+
  theme_classic()+
  xlab("Gap between V and Insert, bp")+
  ylab("Gap between Insert and CH1, bp")
# This is a little bit better, but there are still contigs with negative gaps
# Negative V-insert gap:
all_inserts.1 %>% filter(v.ins.gap < 0) -> all_inserts.wrong
all_inserts.wrong$contig.id
# Wrong assignment of Vcontig.end coordinate
contigs_of_interest <- DNAStringSet(str_sub(all_inserts.wrong$contig.seq, 1, 50))
v_segments <- read.fasta(file = 'd_segments_study/database/v_segments_imgt_simple.fa',
                         seqtype = 'DNA', as.string = TRUE)

new.coord <- c()
for(i in seq(nrow(all_inserts.wrong))){
  psa <- pairwiseAlignment(unlist(v_segments[all_inserts.wrong$Vgene.id[i]]), contigs_of_interest[i],
                           type = "local")
  print(score(psa))
  new.coord <- c(new.coord, end(subject(psa)))
}

all_inserts.wrong$V.contig.coord.end <- new.coord

all_inserts.1 %>% filter(!contig.id %in% all_inserts.wrong$contig.id) %>% bind_rows(all_inserts.wrong) -> all_inserts.1
# let's check the plot again
all_inserts.1 %>% mutate(v.ins.gap = mini.start.insert.coord - V.contig.coord.end,
                         ins.ch1.gap = const.contig.coord.start - max.end.insert.coord) -> all_inserts.1
all_inserts.1 %>% mutate(ch.length = const.contig.coord.end - const.contig.coord.start) -> all_inserts.1
# Let's have a look again
ggplot(data = all_inserts.1)+
  geom_point(aes(x = v.ins.gap, y = ins.ch1.gap))+
  #scale_y_continuous(limits = c(0,100))+
  #scale_x_continuous(limits = c(0,100))+
  theme_classic()+
  xlab("Gap between V and Insert, bp")+
  ylab("Gap between Insert and CH1, bp")

# Good, there are no contigs with negative V-insert gap anymore
all_inserts.1 %>% filter(v.ins.gap < 0) %>% View()
# Well, there are 4 contigs with slight overlapping between V and insert, but the max overlap is 6 bp
# Now let's take a look at the negative insert-ch1 gap
all_inserts.1 %>% filter(ins.ch1.gap < 0) -> all_inserts.wrong
# There are 100 of them
# manually checking some
all_inserts.wrong$contig.seq[1]
# The constant is assigned wrong
contigs_of_interest <- DNAStringSet(all_inserts.wrong$contig.seq)

contigs.const <- c()
for(i in seq(nrow(all_inserts.wrong))){
  contigs.const <- c(contigs.const, detect_segment(contigs_of_interest[i], segments = constant.no.primer))
}

all_inserts.wrong$constants.alignments <- contigs.const

all_inserts.wrong %>% separate(col = "constants.alignments", into = c("const.new", "const.score")) -> all_inserts.wrong

new.coord <- c()
for(i in seq(nrow(all_inserts.wrong))){
  psa <- pairwiseAlignment(unlist(constant.no.primer[all_inserts.wrong$const.new[i]]), contigs_of_interest[i],
                           type = "local")
  print(score(psa))
  new.coord <- c(new.coord, start(subject(psa)))
}

all_inserts.wrong$const.contig.coord.start <- new.coord
# Let's check the contigs with low score of alignment
all_inserts.wrong %>% filter(const.score < 25) -> all_inserts.wrong.low.score
# Just five of them
all_inserts.wrong.low.score$contig.seq[5]
# 679 is wrong, delete
# 769 is wrong, delete
# 883 is wrong, delete
# 1117 is wrong, delete
# 1225 is wrong, delete
# Thus, all of them are the product of unspecific binding
# Let's develop a method to find those contigs
# Their hallmark is an overlap of the suppression primer with the insert
contigs_of_interest <- DNAStringSet(all_inserts.wrong.low.score$contig.seq)
primers <- read.fasta(file = "sup_primers_gs_parts.txt", as.string = TRUE)

prim.anneal <- c()
for(i in seq(nrow(all_inserts.wrong.low.score))){
  prim.anneal <- c(prim.anneal, detect_segment(contigs_of_interest[i], segments = primers))
}

all_inserts.wrong.low.score$prim.anneal <- prim.anneal

all_inserts.wrong.low.score %>% separate(col = "prim.anneal", into = c("prim.anneal", "prim.anneal.score"), sep =";") -> all_inserts.wrong.low.score

new.coord <- c()
for(i in seq(nrow(all_inserts.wrong.low.score))){
  psa <- pairwiseAlignment(unlist(primers[all_inserts.wrong.low.score$prim.anneal[i]]), contigs_of_interest[i],
                           type = "local")
  print(score(psa))
  new.coord <- c(new.coord, start(subject(psa)))
}
all_inserts.wrong.low.score$max.end.insert.coord - new.coord
# If this number is positive, there is an overlap
# let's check whether there are more of such contigs in the data
# But first let's check the plot
all_inserts.wrong %>% mutate(v.ins.gap = mini.start.insert.coord - V.contig.coord.end,
                         ins.ch1.gap = const.contig.coord.start - max.end.insert.coord) -> all_inserts.wrong

all_inserts.wrong %>% select(-c("const.new", "const.score")) -> all_inserts.wrong

all_inserts.1 %>% filter(!contig.id %in% all_inserts.wrong$contig.id) %>% bind_rows(all_inserts.wrong) -> all_inserts.1

all_inserts.1 %>% mutate(v.ins.gap = mini.start.insert.coord - V.contig.coord.end,
                         ins.ch1.gap = const.contig.coord.start - max.end.insert.coord) -> all_inserts.1
all_inserts.1 %>% mutate(ch.length = const.contig.coord.end - const.contig.coord.start) -> all_inserts.1

contigs.to.delete <- c(679, 769, 883, 1117, 1225)
all_inserts.1 %>% filter(!contig.id %in% contigs.to.delete) -> all_inserts.1


# Let's have a look again
ggplot(data = all_inserts.1)+
  geom_point(aes(x = v.ins.gap, y = ins.ch1.gap))+
  #scale_y_continuous(limits = c(0,100))+
  #scale_x_continuous(limits = c(0,100))+
  theme_classic()+
  xlab("Gap between V and Insert, bp")+
  ylab("Gap between Insert and CH1, bp")

# Now it's perfect
# No, not perfect, of course. The gap between Insert and CH1 must be either zero or around the length of J or around the length of DJ
# But for some it's more than 200
# Check the extremes
all_inserts.1 %>% filter(v.ins.gap > 400) %>% View()

write.table(all_inserts.1, file = "all_inserts081019.txt", sep = "\t", quote = FALSE, row.names = FALSE)



################################


contigs_of_interest <- DNAStringSet(all_inserts$contig.seq[is.na(all_inserts$insert.type) & is.na(all_inserts$Jgene.id)])

contigs_j_segments <- c()
for (i in seq(length(contigs_of_interest))){
  contigs_j_segments <- c(contigs_j_segments, detect_segment(contigs_of_interest[i], segments = j_segments))
}

contigs_of_interest <- DNAStringSet(all_inserts$contig.seq[is.na(all_inserts$insert.type) & is.na(all_inserts$Dgene.id)])

contigs_d_segments <- c()
for (i in seq(length(contigs_of_interest))){
  contigs_d_segments <- c(contigs_d_segments, detect_segment(contigs_of_interest[i], segments = d_segments))
}


all_inserts$Jgene.id[is.na(all_inserts$insert.type) & is.na(all_inserts$Jgene.id)] <- contigs_j_segments

all_inserts$Dgene.id[is.na(all_inserts$insert.type) & is.na(all_inserts$Dgene.id)] <- contigs_d_segments

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
