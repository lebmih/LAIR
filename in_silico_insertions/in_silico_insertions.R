#### In silico simulation of the insertions.
#### This script creates a dataset of the insertions with a randomized
#### coordinate to provide a control for the correct characterization of the in
#### vivo occurring inserts.

# Working folder ----------------------------------------------------------
setwd("C:/Users/User/Documents/Lab/LongBCRs/bioinformatics/lair")


# Libraries ---------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(gtable)
library(grid)
library(gridExtra)
library(lebutils)
library(seqinr)
library(Biostrings)
library(DescTools)


# Functions ---------------------------------------------------------------
# This function checks whether the input coordinate is inside the series of
# ranges.
in.ranges <- function(coordinate, starts, ends){
  starts.n <- length(starts)
  ends.n <- length(ends)
  if (starts.n != ends.n) {
    stop("Starts and Ends vectors not of the same length!")
  }
  if(any(starts > ends)) {
    stop("There are inverted ranges!")
  }
  count <- 0
  if(starts.n == 0 | ends.n == 0){
    return(count)
  } else {
    for(i in seq(starts.n)){
      if(coordinate >= starts[i] & coordinate <= ends[i]){
        count <- count + 1
      }
    }
  }

  return(count)
}


# Data import -------------------------------------------------------------
# The contigs table is the output of the main pipeline that characterizes the
# insertions.
contigs <- read.table("contigs_DATE.txt", header = TRUE,
                      sep = "\t", stringsAsFactors = FALSE)
# Assembly report contains information on chromosome length for GRCh38,
# downloaded from https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.26/
assembly.report <- read.table(file = "input/GRCh38_latest_assembly_report.txt",
                              stringsAsFactors = FALSE, header = FALSE,
                              fill = TRUE, sep = "\t")
# Formatting the table.
colnames(assembly.report) <-
  c("seq.name","seq.role","assigned.molecule",
    "assigned.molecule.location.type","genbank.accn",
    "relationship","refseq.accn","assembly.unit",
    "seq.length","UCSC.style.name")
assembly.report <- assembly.report %>%
  dplyr::filter(seq.role == "assembled-molecule")
# Importing the centromeres coordinates.
centromeres <- read.table(file = "hg38_centromeres.txt",
                          stringsAsFactors = FALSE, header = TRUE, fill = TRUE, sep = "\t")
# Picking the most conservative (the widest regions) coordinates for the
# centromeres.
centromeres <- centromeres %>%
  dplyr::select(-c("bin", "name")) %>%
  group_by(chrom) %>%
  summarise(centromer.start = min(chromStart),
            centromer.end = max(chromEnd)) %>%
  mutate(centromer.center = (centromer.start+centromer.end)/2)
# Importing the GRCh38 genes coordinates, downloaded from UCSC Table Browser.
hg38.genes <- read.table("hg38_genes.bed", sep = "\t", stringsAsFactors = FALSE)
colnames(hg38.genes) <- c("chrom", "start", "end", "gene.name",
                          "score", "strand", "start1", "start2",
                          "score1", "exon.number", "exon.length", "exon.start")
# Exon coordinates from UCSC Table Browser.
hg38.exons <- read.table("hg38_exons.bed", sep = "\t", stringsAsFactors = FALSE)
colnames(hg38.exons) <- c("chrom", "start", "end", "exon.name", "score", "strand")
hg38.exons <- hg38.exons %>% mutate(length = end - start + 1)

# Importing the coordinates of the gaps - the unmapped sequences in human
# genome.
gaps <- read.table("hg38_gaps.bed", sep = "\t",
                   stringsAsFactors = FALSE)
colnames(gaps) <- c("chrom", "start", "end")
gaps <- gaps %>% filter(!str_detect(chrom, "_"))


# Generating the set of artificial contigs --------------------------------
if (!"Generated" %in% unique(contigs$ins.type)) {
  contigs.gen <- data.frame()
  for (j in seq(100)){
    contigs1 <- contigs
    contigs1 <- contigs1 %>% dplyr::filter(!is.na(seq.length))
    rand_position <- c()
    for (i in seq(nrow(contigs1))){
      # Extracting gaps in this chromosome to avoid them.
      this.chrom <- contigs1$chrom[i]
      gaps.this.chrom <- gaps %>% dplyr::filter(chrom == this.chrom)
      # Randomizing the coordinate till it is not in the gaps.
      if (nrow(gaps.this.chrom) > 0) {
        ok <- FALSE
        while (!ok) {
          rand_position.tmp <- sample(1:contigs1$seq.length[i],1)
          if (in.ranges(rand_position.tmp, gaps.this.chrom$start,
                        gaps.this.chrom$end) == 0) {
            ok <- TRUE
          }
        }
      } else {
        rand_position.tmp <- sample(1:contigs1$seq.length[i],1)
      }
      
      rand_position <- c(rand_position, rand_position.tmp)
      rm(rand_position.tmp, this.chrom, gaps.this.chrom)
    }
    contigs1$min.g.start <- rand_position - contigs1$insert.length/2
    contigs1$max.g.end <- rand_position + contigs1$insert.length/2
    contigs1$min.g.start[contigs1$min.g.start < 1] <- 1
    contigs1$max.g.end[contigs1$max.g.end > contigs1$seq.length] <-
      contigs1$seq.length[contigs1$max.g.end > contigs1$seq.length]
    contigs.gen <- contigs.gen %>% bind_rows(contigs1)
    if (j %in% seq(0,100,10)) {
      print(paste0(j, "% job done."))
    }
  }
  contigs.gen$population <- "Generated"
  rm(contigs1, i, j, rand_position)
}

contigs.gen[,c("donor","race","cell.type")] <- "Generated"

# Merging two dataframes together to analyze both -------------------------
contigs.gen$ins.type[contigs.gen$ins.type == "J-CH1"] <- "J-CH1.gen"
contigs.gen$ins.type[contigs.gen$ins.type == "VDJ"] <- "VDJ.gen"
contigs.gen$ins.type[contigs.gen$ins.type == "V-CH1"] <- "V-CH1.gen"
# This dataset was used at the first stage of study, the paper includes the
# modification below - this modification accounts for the exon-intron skew in
# the natural data by creating the equal portions of exon-intron generated data.
contigs <- contigs %>% bind_rows(contigs.gen)



# The latest version of in silico contigs generation ----------------------
# This modification splits the natural dataset in the portions depending on the
# exon-intron structure of the insertion origin (i.e. insert class). Then it
# randomizes the insertions taking a random exon for exonic insertions, random
# intron for intronic insertions, etc.
contigs.gen.exon <- contigs.gen %>% filter(ins.class == "exon")
contigs.gen.pexon <- contigs.gen %>% filter(ins.class == "partial_exon")
contigs.gen.intron <- contigs.gen %>% filter(ins.class == "intron")
contigs.gen.inter <- contigs.gen %>% filter(ins.class == "intergenic")
# Randomization for the exonic insertions. The size of exon is retained allowing
# for +- 15 bp.
for(i in seq(nrow(contigs.gen.exon))){
  if(i %in% seq(0,100000,1000)){
    print(i)
  }
  target.len <- contigs.gen.exon$insert.length[i]
  target.chrom <- contigs.gen.exon$chrom[i]
  window <- 15
  ok <- FALSE
  while(ok == FALSE){
    hg38.exons.target <- hg38.exons %>%
      filter(chrom == target.chrom & length >= target.len - window & length <= target.len + window)
    if(nrow(hg38.exons.target) > 0){
      ok <- TRUE
    } else {
      window <- window + 5
    }
  }
  target.exon <- sample(seq(nrow(hg38.exons.target)), 1)
  contigs.gen.exon$min.g.start[i] <- hg38.exons.target$start[target.exon]
  contigs.gen.exon$max.g.end[i] <- hg38.exons.target$end[target.exon]
  
}

contigs.gen.exon <- contigs.gen.exon %>%
  mutate(centr.dist = pmin(abs(min.g.start - centromer.center),
                           abs(centromer.center - max.g.end)),
         tel.dist = pmin(min.g.start, seq.length - max.g.end))


contigs.gen.exon$tel.dist %>% median()

# Randomization for the partial exons.
for(i in seq(nrow(contigs.gen.pexon))){
  if(i %in% seq(0,100000,1000)){
    print(i)
  }
  target.chrom <- contigs.gen.pexon$chrom[i]
  hg38.exons.target <- hg38.exons %>%
    filter(chrom == target.chrom)
  target.exon <- sample(seq(nrow(hg38.exons.target)), 1)
  contigs.gen.pexon$min.g.start[i] <-
    sample(seq(hg38.exons.target$start[target.exon] - contigs.gen.pexon$insert.length[i],
               hg38.exons.target$end[target.exon]), 1)
  contigs.gen.pexon$max.g.end[i] <- contigs.gen.pexon$min.g.start[i] + contigs.gen.pexon$insert.length[i] - 1
  
}

# For the introns I first need to define the introns borders.
hg38.introns <- read.table("hg38_introns.bed", sep = "\t", stringsAsFactors = FALSE)
colnames(hg38.introns) <- c("chrom", "start", "end", "intron.name", "score", "strand")
hg38.introns <- hg38.introns %>% mutate(length = end - start + 1)

for(i in seq(nrow(contigs.gen.intron))){
  if(i %in% seq(0,100000,1000)){
    print(i)
  }
  target.chrom <- contigs.gen.intron$chrom[i]
  target.length <- contigs.gen.intron$insert.length[i]
  
  hg38.introns.target <- hg38.introns %>%
    filter(chrom == target.chrom & length >= target.length)
  target.intron <- sample(seq(nrow(hg38.introns.target)), 1)
  contigs.gen.intron$min.g.start[i] <-
    sample(seq(hg38.introns.target$start[target.intron],
               hg38.introns.target$end[target.intron] - contigs.gen.intron$insert.length[i]), 1)
  contigs.gen.intron$max.g.end[i] <- contigs.gen.intron$min.g.start[i] + contigs.gen.intron$insert.length[i] - 1
  
}

# Same procedure for intergenic insertions.
gaps <- read.table("input/hg38_gaps.bed", sep = "\t",
                   stringsAsFactors = FALSE)
colnames(gaps) <- c("chrom", "start", "end")
gaps <- gaps %>% filter(!str_detect(chrom, "_"))


for (i in seq(nrow(contigs.gen.inter))){
  if(i %in% seq(0,100000,1000)){
    print(i)
  }
  target.chrom <- contigs.gen.inter$chrom[i]
  gaps.this.chrom <- gaps %>% dplyr::filter(chrom == target.chrom)
  hg38.genes.target <- hg38.genes %>% filter(chrom == target.chrom)
  # Randomizing the coordinate till it is not in the gaps and is not in the genes
  if (nrow(gaps.this.chrom) > 0 | nrow(hg38.genes.target) > 0) {
    ok <- FALSE
    while (!ok) {
      rand_position.tmp <- sample(1:contigs.gen.inter$seq.length[i],1)
      if ((in.ranges(rand_position.tmp, gaps.this.chrom$start, gaps.this.chrom$end) == 0) &
          (in.ranges(rand_position.tmp, hg38.genes.target$start, hg38.genes.target$end) == 0)) {
        ok <- TRUE
      }
    }
  } else {
    rand_position.tmp <- sample(1:contigs.gen.inter$seq.length[i],1)
  }
  
  contigs.gen.inter$min.g.start[i] <- rand_position.tmp
  contigs.gen.inter$max.g.end[i] <- contigs.gen.inter$min.g.start[i] + contigs.gen.inter$insert.length[i] - 1
  
  rm(rand_position.tmp, gaps.this.chrom)
}

# Merging the table back together.
contigs.gen.new <- bind_rows(contigs.gen.exon, contigs.gen.pexon,
                             contigs.gen.intron, contigs.gen.inter)
# Recalculating the centromeric and telomeric distances.
contigs.gen.new <- contigs.gen.new %>%
  mutate(centr.dist = pmin(abs(min.g.start - centromer.center),
                           abs(centromer.center - max.g.end)),
         tel.dist = pmin(min.g.start, seq.length - max.g.end))
# Removing the centromeric and telomeric distance for the mtDNA.
contigs.gen$tel.dist[contigs.gen$chrom == "chrM"] <- NA
contigs.gen.new$tel.dist[contigs.gen.new$chrom == "chrM"] <- NA

# Merging the natural contigs with the in silico generated data.
contigs <- contigs.nat %>% bind_rows(contigs.gen.new)

write.table(contigs, file = "contigs_including_insilico_DATE.txt", quote = FALSE,
            sep = "\t", row.names = FALSE)
########## END OF THE SCRIPT
