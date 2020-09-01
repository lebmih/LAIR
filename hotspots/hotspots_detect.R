#### Clustering of the inserts

# Working folder ----------------------------------------------------------
# DELAR11 PC or LEBMIH ASUS
if (Sys.info()["nodename"] == "DELAR11") {
  print("Hi, Mikhail Lebedin!")
  setwd("C:/Users/mlebedi/Documents/bioinformatics/")
} else if (Sys.info()["nodename"] == "LEBMIHASUS") {
  print("Hi, lebmih!")
  setwd("C:/Users/User/Documents/Lab/LongBCRs/bioinformatics/lair")
}

source("wrapper.R")

# Creating report files ---------------------------------------------------
destination.folder <- dirs$reports
dir.create(destination.folder)

# Libraries ---------------------------------------------------------------
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
#library(ggrepel)
#library(reshape2)
#library(Biostrings)
#library(scales)
#library(PtProcess)
#library(gtable)
#library(grid)
#library(gridExtra)
#library(ggpubr)
#library(lebutils)
#library(ggforce)

# Report and constants ----------------------------------------------------
# Set plotbitte to FALSE if you don't want plots, reportbitte to FALSE if you
# don't want reports
constants <- list(plotbitte = TRUE, create.report.and.check = TRUE)
if (!constants$create.report.and.check) {
  print("Alright, no reports.")
} else {
  report.file <- paste(dirs$reports, "/donors_stat_report_", format(Sys.time(),"%d%m%y_%H%M"), ".txt", sep = "")
  write("The report contains the report of donors_general_stat.R script",
        file = report.file, append = TRUE)
  check.file <- paste(destination.folder, "/checkfile_", format(Sys.time(),"%d%m%y_%H%M"), ".txt", sep = "")
  write("This file contains issues that are to be checked manually",
        file = check.file, append = TRUE)
}

# Data import -------------------------------------------------------------
# First mode - used in the paper.
# This mode calculates the distances to the closest neighbour and outputs the
# hotspot if the inserts are closer than 1 Mb to each other.
# Importing the latest table.
contigs <- read.table("contigs_170620.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# Analyzing only the naturally occurring inserts.
contigs <- contigs %>%
  dplyr::filter(ins.type.simple %in% c("VDJ", "Switch", "V-CH1") &
           natural %in% c("primary", "spike-in", "ebv", "exhausted"))
# Keeping the multiple insertion contigs separate.
contigs.mult <- contigs %>% dplyr::filter(is.na(min.g.start) | is.na(max.g.end))
contigs <- contigs %>% dplyr::filter(!is.na(min.g.start) & !is.na(max.g.end))
# Arranging the insertions by chromosome and the coordinate to calculate the
# distance to the closest neighbour.
contigs <- contigs %>% arrange(chrom, min.g.start)
# Assigning the groups based on the distance. The group ID changes when the
# neighbour is further away than 1 Mb.
group.id <- 1
# The group ID is what in the paper will be called "hsX" (X = group ID).
contigs$group.id <- NA
contigs$group.id[1] <- 1
for (i in seq(2, nrow(contigs))) {
  if (is.na(contigs$min.g.start[i]) | is.na(contigs$max.g.end[i])) {
    next
  }
  if (contigs$chrom[i] == contigs$chrom[i-1]){
    if (min(abs(contigs$min.g.start[i] - contigs$max.g.end[i-1]),
            abs(contigs$min.g.start[i-1] - contigs$max.g.end[i]),
            abs(contigs$min.g.start[i] - contigs$min.g.start[i-1]),
            abs(contigs$max.g.end[i] - contigs$max.g.end[i-1])) < 1000000){
      
    } else {
      group.id <- group.id + 1
    }
  } else {
    group.id <- group.id + 1
  }
  contigs$group.id[i] <- group.id
}
# Now group ID is the same for the inserts originating from the same cluster.
# The next step is to summarise the info about the clusters to determine the
# hotspots.
# The summary info:
#   contigs.n - how many contigs bear the insertion from the cluster
#   chrom - the chromosome name
#   min.g.start - genomic coordinate where the cluster starts*
#   max.g.end - genomic coordinate where the cluster ends*
#   genes - genes donating the insertions inside the cluster
#   donors.n - number of donors which bear the insertions from the cluster
#   types - insertions of which types originate from the cluster
#   classes - insertions of which classes originate from the cluster
#   tel.dist - the minimal distance to the telomeres
#   length - the size of the cluster
# * - there is no way to compute the exact border of the actual hotspot, so this
# measure is a heuristic based on the existing data.
clusters <- contigs %>% group_by(group.id) %>%
  summarise(contigs.n = n(),
            chrom = nth(chrom, 1),
            min.g.start = min(min.g.start),
            max.g.end = max(max.g.end),
            genes = paste(unique(gene.name.single), collapse = ";"),
            donors.n = n_distinct(donor),
            types = paste(unique(ins.type.simple), collapse = ";"),
            classes = paste(unique(ins.class), collapse = ";"),
            tel.dist = min(tel.dist)) %>%
  mutate(length = max.g.end - min.g.start)
# Keeping the mitochondrial DNA away, as the 1 Mb threshold was defined for
# nuclear DNA.
clusters <- clusters %>% dplyr::filter(chrom != "chrM")

# Testing for the correlation between the telomeric distance and the number of
# inserts.
cor(clusters$contigs.n, clusters$tel.dist)
# The correlation between telomeric proximity and number of inserts originating
# from the hotspot is -0.16.
ggplot(clusters)+
  geom_point(aes(x = log10(tel.dist+1), y = contigs.n))+
  theme_classic()



# R-loops overlap ---------------------------------------------------------
# Importing the data
rloops <- read.table("rloops/K562_DRIP_peaks_hg38.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(rloops) <- c("chrom", "start", "end")
rloops <- rloops %>% dplyr::filter(!str_detect(chrom, "_"))
# Computing genomic center in order to calculate the center-to-center distance
clusters <- clusters %>% mutate(g.center = (min.g.start + max.g.end)/2)
rloops <- rloops %>% mutate(center = (start + end)/2)

# Calculating the distances and the overlap
dist.to.closest.rloops <- c()
overlap.closest.rloops <- c()

for (i in seq(nrow(clusters))) {
  rloops.this.chrom <- rloops %>% dplyr::filter(chrom == clusters$chrom[i])
  rloops.dist <- min(abs(rloops.this.chrom$center - clusters$g.center[i]))
  rloops.overlap <- count.overlap(rloops.this.chrom$start,
                                  rloops.this.chrom$end,
                                  clusters$min.g.start[i],
                                  clusters$max.g.end[i])
  dist.to.closest.rloops <- c(dist.to.closest.rloops, rloops.dist)
  overlap.closest.rloops <- c(overlap.closest.rloops, rloops.overlap)
  if (i %in% seq(0,nrow(clusters),100)) {
    paste0(round(i*100/nrow(clusters)), "% of the job is done") %>% print()
  }
}
# Control: this number should be 0.
sum(is.infinite(dist.to.closest.rloops))
# Appending the table with the computed values.
clusters$rloops.dist <- dist.to.closest.rloops
clusters$rloops.overlap <- overlap.closest.rloops
# Removing unnecessary vectors and variables.
rm(dist.to.closest.rloops, overlap.closest.rloops, rloops.this.chrom, i, rloops.dist, rloops.overlap)
# Checking the correlation.
cor(clusters$rloops.overlap, clusters$contigs.n)
# There is a substantial correlation - 0.52, though it can be attributed to the
# size of the hotspots - they are more likely to overlap any region by chance
# due to their length.
ggplot(clusters)+
  geom_point(aes(x = contigs.n, y = rloops.overlap), shape = 21, fill = "grey")+
  theme_classic()


# CFS overlap -------------------------------------------------------------
files.path <- "~/Lab/LongBCRs/MethodPaper/CommonFragileSites/fragile_site_gene_bed/"
files <- list.files(path = files.path, pattern = ".bed", full.names = TRUE)
cfs.genes <- data.frame()

for (i in files) {
  cfs.temp <- read.table(file = i,
                         sep = '\t', header = FALSE, stringsAsFactors = FALSE,
                         skip = 3)
  if (is.na(cfs.temp$V2) %>% sum() == length(cfs.temp$V2)){
    cfs.temp %>% dplyr::select(-V2) -> cfs.temp
  }
  colnames(cfs.temp) <- c("chrom","start","end","name","score","strand")
  cfs.genes <- bind_rows(cfs.genes, cfs.temp)
}

rm(cfs.temp, i, files, files.path)

cfs.genes$chrom[cfs.genes$chrom == 'chrx'] <- "chrX"
cfs.genes$chrom[cfs.genes$chrom == 'chry'] <- "chrY"

cfs.genes$start[cfs.genes$name == "NRIp1"] <- 14961235
cfs.genes$end[cfs.genes$name == "NRIp1"] <- 15065936

cfs.genes <- cfs.genes %>% dplyr::filter(!is.na(start))
cfs.genes <- cfs.genes %>% group_by(name) %>%
  summarise(chrom = nth(chrom, 1),
            start = min(start),
            end = max(end),
            score = nth(score, 1),
            strand = nth(strand, 1))
cfs.genes <- cfs.genes %>% ungroup()
cfs.genes$chrom[cfs.genes$chrom == "chromosome2"] <- "chr2"

cfs.genes <- cfs.genes %>% mutate(center = (start + end)/2)

# Calculating the distances
dist.to.closest.cfs <- c()
overlap.closest.cfs <- c()

for (i in seq(nrow(clusters))) {
  cfs.this.chrom <- cfs.genes %>% dplyr::filter(chrom == clusters$chrom[i])
  if(nrow(cfs.this.chrom) == 0){
    dist.to.closest.cfs <- c(dist.to.closest.cfs, NA)
    overlap.closest.cfs <- c(overlap.closest.cfs, NA)
    next
  }
  cfs.dist <- min(abs(cfs.this.chrom$center - clusters$g.center[i]))
  cfs.overlap <- count.overlap(cfs.this.chrom$start,
                               cfs.this.chrom$end,
                               clusters$min.g.start[i],
                               clusters$max.g.end[i])
  dist.to.closest.cfs <- c(dist.to.closest.cfs, cfs.dist)
  overlap.closest.cfs <- c(overlap.closest.cfs, cfs.overlap)
  if (i %in% seq(0,nrow(clusters),100)) {
    paste0(round(i*100/nrow(clusters)), "% of the job is done") %>% print()
  }
}
# Control: this number should be 0
sum(is.infinite(dist.to.closest.cfs))
# Appending the table with the computed values.
clusters$cfs.dist <- dist.to.closest.cfs
clusters$cfs.overlap <- overlap.closest.cfs
# Removing unnecessary vectors and variables.
rm(dist.to.closest.cfs, overlap.closest.cfs, cfs.this.chrom, i, cfs.dist, cfs.overlap)
# Checking the correlation.
cor(clusters$cfs.overlap[!is.na(clusters$cfs.overlap)],
    clusters$contigs.n[!is.na(clusters$cfs.overlap)])
# There is a correlation - 0.24, though it can be attributed to the size of the
# hotspots - they are more likely to overlap any region by chance due to their
# length.
ggplot(clusters)+
  geom_point(aes(x = contigs.n, y = cfs.overlap), shape = 21, fill = "grey")+
  theme_classic()

# ERFS overlap ------------------------------------------------------------
erfs <- read.table("input/erfs_barlow_hg38.bed", header = FALSE, sep = "\t", stringsAsFactors = FALSE)
colnames(erfs) <- c("chrom", "start", "end")
erfs <- erfs %>% mutate(center = (start + end)/2)
# Calculating the distances
dist.to.closest.erfs <- c()
overlap.closest.erfs <- c()

for (i in seq(nrow(clusters))) {
  erfs.this.chrom <- erfs %>% dplyr::filter(chrom == clusters$chrom[i])
  if(nrow(erfs.this.chrom) == 0){
    dist.to.closest.erfs <- c(dist.to.closest.erfs, NA)
    overlap.closest.erfs <- c(overlap.closest.erfs, NA)
    next
  }
  erfs.dist <- min(abs(erfs.this.chrom$center - clusters$g.center[i]))
  erfs.overlap <- count.overlap(erfs.this.chrom$start,
                                erfs.this.chrom$end,
                                clusters$min.g.start[i],
                                clusters$max.g.end[i])
  dist.to.closest.erfs <- c(dist.to.closest.erfs, erfs.dist)
  overlap.closest.erfs <- c(overlap.closest.erfs, erfs.overlap)
  if (i %in% seq(0,nrow(clusters),100)) {
    paste(round(i*100/nrow(clusters)), "% of the job is done", sep = "") %>% print()
  }
}
# Control: this number should be 0
sum(is.infinite(dist.to.closest.erfs))
# Appending the table with the computed values.
clusters$erfs.dist <- dist.to.closest.erfs
clusters$erfs.overlap <- overlap.closest.erfs
# Removing unnecessary vectors and variables.
rm(dist.to.closest.erfs, overlap.closest.erfs, erfs.this.chrom, i,
   erfs.overlap, erfs.dist)
# Checking the correlation.
cor(clusters$erfs.overlap[!is.na(clusters$erfs.overlap)],
    clusters$contigs.n[!is.na(clusters$erfs.overlap)])
# There is a correlation - 0.26, though it can be attributed to the size of the
# hotspots - they are more likely to overlap any region by chance due to their
# length.
ggplot(clusters)+
  geom_point(aes(x = contigs.n, y = log10(erfs.dist+1)), shape = 21, fill = "grey")+
  theme_classic()



# Saving the clusters table -----------------------------------------------
write.table(clusters, file = "clusters_300620.txt", quote = FALSE,
            sep = "\t", row.names = FALSE)
# Extracting only the clusters found in more than one donor
clusters <- clusters %>% dplyr::filter(donors.n > 1)
write.table(clusters, file = "hotspots_300620.txt", quote = FALSE,
            sep = "\t", row.names = FALSE)

#### END OF THE SCRIPT
