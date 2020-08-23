#### Origin of inserted mitochondrial DNA
#### This script visualizes the mapping of the insertions to the human
#### mitochondrial genome. and outputs the Figure 3E.


# Working folder ----------------------------------------------------------
setwd("C:/Users/User/Documents/Lab/LongBCRs/bioinformatics/lair")


# Libraries ---------------------------------------------------------------
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggrepel)
library(reshape2)
library(Biostrings)
library(ggforce)


# Functions ---------------------------------------------------------------
# Counting overlaps between a group of fragments and a single fragment.
count.overlap <- function(subject.starts, subject.ends, query.start, query.end){
  subject.starts <- as.numeric(subject.starts)
  subject.ends <- as.numeric(subject.ends)
  query.start <- as.numeric(query.start)
  query.end <- as.numeric(query.end)
  if(sum(subject.starts > subject.ends, na.rm = TRUE)){
    stop("Some of the subject ranges you provided are flipped")
  }
  if(query.start > query.end){
    stop("Query range you provided is flipped")
  }
  count.overlap <- sum(query.start >= subject.starts & query.start <= subject.ends|
                         query.end >= subject.starts & query.end <= subject.ends|
                         query.start <= subject.starts & query.end >= subject.ends)

  return(count.overlap)
}


# Data import -------------------------------------------------------------
contigs <- read.table(file = "contigs_DATE.txt",
                      sep = '\t', header = TRUE,
                      stringsAsFactors = FALSE)


# Mitochondrial inserts as the segments -----------------------------------
# Focusing on the mitochondrial insertions that occur in vivo and are of defined
# type.
contigs.mt <- contigs %>%
  dplyr::filter(chrom == "chrM" & !ins.type %in% c("Generated", "V-CH1"))
# Assigning an ID to each contig.
contigs.mt$contig.id.mt <- seq(1, nrow(contigs.mt))

# Some contigs (immunoglobulin transcripts) contain inserts from different
# origins, which are mapped independently. In order to do that, the rows are
# split into individual inserts.
contigs.mt.ins <- separate_rows(contigs.mt, "insert.genomic.coord.s.", sep = ",\\W*")
# Mitochondrial genome length from Piovesan et al., 2019.
chrMlen <- 16569
# Extracting the numerical coordinates from the formatted text coordinate.
contigs.mt.ins <- contigs.mt.ins %>%
  separate(insert.genomic.coord.s., c("chrom", "coord"), sep = ":") %>%
  separate(coord, c("start", "end"), sep = "\\..") %>%
  dplyr::filter(chrom == "chrM")
# Transforming the linear coordinates into circular.
contigs.mt.ins <- contigs.mt.ins %>%
  mutate(start = as.integer(start), end = as.integer(end),
         pi.start = start*2*pi/chrMlen, pi.end = end*2*pi/chrMlen,
         y = contig.id.mt*6)
# In some transcripts the mtDNA insertions underwent splicing. It is plotted as
# the dashed lines, coordinates for which are calculated here.
contigs.mt.ins.agg <- contigs.mt.ins %>% group_by(contig.id.mt) %>%
  summarise(min.coord = min(start), max.coord = max(end), ins.count = n(),
            y = nth(y, 1)) %>%
  mutate(pi.start = min.coord*2*pi/chrMlen, pi.end = max.coord*2*pi/chrMlen)
# The segments that span 16569-0 point are handled in a special way to plot them
# correctly.
contigs.mt.ins.agg$cross.zero <- FALSE
contigs.mt.ins.agg$cross.zero[(contigs.mt.ins.agg$max.coord -
                                 contigs.mt.ins.agg$min.coord) > chrMlen/2] <- TRUE
contigs.cross.zero <-
  contigs.mt.ins.agg$contig.id.mt[contigs.mt.ins.agg$cross.zero == TRUE]
contigs.mt.ins.agg.zero <- contigs.mt.ins %>%
  dplyr::filter(contig.id.mt %in% contigs.cross.zero) %>%
  group_by(contig.id.mt) %>%
  summarise(min.coord = max(start), max.coord = min(end), ins.count = n(),
            y = nth(y, 1)) %>%
  mutate(pi.start = (min.coord-chrMlen)*2*pi/chrMlen,
         pi.end = max.coord*2*pi/chrMlen)
contigs.mt.ins.agg <- contigs.mt.ins.agg %>%
  dplyr::filter(!contig.id.mt %in% contigs.cross.zero) %>%
  bind_rows(contigs.mt.ins.agg.zero) %>%
  mutate(contig.id.mt = as.character(contig.id.mt))
contigs.mt.ins <- contigs.mt.ins %>%
  mutate(contig.id.mt = as.character(contig.id.mt))
# The insertions are arranged according to their end coordinate.
contigs.mt.ins.agg <- contigs.mt.ins.agg %>% arrange(desc(max.coord))
# To avoid overlay, the overlap is detected and the segments are separated into
# different levels.
contigs.mt.ins.agg$level <- 0
LVL <- 1

contigs.mt.ins.agg <- contigs.mt.ins.agg %>%
  mutate(min.coord.1 = pmin(min.coord, max.coord),
         max.coord.1 = pmax(min.coord, max.coord),
         min.coord = min.coord.1,
         max.coord = max.coord.1) %>%
  dplyr::select(-c("min.coord.1", "max.coord.1"))

# Sorting into levels.
while (sum(contigs.mt.ins.agg$level == 0) > 0) {
  contigs.mt.ins.agg$level[which(contigs.mt.ins.agg$level == 0)[1]] <- LVL
  for (j in which(contigs.mt.ins.agg$level == 0)) {
    if (count.overlap(subject.starts =
                     contigs.mt.ins.agg$min.coord[contigs.mt.ins.agg$level == LVL],
                     subject.ends =
                     contigs.mt.ins.agg$max.coord[contigs.mt.ins.agg$level == LVL],
                     query.start =
                     contigs.mt.ins.agg$min.coord[j],
                     query.end =
                     contigs.mt.ins.agg$max.coord[j]) == 0) {
      contigs.mt.ins.agg$level[j] <- LVL
    }
  }
  LVL <- LVL + 1
}

contigs.mt.ins <- contigs.mt.ins %>%
  left_join(dplyr::select(contigs.mt.ins.agg, contig.id.mt, level),
            by = "contig.id.mt")


# Plotting parameters.
plot.height <- 1000
plot.width <- 1000
X0 <- 500
Y0 <- 500
circle.rad <- 200

# Insertions plotted as "arcs".
mt.plot <- ggplot()+
  coord_cartesian(xlim = c(0, plot.width), ylim = c(0, plot.height))+
  geom_arc(aes(x0 = X0, y0 = Y0, r = circle.rad,
               start = 0, end = 2*pi), color = "black") +
  geom_arc(data = contigs.mt.ins.agg %>% dplyr::filter(ins.count > 1),
           aes(x0 = X0, y0 = Y0, r = circle.rad+level*4,
               start = pi.start, end = pi.end),
           color = "black",
           alpha = 1.0, size = 0.3, linetype = "dashed")+
  geom_arc(data = contigs.mt.ins,
           aes(x0 = X0, y0 = Y0, r = circle.rad+level*4,
               start = pi.start, end = pi.end),
           color = ins.type,
           alpha = 1.0, size = 0.7, lineend = "square")+
  scale_color_manual(values = c("VDJ" = rgb(0,0,1),
                                "J-CH1" = rgb(1,0,0),
                                "V-CH1" = rgb(0.5,0.5,0.5)))+
  theme_void()+
  theme(legend.position = "none")



# Mitochondrial genes overlay ---------------------------------------------
# Adding the mitochondrial genes on the map.
# Coordinates from Genecards.
mt.genes <- read.table("mt_genes.txt", sep = "\t",
                       stringsAsFactors = FALSE, header = TRUE)
# Plotting the genes on two layers.
mt.genes$id <- rep(c(1,2),ceiling(nrow(mt.genes)/2))[1:nrow(mt.genes)]
# Calculating the circular coordinates.
mt.genes <- mt.genes %>%
  mutate(pi.start = start*2*pi/chrMlen, pi.end = end*2*pi/chrMlen)
# Adding the D-loop feature.
mt.genes$pi.start[mt.genes$feature == "D-loop"] <-
  (mt.genes$start[mt.genes$feature == "D-loop"] - chrMlen)*2*pi/chrMlen

mt.plot <- mt.plot +
  geom_arc(data = mt.genes,
           aes(x0 = X0, y0 = Y0, r = circle.rad+10*id,
               start = pi.start, end = pi.end),
           color = "grey", alpha = 1.0, size = 2, lineend = "square")

# Exporting the full plot into SVG.
ggsave("mt_plot_DATE.svg", plot = mt.plot, width = 10,
       height = 10, device = "svg")

################################################### END OF THE SCRIPT ###---
