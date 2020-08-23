#### Origin of inserted DNA.
#### This script visualizes the mapping of the insertions to the human genome
#### and outputs the Figure 3A.


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
# This function aids the stacking of the triangles to avoid the overplotting.
calculate_offsets <- function(df, method = 'middle'){
  # In "middle" all triangles are plotted pointing to the center of a bin.
  if (method == 'middle'){
    bins_center <- bins + binwidth/2
    df$AX <- bins_center[df$binX]
  }
  # Bins are numbered along a single chromosome, so the chromosome number should
  # be taken into consideration to avoid collisions.
  df <- unite(df, 'binXchr', c('binX', 'chrom'), sep = ";")
  # Counting insert-origins in each bin.
  df <- left_join(df, df %>% group_by(binXchr) %>%
                    summarise(cnt = n()), by = 'binXchr')
  # In "mean" triangle group is pointing to the mean center coordinate of the
  # insert-origins.
  if (method == 'mean'){
    df <- left_join(df, df %>% group_by(binXchr) %>%
                      summarise(meanAX = mean(AX)), by = 'binXchr')
    df$AX <- df$meanAX
    df <- dplyr::select(df, -'meanAX')
  }
  # Inside each bin the insert-origins are arranged by their IDs.
  df <- df %>% group_by(binXchr) %>% mutate(elementrank = min_rank(ins.id))
  # Calculating the position of the individual triangles for them to arrange
  # into an inverted pyramid.
  df$nlevel <- ceiling(-0.5 + 0.5*sqrt(8*df$cnt+1))
  df$level <- -0.5 + 0.5*sqrt(8*df$elementrank+1)
  df$levelr <- ceiling(df$level)
  df <- df %>% group_by(binXchr, levelr) %>%
    mutate(elementposinlevel = min_rank(level))
  df <- unite(df, 'binXchrLevelR', c('binXchr','levelr'), sep = ":")
  dfg <- df %>%
    group_by(binXchrLevelR) %>%
    summarise(elementsonlevel = max(elementposinlevel))
  df <- left_join(df, dfg, by = 'binXchrLevelR')  
  df <- separate(df, 'binXchrLevelR', c('binXchr','levelr'), sep = ":")
  df$levelr <- as.numeric(df$levelr)
  # For each triangle the coordinate is adjusted to fit into the pyramid.
  df$offsetY <- df$levelr - 1
  df$offsetX <- 2*df$elementposinlevel - df$elementsonlevel - 1
  df <- dplyr::select(df, -c('levelr','nlevel','level','elementsonlevel',
                             'elementposinlevel','elementrank','cnt'))
  df <- separate(df, 'binXchr', c('binX', 'chrom'), sep = ";")
  return(df)
}
# Simple unit transformation functions.
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}


# Data import -------------------------------------------------------------
contigs <- read.table(file = "contigs_DATE.txt",
                      sep = '\t', header = TRUE,
                      stringsAsFactors = FALSE)

# Importing the human genome assembly chromosome length from Genome Reference
# Consortium data (GRCh38.p13). Link to the resource:
# https://www.ncbi.nlm.nih.gov/grc/human/data
chrlen <- read.table(file = 'input/GRCh38_chrlength.txt',
                     sep = '\t', header = TRUE,
                     sstringsAsFactors = FALSE)
# Mapping the in vivo occuring insertions, omitting the in silico modelled
# contigs.
contigs <- contigs %>% dplyr::filter(ins.type.simple != "Generated")


# Data processing ---------------------------------------------------------
# Some contigs (immunoglobulin transcripts) contain inserts from different
# origins, which are mapped independently. In order to do that, the rows are
# split into individual inserts.
contigs$contig.id <- seq(nrow(contigs))
columnstosep <- c('insert.id.s.', 'insert.genomic.coord.s.', 'insert.contig.coord.s.','insert.seq.s.')
contigs <- separate_rows(contigs, columnstosep, sep = ',\\W*')
rm(columnstosep)
# Making an ID for each insert.
contigs$ins.id <- seq(nrow(contigs))
# The insertions of an undefined type (in which J segment was not located) are
# not plotted and omitted from the further analysis.
inserts <- inserts %>% dplyr::filter(ins.type.simple != "V-CH1")

inserts <- dplyr::select(inserts, sample.uid, insert.id.s., ins.id, chrom,
                   min.g.start, max.g.end, insert.length, ins.type.simple,
                   IsVandConstInFrame.withoutStop, contig.id, ins.class,
                  population.simple)
# Inserts originating from the undefined scaffolds and the mitochondrial genome
# are not plotted and omitted from the further analysis.
inserts %>% dplyr::filter(!str_detect(chrom, "_") & chrom != "chrM") -> inserts

# To aid the further plotting, the length of the longest chromosome is used as a
# scaling factor.
chrlen <- dplyr::select(chrlen, Chromosome, Length.bp)
chrlen$Length.bp <- as.numeric(str_replace_all(chrlen$Length.bp, ',', ''))
chrlen$relative.length <- chrlen$Length.bp/max(chrlen$Length.bp)
chrlen$chromosome <- str_c('chr', chrlen$Chromosome)

# Using the scaling factor to determine the relative position of the insert on
# the chromosome.
inserts <- inserts %>%
  mutate(ins.gen.center = (min.g.start + max.g.end)/2,
         relative.pos = ins.gen.center/max(chrlen$Length.bp))

# Plotting parameters -----------------------------------------------------
# Canvas parameters
plot.height <- 1500
plot.width <- 1000
x.margin <- 50
y.margin <- 10
# Chromosomes rectangles thickness
thick <- 0 # geom_line used, the thickness is defined inside the ggplot()
# Triangles parameters
tw <- 6 # Width of a triangle.
th <- 10 # Height of a trianlge.
# 
# # Plot characteristics
# plot.params <- list(height = 1500, width = 1000, x.margin = 50,
#                     y.margin = 10, thick = 0, trianglewidth = 6,
#                     triangleheight = 10)

# Dependent Parameters
x.inner <- plot.width - x.margin
y.inner <- plot.height - y.margin
inner.width <- plot.width - 2*x.margin
inner.height <- plot.height - 2*y.margin
N <- nrow(chrlen)
spacing <- inner.height/(N-1)
chrlen$element.i <- seq(N)

# Points coordinates calculation
chrlen$AX <- x.margin
chrlen$AY <- y.inner - spacing*(chrlen$element.i - 1) - thick/2
chrlen$BX <- x.margin + inner.width*chrlen$relative.length
chrlen$BY <- y.inner - spacing*(chrlen$element.i - 1) + thick/2


# Triangles coordinates calculation ---------------------------------------
# Triangle AY coordinate is defined by the chromosome the insert originates
# from.
inserts <- left_join(inserts, dplyr::select(chrlen, chromosome, AY),
                     by = c("chrom" = "chromosome"))
# Triangle AX coordinate is defined by the position on the chromosome.
inserts$AX <- x.margin + inner.width*inserts$relative.pos

# The following conditional block defines the stacking of the triangles. If set
# to TRUE, triangles will stack upon each other. If FALSE, triangle will overlap
# with each other and be on the same level.
dotrianglesrepel = TRUE
if (dotrianglesrepel == TRUE){
  # By default, triangles do not overlap. Adjust this parameter to allow the
  # overlap.
  allowed.overlap = 0
  # Dissecting the chromosome into bins to sort the triangles into them.
  binwidth = tw - allowed.overlap
  bins <- seq(0, inner.width + x.margin, binwidth)
  if (bins[length(bins)] != inner.width + x.margin){
    bins <- c(bins, inner.width + x.margin)
  }
  # Binning.
  inserts$binX <- .bincode(inserts$AX, bins)
  # Offsetting.
  inserts <- calculate_offsets(inserts, method = 'middle')
  inserts$AY <- inserts$AY + th*inserts$offsetY + 2
  inserts$AX <- inserts$AX + (tw/2)*inserts$offsetX
  inserts <- dplyr::select(inserts, -c(offsetY, offsetX))
}

# B and C coordinate of the triangle vertices are calculated based on the A
# coordinate, tw and th.
inserts$BX <- inserts$AX - tw/2
inserts$BY <- inserts$AY + th
inserts$CX <- inserts$AX + tw/2
inserts$CY <- inserts$BY
inserts <- unite(inserts, A, AX, AY, sep = ';', remove = TRUE)
inserts <- unite(inserts, B, BX, BY, sep = ';', remove = TRUE)
inserts <- unite(inserts, C, CX, CY, sep = ';', remove = TRUE)
inserts <- gather(inserts, dot, dot.coord, A:C)
inserts <- separate(inserts, dot.coord, c('dotX', 'dotY'), sep = ';', remove = TRUE)
inserts$dotX <- as.numeric(inserts$dotX)
inserts$dotY <- as.numeric(inserts$dotY)

# Adding centromeres ------------------------------------------------------
# Centromere coordinates are imported from UCSC Table Browser, Gaps section.
gaps <- read.table(file = 'input/hg38_centromeres.txt',
                   sep = '\t',header = TRUE,
                   stringsAsFactors = FALSE)
# The most conservative (the widest) coordinate is used.
gaps <- gaps %>%
  dplyr::select(-c("bin", "name")) %>%
  group_by(chrom) %>%
  summarise(centromer.start = min(chromStart),
            centromer.end = max(chromEnd)) %>%
  mutate(centromer.center = (centromer.start+centromer.end)/2)
centromeres <- gaps
rm(gaps)

centromeres <- dplyr::select(centromeres, chrom, centromer.center)
centromeres$Chromosome <- str_replace_all(centromeres$chrom, 'chr','')
# Merging the chromosome length data with centromeres coordinate.
chrlen <- left_join(chrlen, centromeres, by = "Chromosome", copy = FALSE)
rm(centromeres)
chrlen$chromCenter.relative <- chrlen$centromer.center/max(chrlen$Length.bp)
chrlen$chromCenterY <- chrlen$AY
chrlen$chromCenterX <- x.margin + inner.width*chrlen$chromCenter.relative

chrlen$Length.Mbp <- chrlen$Length.bp / 1000000

chrlen <- dplyr::filter(chrlen, !Length.Mbp < 1)
chrlen$Length.Mbp <- round(chrlen$Length.Mbp, 0)
chrlen$Length.Mbp <- as.character(chrlen$Length.Mbp)
chrlen$Length.Mbp <- paste(chrlen$Length.Mbp, 'Mbp')

color_of_lines <- 'black'

# Plotting ----------------------------------------------------------------
# Main plotting - the chromosomes and the triangles, representing the insert
# origins.
insert.origin.plot <- ggplot(data = chrlen)+
  coord_cartesian(xlim = c(0, plot.width), ylim = c(0, plot.height))+
  geom_segment(aes(x = AX, y = AY, xend = BX, yend = BY), alpha = 1, color = color_of_lines, size = 1)+
  geom_point(aes(x = chromCenterX, y = chromCenterY), size = 2, shape = 21, fill = 'deeppink2')+
  geom_text(aes(x = AX-35, y = AY+3, label = chrom), colour = color_of_lines)+
  geom_text(aes(x = BX + 50, y = AY+3, label = Length.Mbp), size = 3, colour = color_of_lines)+
  geom_polygon(data = inserts,
               aes(x=dotX, y=dotY, group = contig.id, fill = ins.type.simple), alpha = .8)+

  scale_fill_manual(values = c("VDJ" = rgb(0,0,255, maxColorValue = 255),
                                "J-CH1" = rgb(255,0,0, maxColorValue = 255),
                                "V-CH1" = rgb(150,150,150, maxColorValue = 255)))+
  theme_void()


# Adding the frames of the hotspots ---------------------------------------
# Plotting the hotspots on the same canvas.
hotspots <- read.table("hs_DATE.txt",
                       sep = "\t", header = TRUE,
                       stringsAsFactors = FALSE)
# Arranging the hotspots by the number of contigs in them to plot only the
# top-10.
hotspots <- hotspots %>%
  arrange(desc(contigs.n)) %>%
  mutate(AX = x.margin + inner.width*min.g.start/
           max(chrlen$Length.bp),
         BX = x.margin + inner.width*max.g.end/
           max(chrlen$Length.bp)) %>%
  left_join(dplyr::select(chrlen, chromosome, AY),
            by = c("chrom" = "chromosome")) %>%
  mutate(AY = AY - 10, BY = AY + 5)
# Showing only the top-10 hotspots.
insert.origin.plot <- insert.origin.plot +
  geom_rect(data = hotspots[1:10,],
            aes(xmin = AX, ymin = AY, xmax = BX, ymax = BY),
            fill = "violetred")

# Exporting the full plot into SVG.
ggsave("insert_origin_map_DATE.svg", plot = insert.origin.plot,
       device = "svg", width = 15, height = 10)
################################################### END OF THE SCRIPT ###---
