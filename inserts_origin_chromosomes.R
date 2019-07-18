#### Origin of inserted DNA

## Working directory
setwd("~/Lab/LongBCRs/MethodPaper/RunTables")


## Libraries
x <- c('tidyr', 'dplyr', 'ggplot2', 'ggrepel', 'reshape2', 'stringr', 'Biostrings')
lapply(x, require, character.only = TRUE)
rm(x)

## Functions
calculate_offsets <- function(df, method = 'middle'){
  if (method == 'middle'){
    bins_center <- bins + binwidth/2
    df$AX <- bins_center[df$binX]
  }
  
  df <- unite(df, 'binXchr', c('binX', 'chromosome'), sep = ";")
  
  df <- left_join(df, df %>% group_by(binXchr) %>% summarise(cnt = n()), by = 'binXchr')
  
  if (method == 'mean'){
    df <- left_join(df, df %>% group_by(binXchr) %>% summarise(meanAX = mean(AX)), by = 'binXchr')
    df$AX <- df$meanAX
    df <- select(df, -'meanAX')
  }
  
  df <- df %>% group_by(binXchr) %>% mutate(elementrank = min_rank(ins.id))
  
  df$nlevel <- ceiling(-0.5 + 0.5*sqrt(8*df$cnt+1))
  df$level <- -0.5 + 0.5*sqrt(8*df$elementrank+1)
  df$levelr <- ceiling(df$level)
  
  df <- df %>% group_by(binXchr, levelr) %>% mutate(elementposinlevel = min_rank(level))
  
  df <- unite(df, 'binXchrLevelR', c('binXchr','levelr'), sep = ":")
  dfg <- df %>%
    group_by(binXchrLevelR) %>%
    summarise(elementsonlevel = max(elementposinlevel))
  df <- left_join(df, dfg, by = 'binXchrLevelR')  
  df <- separate(df, 'binXchrLevelR', c('binXchr','levelr'), sep = ":")
  df$levelr <- as.numeric(df$levelr)
  
  df$offsetY <- df$levelr - 1
  df$offsetX <- 2*df$elementposinlevel - df$elementsonlevel - 1
  df <- select(df, -c('levelr','nlevel','level','elementsonlevel','elementposinlevel',
                      'elementrank','cnt'))
  df <- separate(df, 'binXchr', c('binX', 'chromosome'), sep = ";")
  return(df)
}


### Body
## Importing the table
all_inserts <- read.table(file = 'inserts_cleanest_06052019.txt',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Widening the table
all_inserts$contig.id <- seq(nrow(all_inserts))
columnstosep <- c('insert.id.s.', 'insert.genomic.coord.s.', 'insert.contig.coord.s.','insert.seq.s.')
all_inserts <- separate_rows(all_inserts, columnstosep, sep = ',\\W*')

# Make an ID for each insert
all_inserts$ins.id <- seq(nrow(all_inserts))

## Filtering out LAIR1 insertions
all_inserts <- filter(all_inserts, !str_detect(all_inserts$Gene.id.s., 'LAIR1'))

## Extracting a piece of data as input for the function
## if you want to drop V-CH1
inserts <- all_inserts[!is.na(all_inserts$insert.type),]
## if you want to show V-CH1
#inserts <- all_inserts

inserts <- select(inserts, run.id, sample.id, insert.id.s., ins.id,
                  insert.genomic.coord.s., insert.length, insert.type,
                  IsVandConstInFrame.withoutStop, contig.id)

inserts <- separate(inserts, insert.genomic.coord.s.,
                     c('chromosome',
                       'insert.genomic.coord.s.'),
                     sep = ':')

inserts <- separate(inserts, insert.genomic.coord.s.,
                     c('insert.genomic.coord.s.start',
                       'insert.genomic.coord.s.end'),
                     sep = '-')

#inserts$chromosome <- str_replace_all(inserts$chromosome, 'chr','')
#inserts$chromosome <- as.numeric(inserts$chromosome)


## Filtering out non-valid chromosomes
filteredchromosomes <- c("chrM","chrUn_gl000220","chrUn_gl000214","chr7_gl000195_random")
inserts <- filter(inserts, !chromosome %in% filteredchromosomes)

## Calculating the inserts center
inserts$insert.genomic.coord.s.start <- as.numeric(inserts$insert.genomic.coord.s.start)
inserts$insert.genomic.coord.s.end <- as.numeric(inserts$insert.genomic.coord.s.end)
inserts$insert.genomic.center <- (inserts$insert.genomic.coord.s.start+inserts$insert.genomic.coord.s.end)/2



## 27.02.2019
# Unbiasing the data
inserts1 <- inserts %>% group_by(contig.id) %>% summarise(run.id = nth(run.id, 1),
                                                          sample.id = nth(sample.id, 1),
                                                          insert.id.s. = nth(insert.id.s., 1),
                                                          ins.id = nth(ins.id, 1),
                                                          chromosome = nth(chromosome, 1),
                                                          insert.genomic.coord.s.start = mean(insert.genomic.coord.s.start),
                                                          insert.genomic.coord.s.end = mean(insert.genomic.coord.s.end),
                                                          insert.genomic.center = mean(insert.genomic.center),
                                                          insert.length = mean(insert.length),
                                                          insert.type = nth(insert.type, 1),
                                                          IsVandConstInFrame.withoutStop = nth(IsVandConstInFrame.withoutStop, 1))

inserts <- inserts1

### Let's look closely to last six inserts in the first chromosome
#inserts <- filter(inserts, chromosome == 'chr2')
#inserts <- inserts[order(-inserts$insert.genomic.center),]
#inserts <- inserts[1:10,]

## Importing the chromosome length data for hg19
chrlen <- read.table(file = 'hg19chrlength.txt', sep = '\t', header = TRUE)
chrlen <- select(chrlen, Chromosome, Length.bp)
chrlen$Length.bp <- as.numeric(str_replace_all(chrlen$Length.bp, ',', ''))
chrlen$relative.length <- chrlen$Length.bp/max(chrlen$Length.bp)
chrlen$chromosome <- str_c('chr', chrlen$Chromosome)

## Calculating the relative position of the insert on the chromosome
inserts$relative.pos <- inserts$insert.genomic.center/max(chrlen$Length.bp)


### PLOT
## Input Parameters
# Plot
plot.height <- 1500
plot.width <- 1000
x.margin <- 50
y.margin <- 10
# Chromosomes rectangles thickness
thick <- 0
# Triangles parameters
tw <- 6
th <- 10

## Dependent Parameters
x.inner <- plot.width - x.margin
y.inner <- plot.height - y.margin
inner.width <- plot.width - 2*x.margin
inner.height <- plot.height - 2*y.margin
N <- nrow(chrlen)
spacing <- inner.height/(N-1)
chrlen$element.i <- seq(N)

## Points coordinates calculation
chrlen$AX <- x.margin
chrlen$AY <- y.inner - spacing*(chrlen$element.i - 1) - thick/2
chrlen$BX <- x.margin + inner.width*chrlen$relative.length
chrlen$BY <- y.inner - spacing*(chrlen$element.i - 1) + thick/2


## Triangles dots
# triangle AY
inserts <- left_join(inserts, select(chrlen, chromosome, AY), by = 'chromosome')
# triangle AX
inserts$AX <- x.margin + inner.width*inserts$relative.pos

# If set to TRUE triangles will stack upon each other
# If FALSE triangle will overlap with each other and be on the same level
dotrianglesrepel = TRUE

if (dotrianglesrepel == TRUE){
  ## Triangles Repel
  # First plots were made with tw/3
  allowed.overlap = 0
  
  # Bins generation
  binwidth = tw - allowed.overlap
  bins <- seq(0, inner.width + x.margin, binwidth)
  if (bins[length(bins)] != inner.width + x.margin){
    bins <- c(bins, inner.width + x.margin)
  }
  
  # Binning
  inserts$binX <- .bincode(inserts$AX, bins)
  
  # Offsetting
  inserts <- calculate_offsets(inserts, method = 'middle')
  inserts$AY <- inserts$AY + th*inserts$offsetY + 2
  inserts$AX <- inserts$AX + (tw/2)*inserts$offsetX
  inserts <- select(inserts, -c(offsetY, offsetX))
}



#inserts$AY <- y.inner - spacing*(inserts$chromosome - 1)
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

#### 12.06.2019
## THIS IS A SHORTCUT!!!
## GET TO THE BOTTOM OF IT!!!
f <- c()
for (i in seq(nrow(inserts))){
  if(inserts$insert.genomic.coord.s.start[i] <
     chrlen$Length.bp[chrlen$Chromosome == str_split(inserts$chromosome[i], "chr")[[1]][2]]){
    f <- c(f, i)
  }
}
inserts <- inserts[f,]

## Rectangle version
## Nahuy takoy graphic
#d <- ggplot(data = chrlen)+
#  coord_cartesian(xlim = c(0, plot.width), ylim = c(0, plot.height))+
#  geom_rect(aes(xmin = AX, ymin = AY, xmax = BX, ymax = BY), alpha = 0.5, color = 'black', fill = 'black')
#d



## Adding centromeres
gaps <- read.table(file = 'hg19chrGaps.txt', sep = '\t', header = TRUE)
centromeres <- filter(gaps, type == 'centromere')
centromeres$chromCenter <- centromeres$chromStart + 3000000/2
centromeres <- select(centromeres, chrom, chromCenter)
centromeres$Chromosome <- str_replace_all(centromeres$chrom, 'chr','')

chrlen <- left_join(chrlen, centromeres, by = "Chromosome", copy = FALSE)
chrlen$chromCenter.relative <- chrlen$chromCenter/max(chrlen$Length.bp)
chrlen$chromCenterY <- chrlen$AY
chrlen$chromCenterX <- x.margin + inner.width*chrlen$chromCenter.relative

inserts$insert.type[is.na(inserts$insert.type)] <- 'V-CH1'

## 19.02.19
chrlen$Length.Mbp <- chrlen$Length.bp / 1000000

chrlen <- filter(chrlen, !Length.Mbp < 1)
chrlen$Length.Mbp <- round(chrlen$Length.Mbp, 0)
chrlen$Length.Mbp <- as.character(chrlen$Length.Mbp)
chrlen$Length.Mbp <- paste(chrlen$Length.Mbp, 'Mbp')

color_of_lines <- 'black'

## Line version
g <- ggplot(data = chrlen)+
  coord_cartesian(xlim = c(0, plot.width), ylim = c(0, plot.height))+
  geom_segment(aes(x = AX, y = AY, xend = BX, yend = BY), alpha = 1, color = color_of_lines, size = 1)+
  geom_point(aes(x = chromCenterX, y = chromCenterY), size = 2, shape = 21, fill = 'deeppink2')+
  geom_text(aes(x = AX-35, y = AY+3, label = chrom), colour = color_of_lines)+
  geom_text(aes(x = BX + 50, y = AY+3, label = Length.Mbp), size = 3, colour = color_of_lines)+
#  geom_polygon(data = inserts, aes(x=dotX, y=dotY, group = ins.id,
#                                   fill = insert.type), alpha = 0.75)+
  geom_polygon(data = inserts, aes(x=dotX, y=dotY, group = ins.id, fill = insert.type,
                                   size = insert.type), alpha = 1.0)+
  
  scale_fill_manual(values = c('V-DJ' = rgb(0,0,255, maxColorValue = 255),
                               'VDJ-CH1' = rgb(255, 0, 0, maxColorValue = 255),
                               'V-CH1' = rgb(255,255,255, maxColorValue = 255)))+
  
  scale_colour_manual(values = c('V-DJ' = rgb(0,0,0, maxColorValue = 255),
                               'VDJ-CH1' = rgb(0,0,0, maxColorValue = 255),
                               'V-CH1' = 'grey45'))+
  scale_size_manual(values = c('V-DJ' = 0.1,
                                   'VDJ-CH1' = 0.1,
                                   'V-CH1' = 0.03))+
  theme_void()

g

g + theme(
  # get rid of panel grids
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  # Change plot and panel background
  plot.background=element_rect(fill = "black"),
  panel.background = element_rect(fill = 'black'),
  # Change legend 
 # legend.position = c(0.6, 0.07),
  #legend.direction = "horizontal",
  legend.background = element_rect(fill = "black", color = NA),
  legend.key = element_rect(color = "gray", fill = "black"),
  legend.title = element_text(color = "white"),
  legend.text = element_text(color = "white")
)

g + theme(legend.position = 'bottom')
g + geom_text(data = genesins,
              aes(x = dotX, y = dotY+ 50, label = gene), size = 3)


genesins <- inserts %>% separate(insert.id.s., c(NA, "gene", NA)) %>% group_by(gene) %>%
  summarise(dotX = nth(dotX,1),
            dotY = nth(dotY,1))
