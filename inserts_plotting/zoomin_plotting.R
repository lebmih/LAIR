#### LIBRARIES
##############
libs <- c('tidyr', 'dplyr', 'ggplot2', 'ggrepel', 'reshape2', 'stringr', 'ggforce')
lapply(libs, require, character.only = TRUE) %>% unlist() %>%
  sum() %>% paste("libraries installed successfully")
rm(libs)
##############

#### FUNCTIONS
#####################
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
#####################

#### INSERTS DATA IMPORT AND PREFILTERING
#########################################
## Importing inserts data
setwd("C:/Users/User/Documents/Lab/LongBCRs/MethodPaper/RunTables/")
all_inserts <- read.table("data_cured_Kathrins.txt", sep = "\t", header = TRUE, na.strings = "NA", 
           stringsAsFactors = FALSE)
## Data prefiltering
# Widening the table - focus on the inserts origins, not on the contigs
all_inserts$contig.id <- seq(nrow(all_inserts))
columnstosep <- c('insert.id.s.', 'insert.genomic.coord.s.', 'insert.contig.coord.s.','insert.seq.s.')
all_inserts <- separate_rows(all_inserts, columnstosep, sep = ',\\W*')
rm(columnstosep)

# Make an ID for each insert
all_inserts$ins.id <- seq(nrow(all_inserts))

## Filtering out LAIR1 insertions
all_inserts <- filter(all_inserts, !str_detect(all_inserts$Gene.id.s., 'LAIR1'))

## Extracting a piece of data as input for the function
## if you want to drop V-CH1 uncomment the following
#inserts <- all_inserts[!is.na(all_inserts$insert.type),]
## if you want to show V-CH1 uncomment the following
inserts <- all_inserts

# Remove the original table not to clutter the environment
rm(all_inserts)

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
                    sep = '\\..')

## Filtering out non-valid chromosomes
filteredchromosomes <- c("chrM","chrUn_gl000220","chrUn_gl000214","chr7_gl000195_random")
inserts <- filter(inserts, !chromosome %in% filteredchromosomes)
rm(filteredchromosomes)

#########################################

#### GENES DATA IMPORT AND PREFILTERING
#######################################
## Here I use NCBI genes database and filter out predicted RNAs (XR and XM)
agni %>% filter(str_detect(name, "NM") | str_detect(name, 'NR')) -> agni_nm
agni_nm %>% arrange(desc(exonCount)) %>% group_by(name2) %>%
  summarise_all(nth, n = 1) %>%
  select(-c('bin','name','cdsStart','cdsEnd','score','cdsStartStat','cdsEndStat','exonFrames')) ->
  agni_nm_unique
agni_nm_unique %>% rename('txStart' = 'start', 'txEnd' = 'end') -> agni_nm_unique

# There are "," on the end of the lines in exonStarts and exonEnds, removing them
agni_nm_unique %>% mutate(exonStarts = str_remove(exonStarts, ',$'),
                          exonEnds = str_remove(exonEnds, ',$')) -> agni_nm_unique
zoominori <- agni_nm_unique
#######################################

#### SPECIFYING PLOT PARAMETERS
###############################
# Plot characteristics
plot.prm <- list(width = 1600, height = 900, x.margin = 50,
                 y.margin = 10, thick = 10, trianglewidth = 6,
                 triangleheight = 10)
plot.prm$inner.width <- plot.prm$width - 2*plot.prm$x.margin
plot.prm$inner.height <- plot.prm$height - 2*plot.prm$y.margin
color_of_lines <- 'black'
arrow.height <- plot.prm$thick * 2
arrow.width <- arrow.height/2
dock.height <- plot.prm$thick
dock.width <- plot.prm$thick/2
ellipse.height <- plot.prm$thick/2
text.offset <- -30
step.between.inserts <- 15
insert.thickness <- 3
###############################

#### SPECIFYING THE ROI
#######################
region_chromosome <- "chr14"
region_start <- 105438000
region_end <- 105499000
plot.prm$bpperpxl <- (region_end - region_start)/plot.prm$inner.width
#######################

#### FILTERING AND PROCESSING THE INSERTS DATA
##############################################
inserts %>% rename('insert.genomic.coord.s.start' = 'start',
                   'insert.genomic.coord.s.end' = 'end') -> inserts
inserts %>% mutate_at(c('start','end'), as.integer) -> inserts

inserts %>% filter((chromosome == region_chromosome) &
                     ((start >= region_start & start <= region_end) |
                        (end >= region_start & end <= region_end))) -> inserts

inserts$start[inserts$start < region_start] <- region_start
inserts$end[inserts$end > region_end] <- region_end
inserts$ins.id <- seq(nrow(inserts))
inserts$startX <- (inserts$start-region_start)/plot.prm$bpperpxl
inserts$endX <- (inserts$end-region_start)/plot.prm$bpperpxl

inserts %>% arrange(startX) -> inserts

inserts$level <- 0
LVL <- 1

while(sum(inserts$level == 0) > 0){
  inserts$level[which(inserts$level == 0)[1]] <- LVL
  for(j in which(inserts$level == 0)){
    if(count.overlap(subject.starts = inserts$start[inserts$level==LVL],
                     subject.ends = inserts$end[inserts$level==LVL],
                     query.start = inserts$start[j],
                     query.end = inserts$end[j]) == 0){
      inserts$level[j] <- LVL
    }
  }
  LVL <- LVL + 1
}

inserts$startY <- plot.prm$height/2 + inserts$level*step.between.inserts

inserts$endY <- inserts$startY
##############################################

#### EXTRACTING GENES EXPRESSION FROM DEGL DATA
###############################################
degl <- read.table("diff_expression_genes_310719.txt",
                   sep = "\t", header = TRUE, na.strings = "NA", 
                   stringsAsFactors = FALSE)

degl$fill <- "grey"
degl$fill[degl$NBCs < quantile(degl$NBCs)[2]] <- "red"
degl$fill[degl$NBCs < quantile(degl$NBCs)[3] & degl$NBCs >= quantile(degl$NBCs)[2]] <- "orange"
degl$fill[degl$NBCs < quantile(degl$NBCs)[4] & degl$NBCs >= quantile(degl$NBCs)[3]] <- "green"
degl$fill[degl$NBCs >= quantile(degl$NBCs)[4]] <- "turquoise1"

degl %>% select(Name, NBCs, start, end, chrom, fill) %>%
  filter((chrom == region_chromosome) &
           ((start >= region_start & start <= region_end) |
              (end >= region_start & end <= region_end))) -> degl

# This is a shortcut, fix that
degl$Name[degl$Name == 'C14orf80'] <- 'TEDC1'

degl %>% rename("Name" = "name2") -> degl
###############################################

#### JOINING GENE DATA WITH GENES EXPRESSION
############################################
zoominori %>% left_join(select(degl, name2, NBCs, fill), by = 'name2') -> zoominori
zoominori$fill[is.na(zoominori$fill)] <- 'grey'
############################################

#### PROCESSING THE GENE DATA
#############################
zoominoriexons <- separate_rows(zoominori, exonStarts, exonEnds, sep = ',\\W*')
zoominoriexons$exon <- "exon"
zoominoriexons$exon[zoominoriexons$exonEnds == zoominoriexons$end] <- "last"
zoominoriexons$shape <- "rect"
zoominoriexons$shape[zoominoriexons$exonEnds == zoominoriexons$end] <- "arrow"

zoominori$exon <- "gene"
zoominori$shape <- "line"
zoominoriexons %>% filter((chrom == region_chromosome) &
                       ((exonStarts >= region_start & exonStarts <= region_end) |
                          (exonEnds >= region_start & exonEnds <= region_end))) -> zoominoriexons
zoominoriexons %>% mutate(start = exonStarts,
                          end = exonEnds) -> zoominoriexons
zoominoriexons %>% select(-c('exonStarts','exonEnds')) -> zoominoriexons
zoominoriexons %>% mutate_at(c('start','end'), as.integer) -> zoominoriexons

zoominori %>% filter((chrom == region_chromosome) &
                       ((start >= region_start & start <= region_end) |
                          (end >= region_start & end <= region_end))) -> zoominori
zoominori %>% select(-c('exonStarts','exonEnds')) -> zoominori

zoominori <- bind_rows(zoominoriexons, zoominori)


zoominori$start[zoominori$start < region_start] <- region_start
zoominori$end[zoominori$end > region_end] <- region_end

zoominori$startX <- (zoominori$start-region_start)/plot.prm$bpperpxl
zoominori$endX <- (zoominori$end-region_start)/plot.prm$bpperpxl
zoominori$startY <- plot.prm$height/2
zoominori$endY <- plot.prm$height/2
zoominori$feat.id <- seq(nrow(zoominori))

zoominori %>% filter(shape %in% c('circle','ellipse')) -> zoominori.circular
if(nrow(zoominori.circular)>0){
  zoominori.circular %>% mutate(x0 = (startX+endX)/2,
                                y0 = startY,
                                a = (endX-startX)/2) %>%
    select(-c('startX','startY','endX','endY')) -> zoominori.circular
}


zoominori %>% filter(shape %in% c('rect','rectangle')) -> zoominori.rectangular
zoominori %>% filter(shape == 'dock') -> zoominori.dock

zoominori.rectangular %>% mutate(AX = startX,
                                 AY = startY + plot.prm$thick/2,
                                 BX = endX,
                                 BY = endY + plot.prm$thick/2,
                                 CX = endX,
                                 CY = endY - plot.prm$thick/2,
                                 DX = startX,
                                 DY = startY - plot.prm$thick/2) %>%
  select(-c('startX','startY','endX','endY')) %>%
  unite('A', c('AX', 'AY'), sep = ";", remove = TRUE) %>%
  unite('B', c('BX', 'BY'), sep = ";", remove = TRUE) %>%
  unite('C', c('CX', 'CY'), sep = ";", remove = TRUE) %>%
  unite('D', c('DX', 'DY'), sep = ";", remove = TRUE) %>%
  gather(dot, dot.coord, A:D) %>%
  separate(dot.coord, c('dotX','dotY'), sep = ';', remove = TRUE) -> zoominori.rectangular


zoominori.dock %>% filter(strand == "-") %>%
  mutate(AX = startX,
         AY = startY + plot.prm$thick/2,
         BX = endX,
         BY = endY + plot.prm$thick/2,
         CX = endX - dock.width,
         CY = endY,
         DX = endX,
         DY = endY - plot.prm$thick/2,
         EX = startX,
         EY = startY - plot.prm$thick/2) -> zoominori.dock.minus
zoominori.dock.minus$CX[zoominori.dock.minus$CX < zoominori.dock.minus$startX] <-
  zoominori.dock.minus$startX[zoominori.dock.minus$CX < zoominori.dock.minus$startX]
zoominori.dock.minus %>%
  select(-c('startX','startY','endX','endY')) %>%
  unite('A', c('AX', 'AY'), sep = ";", remove = TRUE) %>%
  unite('B', c('BX', 'BY'), sep = ";", remove = TRUE) %>%
  unite('C', c('CX', 'CY'), sep = ";", remove = TRUE) %>%
  unite('D', c('DX', 'DY'), sep = ";", remove = TRUE) %>%
  unite('E', c('EX', 'EY'), sep = ";", remove = TRUE) %>%
  gather(dot, dot.coord, A:E) %>%
  separate(dot.coord, c('dotX','dotY'), sep = ';', remove = TRUE) -> zoominori.dock.minus
zoominori.dock %>% filter(strand == "+") %>%
  mutate(AX = startX,
         AY = startY + plot.prm$thick/2,
         BX = endX,
         BY = endY + plot.prm$thick/2,
         CX = endX,
         CY = endY - plot.prm$thick/2,
         DX = startX,
         DY = startY - plot.prm$thick/2,
         EX = startX + dock.width,
         EY = startY) -> zoominori.dock.plus
zoominori.dock.plus$EX[zoominori.dock.plus$EX > zoominori.dock.plus$endX] <-
  zoominori.dock.plus$endX[zoominori.dock.plus$EX > zoominori.dock.plus$endX]
zoominori.dock.plus %>%  
  select(-c('startX','startY','endX','endY')) %>%
  unite('A', c('AX', 'AY'), sep = ";", remove = TRUE) %>%
  unite('B', c('BX', 'BY'), sep = ";", remove = TRUE) %>%
  unite('C', c('CX', 'CY'), sep = ";", remove = TRUE) %>%
  unite('D', c('DX', 'DY'), sep = ";", remove = TRUE) %>%
  unite('E', c('EX', 'EY'), sep = ";", remove = TRUE) %>%
  gather(dot, dot.coord, A:E) %>%
  separate(dot.coord, c('dotX','dotY'), sep = ';', remove = TRUE) -> zoominori.dock.plus
zoominori.dock <- bind_rows(zoominori.dock.plus, zoominori.dock.minus)
rm(zoominori.dock.plus, zoominori.dock.minus)

zoominori %>% filter(shape == 'arrow' & strand == '+') %>%
  mutate(AX = startX,
         AY = startY + plot.prm$thick/2,
         BX = endX - arrow.width,
         BY = startY + plot.prm$thick/2,
         CX = endX - arrow.width,
         CY = startY + arrow.height/2,
         DX = endX,
         DY = endY,
         EX = endX - arrow.width,
         EY = startY - arrow.height/2,
         FX = endX - arrow.width,
         FY = startY - plot.prm$thick/2,
         GX = startX,
         GY = startY - plot.prm$thick/2) -> zoominori.arrow.plus
zoominori.arrow.plus$CX[zoominori.arrow.plus$CX < zoominori.arrow.plus$startX] <-
  zoominori.arrow.plus$startX[zoominori.arrow.plus$CX < zoominori.arrow.plus$startX]
zoominori.arrow.plus$BX[zoominori.arrow.plus$BX < zoominori.arrow.plus$startX] <-
  zoominori.arrow.plus$startX[zoominori.arrow.plus$BX < zoominori.arrow.plus$startX]
zoominori.arrow.plus$EX[zoominori.arrow.plus$EX < zoominori.arrow.plus$startX] <-
  zoominori.arrow.plus$startX[zoominori.arrow.plus$EX < zoominori.arrow.plus$startX]
zoominori.arrow.plus$FX[zoominori.arrow.plus$FX < zoominori.arrow.plus$startX] <-
  zoominori.arrow.plus$startX[zoominori.arrow.plus$FX < zoominori.arrow.plus$startX]
zoominori.arrow.plus %>%
  select(-c('startX','startY','endX','endY')) %>%
  unite('A', c('AX', 'AY'), sep = ";", remove = TRUE) %>%
  unite('B', c('BX', 'BY'), sep = ";", remove = TRUE) %>%
  unite('C', c('CX', 'CY'), sep = ";", remove = TRUE) %>%
  unite('D', c('DX', 'DY'), sep = ";", remove = TRUE) %>%
  unite('E', c('EX', 'EY'), sep = ";", remove = TRUE) %>%
  unite('F', c('FX', 'FY'), sep = ";", remove = TRUE) %>%
  unite('G', c('GX', 'GY'), sep = ";", remove = TRUE) %>%
  gather(dot, dot.coord, A:G) %>%
  separate(dot.coord, c('dotX','dotY'), sep = ';', remove = TRUE) -> zoominori.arrow.plus
zoominori %>% filter(shape == 'arrow' & strand == '-') %>%
  mutate(AX = startX,
         AY = startY,
         BX = startX + arrow.width,
         BY = startY + arrow.height/2,
         CX = startX + arrow.width,
         CY = startY + plot.prm$thick/2,
         DX = endX,
         DY = endY + plot.prm$thick/2,
         EX = endX,
         EY = endY - plot.prm$thick/2,
         FX = startX + arrow.width,
         FY = startY - plot.prm$thick/2,
         GX = startX + arrow.width,
         GY = startY - arrow.height/2) -> zoominori.arrow.minus
zoominori.arrow.minus$BX[zoominori.arrow.minus$BX > zoominori.arrow.minus$endX] <-
  zoominori.arrow.minus$endX[zoominori.arrow.minus$BX > zoominori.arrow.minus$endX]
zoominori.arrow.minus$CX[zoominori.arrow.minus$CX > zoominori.arrow.minus$endX] <-
  zoominori.arrow.minus$endX[zoominori.arrow.minus$CX > zoominori.arrow.minus$endX]
zoominori.arrow.minus$FX[zoominori.arrow.minus$FX > zoominori.arrow.minus$endX] <-
  zoominori.arrow.minus$endX[zoominori.arrow.minus$FX > zoominori.arrow.minus$endX]
zoominori.arrow.minus$GX[zoominori.arrow.minus$GX > zoominori.arrow.minus$endX] <-
  zoominori.arrow.minus$endX[zoominori.arrow.minus$GX > zoominori.arrow.minus$endX]
zoominori.arrow.minus %>%
  select(-c('startX','startY','endX','endY')) %>%
  unite('A', c('AX', 'AY'), sep = ";", remove = TRUE) %>%
  unite('B', c('BX', 'BY'), sep = ";", remove = TRUE) %>%
  unite('C', c('CX', 'CY'), sep = ";", remove = TRUE) %>%
  unite('D', c('DX', 'DY'), sep = ";", remove = TRUE) %>%
  unite('E', c('EX', 'EY'), sep = ";", remove = TRUE) %>%
  unite('F', c('FX', 'FY'), sep = ";", remove = TRUE) %>%
  unite('G', c('GX', 'GY'), sep = ";", remove = TRUE) %>%
  gather(dot, dot.coord, A:G) %>%
  separate(dot.coord, c('dotX','dotY'), sep = ';', remove = TRUE) -> zoominori.arrow.minus
zoominori.arrow <- bind_rows(zoominori.arrow.plus, zoominori.arrow.minus)
rm(zoominori.arrow.plus, zoominori.arrow.minus)

bind_rows(zoominori.rectangular, zoominori.dock, zoominori.arrow) -> zoominori.no.circular
rm(zoominori.rectangular, zoominori.dock, zoominori.arrow)
zoominori.no.circular %>% mutate_at(c('dotX','dotY'), as.numeric) -> zoominori.no.circular
zoominori %>% mutate(centerX = (startX+endX)/2) -> zoominori
#############################

#### EXTRACTING FILL COLOR INFO
###############################
fill_color <- zoominori %>% group_by(feat.id) %>% summarise(fill_color = nth(fill,1)) %>%
  mutate(feat.id = as.character(feat.id))
fill <- fill_color$fill_color
names(fill) <- fill_color$feat.id
rm(fill_color)
###############################

#### PLOTTING
#############
zoominplot <- ggplot(data = zoominori)+
  coord_cartesian(xlim = c(0, plot.prm$width), ylim = c(0, plot.prm$height))+
  geom_segment(aes(x = plot.prm$x.margin,
                   y = plot.prm$height/2,
                   xend = plot.prm$width-plot.prm$x.margin,
                   yend = plot.prm$height/2),
               color = color_of_lines, size = 1, alpha = 0.5)+
  geom_segment(data = zoominori %>% filter(exon == 'gene'),
               aes(x = startX + plot.prm$x.margin,
                   y = startY,
                   xend = endX + plot.prm$x.margin,
                   yend = endY, color = as.character(feat.id)),
               size = 1.5, linetype = 'dashed')+
  
  geom_text(data = zoominori %>% filter(exon == 'gene'),
            aes(x = centerX+plot.prm$x.margin, y = startY + text.offset,
                label = name2), size = 4, colour = color_of_lines, check_overlap = TRUE)+
  geom_polygon(data = zoominori.no.circular,
               aes(x=dotX+plot.prm$x.margin, y=dotY, group = feat.id, fill = as.character(feat.id)),
               color = "black",
               alpha = 1)+
  geom_segment(data = zoominori[1,],
               aes(x = plot.prm$x.margin, y = startY + 2*text.offset,
                   xend = plot.prm$x.margin+1000/plot.prm$bpperpxl,
                   yend = startY + 2*text.offset),
               arrow = arrow(angle = 90, ends = "both", length = unit(10/plot.prm$height, "npc")))+
  geom_text(data = zoominori[1,],
            aes(x = plot.prm$x.margin+500/plot.prm$bpperpxl, y = startY + 2.5*text.offset,
                label = "1 kb"), size = 4, colour = color_of_lines)+
  scale_fill_manual(values = fill)+
  scale_color_manual(values = fill)+
  theme_void()+
  theme(legend.position = "none")

if(nrow(zoominori.circular) != 0){
  zoominplot <- zoominplot + geom_ellipse(data = zoominori.circular,
                                          aes(x0 = x0+x.margin, y0 = y0, a = a, b = ellipse.height, angle = 0),
                                          fill = zoominori.circular$fill)
}

## Adding the inserts to the plot
zoominplot <- zoominplot +
  geom_segment(data = inserts,
               aes(x = startX + plot.prm$x.margin,
                   y = startY,
                   xend = endX + plot.prm$x.margin,
                   yend = endY),
               size = insert.thickness, color = 'red', alpha = 1.0)

## Rendering
zoominplot
#############

library(plotly)
ziplot <- ggplotly(zoominplot)
ziplot
