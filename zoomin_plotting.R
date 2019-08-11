library(ggforce)

#### MY OWN PLOTTING
## DUMMY DATA
setwd("C:/Users/User/Documents/Lab/LongBCRs/Scripts/GenomePlotting/")
zoominori <- read.table("zoominori_dummy.txt", sep = "\t", header = TRUE, na.strings = "NA", 
           stringsAsFactors = FALSE)


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
#######  PLOT PARAMETERS
##########################
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
##########################

### SPECIFYING THE ROI
region_chromosome <- "chr3"
region_start <- 0
region_end <- 11000
zoominori %>% filter((chrom == region_chromosome) &
                       ((start >= region_start & start <= region_end) |
                          (end >= region_start & end <= region_end))) -> zoominori

zoominori$start[zoominori$start < region_start] <- region_start
zoominori$end[zoominori$end > region_end] <- region_end

plot.prm$bpperpxl <- (region_end - region_start)/plot.prm$inner.width 

zoominori$startX <- (zoominori$start-region_start)/plot.prm$bpperpxl
zoominori$endX <- (zoominori$end-region_start)/plot.prm$bpperpxl
zoominori$startY <- plot.prm$height/2
zoominori$endY <- plot.prm$height/2
zoominori$feat.id <- seq(nrow(zoominori))

zoominori %>% filter(shape %in% c('circle','ellipse')) -> zoominori.circular
zoominori.circular %>% mutate(x0 = (startX+endX)/2,
                              y0 = startY,
                              a = (endX-startX)/2) %>%
  select(-c('startX','startY','endX','endY')) -> zoominori.circular

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

fill_color <- zoominori.no.circular %>% group_by(feat.id) %>% summarise(fill_color = nth(fill,1)) %>%
  mutate(feat.id = as.character(feat.id))
fill <- fill_color$fill_color
names(fill) <- fill_color$feat.id
rm(fill_color)

zoominori %>% mutate(centerX = (startX+endX)/2) -> zoominori

zoominplot <- ggplot(data = zoominori)+
  coord_cartesian(xlim = c(0, plot.prm$width), ylim = c(0, plot.prm$height))+
  geom_segment(aes(x = plot.prm$x.margin,
                   y = plot.prm$height/2,
                   xend = plot.prm$width-plot.prm$x.margin,
                   yend = plot.prm$height/2),
               color = color_of_lines, size = 1, alpha = 1)+
  
  geom_text(data = zoominori,
            aes(x = centerX+x.margin, y = startY + text.offset,
                label = text), size = 4, colour = color_of_lines, check_overlap = TRUE)+
  geom_polygon(data = zoominori.no.circular,
               aes(x=dotX+x.margin, y=dotY, group = feat.id, fill = as.character(feat.id)), color = "black",
               alpha = 1)+
  geom_segment(data = zoominori[1,],
               aes(x = x.margin, y = startY + 2*text.offset,
                   xend = x.margin+1000/plot.prm$bpperpxl,
                   yend = startY + 2*text.offset),
               arrow = arrow(angle = 90, ends = "both", length = unit(10/plot.prm$height, "npc")))+
  geom_text(data = zoominori[1,],
            aes(x = x.margin+500/plot.prm$bpperpxl, y = startY + 2.5*text.offset,
                label = "1 kb"), size = 4, colour = color_of_lines)+
  scale_fill_manual(values = fill)+
  theme_void()+
  theme(legend.position = "none")

if(nrow(zoominori.circular) != 0){
  zoominplot <- zoominplot + geom_ellipse(data = zoominori.circular,
                                          aes(x0 = x0+x.margin, y0 = y0, a = a, b = ellipse.height, angle = 0),
                                          fill = zoominori.circular$fill)
}

zoominplot

insertions <- read.table("dummy_insertions.txt", sep = "\t", header = TRUE, na.strings = "NA", 
                         stringsAsFactors = FALSE)

insertions %>% filter((chr == region_chromosome) &
                        ((start >= region_start & start <= region_end) |
                           (end >= region_start & end <= region_end))) -> insertions

insertions$start[insertions$start < region_start] <- region_start
insertions$end[insertions$end > region_end] <- region_end
insertions$ins.id <- seq(nrow(insertions))
insertions$startX <- (insertions$start-region_start)/plot.prm$bpperpxl
insertions$endX <- (insertions$end-region_start)/plot.prm$bpperpxl



insertions %>% arrange(startX) -> insertions

insertions$level <- 0
LVL <- 1

while(sum(insertions$level == 0) > 0){
  insertions$level[which(insertions$level == 0)[1]] <- LVL
  for(j in which(insertions$level == 0)){
    if(count.overlap(subject.starts = insertions$start[insertions$level==LVL],
                     subject.ends = insertions$end[insertions$level==LVL],
                     query.start = insertions$start[j],
                     query.end = insertions$end[j]) == 0){
      insertions$level[j] <- LVL
    }
  }
  LVL <- LVL + 1
}

insertions$startY <- plot.prm$height/2 + insertions$level*step.between.inserts

insertions$endY <- insertions$startY



zoominplot <- zoominplot +
  geom_segment(data = insertions,
               aes(x = startX + x.margin,
                   y = startY,
                   xend = endX + x.margin,
                   yend = endY),
               size = insert.thickness, color = 'red', alpha = 1.0)

zoominplot

library(plotly)
ziplot <- ggplotly(zoominplot)
ziplot
