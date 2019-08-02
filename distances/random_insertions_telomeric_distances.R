## Working directory
setwd("~/Lab/LongBCRs/MethodPaper/RunTables")

## Libraries
x <- c('tidyr', 'dplyr', 'ggplot2', 'ggrepel', 'reshape2', 'stringr', 'Biostrings', 'ggridges',
       'ggpubr', 'cowplot')
lapply(x, require, character.only = TRUE)
rm(x)


### Randomize the insertions
## Plan
## We need to take the insertions table and change the chromosomal coordinate to
## random position inside the same chromosome
## Then we need to calculate overlaps for our true data and for 1000 generations
## of random data

## For the start let's just plot the random insertions
# Replace the insert.genomic.center with random number
# between 1 and chromosome length
degl <- read.table(file = 'degl_expanded_310719.txt',
                      sep = '\t', header = TRUE, stringsAsFactors = FALSE)

inserts <- left_join(inserts, select(chrlen, chromosome, Length.bp), by = 'chromosome')
telomeres.length <- 10000


## Generating the artificial insertions dataset
inserts.generic <- data.frame()

for (j in seq(1:100)){
  inserts.generic.temp <- inserts
  rand_position <- c()
  for (i in seq_along(inserts.generic.temp$chromosome)){
    rand_position <- c(rand_position, sample(1:inserts.generic.temp$Length.bp[i],1))
  }
  inserts.generic.temp$insert.genomic.center <- rand_position
  inserts.generic.temp$insert.genomic.coord.s.start <-
    inserts.generic.temp$insert.genomic.center - inserts.generic.temp$insert.length/2
  inserts.generic.temp$insert.genomic.coord.s.end <-
    inserts.generic.temp$insert.genomic.center + inserts.generic.temp$insert.length/2
  inserts.generic <- bind_rows(inserts.generic, inserts.generic.temp)
  if(j %in% seq(0,100,10)){
    print(paste(j,"out of 100 cycles finished"))
  }
}

rm(i, j, rand_position, inserts.generic.temp)

#inserts.generic$telomeric.distance <- mapply(min,
#                                             inserts.generic$insert.genomic.coord.s.start - telomeres.length,
#                                             inserts.generic$Length.bp - telomeres.length - inserts.generic$insert.genomic.coord.s.end)

#for (j in seq(1:1000)){
#  inserts1 <- inserts
#  rand_position <- c()
#  for (i in seq_along(inserts1$chromosome)){
#    rand_position <- c(rand_position, sample(1:inserts1$Length.bp[i],1))
#  }
#  inserts1$insert.genomic.center <- rand_position
#  inserts1$insert.genomic.coord.s.start <- inserts1$insert.genomic.center - inserts1$insert.length/2
#  inserts1$insert.genomic.coord.s.end <- inserts1$insert.genomic.center + inserts1$insert.length/2
#  inserts1$telomeric.distance <- mapply(min, inserts1$insert.genomic.coord.s.start - telomeres.length,
#                                        inserts1$Length.bp - telomeres.length - inserts1$insert.genomic.coord.s.end)
#  telomeric.distance <- c(telomeric.distance, inserts1$telomeric.distance)
#}

#inserts.vdjch1 <- filter(inserts, insert.type == 'VDJ-CH1')
#inserts.vdj <- filter(inserts, insert.type == 'V-DJ')
#inserts.vch1 <- filter(inserts, insert.type == 'V-CH1')

#teldistanceTotal <- data.frame(inserts = 'Total', telomeric.distance = inserts$telomeric.distance)
#teldistanceVDJCH1 <- data.frame(inserts = 'VDJ-CH1', telomeric.distance = inserts.vdjch1$telomeric.distance)
#teldistanceVDJ <- data.frame(inserts = 'V-DJ', telomeric.distance = inserts.vdj$telomeric.distance)
#teldistanceVCH1 <- data.frame(inserts = 'V-CH1', telomeric.distance = inserts.vch1$telomeric.distance)


#teldistanceGenerated <- data.frame(inserts = 'Generated', telomeric.distance = inserts.generic$telomeric.distance)

#teldistance <- bind_rows(teldistanceVDJCH1, teldistanceVDJ,
#                         teldistanceVCH1, teldistanceGenerated)
## 280719
#teldistance <- bind_rows(teldistanceVDJCH1, teldistanceVDJ, teldistanceGenerated)
inserts.generic$insert.type <- 'Generated'

inserts.all <- bind_rows(inserts, inserts.generic)

inserts.all$insert.type[is.na(inserts.all$insert.type)] <- 'V-CH1'

inserts.all$telomeric.distance <- mapply(min,
                                     inserts.all$insert.genomic.coord.s.start - telomeres.length,
                                     inserts.all$Length.bp - telomeres.length -
                                       inserts.all$insert.genomic.coord.s.end)

inserts.all$telomeric.distance.Mbp <- inserts.all$telomeric.distance/1000000

# H0: Median Telomeric Distance of random insertions equals to that of real insertions
wilcox.test(inserts.all$telomeric.distance[inserts.all$insert.type == 'Generated'],
            inserts.all$telomeric.distance[inserts.all$insert.type == 'VDJ-CH1'],
            mu = 0, alternative = "greater",
            paired = F, correct = F, exact = T)


telomeric_dist_boxplot <- ggplot(data = inserts.all)+
  geom_boxplot(aes(x = insert.type, y = telomeric.distance.Mbp), alpha = 0.2)+
  theme_classic()+
  #xlab("")+
  #scale_x_discrete(labels = c('Random' = 'Random\nInsertions', 'Real' = 'Donors\nInsertions'))+
  ylab("Distance of insertion donor\nto the closest telomere (Mb)")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  xlab("Insert type")+
  annotate('text', x = 2, y = 130, label = 'p = 0.86', size = 4)+
  annotate('text', x = 3, y = 130, label = 'p < 2.2e-16', size = 4)
  #annotate('text', x = 4, y = 130, label = 'p < 2.2e-16', size = 4)
telomeric_dist_boxplot



telomeric_dist_density <- ggplot(data = inserts.all)+
  geom_density(aes(x = telomeric.distance.Mbp, y = ..density.., fill = insert.type), alpha = 0.2)+
  labs(x = "Distance to the closest telomere, Mbp", fill = 'Insert type')+
  theme_classic()
telomeric_dist_density

## Plots saving
ggsave(paste('telomeric_dist_boxplot_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = telomeric_dist_boxplot)


ggsave(paste('telomeric_dist_density_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = telomeric_dist_density)

delete_the_plots_from_R_memory <- TRUE
if(delete_the_plots_from_R_memory == TRUE){rm(telomeric_dist_boxplot, telomeric_dist_boxplot)}


## 20.02.19
gaps <- read.table(file = 'hg19chrGaps.txt', sep = '\t', header = TRUE)


## 21.02.19
## Let's find the distance between rdna and the inserts
rdna <- read.table(file = 'rdna_bed.txt', sep = '\t', header = TRUE)
unique(rdna$feature)
# drop score as we don't need it
rdna <- select(rdna, - score)
                            
## Calculating the center of rDNA regions
rdna$chrom.center <- (rdna$chrom.end + rdna$chrom.start)/2
rdna <- filter(rdna, !str_detect(chromosome, '_'))
rdna <- filter(rdna, !str_detect(chromosome, 'chrM'))
unique(rdna$chromosome)

rm(gaps)

dist.to.closest.rdna <- c()
for (i in seq_along(inserts.all$run.id)){
  rdna.in.this.chrom <- rdna$chrom.center[rdna$chromosome == inserts.all$chromosome[i]]
  rdna.dist <- min(abs(rdna.in.this.chrom - inserts.all$insert.genomic.center[i]))
  dist.to.closest.rdna <- c(dist.to.closest.rdna, rdna.dist)
}

inserts.all$rdna.distance <- dist.to.closest.rdna

inserts.all$rdna.dist.mpb <- inserts.all$rdna.distance / 1000000

rm(rdna.in.this.chrom, dist.to.closest.rdna, i, rdna.dist)

## Testing the statistical significance
inserts.all$insert.type %>% unique() -> unique.insert.types
uni.ins.types.n <- length(unique.insert.types)

wilcox.matrix <- matrix(data = NA, nrow = uni.ins.types.n, ncol = uni.ins.types.n,
                        dimnames = list(unique.insert.types, unique.insert.types))

for(i in colnames(wilcox.matrix)){
  for(j in rownames(wilcox.matrix)){
    wilcox.test(inserts.all$rdna.distance[inserts.all$insert.type == i],
                inserts.all$rdna.distance[inserts.all$insert.type == j],
                mu = 0, alternative = "two.sided",
                paired = F, correct = F, exact = T)$p.value %>%
      format(digits = 2, scientific = TRUE) -> wilcox.matrix[i,j]
  }
}
rm(i,j)



## Plotting the distance to rDNA
rdna_dist_boxplot <- ggplot(data = inserts.all)+
  geom_boxplot(aes(x = insert.type, y = rdna.dist.mpb), alpha = 0.2)+
  theme_classic()+
  #xlab("")+
  #scale_x_discrete(labels = c('Random' = 'Random\nInsertions', 'Real' = 'Donors\nInsertions'))+
  ylim(0, 5)+
  ylab("Distance of insertion donor\nto the closest rDNA (Mb)")+
  xlab("Insert type")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  annotate('text', x = "V-DJ", y = 5,
           label = ifelse(as.numeric(wilcox.matrix["V-DJ","Generated"])<=0.05,
                          wilcox.matrix["V-DJ","Generated"], NA),
           size = 4)+
  annotate('text', x = "VDJ-CH1", y = 5,
           label = ifelse(as.numeric(wilcox.matrix["VDJ-CH1","Generated"])<=0.05,
                          wilcox.matrix["VDJ-CH1","Generated"], NA),
           size = 4)
#annotate('text', x = 4, y = 10, label = 'p = 6.51e-12', size = 5)
rdna_dist_boxplot



rdna_dist_density <- ggplot(data = inserts.all)+
  geom_density_ridges(aes(x = rdna.dist.mpb, y = insert.type, fill = insert.type), alpha = 0.2)+
  labs(x = "Distance to the closest rDNA, Mbp", fill = 'Insert type')+
  theme_classic()+
  xlim(0, 1.2e+01)
rdna_dist_density

## Saving the plots
ggsave(paste('rdna_dist_density_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = rdna_dist_density)

ggsave(paste('rdna_dist_boxplot_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = rdna_dist_boxplot)

if(delete_the_plots_from_R_memory == TRUE){rm(rdna_dist_density, rdna_dist_boxplot)}
rm(wilcox.matrix)



###################################
## 30072019
## Calculating the distance to the closest genes
dist.to.closest.agtotal <- c()
overlap.closest.agtotal <- c()
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



for (i in seq(nrow(inserts.all))){
  agtotal.in.this.chrom.start <- agtotal$start[agtotal$chrom == inserts.all$chromosome[i]]
  agtotal.in.this.chrom.end <- agtotal$end[agtotal$chrom == inserts.all$chromosome[i]]
  agtotal.in.this.chrom.center <- (agtotal.in.this.chrom.start + agtotal.in.this.chrom.end)/2
  agtotal.dist <- min(abs(agtotal.in.this.chrom.center - inserts.all$insert.genomic.center[i]))
  agtotal.overlap <- count.overlap(agtotal.in.this.chrom.start, agtotal.in.this.chrom.end,
                             inserts.all$insert.genomic.coord.s.start[i],
                             inserts.all$insert.genomic.coord.s.end[i])
  dist.to.closest.agtotal <- c(dist.to.closest.agtotal, agtotal.dist)
  overlap.closest.agtotal <- c(overlap.closest.agtotal, agtotal.overlap)
  if(i %in% seq(0,nrow(inserts.all),10000)){
    paste(round(i*100/nrow(inserts.all)), "% of the job is done", sep = "") %>% print()
  }
}

inserts.all$transcribed.dist <- dist.to.closest.agtotal
inserts.all$transcribed.overlap <- overlap.closest.agtotal
inserts.all$transcribed.dist.mpb <- inserts.all$transcribed.dist / 1000000

rm(dist.to.closest.degl, overlap.closest.degl, degl.in.this.chrom.center,
   degl.in.this.chrom.start, degl.in.this.chrom.end, degl.overlap, degl.dist)

# P-VALUES COMPUTATION
wilcox.matrix <- matrix(data = NA, nrow = uni.ins.types.n, ncol = uni.ins.types.n,
                        dimnames = list(unique.insert.types, unique.insert.types))

for(i in colnames(wilcox.matrix)){
  for(j in rownames(wilcox.matrix)){
    wilcox.test(inserts.all$transcribed.dist[inserts.all$insert.type == i],
                inserts.all$transcribed.dist[inserts.all$insert.type == j],
                mu = 0, alternative = "two.sided",
                paired = F, correct = F, exact = T)$p.value %>%
      format(digits = 2, scientific = TRUE) -> wilcox.matrix[i,j]
  }
}
rm(i,j)


melt(wilcox.matrix) %>% mutate(value = as.numeric(levels(value))[value]) %>%
  filter(Var1 == 'Generated', Var2 != 'Generated') -> wilcox.matrix 
wilcox.matrix$significance <- "ns"
wilcox.matrix$significance[wilcox.matrix$value <= 0.05] <- "*"
wilcox.matrix$significance[wilcox.matrix$value <= 0.01] <- "**"
wilcox.matrix$significance[wilcox.matrix$value <= 0.001] <- "***"
wilcox.matrix$significance[wilcox.matrix$value <= 0.0001] <- "****"


genes_dist_boxplot <- ggplot(data = inserts.all)+
  geom_boxplot(aes(x = insert.type, y = transcribed.dist), alpha = 0.2)+
  theme_classic()+
  #xlab("")+
  #scale_x_discrete(labels = c('Random' = 'Random\nInsertions', 'Real' = 'Donors\nInsertions'))+
  ylim(0, 75000)+
  ylab("Distance of insertion donor\nto the closest transcribed locus (bp)")+
  xlab("Insert type")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  annotate('text', x = "V-DJ", y = 75000,
           label = ifelse(wilcox.matrix$significance[wilcox.matrix$Var2 == 'V-DJ'] != 'ns',
                          paste(wilcox.matrix$significance[wilcox.matrix$Var2 == 'V-DJ'],
                                "\np = ",wilcox.matrix$value[wilcox.matrix$Var2 == 'V-DJ'], sep =""), NA),
           size = 4)+
  annotate('text', x = "VDJ-CH1", y = 75000,
           label = ifelse(wilcox.matrix$significance[wilcox.matrix$Var2 == 'VDJ-CH1'] != 'ns',
                          paste(wilcox.matrix$significance[wilcox.matrix$Var2 == 'VDJ-CH1'],
                                "\np = ",wilcox.matrix$value[wilcox.matrix$Var2 == 'VDJ-CH1'], sep =""), NA),
           size = 4)
#annotate('text', x = 4, y = 10, label = 'p = 6.51e-12', size = 5)
genes_dist_boxplot

ggsave(paste('transcribed_dist_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 8, height = 5,
       device = "pdf",
       plot = genes_dist_boxplot)

median(inserts.all$degldist[inserts.all$insert.type == 'Generated'])

mean(inserts.all$degldist[inserts.all$insert.type == 'V-DJ'])
median(inserts.all$degldist[inserts.all$insert.type == 'VDJ-CH1'])


ggplot(data = inserts.all)+
  geom_density_ridges(aes(x = degldist, y = insert.type, fill = insert.type), alpha = 0.3)+
  xlim(0,100000)+
  theme_classic()+
  xlab('Distance to the closest gene, bp')+
  ylab('Insert type')

ggplot(data = inserts.all)+
  geom_violin(aes(y = degldist, x = insert.type, fill = insert.type), alpha = 0.3, adjust = 0.5)+
  ylim(0,100000)+
  theme_classic()+
  ylab('Distance to the closest gene, bp')+
  xlab('Insert type')

closest.degl.names <- c()

degl %>% mutate(center = (start + end)/2) -> degl

for (i in seq(nrow(inserts.all))){
  degl %>% filter(chrom == inserts.all$chromosome[i]) ->
    degl.temp
  
  degl.dist <- min(abs(degl.temp$center - inserts.all$insert.genomic.center[i]))
  
  closest.degl.id <- which(abs(degl.temp$center - inserts.all$insert.genomic.center[i]) ==
                               degl.dist)[1]
  closest.degl.name <- degl.temp$Name[closest.degl.id]
  closest.degl.names <- c(closest.degl.names, closest.degl.name)
  #nbcs.closest.degl <- c(nbcs.closest.degl, arow$NBCs)
  #cbvsnbcs.closest.degl <- c(cbvsnbcs.closest.degl, arow$CBs_vsN)
  if(i %in% seq(0,100000,1000)){
    print(i)
  }
}


inserts.all$closestdegl <- closest.degl.names
inserts.all %>% rename('closestdegl' = 'Name') -> inserts.all

inserts.all %>% left_join(select(degl, -c('start','end','chrom','center')), by = 'Name') ->
  inserts.all

rm(closest.degl.name, closest.degl.names, degl.temp, degl.dist, closest.degl.id, i)


d <- ggplot(data = inserts.all)+
  geom_boxplot(aes(x = insert.type, y = NBCs), alpha = 0.2)+
  #geom_boxplot(aes(x = insert.type, y = CBs), alpha = 0.2)+
  theme_classic()+
  xlab("")+
  ylab("")+
  
  scale_x_discrete(labels = c())+
  ylim(-10, 10)+
  #ylab("Expression difference\nof the closest gene, log2fold")+
  #xlab("Insert type")+
  theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 0.6))
  #annotate('text', x = 2, y = 1, label = paste('p = ',p.GenvsVDJ), size = 4)+
  #annotate('text', x = 3, y = 1, label = paste('p = ',p.GenvsVDJCH1), size = 4)
#annotate('text', x = 4, y = 10, label = 'p = 6.51e-12', size = 5)
d

d14 <- ggplot(data = inserts.all)+
  geom_boxplot(aes(x = insert.type, y = BMPCs_vsN), alpha = 0.2)+
  #geom_boxplot(aes(x = insert.type, y = CBs), alpha = 0.2)+
  theme_classic()+
  #xlab("")+
  #scale_x_discrete(labels = c('Random' = 'Random\nInsertions', 'Real' = 'Donors\nInsertions'))+
  ylim(-10, 10)+
  #ylab("Expression difference\nof the closest gene, log2fold")+
  #xlab("Insert type")+
  #theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))
#annotate('text', x = 2, y = 1, label = paste('p = ',p.GenvsVDJ), size = 4)+
#annotate('text', x = 3, y = 1, label = paste('p = ',p.GenvsVDJCH1), size = 4)
#annotate('text', x = 4, y = 10, label = 'p = 6.51e-12', size = 5)
d14

col_to_test <- c('NBCs','CBs','CCs','MBCs','prePBs','PBs','EPCs','BMPCs',
                 'CBs_vsN','CCs_vsN','MBCs_vsN','prePBs_vsN','PBs_vsN','EPCs_vsN','BMPCs_vsN')

wilcox_genvsvdj <- c()
wilcox_genvsvdjch1 <- c()
wilcox_vdjvsvdjch1 <- c()

for(i in col_to_test){
  wilcox_genvsvdjch1 <-
    c(wilcox_genvsvdjch1,
      wilcox.test(inserts.all[inserts.all$insert.type == 'Generated',i],
                  inserts.all[inserts.all$insert.type == 'VDJ-CH1',i],
                  mu = 0, alternative = "two.sided",
                  paired = F, correct = F, exact = T)$p.value %>%
        format(digits = 2, scientific = TRUE))
  wilcox_genvsvdj <-
    c(wilcox_genvsvdj,
      wilcox.test(inserts.all[inserts.all$insert.type == 'Generated',i],
                  inserts.all[inserts.all$insert.type == 'V-DJ',i],
                  mu = 0, alternative = "two.sided",
                  paired = F, correct = F, exact = T)$p.value %>%
        format(digits = 2, scientific = TRUE))
  wilcox_vdjvsvdjch1 <-
    c(wilcox_vdjvsvdjch1,
      wilcox.test(inserts.all[inserts.all$insert.type == 'V-DJ',i],
                  inserts.all[inserts.all$insert.type == 'VDJ-CH1',i],
                  mu = 0, alternative = "two.sided",
                  paired = F, correct = F, exact = T)$p.value %>%
        format(digits = 2, scientific = TRUE))
}


wilcox.matrix <- data.frame(population = col_to_test,
                             gen.vs.vdj = wilcox_genvsvdj,
                             gen.vs.vdj.sign = "ns", #as.numeric(wilcox_genvsvdj) <= 0.05,
                             gen.vs.vdjch1 = wilcox_genvsvdjch1,
                             gen.vs.vdjch1.sign = "ns", #as.numeric(wilcox_genvsvdjch1) <= 0.05,
                             vdj.vs.vdjch1 = wilcox_vdjvsvdjch1,
                             vdj.vs.vdjch1.sign = "ns",
                            stringsAsFactors = FALSE) #as.numeric(wilcox_vdjvsvdjch1) <= 0.05)

wilcox.matrix$gen.vs.vdjch1 <- as.numeric(as.character(wilcox.matrix$gen.vs.vdjch1))
wilcox.matrix$gen.vs.vdj <- as.numeric(as.character(wilcox.matrix$gen.vs.vdj))
wilcox.matrix$vdj.vs.vdjch1 <- as.numeric(as.character(wilcox.matrix$vdj.vs.vdjch1))


wilcox.matrix$gen.vs.vdjch1.sign[wilcox.matrix$gen.vs.vdjch1 <= 0.05] <- "*"
wilcox.matrix$gen.vs.vdjch1.sign[wilcox.matrix$gen.vs.vdjch1 <= 0.01] <- "**"
wilcox.matrix$gen.vs.vdjch1.sign[wilcox.matrix$gen.vs.vdjch1 <= 0.001] <- "***"
wilcox.matrix$gen.vs.vdjch1.sign[wilcox.matrix$gen.vs.vdjch1 <= 0.0001] <- "****"

wilcox.matrix$gen.vs.vdj.sign[wilcox.matrix$gen.vs.vdj <= 0.05] <- "*"
wilcox.matrix$gen.vs.vdj.sign[wilcox.matrix$gen.vs.vdj <= 0.01] <- "**"
wilcox.matrix$gen.vs.vdj.sign[wilcox.matrix$gen.vs.vdj <= 0.001] <- "***"
wilcox.matrix$gen.vs.vdj.sign[wilcox.matrix$gen.vs.vdj <= 0.0001] <- "****"

wilcox.matrix$vdj.vs.vdjch1.sign[wilcox.matrix$vdj.vs.vdjch1 <= 0.05] <- "*"
wilcox.matrix$vdj.vs.vdjch1.sign[wilcox.matrix$vdj.vs.vdjch1 <= 0.01] <- "**"
wilcox.matrix$vdj.vs.vdjch1.sign[wilcox.matrix$vdj.vs.vdjch1 <= 0.001] <- "***"
wilcox.matrix$vdj.vs.vdjch1.sign[wilcox.matrix$vdj.vs.vdjch1 <= 0.0001] <- "****"


ggarrange(d, d1, d2, d3, d4, d5, d6, d7,
          labels = c('NBCs','CBs','CCs','MBCs','prePBs','PBs','EPCs','BMPCs'),
          ncol = 4, nrow = 2)



inserts.all.melted <- melt(inserts.all, id = c('run.id','sample.id','insert.id.s.',
                                               'ins.id','contig.id','Name', 'insert.type'),
                           measure.vars = col_to_test)
inserts.all.melted %>% rename("variable" = "population", "value" = "expression") -> inserts.all.melted

wilcox.matrix %>% select(-c(vdj.vs.vdjch1,vdj.vs.vdjch1.sign, gen.vs.vdj, gen.vs.vdjch1)) %>%
  rename("gen.vs.vdj.sign" = "V-DJ", "gen.vs.vdjch1.sign" = "VDJ-CH1") %>%
  melt(id = c('population')) %>%
  rename("variable" = "insert.type", "value" = "significance") -> wilcox.matrix.molten


inserts.all.melted1 <- left_join(inserts.all.melted, wilcox.matrix.molten, by = c("insert.type","population"))

inserts.all.melted1 <- transform(inserts.all.melted1,
                                 population=factor(population,
                                                   levels =
                                                     levels(inserts.all.melted$population)))

diff_exp_closest_gene <-
  ggplot(data = inserts.all.melted1)+
  geom_boxplot(aes(x = insert.type, y = expression, fill = insert.type), alpha = 0.2, outlier.shape = NA)+
  geom_text(aes(x = insert.type, y = 4.5,label = significance),
            data = distinct(inserts.all.melted1, population, insert.type, .keep_all = TRUE)%>%
              filter(significance != 'ns'))+
  theme_classic()+
  xlab("")+
  ylim(-5, 5)+
  ylab("Expression difference\nof the closest gene, log2fold")+
  xlab("Insert type")+
  facet_wrap( ~ population, ncol = 4)
diff_exp_closest_gene

ggsave(paste('proximal_expression_boxplots_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 12, height = 6.75,
       device = "pdf",
       plot = diff_exp_closest_gene)

if(delete_the_plots_from_R_memory == TRUE){rm(diff_exp_closest_gene)}
rm(wilcox.matrix)

