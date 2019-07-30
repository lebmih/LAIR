### Randomize the insertions
## Plan
## We need to take the insertions table and change the chromosomal coordinate to
## random position inside the same chromosome
## Then we need to calculate overlaps for our true data and for 1000 generations
## of random data

## For the start let's just plot the random insertions
# Replace the insert.genomic.center with random number
# between 1 and chromosome length
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
}

inserts.generic$telomeric.distance <- mapply(min,
                                             inserts.generic$insert.genomic.coord.s.start - telomeres.length,
                                             inserts.generic$Length.bp - telomeres.length - inserts.generic$insert.genomic.coord.s.end)

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

inserts$telomeric.distance <- mapply(min,
                                     inserts$insert.genomic.coord.s.start - telomeres.length,
                                     inserts$Length.bp - telomeres.length - inserts$insert.genomic.coord.s.end)

inserts$insert.type[is.na(inserts$insert.type)] <- 'V-CH1'

inserts.vdjch1 <- filter(inserts, insert.type == 'VDJ-CH1')
inserts.vdj <- filter(inserts, insert.type == 'V-DJ')
inserts.vch1 <- filter(inserts, insert.type == 'V-CH1')

#teldistanceTotal <- data.frame(inserts = 'Total', telomeric.distance = inserts$telomeric.distance)
teldistanceVDJCH1 <- data.frame(inserts = 'VDJ-CH1', telomeric.distance = inserts.vdjch1$telomeric.distance)
teldistanceVDJ <- data.frame(inserts = 'V-DJ', telomeric.distance = inserts.vdj$telomeric.distance)
teldistanceVCH1 <- data.frame(inserts = 'V-CH1', telomeric.distance = inserts.vch1$telomeric.distance)


teldistanceGenerated <- data.frame(inserts = 'Generated', telomeric.distance = inserts.generic$telomeric.distance)

#teldistance <- bind_rows(teldistanceVDJCH1, teldistanceVDJ,
#                         teldistanceVCH1, teldistanceGenerated)
## 280719
teldistance <- bind_rows(teldistanceVDJCH1, teldistanceVDJ, teldistanceGenerated)


# H0: Median Telomeric Distance of random insertions equals to that of real insertions
wilcox.test(teldistanceVDJCH1$telomeric.distance, teldistanceGenerated$telomeric.distance, mu = 0, alternative = "less",
            paired = F, correct = F, exact = T)


teldistance$telomeric.distance.Mbp <- teldistance$telomeric.distance/1000000

d <- ggplot(data = filter(teldistance, !inserts == 'Total'))+
  geom_boxplot(aes(x = inserts, y = telomeric.distance.Mbp), alpha = 0.2)+
  theme_classic()+
  #xlab("")+
  #scale_x_discrete(labels = c('Random' = 'Random\nInsertions', 'Real' = 'Donors\nInsertions'))+
  ylab("Distance of insertion donor\nto the closest telomere (Mb)")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  xlab("Insert type")+
  annotate('text', x = 2, y = 130, label = 'p = 0.86', size = 4)+
  annotate('text', x = 3, y = 130, label = 'p < 2.2e-16', size = 4)
  #annotate('text', x = 4, y = 130, label = 'p < 2.2e-16', size = 4)
d



dnp <- ggplot(data = filter(teldistance, !inserts == 'Total'))+
  geom_density(aes(x = telomeric.distance.Mbp, y = ..density.., fill = inserts), alpha = 0.2)+
  labs(x = "Distance to the closest telomere, Mbp", fill = 'Insert type')+
  theme_classic()
dnp


## 20.02.19
gaps <- read.table(file = 'hg19chrGaps.txt', sep = '\t', header = TRUE)


## 21.02.19
## Let's find the distance between rdna and the inserts
rdna <- read.table(file = 'rdna_bed.txt', sep = '\t', header = TRUE)
unique(rdna$feature)
# drop score as we don't need it
rdna <- select(rdna, -score)

## Calculating the center of rDNA regions
rdna$chrom.center <- (rdna$chrom.end + rdna$chrom.start)/2
rdna <- filter(rdna, !str_detect(chromosome, '_'))
rdna <- filter(rdna, !str_detect(chromosome, 'chrM'))
unique(rdna$chromosome)

dist.to.closest.rdna <- c()
for (i in seq_along(inserts$run.id)){
  rdna.in.this.chrom <- rdna$chrom.center[rdna$chromosome == inserts$chromosome[i]]
  rdna.dist <- min(abs(rdna.in.this.chrom - inserts$insert.genomic.center[i]))
  dist.to.closest.rdna <- c(dist.to.closest.rdna, rdna.dist)
}

inserts$rdna.distance <- dist.to.closest.rdna

## Calculating the distance to the closest rdna for generated data
dist.to.closest.rdna <- c()
inserts.generic1 <- inserts.generic
inserts.generic <- inserts.generic1[1:530100,]
inserts.generic <- inserts.generic1
for (i in seq_along(inserts.generic$run.id)){
  rdna.in.this.chrom <- rdna$chrom.center[rdna$chromosome == inserts.generic$chromosome[i]]
  rdna.dist <- min(abs(rdna.in.this.chrom - inserts.generic$insert.genomic.center[i]))
  dist.to.closest.rdna <- c(dist.to.closest.rdna, rdna.dist)
}

inserts.generic$rdna.distance <- dist.to.closest.rdna
inserts.generic$insert.type <- 'Generated'

inserts.all <- bind_rows(inserts, inserts.generic)
inserts.all$rdna.dist.mpb <- inserts.all$rdna.distance / 1000000

install.packages("ggridges")
library(ggridges)

rdnaplot <- ggplot(data = inserts.all)+
  geom_density_ridges(aes(x = rdna.dist.mpb, y = insert.type, fill = insert.type), alpha = 0.2)+
  labs(x = "Distance to the closest rDNA, Mbp", fill = 'Insert type')+
  theme_classic()+
  xlim(0, 1.2e+01)
rdnaplot

wilcox.test(inserts.all$rdna.distance[inserts.all$insert.type == 'Generated'],
            inserts.all$rdna.distance[inserts.all$insert.type == 'V-DJ'],
            mu = 0, alternative = "greater",
            paired = F, correct = F, exact = T)

#GeneratedVsVCH1 <- 2.591e-09
GeneratedVsVDJ <- 0.01387
GeneratedVsVDJCH1 <- 0.01386

d <- ggplot(data = inserts.all)+
  geom_boxplot(aes(x = insert.type, y = rdna.dist.mpb), alpha = 0.2)+
  theme_classic()+
  #xlab("")+
  #scale_x_discrete(labels = c('Random' = 'Random\nInsertions', 'Real' = 'Donors\nInsertions'))+
  ylim(0, 10)+
  ylab("Distance of insertion donor\nto the closest rDNA (Mb)")+
  xlab("Insert type")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  annotate('text', x = 2, y = 10, label = 'p = 0.014', size = 4)+
  annotate('text', x = 3, y = 10, label = 'p = 0.014', size = 4)
  #annotate('text', x = 4, y = 10, label = 'p = 6.51e-12', size = 5)
d

##
###################################
## 30072019
## Calculating the distance to the closest 10% overexpressed genes
dist.to.closest.degl <- c()
overlap.closest.degl <- c()
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
  degl.in.this.chrom.start <- degl$start[degl$chromosome == inserts.all$chromosome[i]]
  degl.in.this.chrom.end <- degl$end[degl$chromosome == inserts.all$chromosome[i]]
  degl.in.this.chrom.center <- (degl.in.this.chrom.start + degl.in.this.chrom.end)/2
  degl.dist <- min(abs(degl.in.this.chrom.center - inserts.all$insert.genomic.center[i]))
  degl.overlap <- count.overlap(degl.in.this.chrom.start, degl.in.this.chrom.end,
                             inserts.all$insert.genomic.coord.s.start[i],
                             inserts.all$insert.genomic.coord.s.end[i])
  dist.to.closest.degl <- c(dist.to.closest.degl, degl.dist)
  overlap.closest.degl <- c(overlap.closest.degl, degl.overlap)
}

inserts.all$degldist <- dist.to.closest.degl
inserts.all$degloverlap <- overlap.closest.degl


inserts.all$degldist.mpb <- inserts.all$degldist / 1000000


# P-VALUES COMPUTATION
wilcox.test(inserts.all$degldist.mpb[inserts.all$insert.type == 'V-DJ'],
            inserts.all$degldist.mpb[inserts.all$insert.type == 'VDJ-CH1'],
            mu = 0, alternative = "greater",
            paired = F, correct = F, exact = T)$p.value %>%
  format(digits = 2, scientific = TRUE) -> p.VDJvsVDJCH1
wilcox.test(inserts.all$degldist.mpb[inserts.all$insert.type == 'Generated'],
            inserts.all$degldist.mpb[inserts.all$insert.type == 'V-DJ'],
            mu = 0, alternative = "greater",
            paired = F, correct = F, exact = T)$p.value %>%
  format(digits = 2, scientific = TRUE) -> p.GenvsVDJ
wilcox.test(inserts.all$degldist.mpb[inserts.all$insert.type == 'Generated'],
            inserts.all$degldist.mpb[inserts.all$insert.type == 'VDJ-CH1'],
            mu = 0, alternative = "greater",
            paired = F, correct = F, exact = T)$p.value %>%
  format(digits = 2, scientific = TRUE) -> p.GenvsVDJCH1



d <- ggplot(data = inserts.all)+
  geom_boxplot(aes(x = insert.type, y = degldist.mpb), alpha = 0.2)+
  theme_classic()+
  #xlab("")+
  #scale_x_discrete(labels = c('Random' = 'Random\nInsertions', 'Real' = 'Donors\nInsertions'))+
  ylim(0, 1)+
  ylab("Distance of insertion donor\nto the closest gene (Mb)")+
  xlab("Insert type")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  annotate('text', x = 2, y = 1, label = paste('p = ',p.GenvsVDJ), size = 4)+
  annotate('text', x = 3, y = 1, label = paste('p = ',p.GenvsVDJCH1), size = 4)
#annotate('text', x = 4, y = 10, label = 'p = 6.51e-12', size = 5)
d

median(inserts.all$degldist[inserts.all$insert.type == 'Generated'])

mean(inserts.all$degldist[inserts.all$insert.type == 'V-DJ'])

ggplot(data = inserts.all)+
  geom_density_ridges(aes(x = degloverlap, y = insert.type))+
  xlim(0,2)



nbcs.closest.degl <- c()
cbvsnbcs.closest.degl <- c()
closest.degl.names <- c()

degl %>% mutate(center = (start + end)/2) -> degl
degl.ids.list <- c()

for (i in seq(nrow(inserts.all))){
  degl %>% filter(chromosome == inserts.all$chromosome[i]) ->
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


inserts.all$degldist <- dist.to.closest.degl
inserts.all$degloverlap <- overlap.closest.degl
inserts.all$closestdegl <- closest.degl.names
inserts.all %>% rename('closestdegl' = 'Name') -> inserts.all

inserts.all %>% left_join(select(degl, -c('start','end','chromosome','rel.start',
                                          'rel.end','AY','BY','AX','BX', 'center')), by = 'Name') ->
  inserts.all

library(cowplot)
install.packages('ggpubr')
library(ggpubr)

d <- ggplot(data = inserts.all)+
  geom_boxplot(aes(x = insert.type, y = NBCs), alpha = 0.2)+
  #geom_boxplot(aes(x = insert.type, y = CBs), alpha = 0.2)+
  theme_classic()+
  #xlab("")+
  #scale_x_discrete(labels = c('Random' = 'Random\nInsertions', 'Real' = 'Donors\nInsertions'))+
  ylim(-10, 10)+
  #ylab("Expression difference\nof the closest gene, log2fold")+
  xlab("Insert type")+
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
  xlab("Insert type")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))
#annotate('text', x = 2, y = 1, label = paste('p = ',p.GenvsVDJ), size = 4)+
#annotate('text', x = 3, y = 1, label = paste('p = ',p.GenvsVDJCH1), size = 4)
#annotate('text', x = 4, y = 10, label = 'p = 6.51e-12', size = 5)
d14

col_to_test <- c('NBCs','CBs','CCs','MBCs','prePBs','PBs','EPCs','BMPCs',
                 'CBs_vsN','CCs_vsN','MBCs_vsN','prePBs_vsN','PBs_vsN','EPCs_vsN','BMPCs_vsN')

wilcox_genvsvdjch1 <- c()
wilcox_genvsvdj <- c()
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

as.numeric(wilcox_genvsvdj) < 0.05

wilcox_testing <- data.frame(population = col_to_test,
                             gen.vs.vdj = wilcox_genvsvdj,
                             gen.vs.vdj.sign = as.numeric(wilcox_genvsvdj) < 0.05,
                             gen.vs.vdjch1 = wilcox_genvsvdjch1,
                             gen.vs.vdjch1.sign = as.numeric(wilcox_genvsvdjch1) < 0.05,
                             vdj.vs.vdjch1 = wilcox_vdjvsvdjch1,
                             vdj.vs.vdj.sign = as.numeric(wilcox_vdjvsvdjch1) < 0.05)
t.test(inserts.all[inserts.all$insert.type == 'V-DJ','NBCs'],
       inserts.all[inserts.all$insert.type == 'Generated','NBCs'])

ggarrange(d, d1, d2, d3, d4, d5, d6, d7,
          labels = c('NBCs','CBs','CCs','MBCs','prePBs','PBs','EPCs','BMPCs'),
          ncol = 4, nrow = 2)

library(reshape2)

inserts.all.melted <- melt(inserts.all, id = c('run.id','sample.id','insert.id.s.',
                                               'ins.id','contig.id','Name', 'insert.type'),
                           measure.vars = col_to_test)

d99 <- ggplot(data = inserts.all.melted)+
  geom_boxplot(aes(x = insert.type, y = value), alpha = 0.2)+
  #geom_boxplot(aes(x = insert.type, y = CBs), alpha = 0.2)+
  theme_classic()+
  #xlab("")+
  #scale_x_discrete(labels = c('Random' = 'Random\nInsertions', 'Real' = 'Donors\nInsertions'))+
  ylim(-10, 10)+
  #ylab("Expression difference\nof the closest gene, log2fold")+
  xlab("Insert type")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  facet_wrap( ~ variable, ncol = 4)
#annotate('text', x = 2, y = 1, label = paste('p = ',p.GenvsVDJ), size = 4)+
#annotate('text', x = 3, y = 1, label = paste('p = ',p.GenvsVDJCH1), size = 4)
#annotate('text', x = 4, y = 10, label = 'p = 6.51e-12', size = 5)
d99
