## Centromeric distance

inserts <- left_join(inserts, select(chrlen, chromosome, chromCenter), by = 'chromosome')

## 290719
inserts.all <- left_join(inserts.all, select(chrlen, chromosome, chromCenter), by = 'chromosome')


centromeres.length <- 10000
centromeric.distance <- c()
for (j in seq(1:1000)){
  inserts1 <- inserts
  rand_position <- c()
  for (i in seq_along(inserts1$chromosome)){
    rand_position <- c(rand_position, sample(1:inserts1$Length.bp[i],1))
  }
  inserts1$insert.genomic.center <- rand_position
  inserts1$insert.genomic.coord.s.start <- inserts1$insert.genomic.center - inserts1$insert.length/2
  inserts1$insert.genomic.coord.s.end <- inserts1$insert.genomic.center + inserts1$insert.length/2
  inserts1$centromeric.distance <- mapply(min,
                                          abs(inserts1$chromCenter - centromeres.length/2 - inserts1$insert.genomic.coord.s.end),
                                          abs(inserts1$chromCenter + centromeres.length/2 - inserts1$insert.genomic.coord.s.start))
  centromeric.distance <- c(centromeric.distance, inserts1$centromeric.distance)
}

inserts.all$centromeric.distance <- mapply(min,
                                       abs(inserts.all$chromCenter - centromeres.length/2 - inserts.all$insert.genomic.coord.s.end),
                                       abs(inserts.all$chromCenter + centromeres.length/2 - inserts.all$insert.genomic.coord.s.start))
centrdistanceReal <- data.frame(inserts = 'Real', centromeric.distance = inserts$centromeric.distance)


centrdistanceRandom <- data.frame(inserts = 'Random', centromeric.distance = centromeric.distance)

centrdistance <- bind_rows(centrdistanceReal, centrdistanceRandom)

# H0: Median Telomeric Distance of random insertions equals to that of real insertions
wilcox.test(inserts.all$centromeric.distance[inserts.all$insert.type == 'VDJ-CH1'],
            inserts.all$centromeric.distance[inserts.all$insert.type == 'Generated'],
            mu = 0, alternative = "greater",
            paired = F, correct = F, exact = T)

inserts.all$centromeric.distance <- inserts.all$centromeric.distance/1000000

d <- ggplot(data = inserts.all)+
  geom_boxplot(aes(insert.type, centromeric.distance))+
  theme_classic()+
  xlab("Insert type")+
  scale_x_discrete(labels = c('Generated' = 'Generated', 'V-DJ' = 'V-DJ',
                              'VDJ-CH1' = 'VDJ-CH1'))+
  ylab("Distance of insertion donor\nto the closest centromere (Mb)")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  annotate('text', x = 2, y = 160, label = 'P = 4.3e-08', size = 4)+
  annotate('text', x = 3, y = 160, label = 'P < 2.2e-16', size = 4)
d


