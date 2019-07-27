install.packages('tidyr')
install.packages('dplyr')
install.packages('ggplot2')
install.packages('stringi')
install.packages('ggridges')
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("ShortRead"))
library(ShortRead)
library(tidyr)
library(dplyr)
library(ggplot2)
library(stringi)
library(stringr)
library(ggridges)
#read.table(file ='C:\Users\mlebedi\Documents\minion_FAA\minion_forMisha.tar\minion_forMisha\A1.fastq')
reads <- readFastq("C:/Users/mlebedi/Documents/minion_FAA/minion_forMisha.tar/minion_forMisha/A1.fastq")
browseVignettes("ShortRead")
mean(alphabetScore(readsB2))
install.packages('microseq')
library(microseq)


list.files('C:/Users/mlebedi/Documents/minion_FAA/minion_forMisha.tar/minion_forMisha/') -> files_fastq

total_table <- data.frame()

for(i in files_fastq){
  table <- readFastq(i)
  table$sample <- i
  bind_rows(table, total_table) -> total_table
}

total_table %>% separate(col = sample, into = c('sampleid', 'fastq'), sep = '\\.') %>% select(-fastq) -> total_table1

readsA1 <- readFastq("C:/Users/mlebedi/Documents/minion_FAA/minion_forMisha.tar/minion_forMisha/A1.fastq")
readsA2 <- readFastq("C:/Users/mlebedi/Documents/minion_FAA/minion_forMisha.tar/minion_forMisha/A2.fastq")
readsB1 <- readFastq("C:/Users/mlebedi/Documents/minion_FAA/minion_forMisha.tar/minion_forMisha/B1.fastq")
readsB2 <- readFastq("C:/Users/mlebedi/Documents/minion_FAA/minion_forMisha.tar/minion_forMisha/B2.fastq")


readsA1$sample <- 'A1'
readsA2$sample <- 'A2'
readsB1$sample <- 'B1'
readsB2$sample <- 'B2'

bind_rows(readsA1,readsA2,readsB1,readsB2) -> total_minion_reads

total_minion_reads$Sequence %>% str_length() -> total_minion_reads$read_length

g <- ggplot(data = total_minion_reads)+
  geom_density(aes(read_length, fill = sample), adjust = 0.4, alpha = 0.7)+
  xlab('Read length, bp')+
  ylab('Reads fraction')+
  scale_x_continuous(limits = c(0,1000))+
  theme_classic()+
  scale_fill_discrete(name = 'Sample', labels = c('QSOX2 Background', 'QSOX2 Specific', 'LAIR1 Background', 'LAIR1 Specific'))+
  theme(legend.position = 'bottom')
g

g <- ggplot(data = total_minion_reads%>% filter(sample %in% c('A2','B2')))+
  geom_density_ridges(aes(read_length, y = sample, fill = sample), alpha = 0.7)+
  xlab('Read length, bp')+
  ylab('Reads %')+
  scale_x_continuous(limits = c(0,1000))+
  theme_classic()+
  scale_fill_discrete(name = 'Sample', labels = c('QSOX2 Background', 'QSOX2 Specific', 'LAIR1 Background', 'LAIR1 Specific'))+
  theme(legend.position = 'bottom')
g

g <- ggplot(data = total_minion_reads %>% filter(sample == 'A2'))+
  geom_histogram(aes(read_length, fill = sample),binwidth = 10, alpha = 0.7)+
  xlab('Read length, bp')+
  ylab('Read count')+
  scale_x_continuous(limits = c(0,1000))+
  theme_classic()
g



barcodes <- c('GGTTTTCACCGAGTCCATGGATTCTTTAGCCATTTCCTCAAAGTGCTCAC',
              'GGTTTTCACCGAGTCCATGG',
              'CAGTAAGAGAAGGAAATGCCGGGCTTTATCGCTGCATCTATTATAAGCCC',
              'CAGTAAGAGAAGGAAATGCCGGGCTTTATCGC')
probes <- c('CTGCGCAGAGTGACTGCCTG','GAGTTGCAGGTGAAAGCCCG','AGATGTCCTGGAGTCAAAGT','CCACTGGAGAGGCTCCAATC',
            'CACCTTCGTGTCCCCCCTGT', 'TCGGGGAGTTATTGGCTACT',
            'GCTAGCCCTCCACCAAGGGC', 'CCATCGGTCTTCCCCCTGGC')
probesRC <- probes %>% reverseComplement()

probesOri <- c('AAAGCTCTGGAGGCCCGGAC','CCCGGACTCCCCGGACACAG', 'AGCCCGGCTCCTCAGCTGGA', 'CCCACGCAGAGGCCGTCGGA')
probesOriRC <- probesOri %>% reverseComplement()

adapter <- c('GACTATAGGGCACGCGTGG', 'CCACGCGTGCCCTATAGTC')

for(i in probesOriRC){
  print(readsB2$Sequence %>% str_detect(i) %>% sum())
}

readsA2$Sequence %>% str_detect(adapter[1]) %>% sum()



library(Biostrings)
library(rBLAST)





