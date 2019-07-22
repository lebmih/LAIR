## Timestamp
script_started <- Sys.time()

## Copy-paste the address of your working directory here
setwd("C:/Users/User/Documents/Lab/LongBCRs/MethodPaper/RunTables")

## Libraries
x <- c('tidyr', 'dplyr', 'Hmisc') #'ggplot2', 'ggrepel', 'reshape2', 'stringr', 'Biostrings', 'seqinr')
lapply(x, require, character.only = TRUE)
rm(x) # always clear the environment after such procedures

###########################
###### TIMESTAMP (CODE WRITTEN): 22.07.19

logfile <- c(paste("Starting time:", Sys.time()))

## Importing the last table
all_inserts <- read.table(file = 'Output/inserts_allRuns_22_07_19.txt', sep = '\t',
                          header = TRUE, stringsAsFactors = FALSE)
logfile <- c(logfile,
             paste("Imported table has", nrow(all_inserts), "observations"))

##----------------------------
describe(all_inserts$spike.in)
##----------------------------

## Positive filtering for spike-in contigs
inserts_with_spikeins <- all_inserts[!is.na(all_inserts$spike.in),]

logfile <- c(logfile,
             paste(nrow(inserts_with_spikeins), "observations come from SpikeIn samples"))

## Stupido Excel changes the 1:1000 to '0,736111111', fixing it
logfile <- c(logfile,
             paste(nrow(inserts_with_spikeins[inserts_with_spikeins$spike.in == '0,736111111',]),
                   "instances of Excel changing the spike.in string fixed"))

inserts_with_spikeins$spike.in[inserts_with_spikeins$spike.in == '0,736111111'] <- '1*10-3'

## Also substituting 1:10000 with 1*10-4 to unify the format
inserts_with_spikeins$spike.in[inserts_with_spikeins$spike.in == '1:10000'] <- '1*10-4'

## Logging the types of spike ins identified in the data
logfile <- c(logfile,
             paste("SpikeIn types identified:",
                   paste(unique(inserts_with_spikeins$spike.in), collapse = ", ")))

##----------------------------
describe(inserts_with_spikeins$spike.in)
##----------------------------

## Filtering out everything except LAIR1+ contigs
inserts_with_spikeins_LAIR1only <- inserts_with_spikeins[str_detect(inserts_with_spikeins$Gene.id.s., 'LAIR1'),]

## Logging the number of observations with LAIR1
logfile <- c(logfile,
             paste(nrow(inserts_with_spikeins_LAIR1only), "observations contain LAIR1 as the insertion"))

## Filtering the contigs from run6 because this is the only run with MGO3 and MMJ5 cells spike-ins
## Grouping by sample.id and counting the reads_percentage
inserts_SI_LAIR1only <- inserts_with_spikeins_LAIR1only %>%
  filter(run.id == 'run6') %>%
  group_by(sample.id) %>%
  summarise(run.id = nth(run.id,1),
            donor.id = nth(donor.id,1),
            spike.in = nth(spike.in,1),
            reads = mean(mean.bp.coverage),
            total_reads = mean(PE.reads.number),
            isotype = nth(isotype,1),
            insert.type = nth(insert.type,1)) %>%
  mutate(reads_percentage = reads*100/total_reads)

## Logging the number of samples in the data
logfile <- c(logfile,
             paste("Run 6 is the run of interest as it contains the MGO3 and MMJ5 spike-ins\n",
                   nrow(inserts_SI_LAIR1only), " samples in run6 contain LAIR1 as the insertion",
                   sep = ''))


## OUTPUT
write.table(inserts_SI_LAIR1only,
            file = 'inserts_allRuns_spikein_LAIR1_only_22_07_19.txt',
            sep = '\t',
            row.names = FALSE)

g <- ggplot(inserts_SI_LAIR1only)+
  geom_col(aes(x = sample.id, y = reads_percentage), fill = 'gray', alpha = 1)+
  xlab('Sample')+
  ylab('Reads occupied by LAIR1 spike-in, %')+
  scale_x_discrete(labels = paste(inserts_SI_LAIR1only$sample.id, inserts_SI_LAIR1only$spike.in, sep = '\n'))+
  theme_classic()

ggsave(paste('spike_in_percentage_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       device = "pdf",
       plot = g)

time_elapsed <- Sys.time() - script_started

logfile <- c(logfile,
             paste("Time elapsed: ", round(time_elapsed,1), ' seconds\n',
             "Finish time: ", Sys.time(), sep = ''))

write(logfile, paste('spike_in_percentage_log_',format(Sys.time(), "%d%m%y_%H%M"),'.txt', sep = ''))
