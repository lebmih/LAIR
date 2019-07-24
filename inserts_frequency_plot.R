### Inserts Frequency plotting for In- and Out of Frame
setwd("C:/Users/User/Documents/Lab/LongBCRs/MethodPaper/RunTables")
## Libraries
x <- c('tidyr', 'dplyr', 'ggplot2', 'ggrepel', 'reshape2', 'stringr', 'Biostrings')
lapply(x, require, character.only = TRUE)
rm(x)

all_inserts <- read.table(file = 'data_cured_230719_0048.txt',
                          sep = '\t', header = TRUE, stringsAsFactors = FALSE)

## Determining the inserts frequency for natural donors (excluded nucleofection, activation, etc.)
all_inserts %>% 
  filter(natural == 'true') %>%
  group_by(donor) %>%
  summarise(contig_count = n(),
            cell_count_available = sum(!is.na(cells.number)) > 0,
            cell_count = nth(cells.number, 1),
            in_frame_ins_count = sum(IsVandConstInFrame.withoutStop == TRUE),
            out_frame_ins_count = sum(IsVandConstInFrame.withoutStop == FALSE)) %>%
  mutate(in_frame_freq = in_frame_ins_count/cell_count,
         out_frame_freq = out_frame_ins_count/cell_count) -> in_out_inserts_freq

## Transposition to the longer format
in_out_inserts_freq %>% gather(key = is_in_frame,
                               value = ins_freq,
                               in_frame_freq:out_frame_freq) -> in_out_inserts_freq_long

sum(in_out_inserts_freq$cell_count_available == TRUE)

plotting_upper_cutoff <- 4e-04

outliersInFrame <- sum(in_out_inserts_freq$in_frame_freq > plotting_upper_cutoff, na.rm = TRUE)
outliersOutFrame <- sum(in_out_inserts_freq$out_frame_freq > plotting_upper_cutoff, na.rm = TRUE)

g <- ggplot(in_out_inserts_freq_long)+
  geom_boxplot(aes(x = is_in_frame, y = ins_freq), outlier.shape = NA)+
  geom_dotplot(aes(x = is_in_frame, y = ins_freq), 
               binaxis = "y", stackdir = "center", method = "dotdensity", binwidth = 1e-6,
               alpha = 0.7, dotsize = 7)+
  xlab('Insert type, n = 37')+
  ylab('Insert Frequency')+
  scale_x_discrete(labels = c('In frame', 'Out of frame'))+
  scale_y_continuous(limits = (c(0, plotting_upper_cutoff)))+
  theme_classic()+
  theme(text = element_text(family = 'sans'))

if(outliersInFrame > 0){
  g <- g +
    geom_label(aes(x = 1, y = plotting_upper_cutoff-0.1e-04,
                   label = paste(outliersInFrame,'outliers\nnot shown')),
               family = 'sans', fontface = 'plain', size = 3.8, label.size = 0)
}

if(outliersOutFrame > 0){
  g <- g +
    geom_label(aes(x = 2, y = plotting_upper_cutoff-0.1e-04,
                   label = paste(outliersOutFrame,'outliers\nnot shown')),
               family = 'sans', fontface = 'plain', size = 3.8, label.size = 0)
}

ggsave(paste('insert_frequency_plot_',
             format(Sys.time(), "%d%m%y_%H%M"),
             '.pdf', sep = ''),
       width = 3, height = 4,
       device = "pdf",
       plot = g)

write.table(in_out_inserts_freq, file = 'in_out_inserts_freq.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)
