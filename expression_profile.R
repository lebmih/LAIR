### GENES MAPPER
## Read GTF file from Ensembl
all_genes <- read.table(file = 'GenomicScape/Homo_sapiens.GRCh38.76.gtf',
                          sep = '\t', header = FALSE, stringsAsFactors = FALSE)

colnames(all_genes) <- c('seqname','source','feature','start','end','score','strand','frame','attribute')
all_genes_norm_chr <- all_genes %>%
  filter(seqname %in%unique(all_genes$seqname)[c(1:24,length(unique(all_genes$seqname)))])
all_genes_norm_chr_genes <- all_genes_norm_chr %>% filter(feature == 'gene')

unique(all_genes_norm_chr_genes$source)

all_genes_norm_chr_genes %>% group_by(source) %>% summarise(entries = n()) %>% arrange(desc(entries))

all_genes_norm_chr_genes_only_prot_coding <- all_genes_norm_chr_genes %>% filter(source == 'protein_coding')

# Let's introduce abbreviations
all_genes_norm_chr -> agnc
all_genes_norm_chr_genes -> agnc_g
all_genes_norm_chr_genes_only_prot_coding -> agnc_cds

agnc <- bind_cols(agnc, str_match(agnc$attribute,
                                      'gene_id (.{15,25}); gene_name (.{1,15});')[,2:3] %>%
                      as.data.frame())
colnames(agnc) <- c('seqname','source','feature','start','end',
                      'score','strand','frame','attribute','gene_id',
                      'gene_name')
agnc <- select(agnc, -attribute)

agnc_cds <- bind_cols(agnc_cds, str_match(agnc_cds$attribute,
                                'gene_id (.{15,25}); gene_name (.{1,15});')[,2:3] %>%
                        as.data.frame())
colnames(agnc_cds) <- c('seqname','source','feature','start','end',
                        'score','strand','frame','attribute','gene_id',
                        'gene_name')
agnc_cds <- select(agnc_cds, -attribute)

agnc_g <- bind_cols(agnc_g, str_match(agnc_g$attribute,
                                      'gene_id (.{15,25}); gene_name (.{1,15});')[,2:3] %>%
                      as.data.frame())
colnames(agnc_g) <- c('seqname','source','feature','start','end',
                      'score','strand','frame','attribute','gene_id',
                      'gene_name')
agnc_g <- select(agnc_g, -attribute)

## Let's import the data
diff_express_genes <- read.table(file = 'GenomicScape/differentially_expressed_genes.txt',
                        sep = '\t', header = TRUE, stringsAsFactors = FALSE)

sum(diff_express_genes$Name %in% agnc_cds$gene_name)
# Only 8797 genes out of 9303 are in the table filtered to contain genes only (agnc_cds, 19911 entries)
sum(diff_express_genes$Name %in% agnc_g$gene_name)
# Only 8803 genes out of 9303 are in the table unfiltered for genes (agnc_g, 58603 entries)


## Visual examination of these genes
diff_express_genes %>% filter(!Name %in% agnc_g$gene_name) -> diff_express_genes_notintable

## Those names are written in such format: GENE1 /// GENE2
## Googling, what does it mean is unavailable now
## Let's check whether these genes are present in filtered table
separate_rows(diff_express_genes, Name, sep = "///") -> diff_express_genes_long
diff_express_genes_long$Name %>% str_remove_all(" ") -> diff_express_genes_long$Name

sum(diff_express_genes_long$Name %in% agnc_g$gene_name)
# The number is 9813, meaning that we still miss 399 genes
# Checking manually
sum(agnc_g$gene_name == 'FAM21B', na.rm = TRUE)
# It's not present indeed
# Mostly these genes are LOCnnnnnn
length(unique(diff_express_genes_long$Name))
# It's 9800 meaning we can merge some rows
cols <- colnames(diff_express_genes_long)[-(1:2)]

diff_express_genes_long %>% group_by(Name) %>%
  summarise_at(vars(cols), mean, na.rm = TRUE) ->
  diff_express_genes_long

sum(diff_express_genes_long$Name %in% agnc_g$gene_name)
# The number is 9500 out of 9800, meaning that we still miss 300 genes
# Manually checking them
diff_express_genes_long %>% filter(!Name %in% agnc_g$gene_name) -> diff_express_genes_notintable


# Let's try another table
all_genes_NCBI_notgtf <- read.table(file = 'GenomicScape/genes_NCBI_RefSeq',
                             sep = '\t', header = FALSE, stringsAsFactors = FALSE)

colnames(all_genes_NCBI_notgtf) <- c('bin','name','chrom','strand','txStart','txEnd',
                                     'cdsStart','cdsEnd','exonCount','exonStarts','exonEnds','score',
                                     'name2','cdsStartStat','cdsEndStat','exonFrames')

# ALIASES
all_genes_NCBI_notgtf -> agni

sum(diff_express_genes_long$Name %in% agni$name2)
## 9252 out of 9800! That's better, 94.4%
diff_express_genes_long %>% filter(!Name %in% agni$name2) -> diff_express_genes_notintable

# Let's check whether some of these missing genes are in agnc Ensembl table
sum(diff_express_genes_notintable$Name %in% agnc$gene_name)
# Yes, 332 are found in the ENSEMBL
# Thus, for 9252 genes we can add info from NCBI table, for 332 - from ENSEMBL
# That would be 9584 genes, 97.8%
agni %>% select(name2, chrom, cdsStart, cdsEnd) -> agni_cds

agni %>% group_by(name2) %>% summarise(number_unique_starts = length(unique(cdsStart)),
                                       number_unique_ends = length(unique(cdsEnd)),
                                       number_unique_tx_starts = length(unique(txStart)),
                                       number_unique_tx_ends = length(unique(txEnd)),
                                       number_unique_chr = length(unique(chrom))) ->
  agni_summary

reshape(agni_summary, varying = c('number_unique_starts',
                                  'number_unique_ends',
                                  'number_unique_tx_ends',
                                  'number_unique_tx_starts',
                                  'number_unique_chr'),
        idvar = name2, direction = 'long') -> agni_summary_long

agni_summary %>% melt() -> agni_summary_long


agni_cds %>% filter(!str_detect(chrom,'_')) -> agni_cds

agni_cds$chrom %>% unique() %>% length()

agni_cds %>% group_by(name2) %>% summarise(cdsStart = min(cdsStart),
                                       cdsEnd = max(cdsEnd),
                                       chrom = nth(chrom, 1)) ->
  agni_unique

colnames(agni_unique) <- c('Name','start','end','chrom')

diff_express_genes_long %>% left_join(agni_unique, by = 'Name') -> diff_express_genes_long

## Now with the ENSEMBLE data
sum(is.na(diff_express_genes_long$chrom))

diff_express_genes_long %>% filter(is.na(chrom)) -> degl_no_coords


agnc_g %>% group_by(gene_name) %>% summarise(start = nth(start,1),
                                                       end = nth(end,1),
                                                       chrom = nth(seqname, 1)) ->
  agnc_unique

agnc_unique$chr <- 'chr'

agnc_unique %>% unite(chrom, c('chr','chrom'), sep = '') -> agnc_unique
agnc_unique$chrom[agnc_unique$chrom == 'chrMT'] <- 'chrM'
colnames(agnc_unique) <- c('Name', 'start','end','chrom')

degl_no_coords %>% select(-c('start','end','chrom')) %>%
  left_join(agnc_unique, by = 'Name') -> degl_no_coords

sum(is.na(degl_no_coords$chrom))

diff_express_genes_long %>% filter(!is.na(chrom)) -> diff_express_genes_long
degl_no_coords %>% filter(!is.na(chrom)) -> degl_no_coords
diff_express_genes_long %>% bind_rows(degl_no_coords) -> diff_express_genes_long

## COOL, that's the table
