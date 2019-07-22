setwd("C:/Users/User/Documents/Lab/LongBCRs/MethodPaper/RunTables")

all_inserts <- read.table(file = 'Output/inserts_allRuns_11_04_19.txt', sep = '\t', header = TRUE, stringsAsFactors = FALSE)


install.packages('xlsx')
library(xlsx)
install.packages('rJava')

all_insertsM <- read.xlsx(file = 'Output/inserts_allRuns_11_04_19_MathildesOriginal.xlsx',1)

library(dplyr)
all_insertsM %>% mutate_if(is.factor, as.character) -> all_insertsM


## Generating the table with issues
na_count <- c()
letters_in_cols <- c()
commas_in_cols <- c()
col_class <- c()
for(i in colnames(all_insertsM)){
  na_count <- c(na_count, sum(is.na(all_insertsM[,i])))
  letters_in_cols <- c(letters_in_cols,  sum(str_detect(all_insertsM[,i], "[b-mo-zB-MO-Z]"), na.rm = TRUE))
  commas_in_cols <- c(commas_in_cols,  sum(str_detect(all_insertsM[,i], ","), na.rm = TRUE))
  col_class <- c(col_class, class(all_insertsM[,i]))
}
letters_in_cols <- data.frame(column = colnames(all_insertsM),
                              class = col_class,
                              na_count = na_count,
                              letters_present = letters_in_cols,
                              commas_present = commas_in_cols)
write.table(letters_in_cols, file = 'issue_solving_220719.txt', sep = '\t', row.names = FALSE, quote = FALSE)


## Replacing months with digits in V.contig.coord
all_insertsM$V.contig.coord[str_detect(all_insertsM$V.contig.coord, "[b-mo-zB-MO-Z]")] <-
str_replace_all(all_insertsM$V.contig.coord[str_detect(all_insertsM$V.contig.coord, "[b-mo-zB-MO-Z]")],
                c('Jan'='1','Feb'='2','Mar'='3','Apr'='4','May'='5','Jun'='6','Jul'='7','Aug'='8',
                  'Sep'='9','Oct'='10','Nov'='11','Dec'='12'))

## The same for const.contig.coord
all_insertsM$const.contig.coord[str_detect(all_insertsM$const.contig.coord, "[b-mo-zB-MO-Z]")] <-
  str_replace_all(all_insertsM$const.contig.coord[str_detect(all_insertsM$const.contig.coord, "[b-mo-zB-MO-Z]")],
                  c('Jan'='1','Feb'='2','Mar'='3','Apr'='4','May'='5','Jun'='6','Jul'='7','Aug'='8',
                    'Sep'='9','Oct'='10','Nov'='11','Dec'='12'))

## UNIFYING THE FORMAT OF THE VARIABLES
all_insertsM$spike.in[all_insertsM$spike.in == "1:1000"] <- "1*10-3"
all_insertsM$spike.in[all_insertsM$spike.in == "1:10000"] <- "1*10-4"

## MAKING THE TABLE SAFE TO OPEN AND EDIT IN EXCEL
columns_to_cure <- c('insert.contig.coord.s.', 'insert.genomic.coord.s.',
                     'V.contig.coord', 'D.contig.coord', 'J.contig.coord',
                     'const.contig.coord')
for(i in columns_to_cure){
  all_insertsM[,i] <- str_replace_all(all_insertsM[,i], "-", "..")
}

rm(columns_to_cure)

## NOW THIS TABLE CAN BE OPENED IN EXCEL WITHOUT PROBLEMS
## POLISHING OF TEXT NAs
columns_to_make_real_nas <- c('insert.type', 'spike.in', 'UniProt.info',
                              'Dgene.id','Jgene.id','const.gene','is.exon.complete')
for(i in columns_to_make_real_nas){
  all_insertsM[,i][all_insertsM[,i] == "NA"] <- NA
}

rm(columns_to_make_real_nas)

## POLISHING LOGICAL VARIABLES
logical_columns_to_polish <- c('IsVandConstInFrame.withoutStop', 'is.exon.complete', 'has.exon.exon.junction')

for(i in logical_columns_to_polish){
  all_insertsM[,i][all_insertsM[,i] == "false" | all_insertsM[,i] == "FALSE"] <- FALSE
  all_insertsM[,i][all_insertsM[,i] == "true" | all_insertsM[,i] == "TRUE"] <- TRUE
}

rm(logical_columns_to_polish)


## JOINING THE DATA WITH THE LATEST ADVANCES
# Joining variable
all_insertsM_copy <- all_insertsM
all_insertsM %>% unite(col = 'insert.umi', c('run.id','sample.id','insert.id.s.'), sep = ';') -> all_insertsM
if(length(unique(all_insertsM$insert.umi)) < nrow(all_insertsM)){
  print("INSERT.UMI IS NOT UNIQUE")
}

all_inserts_copy <- all_inserts
all_inserts %>% unite(col = 'insert.umi', c('run.id','sample.id','insert.id.s.'), sep = ';') -> all_inserts
if(length(unique(all_inserts$insert.umi)) < nrow(all_inserts)){
  print("INSERT.UMI IS NOT UNIQUE")
}
# This should yield a 100% intersection
intersect_size <- length(intersect(unique(all_insertsM$insert.umi), unique(all_inserts$insert.umi))) /
  length(unique(all_inserts$insert.umi))
if(intersect_size != 1){
  print('THE TABLES PROVIDED HAVE DIFFERENT CONTIGS')
}
rm(intersect_size)

# Adding donor, race, natural variables to the cured table
all_inserts %>% select(insert.umi,donor,race,natural) -> all_inserts
all_insertsM %>% left_join(all_inserts, by = "insert.umi") -> all_insertsM 

# Ungluing the insert.umi column
all_insertsM %>% separate(col = 'insert.umi', c('run.id','sample.id','insert.id.s.'), sep = ';') -> all_insertsM


# Reordering the columns a bit
colnames(all_insertsM)

desired_cols_order <- c('run.id','sample.id','insert.id.s.','donor','race','natural','spike.in','cells.number',
                        'cell.type','population','isotype','matrix.type','PE.reads.number','mean.bp.coverage',
                        'contig.length','insert.length','insert.type','is.exon.complete','IsVandConstInFrame.withoutStop',
                        'has.exon.exon.junction','Gene.name.s.','Gene.id.s.','UniProt.info','insert.genomic.coord.s.',
                        'insert.contig.coord.s.','mini.start.insert.coord','max.end.insert.coord','Vgene.id','Dgene.id',
                        'Jgene.id','const.gene','V.contig.coord','D.contig.coord','J.contig.coord','const.contig.coord',
                        'contig.seq','contigWithoutInsert.seq','insert.seq.s.','prot.contig.seq.inFrame','donor.id')
length(desired_cols_order)

all_insertsM <- all_insertsM[,desired_cols_order]

rm(desired_cols_order)


## Reversing the contigs
## For all coordinates split the column: 1..100 -> 1  100 (two separate columns)
columns_with_coords_to_separate <- c('const.contig.coord','V.contig.coord','D.contig.coord','J.contig.coord')

for(i in columns_with_coords_to_separate){
  separate(all_insertsM, col = i,
           into = c(paste(i,'start',sep='.'), paste(i,'end',sep='.')),
           sep = '\\..',
           convert = TRUE) -> all_insertsM
}

for(i in columns_with_coords_to_separate){
  all_insertsM[,paste(i,"start","1",sep='.')] <- pmin(all_insertsM[,paste(i,'start',sep='.')],
                                                     all_insertsM[,paste(i,'end',sep='.')])
  all_insertsM[,paste(i,"end","1",sep='.')] <- pmax(all_insertsM[,paste(i,'start',sep='.')],
                                                     all_insertsM[,paste(i,'end',sep='.')])
  all_insertsM[,paste(i,'start',sep='.')] <- all_insertsM[,paste(i,"start","1",sep='.')]
  all_insertsM[,paste(i,'end',sep='.')] <- all_insertsM[,paste(i,"end","1",sep='.')]
  all_insertsM %>% select(-c(paste(i,"start","1",sep='.'), paste(i,"end","1",sep='.'))) -> all_insertsM
}

rm(columns_with_coords_to_separate)


library(Biostrings)
all_insertsM$contig.direction <- 'direct'

for(i in seq(nrow(all_insertsM))){
  if(all_insertsM[i,"V.contig.coord.start"] > all_insertsM[i,"const.contig.coord.end"]){
    all_insertsM$contig.direction[i] <- 'reverse'
    all_insertsM$contig.seq[i] <- as.character(reverseComplement(DNAString(all_insertsM$contig.seq[i])))
    all_insertsM$contigWithoutInsert.seq[i] <-
      as.character(reverseComplement(DNAString(all_insertsM$contigWithoutInsert.seq[i])))
    all_insertsM$insert.seq.s.[i] <-
      paste(
        as.character(
          reverseComplement(
            DNAStringSet(
              unlist(str_split(
                str_remove_all(all_insertsM$insert.seq.s.[i],' '),
                ','))))),
        collapse =',')
  }
}


## Reversing the coordinates

#library(rlang)
#for(i in columns_with_coords_to_invert){
#  all_insertsM[,i][all_insertsM$contig.direction == 'reverse'] <- 
#    all_insertsM %>% 
#    filter(contig.direction == 'reverse') %>%
#    transmute(!!i := contig.length - !!i) %>%
#    select(sym(i)) %>% unlist() %>% as.vector()
#}


all_insertsM$V.contig.coord.start[all_insertsM$contig.direction == 'reverse'] <- all_insertsM %>%
  filter(contig.direction == 'reverse') %>%
  transmute(V.contig.coord.start = contig.length - V.contig.coord.start) %>%
  select(V.contig.coord.start) %>% unlist() %>% as.vector()

all_insertsM$V.contig.coord.end[all_insertsM$contig.direction == 'reverse'] <- 
  all_insertsM %>% 
  filter(contig.direction == 'reverse') %>%
  transmute(V.contig.coord.end = contig.length - V.contig.coord.end) %>%
  select(V.contig.coord.end) %>% unlist() %>% as.vector()

all_insertsM$D.contig.coord.start[all_insertsM$contig.direction == 'reverse'] <- all_insertsM %>%
  filter(contig.direction == 'reverse') %>%
  transmute(D.contig.coord.start = contig.length - D.contig.coord.start) %>%
  select(D.contig.coord.start) %>% unlist() %>% as.vector()

all_insertsM$D.contig.coord.end[all_insertsM$contig.direction == 'reverse'] <- 
  all_insertsM %>% 
  filter(contig.direction == 'reverse') %>%
  transmute(D.contig.coord.end = contig.length - D.contig.coord.end) %>%
  select(D.contig.coord.end) %>% unlist() %>% as.vector()

all_insertsM$J.contig.coord.start[all_insertsM$contig.direction == 'reverse'] <- all_insertsM %>%
  filter(contig.direction == 'reverse') %>%
  transmute(J.contig.coord.start = contig.length - J.contig.coord.start) %>%
  select(J.contig.coord.start) %>% unlist() %>% as.vector()

all_insertsM$J.contig.coord.end[all_insertsM$contig.direction == 'reverse'] <- 
  all_insertsM %>% 
  filter(contig.direction == 'reverse') %>%
  transmute(J.contig.coord.end = contig.length - J.contig.coord.end) %>%
  select(J.contig.coord.end) %>% unlist() %>% as.vector()

all_insertsM$const.contig.coord.start[all_insertsM$contig.direction == 'reverse'] <- all_insertsM %>%
  filter(contig.direction == 'reverse') %>%
  transmute(const.contig.coord.start = contig.length - const.contig.coord.start) %>%
  select(const.contig.coord.start) %>% unlist() %>% as.vector()

all_insertsM$const.contig.coord.end[all_insertsM$contig.direction == 'reverse'] <- 
  all_insertsM %>% 
  filter(contig.direction == 'reverse') %>%
  transmute(const.contig.coord.end = contig.length - const.contig.coord.end) %>%
  select(const.contig.coord.end) %>% unlist() %>% as.vector()


## Aaand finally flip it back
columns_with_coords_to_invert <- c('const.contig.coord','V.contig.coord','D.contig.coord','J.contig.coord')


for(i in columns_with_coords_to_invert){
  all_insertsM[,paste(i,"start","1",sep='.')] <- pmin(all_insertsM[,paste(i,'start',sep='.')],
                                                      all_insertsM[,paste(i,'end',sep='.')])
  all_insertsM[,paste(i,"end","1",sep='.')] <- pmax(all_insertsM[,paste(i,'start',sep='.')],
                                                    all_insertsM[,paste(i,'end',sep='.')])
  all_insertsM[,paste(i,'start',sep='.')] <- all_insertsM[,paste(i,"start","1",sep='.')]
  all_insertsM[,paste(i,'end',sep='.')] <- all_insertsM[,paste(i,"end","1",sep='.')]
  all_insertsM %>% select(-c(paste(i,"start","1",sep='.'), paste(i,"end","1",sep='.'))) -> all_insertsM
}

rm(columns_with_coords_to_invert)
rm(all_insertsM_copy, all_inserts_copy)
rm(i, na_count, commas_in_cols, col_class)
rm(letters_in_cols)

write.table(all_insertsM, file = 'data_cured_230719_0048.txt', sep = '\t',
            row.names = FALSE, quote = FALSE)

