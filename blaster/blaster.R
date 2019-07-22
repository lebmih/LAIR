if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install('Biostrings')
install.packages("devtools")


user_renviron = path.expand(file.path("~", ".Renviron"))
file.edit(user_renviron) # open with another text editor if this fails

## Type this in the file:
#LC_COLLATE  = "English_United States.1252"
#LC_CTYPE    = "English_United States.1252"
#LC_MONETARY = "English_United States.1252"
#LC_NUMERIC  = "English_United States.1252"
#LC_TIME     = "English_United States.1252"

Sys.getlocale()

library(devtools)
install_github("mhahsler/rBLAST")
setwd("C:/Users/User/Documents/Lab/LongBCRs/Blaster/")
install.packages('R.utils')
library(R.utils)
gunzip('hg38.fa.gz')
library(Biostrings)
library(dplyr)
library(tidyr)

## load some test data 
seq <- readDNAStringSet('TOPO10072019.txt')
seq

library(rBLAST)

do.overlap <- function(first.range.start, first.range.end, second.range.start, second.range.end){
  do.overlap <- TRUE
  first.range.start <- as.numeric(first.range.start)
  first.range.end <- as.numeric(first.range.end)
  second.range.start <- as.numeric(second.range.start)
  second.range.end <- as.numeric(second.range.end)
  if(first.range.start > first.range.end | second.range.start > second.range.end ){
    stop("One of the ranges you provided is flipped")
  }
  if(first.range.end < second.range.start | first.range.start > second.range.end){
    do.overlap <- FALSE
  }
  return(do.overlap)
}


bl <- blast(db="C:/Users/User/Documents/Lab/LongBCRs/Blaster/hg38")


## Creating the BLAST databases from Vector Sequences
makeblastdb(file = 'vectorFlanks.fa', dbtype = 'nucl',
            args = '-parse_seqids -title vectorFlanks -out vectorFlanks')
makeblastdb(file = 'vectorWhole.fa', dbtype = 'nucl',
            args = '-parse_seqids -title vectorWhole -out vectorWhole')
vectorFlanks <- blast(db = "C:/Users/User/Documents/Lab/LongBCRs/Blaster/vectorFlanks")
vectorWhole <- blast(db = "C:/Users/User/Documents/Lab/LongBCRs/Blaster/vectorWhole")

## Creating the BLAST databases from Primer Sequences
makeblastdb(file = 'FredAlt_Primers.txt', dbtype = 'nucl',
            args = '-parse_seqids -title faaprimers -out faaprimers')
primers <- blast(db = "C:/Users/User/Documents/Lab/LongBCRs/Blaster/faaprimers")

## Creating the BLAST database from the expected downstream sequence
makeblastdb(file = 'QSOX2_Downstream.fa', dbtype = 'nucl',
            args = '-parse_seqids -title downstream -out downstream')
downstream <- blast(db = "C:/Users/User/Documents/Lab/LongBCRs/Blaster/downstream")


readLength <- c()
vectorLeftStart <- c()
vectorLeftEnd <- c()
vectorRightStart <- c()
vectorRightEnd <- c()
insertedSeqLength <- c()
primersFound <- c()
downstreamFound <- c()
generalFound <- c()

for (i in seq(length(seq))){
  ## Outputing the sequencing read length
  readLength <- c(readLength, length(seq[[i]]))
  ## BLAST against the Vector Flanks
  vectEnds <- predict(vectorFlanks, seq[i], BLAST_args = "-perc_identity 95")
  
  
  ## Finding vector ends can be difficult in case of short (~20 nt) ends
  ## If this is the case, force BLAST to be more sensitive
  ## and then filter the longest alignment for FWD and REV ends:
  if (nrow(vectEnds) < 2){
    vectEnds <- predict(vectorFlanks, seq[i], BLAST_args = "-word_size 7 -perc_identity 95")
    vectEnds1 <- filter(vectEnds, SubjectID == unique(vectEnds$SubjectID)[1]) %>%
      arrange(desc(Alignment.Length))
    vectEnds2 <- filter(vectEnds, SubjectID == unique(vectEnds$SubjectID)[2]) %>%
      arrange(desc(Alignment.Length))
    vectEnds <- bind_rows(vectEnds1[1,], vectEnds2[1,])
    rm(vectEnds1, vectEnds2)
  }
  
  
  ## Check the alignment quality, proceed if only there is at least one aligned piece with less than 3 MM
  if (all(vectEnds$Mismatches < 2) |
      do.overlap(vectEnds$Q.start[1], vectEnds$Q.end[1],
                 vectEnds$Q.start[2], vectEnds$Q.end[2])){
    ## Extract four coordinates and select the 2nd and the 3rd in order - these are the boundaries of the
    ## inserted product. Everything outside these boundaries is considered vector-derived automatically.
    vectCoords <- c(vectEnds$Q.end,
                    vectEnds$Q.start)
    vectCoords <- sort(vectCoords)
    
    ## Assigning the values for the future data frame
    vectorLeftStart <- c(vectorLeftStart, vectCoords[1])
    vectorLeftEnd <- c(vectorLeftEnd, vectCoords[2])
    vectorRightStart <- c(vectorRightStart, vectCoords[3])
    vectorRightEnd <- c(vectorRightEnd, vectCoords[4])
    
    ## Crop the sequence using the calculated boundaries
    insertedSeq <- seq[[i]][vectCoords[2]:vectCoords[3]]
    insertedSeqLength <- c(insertedSeqLength, length(insertedSeq))
    ## Try to find the primers in the sequence
    insertedSeq <- DNAStringSet(insertedSeq)  ## 'predict' function requires DNAStringSet format
    primerBlast <- predict(primers, insertedSeq, BLAST_args = "-word_size 10 -perc_identity 95")
    if (nrow(primerBlast) == 0){
      primersFound <- c(primersFound, "NA")
    } else {
      ## Append the vector of primers with a semicolon separated list of primers aligned
      primerBlast <- unite(primerBlast, 'Query', c('Q.start', 'Q.end'), sep = "-")
      primerBlast <- unite(primerBlast, 'PrimerCoord', c('SubjectID', 'Query'), sep = ":")
      primersFound <- c(primersFound, paste(primerBlast$PrimerCoord, collapse = ";"))
    }
    
    downstreamBlast <- predict(downstream, insertedSeq, BLAST_args = "-perc_identity 95")
    if (nrow(downstreamBlast) == 0){
      downstreamFound <- c(downstreamFound, "NA")
    } else {
      downstreamBlast <- unite(downstreamBlast, 'Query', c('Q.start', 'Q.end'), sep = "-")
      downstreamFound <- c(downstreamFound, paste(downstreamBlast$Query, collapse = ';'))
    }
    
    generalBlast <- predict(bl, insertedSeq, BLAST_args = "-word_size 15 -perc_identity 95")
    if (nrow(generalBlast) == 0){
      generalFound <- c(generalFound, "NA")
    } else {
      ## Append the vector of primers with a semicolon separated list of primers aligned
      generalBlast <- unite(generalBlast, 'Query', c('Q.start', 'Q.end'), sep = "-")
      generalBlast <- unite(generalBlast, 'ChrCoord', c('SubjectID', 'Query'), sep = ":")
      generalFound <- c(generalFound, paste(generalBlast$ChrCoord, collapse = ";"))
    }
    
  } else {
    vectorLeftStart <- c(vectorLeftStart, "NA")
    vectorLeftEnd <- c(vectorLeftEnd, "NA")
    vectorRightStart <- c(vectorRightStart, "NA")
    vectorRightEnd <- c(vectorRightEnd, "NA")
    insertedSeqLength <- c(insertedSeqLength, "NA")
    primersFound <- c(primersFound, "NA")
    downstreamFound <- c(downstreamFound, "NA")
    generalFound <- c(generalFound, "NA")
  }
  
  
}

topoAnnot2 <- data.frame(seqnames = names(seq),
                        readLength,
                        vectorLeftStart,
                        vectorLeftEnd,
                        vectorRightStart,
                        vectorRightEnd,
                        insertedSeqLength,
                        primersFound,
                        downstreamFound,
                        generalFound)

write.table(topoAnnot, file = 'topo100719annot.txt', sep = "\t", row.names = F, col.names = T)
