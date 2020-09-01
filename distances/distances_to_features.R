# Distances to features script.

# Working folder ----------------------------------------------------------
setwd("C:/Users/User/Documents/Lab/LongBCRs/bioinformatics/lair")



# Libraries ---------------------------------------------------------------
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(seqinr)
library(lebutils)
library(stringi)
library(DescTools)


# Functions ---------------------------------------------------------------
in.ranges <- function(coordinate, starts, ends){
  starts.n <- length(starts)
  ends.n <- length(ends)
  if (starts.n != ends.n) {
    stop("Starts and Ends vectors not of the same length!")
  }
  if(any(starts > ends)) {
    stop("There are inverted ranges!")
  }
  count <- 0
  for(i in seq(starts.n)){
    if(coordinate >= starts[i] & coordinate <= ends[i]){
      count <- count + 1
    }
  }
  return(count)
}

# Data import -------------------------------------------------------------
contigs <- read.table("contigs_DATE.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
assembly.report <- read.table(file = "input/GRCh38_latest_assembly_report.txt",
                              stringsAsFactors = FALSE, header = FALSE, fill = TRUE, sep = "\t")

colnames(assembly.report) <- c("seq.name","seq.role","assigned.molecule",
                               "assigned.molecule.location.type","genbank.accn",
                               "relationship","refseq.accn","assembly.unit",
                               "seq.length","UCSC.style.name")
assembly.report <- assembly.report %>% dplyr::filter(seq.role == "assembled-molecule")
assembly.report <- assembly.report %>%
  mutate(chrom = UCSC.style.name, chrom.length = seq.length) %>%
  select(chrom, chrom.length)

centromeres <- read.table(file = "input/hg38_centromeres.txt",
                          stringsAsFactors = FALSE, header = TRUE, fill = TRUE, sep = "\t")
centromeres <- centromeres %>%
  dplyr::select(-c("bin", "name")) %>%
  group_by(chrom) %>%
  summarise(centromer.start = min(chromStart),
            centromer.end = max(chromEnd)) %>%
  mutate(centromer.center = (centromer.start+centromer.end)/2)

# Importing the CFS-associated genes data ---------------------------------
# The same procedure goes for the CFS-associated genes.
files.path <- "~/Lab/LongBCRs/MethodPaper/CommonFragileSites/fragile_site_gene_bed/"
files <- list.files(path = files.path, pattern = ".bed", full.names = TRUE)
cfs.genes <- data.frame()

for (i in files) {
  cfs.temp <- read.table(file = i,
                         sep = '\t', header = FALSE, stringsAsFactors = FALSE,
                         skip = 3)
  if (is.na(cfs.temp$V2) %>% sum() == length(cfs.temp$V2)){
    cfs.temp %>% select(-V2) -> cfs.temp
  }
  colnames(cfs.temp) <- c("chrom","start","end","name","score","strand")
  cfs.genes <- bind_rows(cfs.genes, cfs.temp)
}

rm(cfs.temp, i, files, files.path)

cfs.genes$chrom[cfs.genes$chrom == 'chrx'] <- "chrX"
cfs.genes$chrom[cfs.genes$chrom == 'chry'] <- "chrY"

cfs.genes$start[cfs.genes$name == "NRIp1"] <- 14961235
cfs.genes$end[cfs.genes$name == "NRIp1"] <- 15065936

cfs.genes <- cfs.genes %>% filter(!is.na(start))
cfs.genes <- cfs.genes %>% group_by(name) %>%
  summarise(chrom = nth(chrom, 1),
            start = min(start),
            end = max(end),
            score = nth(score, 1),
            strand = nth(strand, 1))
cfs.genes <- cfs.genes %>% ungroup()
cfs.genes$chrom[cfs.genes$chrom == "chromosome2"] <- "chr2"

# Characterization of CFS themselves --------------------------------------
cfs.genes <- cfs.genes %>% mutate(length = end - start) %>%
  left_join(assembly.report, by = "chrom") %>%
  left_join(centromeres, by = "chrom") %>%
  mutate(centr.dist = pmin(abs(start - centromer.center),
                           abs(centromer.center - end)),
         tel.dist = pmin(start, chrom.length - end))

median(cfs.genes$tel.dist)
# Telomeric distance of CFS - 25.1 Mb
# Centromeric distance of CFS - 30.99 Mb

# Overlap with CFS --------------------------------------------------------
contigs.copy <- contigs %>%
  filter(!is.na(min.g.start) & chrom %in% unique(cfs.genes$chrom)) %>%
  mutate(g.center = (min.g.start + max.g.end)/2)
cfs.genes <- cfs.genes %>% mutate(center = (start + end)/2)

# Calculating the distances
dist.to.closest.cfs <- c()
overlap.closest.cfs <- c()

for (i in seq(nrow(contigs.copy))) {
  cfs.this.chrom <- cfs.genes %>% filter(chrom == contigs.copy$chrom[i])
  cfs.dist <- min(abs(cfs.this.chrom$center - contigs.copy$g.center[i]))
  cfs.overlap <- count.overlap(cfs.this.chrom$start,
                                cfs.this.chrom$end,
                                contigs.copy$min.g.start[i],
                                contigs.copy$max.g.end[i])
  dist.to.closest.cfs <- c(dist.to.closest.cfs, cfs.dist)
  overlap.closest.cfs <- c(overlap.closest.cfs, cfs.overlap)
  if (i %in% seq(0,nrow(contigs.copy),10000)) {
    paste0(round(i*100/nrow(contigs.copy)), "% of the job is done") %>% print()
  }
}
# Control: this number should be 0
sum(is.infinite(dist.to.closest.cfs))

contigs.copy$cfs.dist <- dist.to.closest.cfs
contigs.copy$cfs.overlap <- overlap.closest.cfs

contigs.copy <- contigs.copy %>% select(-g.center)
contigs.nocfs <- contigs %>% filter(is.na(min.g.start) | !chrom %in% unique(cfs.genes$chrom))
contigs <- contigs.copy %>% bind_rows(contigs.nocfs)

rm(dist.to.closest.cfs, overlap.closest.cfs, cfs.this.chrom, i, contigs.copy,
   contigs.nocfs, cfs.dist, cfs.overlap)


# R-loops -----------------------------------------------------------------
rloops <- read.table("rloops/K562_DRIP_peaks_hg38.bed", sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(rloops) <- c("chrom", "start", "end")
rloops <- rloops %>% filter(!str_detect(chrom, "_"))

# Characterization of R-loops themselves ----------------------------------
rloops <- rloops %>% mutate(length = end - start) %>%
  left_join(assembly.report, by = "chrom") %>%
  left_join(centromeres, by = "chrom") %>%
  mutate(centr.dist = pmin(abs(start - centromer.center),
                           abs(centromer.center - end)),
         tel.dist = pmin(start, chrom.length - end))

sum(rloops$tel.dist < 0)
rloops %>% filter(tel.dist < 0) %>% View()
rloops$tel.dist[rloops$chrom == "chrM"] <- NA
median(rloops$centr.dist, na.rm = TRUE)
# Telomeric distance of Rloops - 24.2 Mb
# Centromeric distance of Rloops - 36.6 Mb

# Overlap with R-loops ----------------------------------------------------
contigs.copy <- contigs %>%
  filter(!is.na(min.g.start) & chrom %in% unique(rloops$chrom)) %>%
  mutate(g.center = (min.g.start + max.g.end)/2)
rloops <- rloops %>% mutate(center = (start + end)/2)

# Calculating the distances
dist.to.closest.rloops <- c()
overlap.closest.rloops <- c()

for (i in seq(nrow(contigs.copy))) {
  rloops.this.chrom <- rloops %>% filter(chrom == contigs.copy$chrom[i])
  rloops.dist <- min(abs(rloops.this.chrom$center - contigs.copy$g.center[i]))
  rloops.overlap <- count.overlap(rloops.this.chrom$start,
                               rloops.this.chrom$end,
                               contigs.copy$min.g.start[i],
                               contigs.copy$max.g.end[i])
  dist.to.closest.rloops <- c(dist.to.closest.rloops, rloops.dist)
  overlap.closest.rloops <- c(overlap.closest.rloops, rloops.overlap)
  if (i %in% seq(0,nrow(contigs.copy),10000)) {
    paste0(round(i*100/nrow(contigs.copy)), "% of the job is done") %>% print()
  }
}
# Control: this number should be 0
sum(is.infinite(dist.to.closest.rloops))

contigs.copy$rloops.dist <- dist.to.closest.rloops
contigs.copy$rloops.overlap <- overlap.closest.rloops

contigs.copy <- contigs.copy %>% select(-g.center)
contigs.norloops <- contigs %>% filter(is.na(min.g.start) | !chrom %in% unique(rloops$chrom))
contigs <- contigs.copy %>% bind_rows(contigs.norloops)

rm(dist.to.closest.rloops, overlap.closest.rloops, rloops.this.chrom, i, contigs.copy,
   contigs.norloops, rloops.dist, rloops.overlap)

# Extracting sequences from generated inserts -----------------------------
contigs.gen.bed <- contigs %>% filter(ins.type %in% c("VDJ.gen", "J-CH1.gen", "V-CH1.gen")) %>%
  select(chrom, min.g.start, max.g.end) %>%
  mutate(min.g.start = min.g.start - 1)

contigs.gen.bed.chrm <- contigs.gen.bed %>% filter(chrom == "chrM" & max.g.end > 16569) %>%
  mutate(., id = seq(nrow(.)),
         min.g.start1 = 0,
         max.g.end1 = max.g.end - 16569,
         max.g.end = 16569) %>%
  unite("left", c("min.g.start", "max.g.end")) %>%
  unite("right", c("min.g.start1", "max.g.end1")) %>%
  pivot_longer(cols = c("left", "right"),
               names_to = "type",
               values_to = "coords") %>%
  separate("coords", c("start", "end")) %>%
  select(chrom, start, end, id)

contigs.gen.bed <- contigs.gen.bed %>% filter(!(chrom == "chrM" & max.g.end > 16569))

write.table(contigs.gen.bed, file = "contigs_gen_280420.bed", row.names = FALSE, quote = FALSE,
            col.names = FALSE, sep = "\t")
write.table(contigs.gen.bed.chrm, file = "contigs_gen_mito_280420.bed", row.names = FALSE, quote = FALSE,
            col.names = FALSE, sep = "\t")
# Move the files to the folder with hg38.fa and perform sequence extraction by
# $  bedtools getfasta -fi hg38.fa -tab -bed filename.bed > filename.txt
# For _mito_ add -name argument to the command.
# Then import the .txt files back.
contigs.gen.seq <- read.table("contigs_gen_280420.txt", sep = "\t", stringsAsFactors = FALSE)
contigs.gen.mito.seq <- read.table("contigs_gen_mito_280420.txt", sep = "\t", stringsAsFactors = FALSE)
if(nrow(contigs.gen.seq) == nrow(contigs.gen.bed) &
   nrow(contigs.gen.mito.seq) == nrow(contigs.gen.bed.chrm)){
  print("Alles gut.")
}
rm(contigs.gen.bed, contigs.gen.bed.chrm)

colnames(contigs.gen.seq) <- c("coords", "insert.seq.s")
colnames(contigs.gen.mito.seq) <- c("coords", "insert.seq.s")

contigs.gen.mito.seq <- contigs.gen.mito.seq %>%
  separate(coords, c("id"), sep = ":", remove = TRUE, convert = TRUE,
           extra = "drop") %>% group_by(id) %>%
  summarise(insert.seq.s = paste0(insert.seq.s, collapse = ""))


contigs.gen <- contigs %>% filter(ins.type %in% c("VDJ.gen", "J-CH1.gen", "V-CH1.gen"))
contigs.gen.mito <- contigs.gen %>% filter(chrom == "chrM" & max.g.end > 16569)
contigs.gen <- contigs.gen %>% filter(!(chrom == "chrM" & max.g.end > 16569))
contigs.gen$insert.seq.s. <- str_to_upper(contigs.gen.seq$insert.seq.s)
contigs.gen.mito$insert.seq.s. <- str_to_upper(contigs.gen.mito.seq$insert.seq.s)
rm(contigs.gen.mito.seq, contigs.gen.seq)

contigs.gen <- bind_rows(contigs.gen, contigs.gen.mito)

contigs <- contigs %>% filter(!ins.type %in% c("VDJ.gen", "J-CH1.gen", "V-CH1.gen")) %>%
  bind_rows(contigs.gen)

rm(contigs.gen, contigs.gen.mito)

# Measuring GC content ----------------------------------------------------
gccontent <- c()
for (i in seq(1:nrow(contigs))) {
  gccontent <- c(gccontent, GC(s2c(contigs$insert.seq.s.[i])))
}
contigs$gc.content <- gccontent
rm(i, gccontent)

# Expression --------------------------------------------------------------
rnaseq <- lookfortable(pattern = "rnaseq_\\S*.txt")
# Appending the generated table with expression
contigs.gen <- contigs %>% filter(ins.type %in% c("VDJ.gen", "J-CH1.gen", "V-CH1.gen"))
contigs.gen$contig.id <- seq(max(contigs$contig.id)+1, max(contigs$contig.id) + nrow(contigs.gen))

shared.cols <- intersect(colnames(contigs.gen), colnames(rnaseq))
shared.cols <- shared.cols[-1]

contigs.gen[,shared.cols] <- NA

conv.trx <- c()
conv.trx.genes <- c()
conv.trx.thr <- c()
conv.trx.trx <- c()

for(i in seq(nrow(contigs.gen))){
  if(i %in% seq(0, nrow(contigs.gen), 5000)){
    print(i/nrow(contigs.gen))
  }
  target.chrom <- contigs.gen$chrom[i]
  target.start <- contigs.gen$min.g.start[i]
  target.end <- contigs.gen$max.g.end[i]
  genes.this.chrom <- rnaseq %>% filter(chrom == target.chrom)
  olap.genes <- genes.this.chrom %>% filter((gene.start <= target.start & gene.end >= target.start) |
                                              (gene.start <= target.end & gene.end >= target.end) |
                                              (gene.start >= target.start & gene.end <= target.end))
  if(nrow(olap.genes) == 1){
    contigs.gen[i, shared.cols] <- olap.genes[1, shared.cols]
  } else if(nrow(olap.genes) > 1){
    #whch <- sample(seq(nrow(olap.genes)), 1)
    olap.genes <- olap.genes %>% arrange(trx.above.thr)
    contigs.gen[i, shared.cols] <- olap.genes[1, shared.cols]
    genes <- paste0(olap.genes$gene.name, collapse = ";")
    genes.above.thr <- paste0(olap.genes$trx.above.thr, collapse = ";")
    olap.genes <- olap.genes %>% pivot_longer(cols = shared.cols[-1],
                                              names_to = "test.pop",
                                              values_to = "trx") %>%
      group_by(gene.name) %>% summarise(trx = round(median(trx), 3),
                                        trx.above.thr = nth(trx.above.thr, 1)) %>%
      arrange(desc(trx.above.thr))
    genes.trx <- paste0(olap.genes$trx, collapse = ";")

    conv.trx <- c(conv.trx, contigs.gen$contig.id[i])
    conv.trx.genes <- c(conv.trx.genes, genes)
    conv.trx.thr <- c(conv.trx.thr, genes.above.thr)
    conv.trx.trx <- c(conv.trx.trx, genes.trx)
    rm(genes, genes.above.thr, genes.trx)
  } else {
    next
  }
}


# Doing the same with a copy of real data to compare
contigs.copy <- contigs %>% filter(!ins.type %in% c("VDJ.gen", "J-CH1.gen", "V-CH1.gen"))

shared.cols <- intersect(colnames(contigs.copy), colnames(rnaseq))
shared.cols <- shared.cols[-1]

contigs.copy[,shared.cols] <- NA

conv.trx <- c()
conv.trx.genes <- c()
conv.trx.thr <- c()
conv.trx.trx <- c()

for(i in seq(nrow(contigs.copy))){
  if(i %in% seq(0, nrow(contigs.copy), 100)){
    print(i/nrow(contigs.copy))
  }
  target.chrom <- contigs.copy$chrom[i]
  target.start <- contigs.copy$min.g.start[i]
  target.end <- contigs.copy$max.g.end[i]
  genes.this.chrom <- rnaseq %>% filter(chrom == target.chrom)
  olap.genes <- genes.this.chrom %>% filter((gene.start <= target.start & gene.end >= target.start) |
                                              (gene.start <= target.end & gene.end >= target.end) |
                                              (gene.start >= target.start & gene.end <= target.end))
  if(nrow(olap.genes) == 1){
    contigs.copy[i, shared.cols] <- olap.genes[1, shared.cols]
  } else if(nrow(olap.genes) > 1){
    olap.genes <- olap.genes %>% arrange(desc(trx.above.thr))
    contigs.copy[i, shared.cols] <- olap.genes[1, shared.cols]
    genes <- paste0(olap.genes$gene.name, collapse = ";")
    genes.above.thr <- paste0(olap.genes$trx.above.thr, collapse = ";")
    olap.genes <- olap.genes %>% pivot_longer(cols = shared.cols[-1],
                                              names_to = "test.pop",
                                              values_to = "trx") %>%
      group_by(gene.name) %>% summarise(trx = round(median(trx), 3),
                                        trx.above.thr = nth(trx.above.thr, 1)) %>%
      arrange(desc(trx.above.thr))
    genes.trx <- paste0(olap.genes$trx, collapse = ";")

    conv.trx <- c(conv.trx, contigs.copy$contig.id[i])
    conv.trx.genes <- c(conv.trx.genes, genes)
    conv.trx.thr <- c(conv.trx.thr, genes.above.thr)
    conv.trx.trx <- c(conv.trx.trx, genes.trx)
    rm(genes, genes.above.thr, genes.trx)
  } else {
    next
  }
}

# Checking whether this method yielded different result in comparison to the
# joining the tables by gene names.
contigs.nat <- contigs %>% filter(!ins.type %in% c("VDJ.gen", "J-CH1.gen", "V-CH1.gen"))

for(i in shared.cols){
  setdiff(contigs.copy[,i], contigs.nat[,i]) %>% length() %>% print()
}
# It yields diff result for ~65 observations. Let's take a look at them.

diff.id <- c()
for(i in seq(nrow(contigs.nat))){
  for(j in shared.cols){
    if(is.na(contigs.nat[i,j]) & is.na(contigs.copy[i,j])){
      break
    } else if(is.na(contigs.nat[i,j]) | is.na(contigs.copy[i,j])){
      diff.id <- c(diff.id, contigs.nat$contig.id[i])
    } else {
      if(contigs.nat[i,j] != contigs.copy[i,j]){
        diff.id <- c(diff.id, contigs.nat$contig.id[i])
        break
      }
    }
  }
}

contigs.copy %>% filter(contig.id %in% diff.id) %>% View()

for(j in shared.cols){
  print(contigs.nat[i, j])
}
rnaseq %>% filter(gene.name == "DEFB124") %>% View()
# I checked just a couple manually, the initial way was more robust

contigs <- contigs.nat %>% bind_rows(contigs.gen)

contigs.nat$thecol <- "before"
contigs.copy$thecol <- "after"
contigs.agg <- contigs.nat %>% bind_rows(contigs.copy)

ggplot(contigs.agg %>% filter(!chrom == "chrM") %>% filter(trx.above.thr > 3))+
  geom_violin(aes(x = ins.type, y = log(G4YW_B_naive + 1), fill = thecol, color = population.simple))+
  theme_classic()
# As I see, it also did not change the distribution drastically

pbbc.cols <- shared.cols[str_detect(shared.cols, "naive|SM")]

contigs$contig.id %>% unique() %>% length()
contigs.agg <- contigs %>% select(c("contig.id", pbbc.cols)) %>%
  pivot_longer(cols = pbbc.cols, names_to = "sample",
               values_to = "TPM") %>%
  group_by(contig.id) %>%
  summarise(pbbc.trx = median(TPM, na.rm = TRUE))
contigs <- contigs %>% left_join(contigs.agg, by = "contig.id")
rm(contigs.agg)

convergent <- data.frame(contig.id = conv.trx,
                         conv.genes = conv.trx.genes,
                         conv.trx.above.thr = conv.trx.thr,
                         conv.tpm = conv.trx.trx)
convergent <- convergent %>% filter(!str_detect(conv.genes, "MT-"))

contigs.gen <- contigs.gen %>% select(-c("conv.genes", "conv.trx.above.thr",
                                         "conv.tpm")) %>%
  left_join(convergent, by = "contig.id")

rm(olap.genes, genes.this.chrom, i, j, pbbc.cols, target.chrom, target.end, target.start, diff.id,
   conv.trx, conv.trx.genes, conv.trx.thr, conv.trx.trx, convergent.trx, convergent)

setdiff(contigs$gene.name.single, rnaseq$gene.name)
contigs$gene.name.single %>% unique() %>% length()

unique.genes <- contigs$gene.name.single[!is.na(contigs$gene.name.single)] %>% unique()
sum(unique.genes %in% rnaseq$gene.name)

# Manually fixing incomplete overlap of the gene names --------------------
# gnir - genes not in rnaseq
gnir <- unique.genes[!unique.genes %in% rnaseq$gene.name]
# Looking for the genes in the BioMart
library(biomaRt)
library(EnsDb.Hsapiens.v86)
library(EnsDb.Hsapiens.v79)
mart <- biomaRt::useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
gnir.ids <- getBM(filters = "clone_based_ensembl_gene",
                  attributes = c("clone_based_ensembl_gene", "ensembl_gene_id"),
                  values = gnir, mart = mart)

gnir <- as.data.frame(gnir, stringsAsFactors = FALSE)
gnir <- gnir %>%
  left_join(gnir.ids, by = c("gnir" = "clone_based_ensembl_gene"))
rm(gnir.ids)
# Filling in the unidentified genes manually
gnir$ensembl_gene_id[gnir$gnir == "LINC01811"] <- "ENSG00000226320"
gnir$ensembl_gene_id[gnir$gnir == "LINC01915"] <- "ENSG00000265485"
gnir$ensembl_gene_id[gnir$gnir == "LINC01699"] <- "ENSG00000179452"
gnir$ensembl_gene_id[gnir$gnir == "LINC02196"] <- "ENSG00000250974"
gnir$ensembl_gene_id[gnir$gnir == "LINC02218"] <- "ENSG00000249662"
gnir$ensembl_gene_id[gnir$gnir == "LINC02240"] <- "ENSG00000260192"
gnir$ensembl_gene_id[gnir$gnir == "TEDC1"] <- "ENSG00000185347"
gnir$ensembl_gene_id[gnir$gnir == "LHFPL6"] <- "ENSG00000183722"
gnir$ensembl_gene_id[gnir$gnir == "AFG1L"] <- "ENSG00000135537"
gnir$ensembl_gene_id[gnir$gnir == "AL139099.4"] <- "ENSG00000283029"
gnir$ensembl_gene_id[gnir$gnir == "ELP1"] <- "ENSG00000070061"
gnir$ensembl_gene_id[gnir$gnir == "SLC44A3-AS1"] <- "ENSG00000224081"
gnir$ensembl_gene_id[gnir$gnir == "TGIF2-RAB5IF"] <- "ENSG00000259399"
gnir$ensembl_gene_id[gnir$gnir == "CFAP410"] <- "ENSG00000160226"
gnir$ensembl_gene_id[gnir$gnir == "TENT5D"] <- "ENSG00000174016"
gnir$ensembl_gene_id[gnir$gnir == "RF00100"] <- "ENSG00000202198"
gnir$ensembl_gene_id[gnir$gnir == "AP001092.1"] <- "ENSG00000237410"
# Kicking out ungooglable genes
gnir <- gnir %>% dplyr::filter(!is.na(ensembl_gene_id))
# Now looking for HGNC symbol from EnsDb
geneIDs <- ensembldb::select(EnsDb.Hsapiens.v86, keys = gnir$ensembl_gene_id,
                             keytype = "GENEID",
                             columns = c("SYMBOL","GENEID", "SEQNAME", "SEQSTRAND", "GENESEQSTART", "GENESEQEND"))

gnir <- gnir %>% left_join(geneIDs, by = c("ensembl_gene_id" = "GENEID"))

rnaseq <- rnaseq %>% mutate(gene.length = gene.end - gene.start)


# Appending gene length ---------------------------------------------------
contigs.gnir <- contigs %>% dplyr::filter(gene.name.single %in% gnir$gnir)
contigs.notgnir <- contigs %>% dplyr::filter(!gene.name.single %in% gnir$gnir)
rnaseq.agg <- rnaseq %>% group_by(gene.name) %>%
  summarise(gene.length = nth(gene.length, 1))

contigs.notgnir <- contigs.notgnir %>%
  left_join(dplyr::select(rnaseq.agg, c("gene.name", "gene.length")),
            by = c("gene.name.single" = "gene.name"))

gnir <- gnir %>% mutate(gene.length = GENESEQEND - GENESEQSTART)

contigs.gnir <- contigs.gnir %>%
  left_join(dplyr::select(gnir, c("gnir", "gene.length")),
            by = c("gene.name.single" = "gnir"))
rm(gnir, mart, geneIDs, rnaseq.agg)

contigs <- contigs.notgnir %>% bind_rows(contigs.gnir)
rm(contigs.notgnir, contigs.gnir, contigs.gen, contigs.nat, contigs.copy)
rm(chrom, unique.genes)

write.table(contigs, file = "contigs_280420.txt", quote = FALSE, sep = "\t",
            row.names = FALSE)

# Plotting ----------------------------------------------------------------
contigs <- contigs %>%
  mutate(ins.type = factor(ins.type, levels = c("VDJ.gen", "VDJ",
                                                              "J-CH1.gen", "J-CH1",
                                                              "V-CH1.gen", "V-CH1")),
         population.simple = factor(population.simple, levels = c("naive", "memory_or_plasma",
                                                                  "activated", "ebv", "unsorted",
                                                                  "generated")),
         ins.class = factor(ins.class, levels = c("intergenic", "intron", "partial_exon", "exon")))

rloopsdist <- ggplot(data = contigs %>%
                       dplyr::filter(!ins.type %in% c("V-CH1.gen","V-CH1"),
                                     chrom != "chrM"))+
  geom_boxplot(aes(ins.type, log10(rloops.dist+1), fill = multiple.insertions),
               outlier.shape = 21, outlier.fill = "grey")+
  theme_classic()+
  xlab("Insert type")+
  ylab("Distance of insertion donor\nto the closest DRIP-seq peal, log10(bp+1)")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(0.5,7))
rloopsdist

ggsave("rloops_DATE_mult.svg", plot = rloopsdist, device = "svg", width = 8, height = 6)

pops <- c("naive", "memory_or_plasma", "activated", "ebv", "unsorted")
types <- c("VDJ.gen", "VDJ", "J-CH1.gen", "J-CH1")
classes <- c("intergenic", "intron", "partial_exon", "exon")
mult <- c(F,T)

for(i in types){
  for(j in mult){
    median(contigs$rloops.dist[contigs$ins.type == i & contigs$chrom != "chrM" &
                                 contigs$multiple.insertions == j], na.rm = TRUE) %>% print()
  }

}


for(i in types){
  for(j in mult){
    (sum(contigs$rloops.overlap[contigs$ins.type == i & contigs$chrom != "chrM" & contigs$multiple.insertions == j] > 0, na.rm = TRUE)/
      sum(contigs$ins.type == i & contigs$chrom != "chrM" & contigs$multiple.insertions == j, na.rm = TRUE)) %>% print()
  }

}

wilcox.test(contigs$rloops.dist[contigs$ins.type == "VDJ.gen" & contigs$chrom != "chrM"],
            contigs$rloops.dist[contigs$ins.type == "VDJ" & contigs$chrom != "chrM"])
wilcox.test(contigs$rloops.dist[contigs$ins.type == "J-CH1.gen" & contigs$chrom != "chrM"],
            contigs$rloops.dist[contigs$ins.type == "J-CH1" & contigs$chrom != "chrM"])
wilcox.test(contigs$rloops.dist[contigs$ins.type == "VDJ.gen"],
            contigs$rloops.dist[contigs$ins.type == "VDJ"])

sum(contigs$rloops.overlap[contigs$ins.type == "VDJ.gen"] > 0, na.rm = TRUE)
sum(contigs$rloops.overlap[contigs$ins.type == "VDJ"] > 0, na.rm = TRUE)

sum(contigs$ins.type == "VDJ.gen")
sum(contigs$ins.type == "VDJ")


median(contigs$gene.length[contigs$ins.type == "J-CH1"], na.rm = TRUE)

median(contigs$gc.content[contigs$ins.type == "J-CH1.gen" & contigs$chrom != "chrM"], na.rm = TRUE)


contigs.temp <- contigs %>% dplyr::filter(!ins.type %in% c("V-CH1.gen","V-CH1") &
                                            ins.class != "intergenic" &
                                            trx.above.thr > 2 & chrom != "chrM")

ggplot(data = contigs.temp,
       aes(x = ins.type, y = log2(pbbc.trx+1))) +
  geom_violin(position = position_dodge())+
  geom_boxplot(width=.6, outlier.shape=NA)+
  scale_y_continuous(breaks = c(0,2.0,4.0,6.0,8.0,10.0,12.0,14.0), limits = c(0,15))+
  theme_classic()+
  theme(legend.position = "none")

ggsave("trx_290420_base.svg", device = "svg", width = 8, height = 6)



ggplot(data = rnaseq.agg %>% dplyr::filter(trx.above.thr > 2),
       aes(x = 1, y = log2(pbbc.trx+1)))+
  geom_violin(position = position_dodge())+
  geom_boxplot(width=.6, outlier.shape=NA)+
  scale_y_continuous(breaks = c(0,2.0,4.0,6.0,8.0,10.0,12.0,14.0), limits = c(0,15))+
  theme_classic()+
  theme(legend.position = "none")

rnaseq$chrom %>% unique()
median(rnaseq$G4YW_B_naive[rnaseq$trx.above.thr > 2])

rnaseq <- rnaseq %>%
  left_join(assembly.report, by = "chrom") %>%
  left_join(centromeres, by = "chrom") %>%
  mutate(centr.dist = pmin(abs(gene.start - centromer.center),
                           abs(centromer.center - gene.end)),
         tel.dist = pmin(gene.start, chrom.length - gene.end))


ggplot(data = rnaseq %>% dplyr::filter(trx.above.thr > 3))+
  geom_point(aes(x = log10(gene.length+1), y = log2(G4YW_B_NSM+1)),
             fill = "grey", color = "black", shape = 21)+
  theme_classic()

cor(rnaseq$gene.length, rnaseq$G4YW_B_naive)




contigs %>% dplyr::filter(ins.type %in% c("J-CH1.gen") &
                            ins.class == "exon") %>% nrow()




rnaseq.agg <- rnaseq %>%
  pivot_longer(cols = pbbc.cols, names_to = "sample",
               values_to = "TPM") %>%
  group_by(gene.name) %>%
  summarise(pbbc.trx = median(TPM, na.rm = TRUE),
            chrom = nth(chrom,1),
            trx.above.thr = nth(trx.above.thr,1))

hg38.genes <- read.table("hg38_genes.bed", sep = "\t", stringsAsFactors = FALSE)
colnames(hg38.genes) <- c("chrom", "start", "end", "gene.name",
                          "score", "strand", "start1", "start2",
                          "score1", "exon.number", "exon.length", "exon.start")
hg38.genes <- hg38.genes %>% filter(!str_detect(chrom, "_"))

rnaseq$in.hg38.genes <- NA
for(i in seq(nrow(rnaseq))){
  target.chrom <- rnaseq$chrom[i]
  target.start <- rnaseq$gene.start[i]
  target.end <- rnaseq$gene.end[i]
  hg38.genes.this.chrom <- hg38.genes %>% filter(chrom == target.chrom)
  hg38.genes.this.regions <- hg38.genes.this.chrom %>%
    filter((start <= target.start & end >= target.start) |
             (start <= target.end & end >= target.end) |
             (start >= target.start & end <= target.end))
  if(nrow(hg38.genes.this.regions) > 0){
    rnaseq$in.hg38.genes[i] <- nrow(hg38.genes.this.regions)
  }
  rm(hg38.genes.this.regions)
}
sum(is.na(rnaseq$in.hg38.genes))


contigs.gen.conv <- contigs %>% filter(!is.na(conv.genes))


ggplot(rnaseq)+
  geom_point(aes(x = log2(in.hg38.genes), y = log2(G4YW_B_naive+1)),
             shape = 21, fill = "grey", color = "black")+
  theme_classic()

hg38.exons <- read.table("hg38_exons.bed", sep = "\t", stringsAsFactors = FALSE)
colnames(hg38.exons) <- c("chrom", "start", "end", "exon.name", "score", "strand")
hg38.exons <- hg38.exons %>% mutate(length = end - start + 1)
hg38.exons <- hg38.exons %>% filter(!str_detect(chrom, "_"))

hg38.exons$tpm <- NA

for(i in seq(nrow(hg38.exons))){
  tg.chrom <- hg38.exons$chrom[i]
  tg.start <- hg38.exons$start[i]
  tg.end <- hg38.exons$end[i]
  rnaseq.this.chrom <- rnaseq %>% filter(chrom == tg.chrom)
  rnaseq.this.regions <- rnaseq.this.chrom %>%
    filter((gene.start <= tg.start & gene.end >= tg.start) |
             (gene.start <= tg.end & gene.end >= tg.end) |
             (gene.start >= tg.start & gene.end <= tg.end))
  if(nrow(rnaseq.this.regions) > 0){
    hg38.exons$tpm[i] <- median(rnaseq.this.regions$G4YW_B_naive, na.rm = TRUE)
  }
  
}

ggplot(hg38.exons)+
  geom_point(aes(x = log2(length), y = log2(tpm+1)),
             shape = 21, fill = "grey", color = "black")+
  theme_classic()

contigs.new <- contigs %>% filter(!ins.type %in% c("VDJ.gen", "V-CH1.gen", "J-CH1.gen")) %>%
  bind_rows(contigs.gen)

ggplot(contigs.new %>% filter(!ins.type %in% c("V-CH1", "V-CH1.gen"),
                          chrom != "chrM"))+
  geom_boxplot(aes(x = ins.type, y = log10(gene.length)))+
  theme_classic()+
  scale_y_continuous(limits = c(2.5,6.5))

contigs.nat <- contigs %>% filter(ins.type %in% c("VDJ", "V-CH1", "J-CH1"))

write.table(contigs.nat, file = "contigs_nat_010520.txt", sep = "\t", row.names = FALSE, quote = FALSE)



contigs.nat %>% filter(chrom == "chrM" & ((min.g.start >= 14747 & min.g.start <= 15887) |
                         (max.g.end >= 14747 & max.g.end <= 15887)))  %>% nrow()


contigs.nat %>% filter(chrom == "chrM" & ((min.g.start >= 10760 & min.g.start <= 12137) |
                                            (max.g.end >= 10760 & max.g.end <= 12137)))  %>% nrow()

wilcox.test(contigs$cfs.dist[contigs$ins.type == "VDJ.gen"],
            contigs$cfs.dist[contigs$ins.type == "VDJ"])
wilcox.test(contigs$cfs.dist[contigs$ins.type == "J-CH1.gen"],
            contigs$cfs.dist[contigs$ins.type == "J-CH1"])
median(contigs$cfs.dist[contigs$ins.type == "VDJ.gen"], na.rm = TRUE)


# 08.05.20

rloopsdist <- ggplot(data = contigs %>%
                       dplyr::filter(!ins.type %in% c("V-CH1.gen","V-CH1"),
                                     chrom != "chrM",
                                     population.simple %in% c("naive", "memory_or_plasma", "activated")))+
  geom_boxplot(aes(ins.type, log10(aid.offt.olap.dist+1), fill = population.simple),
               outlier.shape = 21, outlier.fill = "grey")+
  theme_classic()+
  xlab("Insert type")+
  ylab("Distance of insertion donor\nto the closest DRIP-seq peal, log10(bp+1)")+
  theme(axis.text.x = element_text(size = 16, angle = 0, vjust = 0.6))+
  theme(legend.position = "none")+
  scale_y_continuous(limits = c(1,8), breaks = seq(1,8,1))
rloopsdist

ggsave("aid_080520_pops1.svg", device = "svg", width = 7, height = 5)


wilcox.test(contigs$aid.offt.olap.dist[contigs$ins.type == "J-CH1" 
                                       &
                                         contigs$chrom != "chrM" & contigs$population.simple == "naive"],
            contigs$aid.offt.olap.dist[contigs$ins.type == "J-CH1" &
                                         contigs$chrom != "chrM" & contigs$population.simple == "memory_or_plasma"])


wilcox.test(contigs$aid.offt.olap.dist[contigs$ins.type == "J-CH1.gen" & contigs$chrom != "chrM"],
            contigs$aid.offt.olap.dist[contigs$ins.type == "J-CH1" & contigs$chrom != "chrM"])
wilcox.test(contigs$aid.offt.olap.dist[contigs$ins.type == "VDJ.gen" & contigs$chrom != "chrM"],
            contigs$aid.offt.olap.dist[contigs$ins.type == "VDJ" & contigs$chrom != "chrM"])
