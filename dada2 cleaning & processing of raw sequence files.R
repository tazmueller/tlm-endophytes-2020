
# dada2 processing and taxonomy assignment

#In terminal: 

#load modules here
module load cutadapt

# remove adapters
mkdir 01_adapter
for i in *_R1_001.fastq .gz; do echo "~/.local/bin/cutadapt --cores 8 --minimum-length 150 -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT -o 01_adapter/${i//_R1_001.fastq.gz/-cut_R1_001.fastq.gz} -p 01_adapter/${i//_R1_001.fastq.gz/-cut_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" >> run_cutadapt.sh; done
chmod +x run_cutadapt.sh
./run_cutadapt.sh
cd 01_adapter    
mkdir 01_logs
mv *log.txt 01_logs

# remove primers
mkdir ../02_filtered  
for i in *_R1_001.fastq.gz; do echo "~/.local/bin/cutadapt --cores 8 --minimum-length 100 --discard-untrimmed -g GTGCCAGCMGCCGCGGTAA -G GGACTACHVGGGTWTCTAAT --discard-untrimmed -o ../02_filtered/${i//-cut_R1_001.fastq.gz/-trimmed_R1_001.fastq.gz} -p ../02_filtered/${i//-cut_R1_001.fastq.gz/-trimmed_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../02_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.sh; done
chmod + x run_cutadapt2.sh
./run_cutadapt2.sh
cd ../02_filtered/
  mkdir 02_logs
mv *log.txt 02_logs


# DADA2 in R
#R

library(dada2)
library(ggplot2)
library(patchwork)
library(DECIPHER)
setwd("/fastq_files/")
path <- (".")
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# quality plots
pdf("quality_forwards.pdf")
plotQualityProfile(fnFs)
dev.off()
pdf("quality_reverse.pdf")
plotQualityProfile(fnRs)
dev.off()

# V4
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(200,160), maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=4)
head(out)

# error models
errF <- learnErrors(filtFs, multithread=4)
errR <- learnErrors(filtRs, multithread=4)
# plots
pdf("error_model_forwards.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()
pdf("error_model_reverse.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

#dereplicate reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

dadaFs <- dada(derep_forward, err=errF, multithread=4, pool="pseudo")
dadaRs <- dada(derep_reverse, err=errR, multithread=4, pool="pseudo")

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, trimOverhang=TRUE, minOverlap=50)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, "seqtab.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, "seqtab_nochim.rds")

# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                          filtered=out[,2], dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/out[,1]*100, 1))
write.table(summary_tab, file = "sequence_process_summary.txt", sep = "\t", quote=FALSE)

seqtab.nochim <- readRDS("seqtab_nochim.rds")
dna <- DNAStringSet(getSequences(seqtab.nochim)) # Create a DNAStringSet from the ASVs
load("/home/umii/goul0109/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET
ids <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE) # use all processors
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species") # ranks of interest
# Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
taxid <- t(sapply(ids, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))
colnames(taxid) <- ranks; rownames(taxid) <- getSequences(seqtab.nochim)
saveRDS(taxid, file = "taxID.rds")
quit("no") 

############
#ITS
############

#load modules here
module load cutadapt

# remove adapters
mkdir 01_adapter
for i in *_R1_001.fastq.gz; do echo "~/.local/bin/cutadapt --cores 8 --minimum-length 150 -a TATGGTAATTGTGTGNCAGCNGCCGCGGTAA -g ATTAGANACCCNNGTAGTCCGGCTGGCTGACT -A AGTCAGCCAGCCGGACTACNVGGGTNTCTAAT -o 01_adapter/${i//_R1_001.fastq.gz/-cut_R1_001.fastq.gz} -p 01_adapter/${i//_R1_001.fastq.gz/-cut_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > 01_adapter/cutadapt.${i//_R1_001.fastq.gz/.log.txt}" >> run_cutadapt.sh; done
chmod +x run_cutadapt.sh
./run_cutadapt.sh
cd 01_adapter    
mkdir 01_logs
mv *log.txt 01_logs

R

library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")

path <- "."  ## CHANGE ME to the directory containing the fastq files.
list.files(path)
fnFs <- sort(list.files(path, pattern = "R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "R2_001.fastq.gz", full.names = TRUE))
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GCTGCGTTCTTCATCGATGC"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients




list.files(path, pattern = "*.gz")
fnFs.filtN <- file.path(path, "../02_N", basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, "../02_N", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
quit("no")

mkdir ../03_filtered  
for i in *_R1_001.fastq.gz; do echo "~/.local/bin/cutadapt --cores 8 --minimum-length 75 --discard-untrimmed -g CTTGGTCATTTAGAGGAAGTAA -a GCATCGATGAAGAACGCAGC -G GCTGCGTTCTTCATCGATGC -A TTACTTCCTCTAAATGACCAAG -n 2 -o ../03_filtered/${i//-cut_R1_001.fastq.gz/-trimmed_R1_001.fastq.gz} -p ../03_filtered/${i//-cut_R1_001.fastq.gz/-trimmed_R2_001.fastq.gz} ${i} ${i//_R1_/_R2_} > ../03_filtered/cutadapt.${i//_R1_001.fastq.gz/.adapter.log.txt}" >> run_cutadapt2.sh; done
chmod +x run_cutadapt2.sh
./run_cutadapt2.sh
cd ../03_filtered
mkdir 03_logs
mv *log.txt 03_logs

R

#####
#####
#####
library(ggplot2)
packageVersion("ggplot2")
library(dada2)
packageVersion("dada2")
library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
library(patchwork)
packageVersion("patchwork")
library(DECIPHER)
packageVersion("DECIPHER")

path <- (".")
list.files(path)
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


## checking primers ##
FWD <- "CTTGGTCATTTAGAGGAAGTAA"  ## CHANGE ME to your forward primer sequence
REV <- "GCTGCGTTCTTCATCGATGC"  ## CHANGE ME...

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

fnFs.filtN <- file.path(path, basename(fnFs)) # Put N-filterd files in filtN/ subdirectory
fnRs.filtN <- file.path(path, basename(fnRs))
primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))
#####################

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "04_dada2_filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "04_dada2_filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# quality plots
pdf("quality_forwards.pdf")
plotQualityProfile(fnFs)
dev.off()
pdf("quality_reverse.pdf")
plotQualityProfile(fnRs)
dev.off()

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN = 0, maxEE = c(4, 4), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = TRUE)


# error models
errF <- learnErrors(filtFs, multithread=8)
errR <- learnErrors(filtRs, multithread=8)
# plots
pdf("error_model_forwards.pdf")
plotErrors(errF, nominalQ=TRUE)
dev.off()
pdf("error_model_reverse.pdf")
plotErrors(errR, nominalQ=TRUE)
dev.off()

#dereplicate reads
derep_forward <- derepFastq(filtFs, verbose=TRUE)
names(derep_forward) <- sample.names # the sample names in these objects are initially the file names of the samples, this sets them to the sample names for the rest of the workflow
derep_reverse <- derepFastq(filtRs, verbose=TRUE)
names(derep_reverse) <- sample.names

dadaFs <- dada(derep_forward, err=errF, multithread=8, pool="pseudo")
dadaRs <- dada(derep_reverse, err=errR, multithread=8, pool="pseudo")

merged_amplicons <- mergePairs(dadaFs, derep_forward, dadaRs, derep_reverse, trimOverhang=TRUE, minOverlap=50)

seqtab <- makeSequenceTable(merged_amplicons)
dim(seqtab)
saveRDS(seqtab, "../05_seqtab.rds")

seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, "../seqtab_nochim.rds")

# set a little function
getN <- function(x) sum(getUniques(x))

# making a little table
summary_tab <- data.frame(row.names=sample.names, dada2_input=out[,1],
                          filtered=out[,2], dada_f=sapply(dadaFs, getN),
                          dada_r=sapply(dadaRs, getN), merged=sapply(merged_amplicons, getN),
                          nonchim=rowSums(seqtab.nochim),
                          final_perc_reads_retained=round(rowSums(seqtab.nochim)/out[,1]*100, 1))
write.table(summary_tab, file = "../sequence_process_summary_04_05_06.txt", sep = "\t", quote=FALSE)

# taxonomy ID
unite.ref <- "/home/umii/public/gopher-data/ITS_UNITE/sh_general_release_dynamic_02.02.2019_dev.fasta"  # CHANGE ME to location on your machine
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
saveRDS(taxa, file = "unite_all_taxa.rds")
#####