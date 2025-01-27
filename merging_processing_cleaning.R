###########################################
# Loading packages
###########################################

library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
install.packages("tidyverse")
library(tidyverse)
install.packages("vegan")
library(vegan)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")


###########################################
# Reading & converting data
###########################################

#input sequences & remove characters from rownames
inputtable.full <- readRDS('seqtab_nochim.rds')
row.names(inputtable.full) <- sub("-.*", "", row.names(inputtable.full))

#input taxonomy (from dada2 rds)
taxonomy <- as.data.frame(readRDS('moeller_unite_all_taxa.rds'))
write.csv(taxonomy, "euk_taxonomy.csv")

#input metadata
metadata <- read.csv('endo_serv_metadata.csv', head = TRUE)

#assign rownames
row.names(metadata) <- metadata$sampleID


#convert host & sites to factors
metadata$host <- as.factor(metadata$host)
metadata$host <- factor(metadata$host, levels = c("xantiana", "unguiculata", "marah", "none"))
metadata$siteType <- factor(metadata$siteType, levels = c("xan_only", "ung_incl", "control"))

#check rownames
row.names(inputtable.full)


###########################################
# Making phyloseq object
###########################################
#convert otu table to matrix
full.otu <- data.matrix(inputtable.full)

#check class type
class(taxonomy)
class(full.otu)

#Inverting dataset (if you get "NULL" run this again line by line)

full.otu <- t(full.otu) #BE CAREFUL NOT TO RUN THIS MULTIPLE TIMES OR IT WILL INVERT THE DATAFRAME

#converting taxa table to matrix
taxonomy <- as.matrix(taxonomy)
taxonomy

#creating otu table object
otu.table <- otu_table(as.matrix(full.otu), taxa_are_rows = TRUE)

#creating taxonomy table object
taxa <- tax_table(as.matrix(taxonomy))

#creating sample data object
sdata <- sample_data(metadata)

#ensure rownames between all are the same
all(rownames(full.otu) == rownames(taxonomy))

#convert to phyloseq
ps.full <- phyloseq(otu.table, taxa, sdata)

#check sample names
sample_names(otu.table)
sample_names(taxa)

#check variables
get_variable(ps.full)
ps.full

#check OTU tables
otu_table(ps.full)

# how many reads per sample?
sort(sample_sums(ps.full))


###########################################
# Cleaning & processing
###########################################

# remove any taxa with less than 10 total reads across all samples
ps.full.t <- prune_taxa(taxa_sums(ps.full) >= 10, ps.full)
#3015 taxa, 124 samples (including controls)


# remove taxa with NA for Kingdom (Domain), or Kingdom other than Fungi
ps.full.id <- subset_taxa(ps.full.t, Kingdom == "k__Fungi")
ps.full.id

check <- tax_table(ps.full.id)
check <- as.data.frame(check)
##################################### read count: see how the cleaning affected taxa 
nrow(tax_table(ps.full)) ########   3450 (Original dataset)
nrow(tax_table(ps.full.t)) ######   3015 (Removing taxa with less than 10 total reads)
nrow(tax_table(ps.full.id)) #####   2533 (Anything non-fungal removed Kingdom Fungi)
sort(sample_sums(ps.full.t))

# how did cleaning affect read numbers?
sort(sample_sums(ps.full.id))




###########################################
### Removing rinse controls 
###########################################

#Get CTRL sample names:
#sample_names(ps.full.id)
#List of CTRL samples:
#CTR39
#CTR93
#CTR94
#CTR95
#CTR96 [this sample is the one that is pooled rinsewater, so will be the one subtracted]




#subset negative controls

contam.96 <- subset_samples(ps.full.id, sampleID == "CTR96")


# prune taxa not in controls
contam.96 <- prune_taxa(taxa_sums(contam.96) > 0, contam.96)

contam.96 #84 taxa


# get list of sequences in NegCtrls

contam.96.df <- as.data.frame(as(otu_table(contam.96), "matrix"))


# get list of taxa in NegCtrls
contam.96.df <- as.data.frame(as(tax_table(contam.96), "matrix"))

# how abundant were these contaminants for each control?
mean(taxa_sums(contam.96))
taxa_sums(contam.96)

#create full dataframe
samples.full.df <- as.data.frame(as(otu_table(ps.full.id), "matrix"))


#create working subset dataframe
samples.subtract.full.df <- samples.full.df


#unspreading to work with rownames code
samples.full.df.s <- t(samples.full.df)
samples.subtract.full.df.s <- as.data.frame(samples.full.df.s)
class(samples.subtract.full.df.s)


# get row numbers for NegCtrls
NegCtrl.row.96 <- which(rownames(samples.subtract.full.df.s) == "CTR96")

samples.subtract.full.df <- as.data.frame(apply(samples.subtract.full.df.s, 2, function(x) x-(x[NegCtrl.row.96])))



# make sure the subtraction looks right
samples.full.df.s[, 50:55]
samples.subtract.full.df[, 50:55]

write.csv (samples.subtract.full.df, "subtract.csv")
# replace negative abundance values with 0
samples.subtract.full.df <- as.data.frame(apply(samples.subtract.full.df, 2, function(x) ifelse(x < 0, 0, x)))


# check to make sure it looks right again
samples.subtract.full.df[, 15:25]


# make new phyloseq object with new OTU table
ps.full.adj <- ps.full.id
otu_table(ps.full.adj) <- otu_table(as.matrix(samples.subtract.full.df), taxa_are_rows = F)
sort(sample_sums(ps.full.adj))


#subset out taxa that are not 0
ps.full.adj <- prune_taxa(taxa_sums(ps.full.adj) > 0, taz.full.adj)
ps.full.adj

# final total = 2527 taxa



###########################################
###  Final clean & assemble
###########################################

#subset out all controls from sample data
ps.full.nc <- subset_samples(ps.full.adj, sampleID != "CTR39" & sampleID != "CTR93" & sampleID != "CTR94" & sampleID != "CTR95" & sampleID != "CTR96" )


#check sample read #s
sort(sample_sums(ps.full.id))
sort(sample_sums(ps.full.adj))
sort(sample_sums(ps.full1))

# remove any taxa with no reads
ps.full.f <- prune_taxa(taxa_sums(ps.full.nc) > 0, ps.full.nc) 
ps.full.f
# final total = 2524 taxa

# get mean read abundance per sample
mean(sample_sums(ps.full.f)) # 23835
median(sample_sums(ps.full.f)) # 23705
sort(sample_sums(ps.full.f))

samplesums <- as.data.frame(sort(sample_sums(ps.full.f)))
write.csv(samplesums, "post-clean sample read counts.csv")

# use short names for ASVs
dna.its.full <- Biostrings::DNAStringSet(taxa_names(ps.full.f))
names(dna.its.full) <- taxa_names(ps.full.f)
ps.full.f <- merge_phyloseq(ps.full.f, dna.its.full)
taxa_names(ps.full.f) <- paste0("ASV", seq(ntaxa(ps.full.f)))
ps.full.f
tax_table(ps.full.f)
taxa_table <- tax_table(ps.full.f)

#write sample data to CSV
ps.full.sampledata <- as.data.frame(as(sample_data(ps.full.f), "matrix"))
write.csv(ps.full.sampledata, "FinalSampleData.csv")

#save cleaned RDS
saveRDS(ps.full.f, "ps.full.f.rds")
