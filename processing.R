###########################################
# Loading packages
###########################################

library(ShortRead)
packageVersion("ShortRead")
library(Biostrings)
packageVersion("Biostrings")
install.packages("tidyverse")
library(tidyverse)
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
metadata <- read.csv('endo_surv_metadata.csv', head = TRUE)

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
ps.full.adj <- prune_taxa(taxa_sums(ps.full.adj) > 0, ps.full.adj)
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
sort(sample_sums(ps.full))

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



###########################################
# Pruning
###########################################
library("vegan")
library(phyloseq)

# remove samples with < 7470 reads
ps.fungi <- prune_samples(sample_sums(ps.fungi) > 7470, ps.fungi) # remove sample with < 7470 reads

#113 samples remaining

ps.fungi <- prune_taxa(taxa_sums(ps.fungi) > 0, ps.fungi) # remove any taxa with no reads after previous pruning

ps.fungi
#2473 taxa

sample_names(ps.fungi)
sort(sample_sums(ps.fungi))
#Samples 18M, OSRM, MOM, CHFM, GCX, and LOX were removed from the dataset

mean(sample_sums(ps.fungi)) # 24859
median(sample_sums(ps.fungi)) # 24173

#save PRUNED RDS
#saveRDS(ps.fungi, "ps.fungi.pruned.rds")



###########################################
# Rarefying
###########################################

#setwd ("~/C1 Laptop Local/C1 LAPTOP WORKING") #change to your directory
#ps.fungi <- readRDS("ps.fungi.pruned.rds") 

#organize otu table to be readable by vegan package
sample_data(ps.fungi)
r.tpf<- otu_table(ps.fungi)
r.tpft <- t(r.tpf)
class(r.tpft) <- "matrix"
tax_table(ps.fungi)

#look at rarefaction curve of original data

rare.curve <- rarecurve(r.tpft, step=50, cex=0.5, ylab="ASVs", label=T,
                        main="Rarefaction Curve for all samples")


#rarefy- seed = set.seed(123)
ps.fungi.r = rarefy_even_depth(ps.fungi, rngseed=123, sample.size=min(sample_sums(ps.fungi)), replace=F, trimOTUs = T, verbose = T)

#2462 taxa

#save rarefied dataset
save(ps.fungi.r, file="ps.fungi.rare123.RData")
saveRDS(ps.fungi.r, file="ps.fungi.rare123.rds")




###########################################
# Filling in NA taxa
###########################################


require("tidyverse")
require("phyloseq")

name_na_taxa <- function(ps_obj, include_rank = T, na_label = "Unknown <tax> (<rank>)"){
  
  # Check arguments
  if(!grepl("<tax>", na_label)){
    stop("Error: include '<tax>' in the na_label")
  }
  if (include_rank){
    if(!grepl("<rank>", na_label)){
      stop("Error: include_rank = TRUE; include '<rank>' in the na_label")
    }
  } else {
    if(grepl("<rank>", na_label)){
      stop("Error: include_rank = FALSE; remove '<rank>' from the na_label")
    }
  }
  
  # Convert to long data
  taxa_long <- tax_table(ps_obj) %>%
    data.frame(row_name = row.names(.)) %>%
    pivot_longer(!row_name,
                 names_to = "rank",
                 values_to = "tax")
  
  # Fill in NAs using the value above
  taxa_long <- taxa_long %>%
    mutate(na = is.na(tax)) %>%
    group_by(row_name) %>%
    fill(tax)
  
  # Create na_labels
  taxa_long <- taxa_long %>%
    mutate(expr = ifelse(na,
                         na_label,
                         tax),
           na_label = str_replace(expr, "<tax>", tax))
  
  # Add the last annotated rank
  if (include_rank){
    taxa_long <- taxa_long %>%
      mutate(last_rank = ifelse(na,
                                NA,
                                rank)) %>%
      fill(last_rank) %>%
      mutate(na_label = str_replace(na_label, "<rank>", last_rank))
  }
  
  # Convert back to tax_table
  taxa_mat <- taxa_long %>%
    select(row_name, rank, na_label) %>%
    pivot_wider(names_from = rank, values_from = na_label) %>%
    as.matrix()
  row.names(taxa_mat) <- taxa_mat[,"row_name"]
  taxa_mat <- taxa_mat[,colnames(taxa_mat) != "row_name"]
  tax_table(ps_obj) <- taxa_mat
  
  return(ps_obj)
}


ps.fungi <- name_na_taxa(ps.fungi, include_rank = F, na_label = "Unknown <tax>")
tax_table(ps.fungi)
ps.fungi.r <- name_na_taxa(ps.fungi.r, include_rank = F, na_label = "Unknown <tax>")
tax_table(ps.fungi.r)
tpfr.df.taxa <- as.data.frame(tax_table(ps.fungi.r))
tpfr.df.taxa <- as.data.frame(tax_table(ps.fungi.r))
class(tpfr.df.taxa)

tpfr.df.taxa$ASV <- rownames(tpfr.df.taxa)







###########################################
# Calculating alpha diversity
###########################################

adiv <- estimate_richness(ps.fungi.r, measures = c("Observed", "Shannon", "Simpson", "Chao1", "InvSimpson"))
adiv$sampleID <- rownames(adiv)
adiv$sampleID <- sub("^X", "", adiv$sampleID)
rownames(adiv) <- NULL

metadata.r <- data.frame(sample_data(ps.fungi.r))

meta_div_full <- dplyr::full_join(metadata.r, adiv, by = "sampleID")





###########################################
# Subsetting metadata and phyloseq objects
###########################################


### Marah only
ps.marah <- subset_samples(ps.fungi.r, host == "marah")
ps.marah <- prune_taxa(taxa_sums(ps.marah) > 0, ps.marah) # remove any taxa with no reads
ps.marah


### Unguiculata only
ps.ung <- subset_samples(ps.fungi.r, host == "unguiculata")
ps.ung <- prune_taxa(taxa_sums(ps.ung) > 0, ps.ung) # remove any taxa with no reads
ps.ung

### Xantiana only
ps.xan <- subset_samples(ps.fungi.r, host == "xantiana")
ps.xan <- prune_taxa(taxa_sums(ps.xan) > 0, ps.xan) # remove any taxa with no reads
ps.xan


### UMX sites only (hereafter will be referenced descriptively in the objects as "ungsites" or "ungs" because unguiculata occurs at those sites)
ps.ungsites <- subset_samples(ps.fungi.r, siteType == "ung_incl")
ps.ungsites <- prune_taxa(taxa_sums(ps.ungsites) > 0, ps.ungsites) # remove any taxa with no reads
ps.ungsites

#xantiana at UMX sites only
ps.xan.ungsites <- subset_samples(ps.ungsites, host == "xantiana")
ps.xan.ungsites <- prune_taxa(taxa_sums(ps.xan.ungsites) > 0, ps.xan.ungsites)
ps.xan.ungsites


#marah at UMX sites only
ps.mar.ungsites <- subset_samples(ps.marah, siteType == "ung_incl")
ps.mar.ungsites <- prune_taxa(taxa_sums(ps.mar.ungsites) > 0, ps.mar.ungsites)
ps.mar.ungsites

### Xan sites only
ps.xansites <- subset_samples(ps.fungi, siteType == "xan_only")
ps.xansites <- prune_taxa(taxa_sums(ps.xansites) > 0, ps.xansites)
ps.xansites

#Xan at xantiana sites only
ps.xan.xansites <- subset_samples(ps.xansites, host == "xantiana")
ps.xan.xansites <- prune_taxa(taxa_sums(ps.xan.xansites) > 0, ps.xan.xansites)
ps.xan.xansites

#Marah at xantiana sites only
ps.mar.xansites <- subset_samples(ps.xansites, host == "marah")
ps.mar.xansites <- prune_taxa(taxa_sums(ps.mar.xansites) > 0, ps.mar.xansites)
ps.mar.xansites

#marah at UMX sites only
ps.mar.ungsites <- subset_samples(ps.marah, siteType == "ung_incl")
ps.mar.ungsites <- prune_taxa(taxa_sums(ps.mar.ungsites) > 0, ps.mar.ungsites)
ps.mar.ungsites
sample_names (ps.mar.ungsites)

### Xan and Marah at XM sites
ps.xanmar <- subset_samples(ps.fungi.r, host != "unguiculata")
ps.xanmar <- prune_taxa(taxa_sums(ps.xanmar) > 0, ps.xanmar) # remove any taxa with no reads
ps.xanmar
ps.xanmar

meta_xan <- meta_div_full %>%
  filter(host == "xantiana")

meta_mar <- meta_div_full %>%
  filter(host == "marah")

meta_xan_ungs <- meta_div_full %>%
  filter(host == "xantiana") %>%
  filter(siteType == "ung_incl")

meta_mar_ungs <- meta_div_full %>%
  filter(host == "marah") %>%
  filter(siteType == "ung_incl")

meta_ung_ungs <- meta_div_full %>%
  filter(host == "unguiculata")

meta_xan_subs <- meta_div_full %>%
  filter(host == "xantiana") %>%
  filter(siteType == "xan_only")


meta_mar_subs <- meta_div_full %>%
  filter(host == "marah") %>%
  filter(siteType == "xan_only")

meta_ungs <- meta_div_full %>%
  filter(siteType == "ung_incl")

meta_xanmar <- meta_div_full %>%
  filter(host != "unguiculata")


###########################################
# subsetting metadata/phyloseq objects- non-rarefied
###########################################
ps.fungi.nr <- ps.fungi
metadata.nr <- data.frame(sample_data(ps.fungi))

### Marah only
ps.marah.nr <- subset_samples(ps.fungi.nr, host == "marah")
ps.marah.nr <- prune_taxa(taxa_sums(ps.marah.nr) > 0, ps.marah.nr) # remove any taxa with no reads
ps.marah.nr


### Unguiculata only
ps.ung.nr <- subset_samples(ps.fungi.nr, host == "unguiculata")
ps.ung.nr <- prune_taxa(taxa_sums(ps.ung.nr) > 0, ps.ung.nr) # remove any taxa with no reads
ps.ung.nr

### Xantiana only
ps.xan.nr <- subset_samples(ps.fungi.nr, host == "xantiana")
ps.xan.nr <- prune_taxa(taxa_sums(ps.xan.nr) > 0, ps.xan.nr) # remove any taxa with no reads
ps.xan.nr


### UMX sites only (hereafter will be referenced descriptively in the objects as "ungsites" or "ungs" because unguiculata occurs at those sites)
ps.ungsites.nr <- subset_samples(ps.fungi.nr, siteType == "ung_incl")
ps.ungsites.nr <- prune_taxa(taxa_sums(ps.ungsites.nr) > 0, ps.ungsites.nr) # remove any taxa with no reads
ps.ungsites.nr

#xantiana at UMX sites only
ps.xan.ungsites.nr <- subset_samples(ps.ungsites.nr, host == "xantiana")
ps.xan.ungsites.nr <- prune_taxa(taxa_sums(ps.xan.ungsites.nr) > 0, ps.xan.ungsites.nr)
ps.xan.ungsites.nr


#marah at UMX sites only
ps.mar.ungsites.nr <- subset_samples(ps.marah.nr, siteType == "ung_incl")
ps.mar.ungsites.nr <- prune_taxa(taxa_sums(ps.mar.ungsites.nr) > 0, ps.mar.ungsites.nr)
ps.mar.ungsites.nr

### Xan sites only
ps.xansites.nr <- subset_samples(ps.fungi.nr, siteType == "xan_only")
ps.xansites.nr <- prune_taxa(taxa_sums(ps.xansites.nr) > 0, ps.xansites.nr)
ps.xansites.nr

#Xan at xantiana sites only
ps.xan.xansites.nr <- subset_samples(ps.xansites.nr, host == "xantiana")
ps.xan.xansites.nr <- prune_taxa(taxa_sums(ps.xan.xansites.nr) > 0, ps.xan.xansites.nr)
ps.xan.xansites.nr

#Marah at xantiana sites only
ps.mar.xansites.nr <- subset_samples(ps.xansites.nr, host == "marah")
ps.mar.xansites.nr <- prune_taxa(taxa_sums(ps.mar.xansites.nr) > 0, ps.mar.xansites.nr)
ps.mar.xansites.nr

#marah at UMX sites only
ps.mar.ungsites.nr <- subset_samples(ps.marah.nr, siteType == "ung_incl")
ps.mar.ungsites.nr <- prune_taxa(taxa_sums(ps.mar.ungsites.nr) > 0, ps.mar.ungsites.nr)
ps.mar.ungsites.nr
sample_names (ps.mar.ungsites.nr)

### Xan and Marah at XM sites
ps.xanmar.nr <- subset_samples(ps.fungi.nr, host != "unguiculata")
ps.xanmar.nr <- prune_taxa(taxa_sums(ps.xanmar.nr) > 0, ps.xanmar.nr) # remove any taxa with no reads
ps.xanmar.nr

meta_xan_nr <- metadata.nr %>%
  filter(host == "xantiana")

meta_mar_nr <- metadata.nr %>%
  filter(host == "marah")

meta_xan_ungs_nr <- metadata.nr %>%
  filter(host == "xantiana") %>%
  filter(siteType == "ung_incl")

meta_mar_ungs_nr <- metadata.nr %>%
  filter(host == "marah") %>%
  filter(siteType == "ung_incl")

meta_ung_ungs_nr <- metadata.nr %>%
  filter(host == "unguiculata")

meta_xan_subs_nr <- metadata.nr %>%
  filter(host == "xantiana") %>%
  filter(siteType == "xan_only")


meta_mar_subs_nr <- metadata.nr %>%
  filter(host == "marah") %>%
  filter(siteType == "xan_only")

meta_ungs_nr <- metadata.nr %>%
  filter(siteType == "ung_incl")

meta_xanmar_nr <- metadata.nr %>%
  filter(host != "unguiculata")

