# tlm-endophytes-2020
Data and code for Mueller & Moeller study on foliar fungal endophyte communities in the Kern River Valley.


FILES IN THIS REPOSITORY:

endo_surv_metadata.csv
dada2 cleaning & processing of raw sequence files.R
processing.R
analysis.R
modis_EVI_extraction.txt


"endo_surv_metadata.csv" contains metadata & environmental data for samples at each site.
Variables:
  sampleID - unique ID for each sequenced sample. 
  siteID - abbreviation code for collection site
  siteName - full name for collection site
  siteType - whether site was in XM or UXM categories (i.e. if C. unguiculata was present/collected)
  host - species of plant host was collected from
  x - longitude
  y - latitude
  elevation - elevation in meters
  soilType - soil parent material type
  ppt - 30-year normals for annual total precipitation (mm) (sourced from PRISM at 800m grid cell resolution)
  solslope - 30-year normals for average daily total shortwave radiation on a sloped ground surface (MJ per square m per day)(sourced from PRISM at 800m grid cell resolution)
  annualtmean - 30-year normals for annual mean temperature (celsius) (sourced from PRISM at 800m grid cell resolution)
  avg_evi - Average Enhanced Vegetation Index (sourced from satellite data at the 250m resolution- see modis_evi_extraction.txt, calculated by averaging all values taken every 2 weeks between 2000 - 2019)


"dada2 cleaning & processing of raw sequence files.R" is an R script that contains code (some of which must be run in R and some of which must be run in zsh (as commented in the script)  to process  fastq files, align, clean and demultiplex sequences, remove primers, create an OTU table, and assign taxonomy.
   
"processing.R" is an R script that contains code to assemble OTU table, sample metadata, and taxonomy table into a phyloseq object, clean & process sequences to filter for only fungal taxa, remove rinsewater controls, pruning, rarefying, and subsetting dataset, and calculating alpha diversity.
   
"analysis.R" is an R script that contains code used for statistical analyses, including: calculating summary statistics, running alpha diversity linear models and Tukey tests, making spatial & environmental distance matrices, mantel tests, differential abundance analysis, and linear mixed-effect models of shared & unique taxa.

"modis_EVI_extraction.txt" is a textfile that contains code used to retrieve and access MODIS satellite data, to be run in Google Earth Engine's code editor, from (https://developers.google.com/earth-engine/datasets/catalog/MODIS_061_MOD13Q1)

   

