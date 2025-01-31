

###################################################################################
################################   ANALYSIS   #####################################
###################################################################################

###########################################
#Summary statistics
###########################################
library(plotrix)
mean(sample_sums(ps.xan.ungsites.nr)) 
std.error(sample_sums(ps.xan.ungsites.nr)) 

summary(ps.mar.ungsites.nr)
mean(sample_sums(ps.mar.ungsites.nr)) 
std.error(sample_sums(ps.mar.ungsites.nr)) 

summary(ps.xan.nr)
mean(sample_sums(ps.xa.nr)) 
std.error(sample_sums(ps.xan.nr)) 

summary(ps.marah.nr)
mean(sample_sums(ps.marah.nr)) 
std.error(sample_sums(ps.marah.nr)) 

summary(ps.ung.nr)
mean(sample_sums(ps.ung.nr)) 
std.error(sample_sums(ps.ung.nr)) 


###########################################
#Alpha diversity linear models and Tukey tests
###########################################
library(car)

mar_xan_obs_lm <- lm(Observed ~ host + PC1 + PC2 + host*PC1 + host*PC2 + PC1*PC2, data = meta_xanmar)
Anova(mar_xan_obs_lm)
summary(mar_xan_obs_lm)

mar_xan_shan_lm <- lm(Shannon ~ host + PC1 + PC2 + host*PC1 + host*PC2 + PC1*PC2, data = meta_xanmar)
Anova(mar_xan_shan_lm)
summary(mar_xan_shan_lm)


ungs_obs_lm <- lm(Observed ~ host + PC1 + PC2 + host:PC1 + host:PC2 + PC1:PC2, data = meta_ungs)
Anova(ungs_obs_lm)
summary(ungs_obs_lm)

ungs_lm_shan <- lm(Shannon ~ host + PC1 + PC2 + host:PC1 + host:PC2 + PC1:PC2, data = meta_ungs)
Anova(ungs_lm_shan)
summary(ungs_lm_shan)

library(tidyverse)
library(janitor)
library(emmeans)

#plot (ungs_obs_lm, host ~ PC1)
emmip(ungs_obs_lm, host ~ PC1*PC2, cov.reduce = range)
emmip(ungs_lm_shan, host ~ PC1 * PC2, cov.reduce = range)



ungs_obs_lm_test <- lm(Observed ~ host*PC1*PC2, data = meta_ungs)
Anova(ungs_obs_lm_test)
Anova(ungs_obs_lm)
summary(ungs_obs_lm_test)

ungs_obs_emm <- emmeans(ungs_obs_lm, ~ host * PC2 * PC1)


# posthoc test of species
categoricalEM <- emmeans(ungs_obs_lm_test, "host")
pairs(categoricalEM)

categoricalEM_host_obs <- emmeans(ungs_obs_lm, "host")
pairs(categoricalEM_host_obs)

categoricalEM_host_shan <- emmeans(ungs_lm_shan, "host")
pairs(categoricalEM_host_shan)

categoricalEM_host <- emmeans(ungs_obs_lm, "host")
pairs(categoricalEM_host)

# posthoc test of interaction

interactionEM_host_PC2_obs <- emtrends(ungs_obs_lm, "host", var="PC2")
pairs(interactionEM_host_PC2_obs)


interactionEM_host_PC2_shan <- emtrends(ungs_lm_shan, "host", var="PC2")
pairs(interactionEM_host_PC2_shan)









###########################################
# Making geographic distance matrix
###########################################

#loading geodist packages (site_dist) + sitelist in geographic north-south order for reference
library(geodist)
site_dist <- metadata

site_dist <- site_dist %>%
  dplyr::filter(siteName != "CTR")

#reference vector of sites in geographic order: 
ref_geo_ordered_sites <- c("site 6", "salmon creek falls", "gold ledge", "corral creek", "site 22", "chico flat", 
                           "camp 3", "headquarters", "site 63", "golf course", "highway 155", "old state road", "lower sawmill", 
                           "upper sawmill", "keyesville", "squirrel mountain", "green rock west", "site 8", "borel road", "black gulch", 
                           "bodfish crest", "hobo", "freeway ridge", "old kern canyon road east", "raft", "old kern canyon road west",
                           "delonegha east", "delonegha west", "mill creek trailhead", "democrat", "china gardens", "cow flat", "cattle pens", 
                           "mimulus pictus", "lucas creek", "doctor", "upper richbar north", "live oak", "drogon", "good place", "saddle mountain road",
                           "site 49", "49B", "site 18", "muddy oaks", "site 17", "site 16", "site 14", "site 85", "site 84")

site_dist_unique = subset(site_dist, select = c("siteID","x", "y"))
site_dist_unique <- unique(site_dist_unique)

unique_distance_matrix <- geodist(site_dist_unique, site_dist_unique, measure = 'geodesic' )/1000
colnames(unique_distance_matrix) <- site_dist_unique$siteID
rownames(unique_distance_matrix) <- site_dist_unique$siteID
unique_distance_matrix




###########################################
# Environmental distance correlation tests
###########################################

library(factoextra)
library(dplyr)
library(ggrepel)
library(ggplot2)

#Making the PCA plots

# Assuming meta_div_full is already loaded
# Prepare the data
sites_unique <- meta_div_full %>%
  select(-c(host, sampleID, Observed, Chao1, se.chao1, Shannon, Simpson, InvSimpson)) %>%
  unique()

site_dist_xan.clean <- sites_unique %>%
  mutate(label = siteName)

unique_pca_env <- site_dist_xan.clean %>% 
  select(siteID, siteType, label, ppt, solslope, annualtmean, avg_evi) %>%
  tidyr::drop_na() # drop any rows with NAs, which can't be included in PCA

xan_pca_env <- unique_pca_env
xan_epca <- prcomp(xan_pca_env[,4:ncol(xan_pca_env)], center=TRUE, scale. = TRUE)

# Calculate the percentage of variance explained by each principal component
var_explained <- xan_epca$sdev^2 / sum(xan_epca$sdev^2) * 100

# Create the site plot
pca_xan_env_plot <- fviz_pca_ind(xan_epca, geom = "point") +
  geom_point(aes(color = xan_pca_env$siteType), size = 2) +
  geom_text_repel(aes(label = xan_pca_env$siteID), max.overlaps = Inf, size = 2.6) +
  scale_color_manual(values = c("xan_only" = "deeppink", "ung_incl" = "slateblue4"),
                     name = "Site Type",
                     labels = c("xan_only" = "MX site", "ung_incl" = "UMX site")) +
  labs(x = paste0("Environment PCA axis 1 (", round(var_explained[1], 1), "%)"),
       y = paste0("Environment PCA axis 2 (", round(var_explained[2], 1), "%)"),
       title = "All sites") +
  theme_minimal() +
  theme(legend.position = "bottom")

pca_xan_env_plot

# Create the biplot with custom labels
custom_labels <- c("ppt" = "Precipitation", 
                   "solslope" = "Solar Radiation", 
                   "annualtmean" = "Temperature", 
                   "avg_evi" = "EVI")

biplot_xan_env <- fviz_pca_var(xan_epca, geom = "arrow") +
  geom_text_repel(aes(x = xan_epca$rotation[,1]*1.3, 
                      y = xan_epca$rotation[,2]*1.35, 
                      label = custom_labels[rownames(xan_epca$rotation)]),
                  size = 4,
                  point.padding = 0.5,
  ) +
  labs(x = paste0("Environment PCA axis 1 (", round(var_explained[1], 1), "%)"),
       y = paste0("Environment PCA axis 2 (", round(var_explained[2], 1), "%)"),
       title = "Environmental parameters") +
  theme_minimal() +
  xlim(-1.2, 1) +  # Adjust these limits as needed
  ylim(-1.2, 1)    # Adjust these limits as needed

# Display the plots
print(biplot_xan_env)
print(pca_xan_env_plot)
print(biplot_xan_env)

# Save the plots
ggsave("pca_all_sites.png", pca_xan_env_plot, width = 10, height = 8, dpi = 300)
ggsave("pca_biplot.png", biplot_xan_env, width = 10, height = 8, dpi = 300)


names_xan <- xan_pca_env$sampleID

#get the env matrix for all sites

# extract axis scores
scores.xan_epca<-cbind(xan_pca_env[,c('siteID')], xan_epca$x[,1:4])
rownames(scores.xan_epca) <- xan_pca_env$siteID

write.csv(scores.xan_epca, "Sitewide Unique PCA scores.cl.up.csv")

# extract axis loading
loading.xan_epca<-cbind(rownames(xan_epca$rotation), xan_epca$rotation[,1:4])

write.csv(loading.xan_epca, "unique sites pca axis loading.csv")

#Calculating environmental distance from multidimensional PCA axis space
xan_epca<-as.data.frame(
  as.matrix(
    dist(
      scores.xan_epca[,3:ncol(
        scores.xan_epca)], method="euclidean", diag=T, upper=T
    ), nrow=length(names_xan), ncol=length(names_xan)
  )
)
xan_epca
all_envdist_multi_pca <- xan_epca

#graph the correlation for all sites
unique_distance_matrix
all_envdist_multi_pca <- as.matrix(all_envdist_multi_pca)
class(all_envdist_multi_pca)
ud = as.vector(unique_distance_matrix)
aep = as.vector(all_envdist_multi_pca)

#new data frame with vectorized distance matrices
mat = data.frame(ud,aep)

#env dist vs geographic distance
edist_geodist = ggplot(mat, aes(y = aep, x = ud)) + 
  geom_point(size = 1, alpha = 0.5) + 
  labs(x = "Euclidean PCA distance", y = "Geographic distance (km)") + 
  theme( axis.text.x = element_text(face = "bold",colour = "black", size = 12), 
         axis.text.y = element_text(face = "bold", size = 11, colour = "black"), 
         axis.title= element_text(face = "bold", size = 14, colour = "black"), 
         panel.background = element_blank(), 
         panel.border = element_rect(fill = NA, colour = "black"))
edist_geodist

#for UMX sites:

# Assuming xan_pca_env, scores.xan_epca, and unique_distance_matrix are already in the environment

# Convert scores.xan_epca to a data frame if it's not already one
if(!is.data.frame(scores.xan_epca)) {
  scores.xan_epca <- as.data.frame(scores.xan_epca)
  colnames(scores.xan_epca) <- c("siteID", paste0("PC", 1:(ncol(scores.xan_epca)-1)))
}

# Create a subset of 19 "ung_incl" sites
umx_subset <- xan_pca_env %>%
  filter(siteType == "ung_incl") %>%
  slice_sample(n = 19)

# Extract PCA scores for the subset
scores.umx_epca <- scores.xan_epca[scores.xan_epca$siteID %in% umx_subset$siteID, ]

# Calculate environmental distance matrix for UMX subset
umx_epca <- as.data.frame(
  as.matrix(
    dist(
      scores.umx_epca[, -1],  # Exclude the siteID column
      method = "euclidean", 
      diag = TRUE, 
      upper = TRUE
    )
  )
)

# Extract geographic distance for UMX subset
umx_geodist <- unique_distance_matrix[umx_subset$siteID, umx_subset$siteID]

# Convert distance matrices to vectors
umx_env_dist <- as.vector(as.matrix(umx_epca))
umx_geo_dist <- as.vector(as.matrix(umx_geodist))

# Create a new data frame with vectorized distance matrices
umx_mat <- data.frame(
  env_dist = umx_env_dist,
  geo_dist = umx_geo_dist
)

# Create the plot
umx_edist_geodist <- ggplot(umx_mat, aes(x = env_dist, y = geo_dist)) +
  geom_point(size = 1, alpha = 0.5) +
  labs(
    x = "Environmental PCA Distance",
    y = "Geographic Distance (km)",
    title = "UMX Sites: Environmental vs Geographic Distance"
  ) +
  theme(
    axis.text.x = element_text(face = "bold", colour = "black", size = 12),
    axis.text.y = element_text(face = "bold", size = 11, colour = "black"),
    axis.title = element_text(face = "bold", size = 14, colour = "black"),
    panel.background = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black"),
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5)
  )

# Display the plot
print(umx_edist_geodist)

# Save the plot
ggsave("umx_env_vs_geo_distance.png", umx_edist_geodist, width = 10, height = 8, dpi = 300)

# Perform Mantel test on UMX subset
umx_mantel <- vegan::mantel(umx_epca, umx_geodist, method = "spearman", permutations = 999)

# Print Mantel test results
print(umx_mantel)





###########################################
# Mantel tests of environmental & spatial distance
###########################################
      #####################################
      # for xantiana at XM sites
      ####################################

#making geo distance matrix for xan at xm sites

site_dist_xans.xan <- site_dist %>%
  #filter(siteType == "xan_only") %>%
  dplyr::filter(host == "xantiana") %>%
  dplyr::filter( sampleID != "GCX") %>%
  dplyr::filter( sampleID != "LOX") %>%
  dplyr::select("sampleID","x", "y")

site_dist_xans.xan
xan.xan_distance_matrix <- geodist(site_dist_xans.xan, site_dist_xans.xan, measure = 'geodesic' )/1000
colnames(xan.xan_distance_matrix) <- site_dist_xans.xan$sampleID
rownames(xan.xan_distance_matrix) <- site_dist_xans.xan$sampleID

xan.xan_distance_matrix

write.csv(xan.xan_distance_matrix, "xan_xansites_cl_up.csv")

#making environmental PCA distance matrix 

meta_div_full


#XAN at XM sites
site_dist_XAN.clean <- meta_div_full %>%
  filter(host == "marah") %>%
  filter( sampleID != "GCX") %>%
  filter( sampleID != "LOX") %>%
  mutate(label = siteName)

xan_pca_env<- site_dist_xan.clean %>% 
  select(sampleID, PC1, PC2, PC3, PC4) %>%
  tidyr::drop_na() # drop any rows with NAs, which can't be included in PCA

scores.xan_epca<-xan_pca_env

rownames(scores.xan_epca) <- xan_pca_env$sampleID

names_xan <- xan_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
names_xan <- xan_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
xan_envdist_multi_pca<-as.data.frame(
  as.matrix(
    dist(
      scores.xan_epca[,3:ncol(
        scores.xan_epca)], method="manhattan", diag=T, upper=T
    ), nrow=length(names_xan), ncol=length(names_xan)
  )
)

#creating species richness & alpha div matrices for xan at all sites

#### environment x diversity questions

### Observed richness

obs_mat_xan2 <- meta_xan %>%
  select(sampleID, Observed)

#create observed matrix for each species
obs_matdiff_xan2 <- outer(obs_mat_xan2$Observed, obs_mat_xan2$Observed, "-")
colnames(obs_matdiff_xan2) <- obs_mat_xan2$sampleID
rownames(obs_matdiff_xan2) <- obs_mat_xan2$sampleID
obs_matdist_xan2 <- as.dist(obs_matdiff_xan2)
obs_matdist_xan2 <- abs(obs_matdist_xan2)
obs_matdist_xan2

### Shannon diversity

shan_mat_xan2 <- meta_xan %>%
  select(sampleID, Shannon)

#create Shannon matrix
shan_matdiff_xan2 <- outer(shan_mat_xan2$Shannon, shan_mat_xan2$Shannon, "-")
colnames(shan_matdiff_xan2) <- shan_mat_xan2$sampleID
rownames(shan_matdiff_xan2) <- shan_mat_xan2$sampleID
shan_matdist_xan2 <- as.dist(shan_matdiff_xan2)
shan_matdist_xan2 <- abs(shan_matdist_xan2)
shan_matdist_xan2


#mantel tests for env distance, geo distance, & community distance for xantiana at UMX sites

xan_edist <- xan_envdist_multi_pca

xan_geodist_man <- as.dist(xan.xan_distance_matrix)


xan_edist
write.csv(xan_edist, "xan_pca4_envdist.cl.up.csv")

xan_edist_man <- read.csv("xan_pca4_envdist.cl.up.csv", header = TRUE, row.names = 1)
xan_edist_man <- as.dist(xan_edist_man)


##### mantels
xan_env_geo_mantel <- vegan::mantel(xan_geodist_man, xan_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_env_geo_mantel

ps.xan.bray <- phyloseq::distance(ps.xan, method = "bray")
ps.xan.jaccard <- phyloseq::distance(ps.xan, method = "jaccard")

plot(xan_edist_man,ps.xan.bray,pch=1,cex=0.5,col="black",bty="l")
plot(xan_edist_man,xan_geodist_man,pch=1,cex=0.5,col="black",bty="l")

xan_geodist_man
xan_edist_man

#env mantel
xan_env_div_mantel <- vegan::mantel(ps.xan.bray, xan_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_env_div_mantel

#geo mantel
xan_geo_div_mantel <- vegan::mantel(ps.xan.bray, xan_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_geo_div_mantel


#partial mantel of env controlling for geo
xan_env_part <- mantel.partial(ps.xan.bray, xan_edist_man, xan_geodist_man, method="pearson", permutations=9999)
xan_env_part

#partial mantel of geo controlling for env
xan_geo_part <- mantel.partial(ps.xan.bray, xan_geodist_man, xan_edist_man, method="pearson", permutations=9999)
xan_geo_part

##### jaccard versions
#env mantel
xan_env_div_mantel.j <- vegan::mantel(ps.xan.jaccard, xan_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_env_div_mantel.j

#geo mantel
xan_geo_div_mantel.j <- vegan::mantel(ps.xan.jaccard, xan_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_geo_div_mantel.j


#partial mantel of env controlling for geo
xan_env_part.j <- mantel.partial(ps.xan.jaccard, xan_edist_man, xan_geodist_man, method="pearson", permutations=9999)
xan_env_part.j

#partial mantel of geo controlling for env
xan_geo_part.j <- mantel.partial(ps.xan.jaccard, xan_geodist_man, xan_edist_man, method="pearson", permutations=9999)
xan_geo_part.j





#### Observed: 

#env mantel
xan_env_mantel_obs <- vegan::mantel(obs_matdist_xan2, xan_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_env_mantel_obs

#geo mantel
xan_geo_mantel_obs <- vegan::mantel(obs_matdist_xan2, xan_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_geo_mantel_obs
#.

#partial mantel of env controlling for geo
xan_env_part_obs <- mantel.partial(obs_matdist_xan2, xan_edist_man, xan_geodist_man, method="pearson", permutations=9999)
xan_env_part_obs

xan_geo_part_obs <- mantel.partial(obs_matdist_xan2, xan_geodist_man, xan_edist_man, method="pearson", permutations=9999)
xan_geo_part_obs
#.


#### Shannon: 

#env mantel
xan_env_mantel_shan <- vegan::mantel(shan_matdist_xan2, xan_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_env_mantel_shan

#geo mantel
xan_geo_mantel_shan <- vegan::mantel(shan_matdist_xan2, xan_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_geo_mantel_shan
#**

#partial mantel of env controlling for geo
xan_env_part_shan <- mantel.partial(shan_matdist_xan2, xan_edist_man, xan_geodist_man, method="pearson", permutations=9999)
xan_env_part_shan

xan_geo_part_shan <- mantel.partial(shan_matdist_xan2, xan_geodist_man, xan_edist_man, method="pearson", permutations=9999)
xan_geo_part_shan
#**


library(vegan)
xan.xans.corlog <- mantel.correlog(ps.xan.bray, D.geo=xan_geodist_man, n.class=0, break.pts=NULL, 
                                   cutoff=FALSE, r.type="pearson", nperm=999, mult="holm", progressive=TRUE)

summary(xan.xans.corlog)
plot(xan.xans.corlog)



      #####################################
      # for marah at XM sites
      ####################################
#geodist for marah at all sites

site_dist_mar <- site_dist %>%
  #filter(siteType == "xan_only") %>%
  filter(host == "marah") %>%
  filter(sampleID != "18M") %>%
  filter(sampleID != "MOM") %>%
  filter(sampleID != "FRM") %>%
  filter(sampleID != "OSRM") %>%
  select("sampleID","x", "y")

site_dist_mar
mar_distance_matrix <- geodist(site_dist_mar, site_dist_mar, measure = 'geodesic' )/1000
colnames(mar_distance_matrix) <- site_dist_mar$sampleID
rownames(mar_distance_matrix) <- site_dist_mar$sampleID

mar_distance_matrix

write.csv(mar_distance_matrix, "mar_allsites_cl_up.csv")

#making environmental PCA distance matrix

meta_div_full


#MARAH ALL SITES
site_dist_mar.clean <- meta_div_full %>%
  #filter(siteType == "xan_only") %>%
  filter(host == "marah") %>%
  filter( sampleID != "18M") %>%
  filter( sampleID != "MOM") %>%
  filter( sampleID != "CHFM") %>%
  filter( sampleID != "OSRM") %>%
  mutate(label = siteName)

mar_pca_env<- site_dist_mar.clean %>% 
  select(sampleID, PC1, PC2, PC3, PC4) %>%
  tidyr::drop_na() # drop any rows with NAs, which can't be included in PCA

scores.mar_epca<-mar_pca_env

rownames(scores.mar_epca) <- mar_pca_env$sampleID

names_mar <- mar_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
names_mar <- mar_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
mar_envdist_multi_pca<-as.data.frame(
  as.matrix(
    dist(
      scores.mar_epca[,3:ncol(
        scores.mar_epca)], method="manhattan", diag=T, upper=T
    ), nrow=length(names_mar), ncol=length(names_mar)
  )
)

# creating species richness & alpha div matrices for marah at all sites

#### environment x diversity questions

### Observed richness

obs_mat_mar2 <- meta_mar2 %>%
  select(sampleID, Observed)

#create observed matrix for each species
obs_matdiff_mar2 <- outer(obs_mat_mar2$Observed, obs_mat_mar2$Observed, "-")
colnames(obs_matdiff_mar2) <- obs_mat_mar2$sampleID
rownames(obs_matdiff_mar2) <- obs_mat_mar2$sampleID
obs_matdist_mar2 <- as.dist(obs_matdiff_mar2)
obs_matdist_mar2 <- abs(obs_matdist_mar2)
obs_matdist_mar2

### Shannon diversity

shan_mat_mar2 <- meta_mar2 %>%
  select(sampleID, Shannon)

#create observed matrix for each species
shan_matdiff_mar2 <- outer(shan_mat_mar2$Shannon, shan_mat_mar2$Shannon, "-")
colnames(shan_matdiff_mar2) <- shan_mat_mar2$sampleID
rownames(shan_matdiff_mar2) <- shan_mat_mar2$sampleID
shan_matdist_mar2 <- as.dist(shan_matdiff_mar2)
shan_matdist_mar2 <- abs(shan_matdist_mar2)
shan_matdist_mar2


#mantel tests for envd + geod + comd, impt env dist and geo dist- mar sites

#mar_geodist
mar_edist <- mar_envdist_multi_pca

#.geo.dist.mar <- read.csv("mar_sample_geodist.csv", header = TRUE, row.names = 1)

mar_geodist_man <- as.dist(mar_distance_matrix)
mar_geodist <- mar_geodist_man


mar_edist
write.csv(mar_edist, "mar_pca4_envdist.csv")

mar_edist_man <- read.csv("mar_pca4_envdist.csv", header = TRUE, row.names = 1)
mar_edist_man <- as.dist(mar_edist_man)

sample_names(ps.marah)
##### mantels
mar_env_geo_mantel <- vegan::mantel(mar_geodist_man, mar_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_env_geo_mantel

ps.mar.bray <- phyloseq::distance(ps.marah, method = "bray")
ps.mar.jaccard <- phyloseq::distance(ps.marah, method = "jaccard")

ps.marah

plot(mar_edist_man,ps.mar.bray,pch=1,cex=0.5,col="black",bty="l")
plot(mar_edist_man, mar_geodist_man, pch=1,cex=0.5,col="black",bty="l")
plot(mar_geodist_man,ps.mar.bray,pch=1,cex=0.5,col="black",bty="l")

mar_geodist_man
mar_edist_man
#env mantel
mar_env_div_mantel <- vegan::mantel(ps.mar.bray, mar_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_env_div_mantel

#geo mantel
mar_geo_div_mantel <- vegan::mantel(ps.mar.bray, mar_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_geo_div_mantel

#partial mantel of env controlling for geo
mar_env_part <- mantel.partial(ps.mar.bray, mar_edist_man, mar_geodist_man, method="pearson", permutations=9999)
mar_env_part

#partial mantel of geo controlling for env
mar_geo_part <- mantel.partial(ps.mar.bray, mar_geodist_man, mar_edist_man, method="pearson", permutations=9999)
mar_geo_part

##### jaccard versions
#env mantel
mar_env_div_mantel.j <- vegan::mantel(ps.mar.jaccard, mar_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_env_div_mantel.j

#geo mantel
mar_geo_div_mantel.j <- vegan::mantel(ps.mar.jaccard, mar_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_geo_div_mantel.j


#partial mantel of env controlling for geo
mar_env_part.j <- mantel.partial(ps.mar.jaccard, mar_edist_man, mar_geodist_man, method="pearson", permutations=9999)
mar_env_part.j

#partial mantel of geo controlling for env
mar_geo_part.j <- mantel.partial(ps.mar.jaccard, mar_geodist_man, mar_edist_man, method="pearson", permutations=9999)
mar_geo_part.j

#### Observed: 

#env mantel
mar_env_mantel_obs <- vegan::mantel(obs_matdist_mar2, mar_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_env_mantel_obs

#geo mantel
mar_geo_mantel_obs <- vegan::mantel(obs_matdist_mar2, mar_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_geo_mantel_obs

#partial mantel of env controlling for geo
mar_env_part_obs <- mantel.partial(obs_matdist_mar2, mar_edist_man, mar_geodist_man, method="pearson", permutations=9999)
mar_env_part_obs

mar_geo_part_obs <- mantel.partial(obs_matdist_mar2, mar_geodist_man, mar_edist_man, method="pearson", permutations=9999)
mar_geo_part_obs


#### Shannon: 

#env mantel
mar_env_mantel_shan <- vegan::mantel(shan_matdist_mar2, mar_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_env_mantel_shan

#geo mantel
mar_geo_mantel_shan <- vegan::mantel(shan_matdist_mar2, mar_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_geo_mantel_shan

#partial mantel of env controlling for geo
mar_env_part_shan <- mantel.partial(shan_matdist_mar2, mar_edist_man, mar_geodist_man, method="pearson", permutations=9999)
mar_env_part_shan

mar_geo_part_shan <- mantel.partial(shan_matdist_mar2, mar_geodist_man, mar_edist_man, method="pearson", permutations=9999)
mar_geo_part_shan



      #####################################
      # for marah at UMX sites
      #####################################

MAR AT UNG SITES 
#geodist for marah at unguic sites

site_dist
site_dist_mar_ungs <- site_dist %>%
  filter(siteType == "ung_incl")%>%
  filter(host == "marah") %>%
  filter(sampleID != "18M") %>%
  filter(sampleID != "MOM") %>%
  filter(sampleID != "CHFM") %>%
  filter(sampleID != "OSRM") %>%
  select("sampleID","x", "y")

site_dist_mar_ungs
mar_ungs_distance_matrix <- geodist(site_dist_mar_ungs, site_dist_mar_ungs, measure = 'geodesic' )/1000
colnames(mar_ungs_distance_matrix) <- site_dist_mar_ungs$sampleID
rownames(mar_ungs_distance_matrix) <- site_dist_mar_ungs$sampleID

mar_ungs_distance_matrix

write.csv(mar_ungs_distance_matrix, "mar_ungsites_cl_up.csv")


#making environmental PCA distance matrix 

meta_div_full


#MARAH UNG SITES
site_dist_mar_ungs.clean <- meta_div_full %>%
  filter(siteType == "ung_incl")%>%
  filter(host == "marah") %>%
  filter(sampleID != "18M") %>%
  filter(sampleID != "MOM") %>%
  filter(sampleID != "CHFM") %>%
  filter(sampleID != "OSRM") %>%
  mutate(label = siteName)

mar_ungs_pca_env<- site_dist_mar_ungs.clean %>% 
  select(sampleID, PC1, PC2, PC3, PC4) %>%
  tidyr::drop_na() # drop any rows with NAs, which can't be included in PCA

scores.mar_ungs_epca<-mar_ungs_pca_env

rownames(scores.mar_ungs_epca) <- mar_ungs_pca_env$sampleID

names_mar_ungs <- mar_ungs_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
names_mar_ungs <- mar_ungs_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
mar_ungs_envdist_multi_pca<-as.data.frame(
  as.matrix(
    dist(
      scores.mar_ungs_epca[,3:ncol(
        scores.mar_ungs_epca)], method="manhattan", diag=T, upper=T
    ), nrow=length(names_mar_ungs), ncol=length(names_mar_ungs)
  )
)

#creating species richness & alpha div matrices for marah at UNG sites

#### environment x diversity questions

### Observed richness

obs_mat_mar_ungs2 <- meta_mar_ungs2 %>%
  select(sampleID, Observed)

#create observed matrix for each species
obs_matdiff_mar_ungs2 <- outer(obs_mat_mar_ungs2$Observed, obs_mat_mar_ungs2$Observed, "-")
colnames(obs_matdiff_mar_ungs2) <- obs_mat_mar_ungs2$sampleID
rownames(obs_matdiff_mar_ungs2) <- obs_mat_mar_ungs2$sampleID
obs_matdist_mar_ungs2 <- as.dist(obs_matdiff_mar_ungs2)
obs_matdist_mar_ungs2 <- abs(obs_matdist_mar_ungs2)
obs_matdist_mar_ungs2

### Shannon diversity

shan_mat_mar_ungs2 <- meta_mar_ungs2 %>%
  select(sampleID, Shannon)

#create observed matrix for each species
shan_matdiff_mar_ungs2 <- outer(shan_mat_mar_ungs2$Shannon, shan_mat_mar_ungs2$Shannon, "-")
colnames(shan_matdiff_mar_ungs2) <- shan_mat_mar_ungs2$sampleID
rownames(shan_matdiff_mar_ungs2) <- shan_mat_mar_ungs2$sampleID
shan_matdist_mar_ungs2 <- as.dist(shan_matdiff_mar_ungs2)
shan_matdist_mar_ungs2 <- abs(shan_matdist_mar_ungs2)
shan_matdist_mar_ungs2

#mantel tests for envd + geod + comd, impt env and geo dist- mar AT UNG sites
#mar_ungs_geodist
mar_ungs_edist <- mar_ungs_envdist_multi_pca
mar_ungs_edist_man <- mar_ungs_edist
#geo.dist.mar <- read.csv("mar_sample_geodist.csv", header = TRUE, row.names = 1)

mar_ungs_geodist_man <- as.dist(mar_ungs_distance_matrix)
mar_ungs_geodist <- mar_ungs_geodist_man

mar_ungs_edist
write.csv(mar_ungs_edist, "mar_ungs_pca4_envdist.csv")

mar_ungs_edist_man <- read.csv("mar_ungs_pca4_envdist.csv", header = TRUE, row.names = 1)
mar_ungs_edist_man <- as.dist(mar_ungs_edist_man)


##### mantels
mar_ungs_env_geo_mantel <- vegan::mantel(mar_ungs_geodist_man, mar_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_ungs_env_geo_mantel

ps.mar.ungs.bray <- phyloseq::distance(ps.mar.ungsites, method = "bray")
ps.mar.ungs.jaccard <- phyloseq::distance(ps.mar.ungsites, method = "jaccard")

ps.marah

plot(mar_ungs_edist_man,ps.mar.ungs.bray,pch=1,cex=0.5,col="black",bty="l")
plot(mar_ungs_edist_man, mar_ungs_geodist, pch=1,cex=0.5,col="black",bty="l")
plot(mar_ungs_geodist,ps.mar.ungs.bray,pch=1,cex=0.5,col="black",bty="l")

mar_ungs_geodist_man
mar_ungs_edist_man
#env mantel
mar_ungs_env_div_mantel <- vegan::mantel(ps.mar.ungs.bray, mar_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_ungs_env_div_mantel

#geo mantel
mar_ungs_geo_div_mantel <- vegan::mantel(ps.mar.ungs.bray, mar_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_ungs_geo_div_mantel

#partial mantel of env controlling for geo
mar_ungs_env_part <- mantel.partial(ps.mar.ungs.bray, mar_ungs_edist_man, mar_ungs_geodist_man, method="pearson", permutations=9999)
mar_ungs_env_part

#partial mantel of geo controlling for env
mar_ungs_geo_part <- mantel.partial(ps.mar.ungs.bray, mar_ungs_geodist_man, mar_ungs_edist_man, method="pearson", permutations=9999)
mar_ungs_geo_part

##### jaccard versions
#env mantel
mar_ungs_env_div_mantel.j <- vegan::mantel(ps.mar.ungs.jaccard, mar_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_ungs_env_div_mantel.j

#geo mantel
mar_ungs_geo_div_mantel.j <- vegan::mantel(ps.mar.ungs.jaccard, mar_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_ungs_geo_div_mantel.j

#partial mantel of env controlling for geo
mar_ungs_env_part.j <- mantel.partial(ps.mar.ungs.jaccard, mar_ungs_edist_man, mar_ungs_geodist_man, method="pearson", permutations=9999)
mar_ungs_env_part.j

#partial mantel of geo controlling for env
mar_ungs_geo_part.j <- mantel.partial(ps.mar.ungs.jaccard, mar_ungs_geodist_man, mar_ungs_edist_man, method="pearson", permutations=9999)
mar_ungs_geo_part.j




#### Observed: 

#env mantel
mar_ungs_env_mantel_obs <- vegan::mantel(obs_matdist_mar_ungs2, mar_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_ungs_env_mantel_obs

#geo mantel
mar_ungs_geo_mantel_obs <- vegan::mantel(obs_matdist_mar_ungs2, mar_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_ungs_geo_mantel_obs

#partial mantel of env controlling for geo
mar_ungs_env_part_obs <- mantel.partial(obs_matdist_mar_ungs2, mar_ungs_edist_man, mar_ungs_geodist_man, method="pearson", permutations=9999)
mar_ungs_env_part_obs

mar_ungs_geo_part_obs <- mantel.partial(obs_matdist_mar_ungs2, mar_ungs_geodist_man, mar_ungs_edist_man, method="pearson", permutations=9999)
mar_ungs_geo_part_obs

#### Shannon: 

#env mantel
mar_ungs_env_mantel_shan <- vegan::mantel(shan_matdist_mar_ungs2, mar_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_ungs_env_mantel_shan

#geo mantel
mar_ungs_geo_mantel_shan <- vegan::mantel(shan_matdist_mar_ungs2, mar_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
mar_ungs_geo_mantel_shan
#*

#partial mantel of env controlling for geo
mar_ungs_env_part_shan <- mantel.partial(shan_matdist_mar_ungs2, mar_ungs_edist_man, mar_ungs_geodist_man, method="pearson", permutations=9999)
mar_ungs_env_part_shan
#*

mar_ungs_geo_part_shan <- mantel.partial(shan_matdist_mar_ungs2, mar_ungs_geodist_man, mar_ungs_edist_man, method="pearson", permutations=9999)
mar_ungs_geo_part_shan

      #####################################
      # for xantiana at UMX sites
      #####################################
#geodist for XAN at unguic sites

site_dist
site_dist_xan_ungs <- site_dist %>%
  filter(siteType == "ung_incl")%>%
  filter(host == "xantiana") %>%
  filter(sampleID != "GCX") %>%
  filter(sampleID != "LOX") %>%
  select("sampleID","x", "y")

site_dist_xan_ungs
xan_ungs_distance_matrix <- geodist(site_dist_xan_ungs, site_dist_xan_ungs, measure = 'geodesic' )/1000
colnames(xan_ungs_distance_matrix) <- site_dist_xan_ungs$sampleID
rownames(xan_ungs_distance_matrix) <- site_dist_xan_ungs$sampleID

xan_ungs_distance_matrix

write.csv(xan_ungs_distance_matrix, "xan_ungsites_cl_up.csv")

#making environmental PCA distance matrix 

meta_div_full


#XAN UNG SITES
site_dist_xan_ungs.clean <- site_dist %>%
  filter(siteType == "ung_incl")%>%
  filter(host == "xantiana") %>%
  filter(sampleID != "GCX") %>%
  filter(sampleID != "LOX") %>%
  mutate(label = siteName)

xan_ungs_pca_env<- site_dist_xan_ungs.clean %>% 
  select(sampleID, PC1, PC2, PC3, PC4) %>%
  tidyr::drop_na() # drop any rows with NAs, which can't be included in PCA

scores.xan_ungs_epca<-xan_ungs_pca_env

rownames(scores.xan_ungs_epca) <- xan_ungs_pca_env$sampleID

names_xan_ungs <- xan_ungs_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
names_xan_ungs <- xan_ungs_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
xan_ungs_envdist_multi_pca<-as.data.frame(
  as.matrix(
    dist(
      scores.xan_ungs_epca[,3:ncol(
        scores.xan_ungs_epca)], method="manhattan", diag=T, upper=T
    ), nrow=length(names_xan_ungs), ncol=length(names_xan_ungs)
  )
)
#creating species richness & alpha div matrices for xantiana at UNG sites}

#### environment x diversity questions

### Observed richness

obs_mat_xan_ungs2 <- meta_xan_ungs2 %>%
  select(sampleID, Observed)

#create observed matrix for each species
obs_matdiff_xan_ungs2 <- outer(obs_mat_xan_ungs2$Observed, obs_mat_xan_ungs2$Observed, "-")
colnames(obs_matdiff_xan_ungs2) <- obs_mat_xan_ungs2$sampleID
rownames(obs_matdiff_xan_ungs2) <- obs_mat_xan_ungs2$sampleID
obs_matdist_xan_ungs2 <- as.dist(obs_matdiff_xan_ungs2)
obs_matdist_xan_ungs2 <- abs(obs_matdist_xan_ungs2)
obs_matdist_xan_ungs2

### Shannon diversity

shan_mat_xan_ungs2 <- meta_xan_ungs2 %>%
  select(sampleID, Shannon)

#create observed matrix for each species
shan_matdiff_xan_ungs2 <- outer(shan_mat_xan_ungs2$Shannon, shan_mat_xan_ungs2$Shannon, "-")
colnames(shan_matdiff_xan_ungs2) <- shan_mat_xan_ungs2$sampleID
rownames(shan_matdiff_xan_ungs2) <- shan_mat_xan_ungs2$sampleID
shan_matdist_xan_ungs2 <- as.dist(shan_matdiff_xan_ungs2)
shan_matdist_xan_ungs2 <- abs(shan_matdist_xan_ungs2)
shan_matdist_xan_ungs2

# mantel tests for envd + geod + comd, impt env and geo dist- xan AT UNG sites
#xan_ungs_geodist
xan_ungs_edist <- xan_ungs_envdist_multi_pca


xan_ungs_geodist_man <- as.dist(xan_ungs_distance_matrix)
xan_ungs_geodist <- xan_ungs_geodist_man


xan_ungs_edist
write.csv(xan_ungs_edist, "xan_ungs_pca4_envdist.csv")

xan_ungs_edist_man <- read.csv("xan_ungs_pca4_envdist.csv", header = TRUE, row.names = 1)
xan_ungs_edist_man <- as.dist(xan_ungs_edist_man)


##### mantels
xan_ungs_env_geo_mantel <- vegan::mantel(xan_ungs_geodist_man, xan_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_ungs_env_geo_mantel

ps.xan.ungs.bray <- phyloseq::distance(ps.xan.ungsites, method = "bray")
ps.xan.ungs.jaccard <- phyloseq::distance(ps.xan.ungsites, method = "jaccard")


plot(xan_ungs_edist_man,ps.xan.ungs.bray,pch=1,cex=0.5,col="black",bty="l")
plot(xan_ungs_edist_man, xan_ungs_geodist, pch=1,cex=0.5,col="black",bty="l")
plot(xan_ungs_geodist,ps.xan.ungs.bray,pch=1,cex=0.5,col="black",bty="l")

xan_ungs_geodist
xan_ungs_edist_man
#env mantel
xan_ungs_env_div_mantel <- vegan::mantel(ps.xan.ungs.bray, xan_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_ungs_env_div_mantel

#geo mantel
xan_ungs_geo_div_mantel <- vegan::mantel(ps.xan.ungs.bray, xan_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_ungs_geo_div_mantel

#partial mantel of env controlling for geo
xan_ungs_env_part <- mantel.partial(ps.xan.ungs.bray, xan_ungs_edist_man, xan_ungs_geodist_man, method="pearson", permutations=9999)
xan_ungs_env_part

#partial mantel of geo controlling for env
xan_ungs_geo_part <- mantel.partial(ps.xan.ungs.bray, xan_ungs_geodist_man, xan_ungs_edist_man, method="pearson", permutations=9999)
xan_ungs_geo_part

##### jaccard versions
#env mantel
xan_ungs_env_div_mantel.j <- vegan::mantel(ps.xan.ungs.jaccard, xan_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_ungs_env_div_mantel.j

#geo mantel
xan_ungs_geo_div_mantel.j <- vegan::mantel(ps.xan.ungs.jaccard, xan_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_ungs_geo_div_mantel.j

#partial mantel of env controlling for geo
xan_ungs_env_part.j <- mantel.partial(ps.xan.ungs.jaccard, xan_ungs_edist_man, xan_ungs_geodist_man, method="pearson", permutations=9999)
xan_ungs_env_part.j

#partial mantel of geo controlling for env
xan_ungs_geo_part.j <- mantel.partial(ps.xan.ungs.jaccard, xan_ungs_geodist_man, xan_ungs_edist_man, method="pearson", permutations=9999)
xan_ungs_geo_part.j

#### Observed: 

#env mantel
xan_ungs_env_mantel_obs <- vegan::mantel(obs_matdist_xan_ungs2, xan_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_ungs_env_mantel_obs

#geo mantel
xan_ungs_geo_mantel_obs <- vegan::mantel(obs_matdist_xan_ungs2, xan_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_ungs_geo_mantel_obs

#partial mantel of env controlling for geo
xan_ungs_env_part_obs <- mantel.partial(obs_matdist_xan_ungs2, xan_ungs_edist_man, xan_ungs_geodist_man, method="pearson", permutations=9999)
xan_ungs_env_part_obs

xan_ungs_geo_part_obs <- mantel.partial(obs_matdist_xan_ungs2, xan_ungs_geodist_man, xan_ungs_edist_man, method="pearson", permutations=9999)
xan_ungs_geo_part_obs

#### Shannon: 

#env mantel
xan_ungs_env_mantel_shan <- vegan::mantel(shan_matdist_xan_ungs2, xan_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_ungs_env_mantel_shan

#geo mantel
xan_ungs_geo_mantel_shan <- vegan::mantel(shan_matdist_xan_ungs2, xan_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
xan_ungs_geo_mantel_shan
#**

#partial mantel of env controlling for geo
xan_ungs_env_part_shan <- mantel.partial(shan_matdist_xan_ungs2, xan_ungs_edist_man, xan_ungs_geodist_man, method="pearson", permutations=9999)
xan_ungs_env_part_shan

xan_ungs_geo_part_shan <- mantel.partial(shan_matdist_xan_ungs2, xan_ungs_geodist_man, xan_ungs_edist_man, method="pearson", permutations=9999)
xan_ungs_geo_part_shan
#*


      #####################################
      # for unguiculata at UMX sites
      #####################################
#geodist for UNG at unguic sites

site_dist
site_dist_ung_ungs <- site_dist %>%
  filter(siteType == "ung_incl")%>%
  filter(host == "unguiculata") %>%
  select("sampleID","x", "y")

site_dist_ung_ungs
ung_ungs_distance_matrix <- geodist(site_dist_ung_ungs, site_dist_ung_ungs, measure = 'geodesic' )/1000
colnames(ung_ungs_distance_matrix) <- site_dist_ung_ungs$sampleID
rownames(ung_ungs_distance_matrix) <- site_dist_ung_ungs$sampleID

ung_ungs_distance_matrix

write.csv(ung_ungs_distance_matrix, "ung_ungsites_cl_up.csv")


# making environmental PCA distance matrix 

meta_div_full


#UNG UNG SITES
site_dist_ung_ungs.clean <- site_dist %>%
  filter(siteType == "ung_incl")%>%
  filter(host == "unguiculata") %>%
  mutate(label = siteName)

ung_ungs_pca_env<- site_dist_ung_ungs.clean %>% 
  select(sampleID, PC1, PC2, PC3, PC4) %>%
  tidyr::drop_na() # drop any rows with NAs, which can't be included in PCA

scores.ung_ungs_epca<-ung_ungs_pca_env

rownames(scores.ung_ungs_epca) <- ung_ungs_pca_env$sampleID

names_ung_ungs <- ung_ungs_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
names_ung_ungs <- ung_ungs_pca_env$sampleID

#Calculating environmental distance from multidimensional PCA axis space
ung_ungs_envdist_multi_pca<-as.data.frame(
  as.matrix(
    dist(
      scores.ung_ungs_epca[,3:ncol(
        scores.ung_ungs_epca)], method="manhattan", diag=T, upper=T
    ), nrow=length(names_ung_ungs), ncol=length(names_ung_ungs)
  )
)

#creating species richness & alpha div matrices for unguiculata at UNG sites

#### environment x diversity questions

### Observed richness

obs_mat_ung_ungs2 <- meta_ung_ungs2 %>%
  select(sampleID, Observed)

#create observed matrix for each species
obs_matdiff_ung_ungs2 <- outer(obs_mat_ung_ungs2$Observed, obs_mat_ung_ungs2$Observed, "-")
colnames(obs_matdiff_ung_ungs2) <- obs_mat_ung_ungs2$sampleID
rownames(obs_matdiff_ung_ungs2) <- obs_mat_ung_ungs2$sampleID
obs_matdist_ung_ungs2 <- as.dist(obs_matdiff_ung_ungs2)
obs_matdist_ung_ungs2 <- abs(obs_matdist_ung_ungs2)
obs_matdist_ung_ungs2

### Shannon diversity

shan_mat_ung_ungs2 <- meta_ung_ungs2 %>%
  select(sampleID, Shannon)

#create observed matrix for each species
shan_matdiff_ung_ungs2 <- outer(shan_mat_ung_ungs2$Shannon, shan_mat_ung_ungs2$Shannon, "-")
colnames(shan_matdiff_ung_ungs2) <- shan_mat_ung_ungs2$sampleID
rownames(shan_matdiff_ung_ungs2) <- shan_mat_ung_ungs2$sampleID
shan_matdist_ung_ungs2 <- as.dist(shan_matdiff_ung_ungs2)
shan_matdist_ung_ungs2 <- abs(shan_matdist_ung_ungs2)
shan_matdist_ung_ungs2

#mantel tests for envd + geod + comd, impt env/geo dist- ung AT UNG sites
#ung_ungs_geodist
ung_ungs_edist <- ung_ungs_envdist_multi_pca


ung_ungs_geodist_man <- as.dist(ung_ungs_distance_matrix)


ung_ungs_edist
#write.csv(ung_ungs_edist, "ung_ungs_pca4_envdist.csv")

#ung_ungs_edist_man <- read.csv("ung_ungs_pca4_envdist.csv", header = TRUE, row.names = 1)
ung_ungs_edist_man <- as.dist(ung_ungs_edist)


##### mantels
ung_ungs_env_geo_mantel <- vegan::mantel(ung_ungs_geodist_man, ung_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
ung_ungs_env_geo_mantel

ps.ung.ungs.bray <- phyloseq::distance(ps.ung, method = "bray")
ps.ung.ungs.jaccard <- phyloseq::distance(ps.ung, method = "jaccard")


plot(ung_ungs_edist_man,ps.ung.ungs.bray,pch=1,cex=0.5,col="black",bty="l")
plot(ung_ungs_edist_man, ung_ungs_geodist_man, pch=1,cex=0.5,col="black",bty="l")
plot(ung_ungs_geodist_man,ps.ung.ungs.bray,pch=1,cex=0.5,col="black",bty="l")

ung_ungs_geodist_man
ung_ungs_edist_man
#env mantel
ung_ungs_env_div_mantel <- vegan::mantel(ps.ung.ungs.bray, ung_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
ung_ungs_env_div_mantel

#geo mantel
ung_ungs_geo_div_mantel <- vegan::mantel(ps.ung.ungs.bray, ung_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
ung_ungs_geo_div_mantel

#partial mantel of env controlling for geo
ung_ungs_env_part <- mantel.partial(ps.ung.ungs.bray, ung_ungs_edist_man, ung_ungs_geodist_man, method="pearson", permutations=9999)
ung_ungs_env_part

#partial mantel of geo controlling for env
ung_ungs_geo_part <- mantel.partial(ps.ung.ungs.bray, ung_ungs_geodist_man, ung_ungs_edist_man, method="pearson", permutations=9999)
ung_ungs_geo_part

##### jaccard versions
#env mantel
ung_ungs_env_div_mantel.j <- vegan::mantel(ps.ung.ungs.jaccard, ung_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
ung_ungs_env_div_mantel.j

#geo mantel
ung_ungs_geo_div_mantel.j <- vegan::mantel(ps.ung.ungs.jaccard, ung_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
ung_ungs_geo_div_mantel.j

#partial mantel of env controlling for geo
ung_ungs_env_part.j <- mantel.partial(ps.ung.ungs.jaccard, ung_ungs_edist_man, ung_ungs_geodist_man, method="pearson", permutations=9999)
ung_ungs_env_part.j

#partial mantel of geo controlling for env
ung_ungs_geo_part.j <- mantel.partial(ps.ung.ungs.jaccard, ung_ungs_geodist_man, ung_ungs_edist_man, method="pearson", permutations=9999)
ung_ungs_geo_part.j

#### Observed: 

#env mantel
ung_ungs_env_mantel_obs <- vegan::mantel(obs_matdist_ung_ungs2, ung_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
ung_ungs_env_mantel_obs

#geo mantel
ung_ungs_geo_mantel_obs <- vegan::mantel(obs_matdist_ung_ungs2, ung_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
ung_ungs_geo_mantel_obs

#partial mantel of env controlling for geo
ung_ungs_env_part_obs <- mantel.partial(obs_matdist_ung_ungs2, ung_ungs_edist_man, ung_ungs_geodist_man, method="pearson", permutations=9999)
ung_ungs_env_part_obs

ung_ungs_geo_part_obs <- mantel.partial(obs_matdist_ung_ungs2, ung_ungs_geodist_man, ung_ungs_edist_man, method="pearson", permutations=9999)
ung_ungs_geo_part_obs

#### Shannon: 

#env mantel
ung_ungs_env_mantel_shan <- vegan::mantel(shan_matdist_ung_ungs2, ung_ungs_edist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
ung_ungs_env_mantel_shan
#*

#geo mantel
ung_ungs_geo_mantel_shan <- vegan::mantel(shan_matdist_ung_ungs2, ung_ungs_geodist_man, method = "spearman", permutations = 9999, na.rm = TRUE)
ung_ungs_geo_mantel_shan
#.

#partial mantel of env controlling for geo
ung_ungs_env_part_shan <- mantel.partial(shan_matdist_ung_ungs2, ung_ungs_edist_man, ung_ungs_geodist_man, method="pearson", permutations=9999)
ung_ungs_env_part_shan
#*

ung_ungs_geo_part_shan <- mantel.partial(shan_matdist_ung_ungs2, ung_ungs_geodist_man, ung_ungs_edist_man, method="pearson", permutations=9999)
ung_ungs_geo_part_shan
#.


###########################################
# Differential abundance analysis
###########################################


library(ANCOMBC)
library(DT)
library(phyloseq)
library(scales)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC", force = TRUE)
      #####################################
      # UMX sites
      #####################################

#ancom-bc2 with UMX sites, including structural zeros- GENUS}

set.seed(123)
# It should be noted that we have set the number of bootstrap samples (B) equal 
# to 10 in the 'trend_control' function for computational expediency. 
# However, it is recommended that users utilize the default value of B, 
# which is 100, or larger values for optimal performance.
output.gen = ancombc2(data = ps.ungsites, assay_name = "counts", tax_level = "Genus",
                      fix_formula = "host + PC1 + PC2", rand_formula = NULL,
                      p_adj_method = "holm", pseudo_sens = TRUE,
                      prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                      group = "host", struc_zero = TRUE, neg_lb = FALSE,
                      alpha = 0.05, n_cl = 2, verbose = TRUE,
                      global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                      iter_control = list(tol = 1e-2, max_iter = 20, 
                                          verbose = TRUE),
                      em_control = list(tol = 1e-5, max_iter = 100),
                      lme_control = lme4::lmerControl(),
                      mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                      trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                  nrow = 2, 
                                                                  byrow = TRUE),
                                                           matrix(c(-1, 0, 1, -1),
                                                                  nrow = 2, 
                                                                  byrow = TRUE),
                                                           matrix(c(1, 0, 1, -1),
                                                                  nrow = 2, 
                                                                  byrow = TRUE)),
                                           node = list(2, 2, 1),
                                           solver = "ECOS",
                                           B = 100))

res_prim_gen = output.gen$res
res_pair_gen = output.gen$res_pair
res_global_gen = output.gen$global

tab_zero = output.gen$zero_ind
tab_zero %>%
  datatable(caption = "The detection of structural zeros")


#output coding
res_prim_gen = output.gen$res
res_pair_gen = output.gen$res_pair
res_global_gen = output.gen$res_global
tab_zero_gen = output.gen$zero_ind
library(dplyr)
#tab_zero_gen %>%
#   datatable(caption = "The detection of structural zeros")


#ungsites primary analysis of covariates- SS filter, genus- PC1


colnames(res_prim_gen)
df_PC1 = res_prim_gen %>%
  dplyr::select(taxon, ends_with("PC1")) 

df_fig_PC1 = df_PC1 %>%
  dplyr::filter(diff_PC1 == 1) %>% 
  dplyr::arrange(desc(lfc_PC1)) %>%
  dplyr::mutate(direct = ifelse(lfc_PC1 > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_PC1 == 1, "purple", "black"))

df_fig_PC1$taxon = factor(df_fig_PC1$taxon, levels = df_fig_PC1$taxon)
df_fig_PC1$direct = factor(df_fig_PC1$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

library(ggplot2)
fig_PC1 = df_fig_PC1 %>%
  ggplot(aes(x = taxon, y = lfc_PC1, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_PC1 - se_PC1, ymax = lfc_PC1 + se_PC1), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = "Endophyte genus", y = "Abundance LFC", 
       title = "Log fold change in abundance along PC1") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  scale_x_discrete(labels = label_wrap(8))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig_PC1$color), legend.position = "bottom")
fig_PC1

fig_PC1_vert = df_fig_PC1 %>%
  ggplot(aes(x = taxon, y = lfc_PC1, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_PC1 - se_PC1, ymax = lfc_PC1 + se_PC1), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = "Endophyte genus", y = "Abundance LFC", 
       title = "Abundance LFC along PC1") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  scale_x_discrete(labels = label_wrap(8))+
  theme_bw() + 
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(), axis.text.y = element_text(hjust = 1, color = df_fig_PC1$color),legend.position = "right")
fig_PC1_vert


#umx sites primary analysis of covariates- SS filter, genus- PC2


colnames(res_prim_gen)
df_PC2 = res_prim_gen %>%
  dplyr::select(taxon, ends_with("PC2")) 

df_fig_PC2 = df_PC2 %>%
  dplyr::filter(diff_PC2 == 1) %>% 
  dplyr::arrange(desc(lfc_PC2)) %>%
  dplyr::mutate(direct = ifelse(lfc_PC2 > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_PC2 == 1, "purple", "black"))

df_fig_PC2$taxon = factor(df_fig_PC2$taxon, levels = df_fig_PC2$taxon)
df_fig_PC2$direct = factor(df_fig_PC2$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

library(ggplot2)
fig_PC2 = df_fig_PC2 %>%
  ggplot(aes(x = taxon, y = lfc_PC2, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_PC2 - se_PC2, ymax = lfc_PC2 + se_PC2), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = "Endophyte genus", y = "Abundance LFC", 
       title = "LFC in abundance along PC2") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  scale_x_discrete(labels = label_wrap(8))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig_PC2$color), legend.position = "bottom")
fig_PC2

fig_PC2_vert = df_fig_PC2 %>%
  ggplot(aes(x = taxon, y = lfc_PC2, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_PC2 - se_PC2, ymax = lfc_PC2 + se_PC2), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = "Endophyte genus", y = "Abundance LFC", 
       title = "Abundance LFC along PC2") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  scale_x_discrete(labels = label_wrap(8))+
  theme_bw() + 
  coord_flip() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(), axis.text.y = element_text(hjust = 1, color = df_fig_PC2$color),legend.position = "right")
fig_PC2_vert


# ss filter figure pairwise genus



colnames(res_pair_gen)
df_host = res_prim_gen %>%
  dplyr::select(taxon, contains("host")) 

df_fig_pair1 = res_pair_gen %>%
  dplyr::filter(diff_hostunguiculata == 1 |
                  diff_hostmarah == 1 | 
                  diff_hostmarah_hostunguiculata == 1) %>%
  dplyr::mutate(lfc1 = ifelse(diff_hostunguiculata == 1, 
                              round(lfc_hostunguiculata, 2), 0),
                lfc2 = ifelse(diff_hostmarah == 1, 
                              round(lfc_hostmarah, 2), 0),
                lfc3 = ifelse(diff_hostmarah_hostunguiculata == 1, 
                              round(lfc_hostmarah_hostunguiculata, 2), 0)) %>%
  tidyr::pivot_longer(cols = lfc1:lfc3, 
                      names_to = "host", values_to = "value") %>%
  dplyr::arrange(taxon)

df_fig_pair2 = res_pair_gen %>%
  dplyr::filter(diff_hostunguiculata == 1 |
                  diff_hostmarah == 1 | 
                  diff_hostmarah_hostunguiculata == 1) %>%
  dplyr::mutate(lfc1 = ifelse(passed_ss_hostunguiculata == 1 & diff_hostunguiculata == 1, 
                              "purple", "black"),
                              lfc2 = ifelse(passed_ss_hostmarah == 1 & diff_hostmarah == 1, 
                                            "purple", "black"),
                                            lfc3 = ifelse(passed_ss_hostmarah_hostunguiculata == 1 & diff_hostmarah_hostunguiculata == 1, 
                                                          "purple", "black")) %>%
                                                            tidyr::pivot_longer(cols = lfc1:lfc3, 
                                                                                names_to = "host", values_to = "color") %>%
  dplyr::arrange(taxon)

df_fig_pair = df_fig_pair1 %>%
  dplyr::left_join(df_fig_pair2, by = c("taxon", "host"))

df_fig_pair$host = recode(df_fig_pair$host, 
                          `lfc1` = "unguiculata - xantiana",
                          `lfc2` = "Marah - xantiana",
                          `lfc3` = "Marah - unguiculata")

df_fig_pair$host = factor(df_fig_pair$host, 
                          levels = c("unguiculata - xantiana",
                                     "Marah - xantiana", 
                                     "Marah - unguiculata"))
library(scales)
lo = floor(min(df_fig_pair$value))
up = ceiling(max(df_fig_pair$value))
mid = (lo + up)/2
fig_pair = df_fig_pair %>%
  ggplot(aes(x = host, y = taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradientn(colors = c("#00bfc4","white","#f8766d"), 
                       values = rescale(c(-4,0,6.5)),
                       guide = "colorbar", limits=c(-4,6.5))+
  geom_text(aes(host, taxon, label = value, color = color), size = 4) +
  scale_color_identity(guide = FALSE) +
  scale_y_discrete(labels = label_wrap(10)) +
  labs(x = "Host pair", y = "Endophyte genus", title = "Abundance LFC between hosts") +   
  theme_minimal()+
  theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 60, hjust = 1))


fig_pair


#putting together figures- Figure 5***}
fig_PC2
fig_PC1
fig_pair

fig_PC2_vert
fig_PC1_vert

library(gridExtra)
library(multcompView)
library(scales)

lay_hor <- rbind(c(1,2),
                 c(1,3))
figure_5_hor <- grid.arrange(fig_pair, fig_PC1, fig_PC2, layout_matrix = lay_hor)

lay_vert <- rbind(c(1,2),
                  c(1,3))
figure_5_vert <- grid.arrange(fig_pair, fig_PC1_vert, fig_PC2_vert, layout_matrix = lay_vert)


ggsave("combined_fig5_hor.png", figure_5_hor, width = 10, height = 8, dpi = 600)
ggsave("combined_fig5_hor.svg", figure_5_hor, width = 10, height = 8, dpi = 600)

ggsave("combined_fig5_vert.png", figure_5_vert, width = 10, height = 8, dpi = 600)
ggsave("combined_fig5_vert.svg", figure_5_vert, width = 10, height = 8, dpi = 600)


      #####################################
      # MX sites
      #####################################

#ancom-bc2 with MX sites including structural zeros- genus

set.seed(123)
# It should be noted that we have set the number of bootstrap samples (B) equal 
# to 10 in the 'trend_control' function for computational expediency. 
# However, it is recommended that users utilize the default value of B, 
# which is 100, or larger values for optimal performance.
output_mx_gen = ancombc2(data = ps.mx, assay_name = "counts", tax_level = "Genus",
                             fix_formula = "host + PC1 + PC2", rand_formula = NULL,
                             p_adj_method = "holm", pseudo_sens = TRUE,
                             prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05,
                             group = "host", struc_zero = TRUE, neg_lb = FALSE,
                             alpha = 0.05, n_cl = 2, verbose = TRUE,
                             global = TRUE, pairwise = TRUE, dunnet = TRUE, trend = TRUE,
                             iter_control = list(tol = 1e-2, max_iter = 20, 
                                                 verbose = TRUE),
                             em_control = list(tol = 1e-5, max_iter = 100),
                             lme_control = lme4::lmerControl(),
                             mdfdr_control = list(fwer_ctrl_method = "holm", B = 100),
                             trend_control = list(contrast = list(matrix(c(1, 0, -1, 1),
                                                                         nrow = 2, 
                                                                         byrow = TRUE),
                                                                  matrix(c(-1, 0, 1, -1),
                                                                         nrow = 2, 
                                                                         byrow = TRUE),
                                                                  matrix(c(1, 0, 1, -1),
                                                                         nrow = 2, 
                                                                         byrow = TRUE)),
                                                  node = list(2, 2, 1),
                                                  solver = "ECOS",
                                                  B = 100))

#output coding
res_prim_xm_gen = output_mx_gen$res
res_pair_xm_gen = output_mx_gen$res_pair

tab_zero_xm_gen = output_mx_gen$zero_ind
tab_zero_xm_gen %>%
  datatable(caption = "The detection of structural zeros")

#XM differential abundance between hosts global family}

res_prim_xm_sig_gen <- res_prim_xm_gen %>%
  dplyr::filter(diff_hostmarah == TRUE)
res_prim_xm_sig_gen

# XM primary analysis of covariates- SS filter, family- PC2}


colnames(res_prim_xm_gen)

df_PC2_xm_gen = res_prim_xm_gen %>%
  dplyr::select(taxon, ends_with("PC2")) 

df_fig_PC2_xm_gen = df_PC2_xm_gen %>%
  dplyr::filter(diff_PC2 == 1) %>% 
  dplyr::arrange(desc(lfc_PC2)) %>%
  dplyr::mutate(direct = ifelse(lfc_PC2 > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_PC2 == 1, "purple", "black"))

df_fig_PC2_xm_gen$taxon = factor(df_fig_PC2_xm_gen$taxon, levels = df_fig_PC2_xm_gen$taxon)
df_fig_PC2_xm_gen$direct = factor(df_fig_PC2_xm_gen$direct, 
                                  levels = c("Positive LFC", "Negative LFC"))

library(ggplot2)

fig_PC2_xm_gen = df_fig_PC2_xm_gen %>%
  ggplot(aes(x = taxon, y = lfc_PC2, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_PC2 - se_PC2, ymax = lfc_PC2 + se_PC2), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "LFC of abundance with increase of 1 for PC2") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  scale_x_discrete(labels = label_wrap(8))+
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig_PC2_xm_gen$color), legend.position = "bottom")

fig_PC2_xm_gen


# mx primary analysis of covariates- SS filter, family- PC1}


colnames(res_prim_xm_gen)

df_PC1_xm_gen = res_prim_xm_gen %>%
  dplyr::select(taxon, ends_with("PC1")) 

df_fig_PC1_xm_gen = df_PC1_xm_gen %>%
  dplyr::filter(diff_PC1 == 1) %>% 
  dplyr::arrange(desc(lfc_PC1)) %>%
  dplyr::mutate(direct = ifelse(lfc_PC1 > 0, "Positive LFC", "Negative LFC"),
                color = ifelse(passed_ss_PC1 == 1, "purple", "black"))

df_fig_PC1_xm_gen$taxon = factor(df_fig_PC1_xm_gen$taxon, levels = df_fig_PC1_xm_gen$taxon)
df_fig_PC1_xm_gen$direct = factor(df_fig_PC1_xm_gen$direct,
                                  levels = c("Positive LFC", "Negative LFC"))

library(ggplot2)
fig_PC1_xm_gen = df_fig_PC1_xm_gen %>%
  ggplot(aes(x = taxon, y = lfc_PC1, fill = direct)) + 
  geom_bar(stat = "identity", width = 0.7, color = "black", 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = lfc_PC1 - se_PC1, ymax = lfc_PC1 + se_PC1), 
                width = 0.2, position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "LFC of abundance with increase of 1 for PC1") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1,
                                   color = df_fig_PC1_xm_gen$color), legend.position = "bottom")
fig_PC1_xm_gen



###########################################
# Analysis of shared & unique taxa
###########################################


# Assuming you have a phyloseq object named 'ps'
# Convert the OTU table from the phyloseq object into a data frame, appending  to each object name
otu_table <- as.data.frame(otu_table(ps.fungi.r))
sample_data <- meta_div_full


# Include the SampleID from the row names in both the OTU table and sample data
otu_table$sampleID <- rownames(otu_table)
#sample_data$sampleID <- rownames(sample_data)

# Melt the OTU table to long format, suitable for merging
otu_long <- pivot_longer(otu_table, cols = -sampleID, names_to = "Taxa", values_to = "Count")
otu_long <- otu_long %>% filter(Count > 0)  # Removing zero counts to simplify the data


# Joining the sample metadata with the OTU data, both appended with 
full_data <- left_join(otu_long, sample_data, by = "sampleID")

# Displaying some of the data to ensure everything is correct
print(head(full_data))


# Ensure host_list is correct and exists
host_list <- unique(full_data$host)

# First, pivot the data to create individual columns for each host species based on counts
full_data_wide <- full_data %>%
  pivot_wider(names_from = host, values_from = Count, values_fill = list(Count = 0))  # Fill missing counts with 0

taxa_summary <- full_data_wide %>%
  group_by(siteName, siteType, siteID, x, y, Taxa) %>%  # Include 'siteID' in the grouping
  summarise(
    Marah = sum(marah, na.rm = TRUE),
    Xantiana = sum(xantiana, na.rm = TRUE),
    Unguiculata = sum(unguiculata, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  mutate(
    Shared_All = if_else(Marah > 0 & Xantiana > 0 & Unguiculata > 0, pmin(Marah, Xantiana, Unguiculata, na.rm = TRUE), 0),
    Shared_Marah_Xantiana = if_else(Marah > 0 & Xantiana > 0 & Unguiculata == 0, pmin(Marah, Xantiana, na.rm = TRUE), 0),
    Shared_Marah_Unguiculata = if_else(Marah > 0 & Unguiculata > 0 & Xantiana == 0, pmin(Marah, Unguiculata, na.rm = TRUE), 0),
    Shared_Xantiana_Unguiculata = if_else(Xantiana > 0 & Unguiculata > 0 & Marah == 0, pmin(Xantiana, Unguiculata, na.rm = TRUE), 0),
    Unique_Marah = if_else(Marah > 0 & Xantiana == 0 & Unguiculata == 0, Marah, 0),
    Unique_Xantiana = if_else(Xantiana > 0 & Marah == 0 & Unguiculata == 0, Xantiana, 0),
    Unique_Unguiculata = if_else(Unguiculata > 0 & Marah == 0 & Xantiana == 0, Unguiculata, 0)
  )

#rename columns
taxa_summary <- taxa_summary %>%
  rename(longitude = x, latitude = y)

# Check for NA values and summarize data to understand the range and content
summary(taxa_summary)

# Assuming 'longitude', 'latitude', 'Unique_Marah', 'Shared_Marah_Xantiana', 'Shared_All' need checking
taxa_summary <- taxa_summary %>%
  filter(!is.na(longitude) & !is.na(latitude) & !is.na(Unique_Marah) & 
           !is.na(Shared_Marah_Xantiana) & !is.na(Shared_All))

colnames(taxa_summary)

# Recalculate summary after filtering to ensure correctness
summary(taxa_summary)



taxa_summary_long <- pivot_longer(taxa_summary, 
                                      cols = c(Unique_Marah, Unique_Xantiana, Unique_Unguiculata,
                                               Shared_Marah_Xantiana, Shared_Marah_Unguiculata,
                                               Shared_Xantiana_Unguiculata, Shared_All),
                                      names_to = "Category", values_to = "Value")


taxa_summary_long_trimmed <- taxa_summary_long %>%
  filter(Value != 0)


taxa_summary_boolean <- taxa_summary %>%
  mutate(Shared_All = ifelse(Shared_All == "0",FALSE,TRUE)) %>%
  mutate(Shared_Marah_Xantiana = ifelse(Shared_Marah_Xantiana == "0",FALSE,TRUE)) %>%
  mutate(Shared_Marah_Unguiculata = ifelse(Shared_Marah_Unguiculata == "0",FALSE,TRUE)) %>%
  mutate(Shared_Xantiana_Unguiculata = ifelse(Shared_Xantiana_Unguiculata == "0",FALSE,TRUE)) %>%
  mutate(Unique_Marah = ifelse(Unique_Marah == "0",FALSE,TRUE)) %>%
  mutate(Unique_Xantiana = ifelse(Unique_Xantiana == "0",FALSE,TRUE)) %>%
  mutate(Unique_Unguiculata = ifelse(Unique_Unguiculata == "0",FALSE,TRUE))

taxa_summary_binary <- taxa_summary %>%
  mutate(Marah = ifelse(Marah == "0",0,1)) %>%
  mutate(Xantiana = ifelse(Xantiana == "0",0,1)) %>%
  mutate(Unguiculata = ifelse(Unguiculata == "0",0,1)) %>%
  mutate(Shared_All = ifelse(Shared_All == "0",0,1)) %>%
  mutate(Shared_Marah_Xantiana = ifelse(Shared_Marah_Xantiana == "0",0,1)) %>%
  mutate(Shared_Marah_Unguiculata = ifelse(Shared_Marah_Unguiculata == "0",0,1)) %>%
  mutate(Shared_Xantiana_Unguiculata = ifelse(Shared_Xantiana_Unguiculata == "0",0,1)) %>%
  mutate(Unique_Marah = ifelse(Unique_Marah == "0",0,1)) %>%
  mutate(Unique_Xantiana = ifelse(Unique_Xantiana == "0",0,1)) %>%
  mutate(Unique_Unguiculata = ifelse(Unique_Unguiculata == "0",0,1))



taxa_summary_sum <- taxa_summary_binary %>%
  group_by(siteID) %>%
  summarise(Shared_All_Sum = sum(Shared_All), Shared_Marah_Xantiana_Sum = sum(Shared_Marah_Xantiana), Shared_Marah_Unguiculata_Sum = sum(Shared_Marah_Unguiculata), Shared_Xantiana_Unguiculata_Sum = sum(Shared_Xantiana_Unguiculata), Unique_Marah_Sum = sum(Unique_Marah), Unique_Xantiana_Sum = sum(Unique_Xantiana),Unique_Unguiculata_Sum = sum(Unique_Unguiculata))


taxa_summary_sum_ungs <- taxa_summary_binary %>%
  filter(siteType == "ung_incl") %>%
  group_by(siteID) %>%
  summarise(Shared_All_Sum = sum(Shared_All), Shared_Marah_Xantiana_Sum = sum(Shared_Marah_Xantiana), Shared_Marah_Unguiculata_Sum = sum(Shared_Marah_Unguiculata), Shared_Xantiana_Unguiculata_Sum = sum(Shared_Xantiana_Unguiculata), Unique_Marah_Sum = sum(Unique_Marah), Unique_Xantiana_Sum = sum(Unique_Xantiana),Unique_Unguiculata_Sum = sum(Unique_Unguiculata))

taxa_summary_sum_grouped <- taxa_summary_sum %>%
  mutate(Shared2 = Shared_Marah_Xantiana_Sum + Shared_Marah_Unguiculata_Sum + Shared_Xantiana_Unguiculata_Sum) %>%
  mutate(UniquesAll = Unique_Marah_Sum + Unique_Xantiana_Sum + Unique_Unguiculata_Sum)


taxa_summary_boolean <- taxa_summary_boolean[, c("Marah", "Xantiana", "Unguiculata", "Shared_All", "Shared_Marah_Xantiana", "Shared_Marah_Unguiculata", "Shared_Xantiana_Unguiculata", "Unique_Marah", "Unique_Xantiana", "Unique_Unguiculata", "siteID", "Taxa","longitude", "latitude")]

ung_sites <- c("corral creek", "site 22", "headquarters", "highway 155", "green rock west", "hobo", "delonegha east", "mill creek trailhead", "democrat", "china gardens", "cow flat", "cattle pens", "mimulus pictus", "upper richbar north", "live oak", "drogon", "good place", "muddy oaks", "site 16")



taxa_summary_sum_long <- read.csv("taxa_summary_sum_long.csv")

taxa_summary_sum_long_ungs <- taxa_summary_sum_long %>%
  filter(siteType == 2)
write.csv(taxa_summary_sum_long_ungs, "ung_sites_occurrence_categories.csv")

taxa_summary_sum_long_ungs$Category <- recode(taxa_summary_sum_long_ungs$Category, 
                                                  'Shared_All_Sum' = 'All Hosts', 
                                                  'Shared_Marah_Unguiculata_Sum' = 'Marah & Unguiculata', 
                                                  'Shared_Marah_Xantiana_Sum' = 'Xantiana & Marah', 
                                                  'Shared_Xantiana_Unguiculata_Sum' = 'Xantiana & Unguiculata', 
                                                  'Unique_Marah_Sum' = 'Marah', 
                                                  'Unique_Unguiculata_Sum' = 'Unguiculata', 
                                                  'Unique_Xantiana_Sum' = 'Xantiana')


library(dplyr)
library(lme4)
setwd("~/C1 Laptop Local/C1 LAPTOP WORKING")
USOC <- read.csv("ung_sites_occurrence_categories.csv")
mx.ungs <- lmer(Value ~ Category + (1 | siteID), data = USOC)
saveRDS(mx.ungs, "mx.ungs.rds")
XMOC <- read.csv("mx_categories.csv")
mxed.mx <- lmer(Value ~ Category + (1 | siteID), data = XMOC)
saveRDS(mxed.mx, "mxed.mx.rds")
library(emmeans)
pairwise_ungs_mixed <- emmeans(model, list(pairwise ~ Category), adjust = "tukey")
saveRDS(pairwise_ungs_mixed, "pairwise_ungs_mixed.rds")

mx.ungs <- readRDS("mx.ungs.rds")
mxed.mx <- readRDS("mxed.xanmar.rds")
pairwise_ungs_mixed <- readRDS("pairwise_ungs_mixed.rds")
pairwise_ungs_mixed
# Load necessary libraries
library(ggplot2)
library(dplyr)
library(emmeans)
library(multcompView)
library(scales)

taxa_summary_sum_long_ungs2 <- taxa_summary_sum_long_ungs %>%
  mutate(Category = case_when(
    Category == "Shared_All_Sum" ~ "All Hosts",
    Category == "Shared_Marah_Unguiculata_Sum" ~ "Marah & C. unguiculata",
    Category == "Shared_Marah_Xantiana_Sum" ~ "C. xantiana & Marah",
    Category == "Shared_Xantiana_Unguiculata_Sum" ~ "C. xantiana & C. unguiculata",
    Category == "Unique_Marah_Sum" ~ "Marah",
    Category == "Unique_Unguiculata_Sum" ~ "C. unguiculata",
    Category == "Unique_Xantiana_Sum" ~ "C. xantiana",
    TRUE ~ Category
  ))
taxa_summary_sum_long_ungs2$Category = factor(taxa_summary_sum_long_ungs2$Category, levels = c("Marah", "C. xantiana", "C. unguiculata", "Marah & C. unguiculata", "C. xantiana & C. unguiculata", "C. xantiana & Marah", "All Hosts"), ordered = TRUE)

taxa_means <- taxa_summary_sum_long_ungs2 %>%
  group_by(Category) %>%
  summarize(Value = mean(Value))

# Create the modified barplot
ung.mxm.barplot <- ggplot(data = taxa_means, 
                          aes(x = Category, y = Value, fill = Category)) +
  geom_col(position = "dodge", width = 0.7, orientation = "x") +
  ylab("Average Number of Taxa Across Sites") +
  xlab("Co-Occurrence Category") +
  scale_fill_manual(values = c("#7CAE00", "#FF61CC", "#00B8E7", "#00C19A", "#C77CFF", "#CD9600", "gray20")) +
  scale_x_discrete(labels = label_wrap(10)) +
  theme(legend.position = "none")

# Display the plot
ung.mxm.barplot

# Display the plot
ung.mxm.barplot

# Display the plot
ung.mxm.barplot

ung.mxm.barplot <- ggplot(data = taxa_summary_sum_long_ungs2, 
                          aes(x=Category, y=Value, fill = Category)) +
  geom_col(position = "dodge", width=0.7, orientation = "x") +
  ylab("Average Number of Taxa") +
  xlab("Occurrence Category")+
  scale_fill_manual(values = c("#7CAE00", "#FF61CC", "#00B8E7", "#00C19A", "#C77CFF", "#CD9600", "gray20"))+
  scale_x_discrete(labels = label_wrap(10)) +
  theme(legend.position = "none")
ung.mxm.barplot