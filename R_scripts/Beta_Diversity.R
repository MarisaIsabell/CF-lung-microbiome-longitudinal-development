library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)

# load functions 
source("~/R_scripts/Analysis_functions.R")

## Beta Diversity Boxplot for Psa clusters and Commensals

calculate the beta-diversity in patients that blong to the cluster Psae compare to the cluster with commensals 

ps_cluster <- readRDS("~/data/ps_cluster.rds")
metadata <- meta(ps_cluster)

# make the distance matrix 
dist <- phyloseq::distance(ps_cluster, method = "horn", type = "sample")
# melt the matrix to a df
library(reshape2)
library(dplyr)
dist_df <- melt(as.matrix(dist), varnames = c("sample1", "sample2"))

# select only the data we need from metadata 
patients <- metadata[, c( "Idf..Nummer" , "clinic_id", "Hclust", "Shannon", "beta.mean.horn")]
# assign "Pseudomonas" to cluster 3,20,11 and commensals to the other clusters
patients$cluster <- "Commensals"
pse_cluster <-c(3,20,11)
patients$cluster[patients$Hclust %in% pse_cluster] <- "Pseudomonas"

# assign the patient and cluster information to each comparison 
dist_df <- left_join(dist_df, patients, by = c("sample1" = "Idf..Nummer"))
colnames(dist_df)[colnames(dist_df) == "clinic_id"] <- "patient1"
dist_df <- left_join(dist_df, patients, by = c("sample2" = "Idf..Nummer"))
colnames(dist_df)[colnames(dist_df) == "clinic_id"] <- "patient2"

# only keep comparisons within the sampe cluster 
dist_df$keep <- apply(dist_df, 1, function(x)
  {if (x["cluster.x"] == x["cluster.y"]){
    dist_df$keep <- "yes"}
  else{
    dist_df$keep <- "no"
  }
  })
dist_df <- dplyr::filter(dist_df, keep == "yes" )


# visualize the boxplot 
library(ggsignif)
plot <- ggplot(dist_df, aes(x = cluster.x, y = value, fill = cluster.y))
plot+
  geom_boxplot()+
  ylab("beta Diversity | Horn")+
  xlab("")+
  scale_x_discrete(limits=c("Pseudomonas", "Commensals"))+
  theme(legend.position = "None")+
  geom_signif(comparison = list(c("Pseudomonas", "Commensals")), 
              map_signif_level = TRUE,
              test = "wilcox.test",
              textsize = 4,
              y_position = 1.05,
              vjust = 0.5)
