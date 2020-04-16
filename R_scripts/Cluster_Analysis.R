library(ggplot2)
library(phyloseq)
library(microbiome)
library(dplyr)
library(cluster)

# load functions 
source("C:~/R_scripts/Analysis_functions.R")

# load the ps object
ps <- readRDS(file= "~/data/pspaper.rds")


# Hierarchical Clustering ####

dist <- distance(ps, "horn")

# making a tree based on the dist 
tree_ward.D <- hclust(dist, method = "ward.D")
tree_ward.D2 <- hclust(dist, method = "ward.D2")
tree_single <- hclust(dist, method = "single")
tree_complete <- hclust(dist, method = "complete")
tree_average <- hclust(dist, method = "average")
tree_mcquitty <- hclust(dist, method = "mcquitty")
tree_median <- hclust(dist, method = "median")
tree_centroid <- hclust(dist, method = "centroid")

# for each tree we calculate the cophenetic distance matrix 
coph_ward.D <- cophenetic(tree_ward.D)
coph_ward.D2 <- cophenetic(tree_ward.D2)
coph_single <- cophenetic(tree_single)
coph_complete <- cophenetic(tree_complete)
coph_average <- cophenetic(tree_average)
coph_mcquitty <- cophenetic(tree_mcquitty)
coph_median <- cophenetic(tree_median)
coph_centroid <- cophenetic(tree_centroid)

# determine the difference of each cophenetic distance matrix to the dist 
cor(dist, coph_ward.D, method = "spearman") # 0.60 
cor(dist, coph_ward.D2, method = "spearman") # 0.62
cor(dist, coph_average, method = "spearman") # 0.64
cor(dist, coph_median, method = "spearman") # 0.27
cor(dist, coph_complete, method = "spearman") # 0.63

# make the clusters with the cutree function
mycluster <- function(x, k) list(cluster=cutree(tree_ward.D, k=k))
# the next step will take a while. 
gshc <- clusGap(t(otu_table(ps)), FUN = mycluster, K.max = 70, B = 100)
# plot the gshc 
plot_clusgap(gshc)
# saves the cluster information to each sample name in a dataframe
cl.split <- as.data.frame(cutree(tree_ward.D, k = 24))

# Assign the cluster information to the samples
rownames(cluster) <- cluster$X
names(cluster)[names(cluster) == "cutree.tree_ward.D..k...24."] <- "Hclust"

# add the cluster information to the metadata 
new_metadata <- meta(ps)
new_metadata <- left_join(new_metadata, cluster, by = c("Idf..Nummer" = "X"))

new_metadata$Hclust<- as.numeric(new_metadata$Hclust)
rownames(new_metadata) <- new_metadata$Idf..Nummer # without having the right rownames, the merge will introduce NA for the Hclust and Kclust column. 

# join the metadata back to the phyloseq object 
ps_cluster <- merge_phyloseq(ps, sample_data(new_metadata))


# make otu table, abundant >5% 
#we decided to plot all rsv which occure in at least one sample above 5% abundance 

ps_top <- Taxa_abovepercent_and_others(ps_cluster, 0.05)

# extract the otu table 
otutable <- psmelt(otu_table(ps_top, taxa_are_rows = TRUE))
# merge with the taxtable_top to assign the taxonomy to the otu-name 
taxtable <- tax_table(ps_cluster) %>% 
  data.frame(stringsAsFactors = FALSE)
taxtable$seq <- rownames(taxtable)
otutable26 <- dplyr::left_join(otutable, taxtable, by = c("OTU" = "seq"))
# merge with the cluster information 
otutable26C <- dplyr::left_join(otutable26, cluster, by = c("Sample" = "X"))
#rearrange the otutable for the plotting in cluster 
otutable26C <- plyr::mutate(otutable26C,Sample = forcats::fct_reorder(Sample, as.numeric(Hclust), fun=median))
levels <- levels(as.factor(otutable26C$Genus))
levels[length(levels) + 1] <- "Others"
otutable26C$Genus <- factor(otutable26C$Genus, levels = levels)
otutable26C$Genus[is.na(otutable26C$Genus)] <- "Others"

# reorder the df 
otutable26C <- otutable26C[order(otutable26C$Genus),]


# add clinic_id to the cluster table (for labeling of the bar plot ) 
metadata <- meta(ps)
sp <- dplyr::select(metadata, Idf..Nummer, clinic_id)

cluster <- dplyr::left_join(cluster, sp, by = c("X" = "Idf..Nummer"))

### plot the Bar plot  each samples, arranged by Clusters 

# define colors
library(RColorBrewer)

genusList <- as.data.frame(unique(otutable26C$Genus))
genus <- genusList[[1]]
colorsPalette <-  c(1:36)
names(colorsPalette) <- genus
 colorsPalette["Others"] <- "gray"
    colorsPalette["Actinomyces"] <- "turquoise4"
    colorsPalette["Atopobium"] <- "steelblue"
    colorsPalette["Rothia"] <- "mediumblue"
    colorsPalette["Alloprevotella"] <- "green3"
    colorsPalette["Bergeyella"] <- "olivedrab"
    colorsPalette["Capnocytophaga"] <- "forestgreen"
    colorsPalette["F0058"] <- "greenyellow"
    colorsPalette["Porphyromonas"] <- "mediumspringgreen"
    colorsPalette["Prevotella"] <- "mediumaquamarine"
    colorsPalette["Campylobacter"] <- "lavender"
    colorsPalette["Abiotrophia"] <- "lightgoldenrod1"
    colorsPalette["Gemella"] <-"chocolate1"
    colorsPalette["Granulicatella"] <- "darkgoldenrod1"
    colorsPalette["Megasphaera"] <- "lightsalmon"
    colorsPalette["Oribacterium"] <- "burlywood1"
    colorsPalette["Parvimonas"] <- "gold2"
    colorsPalette["Staphylococcus"] <- "yellow"
    colorsPalette["Streptococcus"] <- "orange2"
    colorsPalette["Veillonella"] <- "yellow3"
    colorsPalette["Fusobacterium"] <- "turquoise3"
    colorsPalette["Leptotrichia"] <- "deepskyblue"
    colorsPalette["Candidatus_Saccharimonas"] <- "lavenderblush"
    colorsPalette["Achromobacter"] <- "coral1"
    colorsPalette["Acinetobacter"] <- "mediumpurple1"
    colorsPalette["Actinobacillus"] <- "palevioletred"
    colorsPalette["Aggregatibacter"] <- "brown4"
    colorsPalette["Bordetella"] <- "hotpink"
    colorsPalette["Eikenella"] <- "indianred1"
    colorsPalette["Escherichia/Shigella"] <-"lightpink"
    colorsPalette["Haemophilus"] <- "mediumorchid3"
    colorsPalette["Kingella"] <- "maroon3"
    colorsPalette["Lautropia"] <- "salmon3"
    colorsPalette["Moraxella"] <- "plum3"
    colorsPalette["Neisseria"] <- "rosybrown3"
    colorsPalette["Pseudomonas"] <- "red2"
    colorsPalette["Stenotrophomonas"] <- "violetred2"
    colorsPalette["Treponema"] <- "lightcyan"
    colorsPalette["Mycoplasma"] <- "slategray1"

    
# rearrange genus 
otutable26C$Genus <- factor(otutable26C$Genus, levels = c("Actinomyces", "Rothia", "Alloprevotella", "Bergeyella", "Capnocytophaga", "F0058", "Porphyromonas", "Prevotella", "Campylobacter", "Abiotrophia", "Gemella", "Granulicatella", "Megasphaera", "Oribacterium", "Parvimonas", "Staphylococcus", "Streptococcus", "Veillonella", "Fusobacterium", "Leptotrichia", "Candidatus_Saccharimonas", "Achromobacter", "Acinetobacter", "Actinobacillus", "Aggregatibacter", "Bordetella", "Eikenella", "Escherichia/Shigella", "Haemophilus", "Kingella", "Lautropia", "Moraxella", "Neisseria", "Pseudomonas", "Stenotrophomonas", "Treponema", "Mycoplasma", "Others"))


# make distance matrix, calculate the tree and save it as dendogran
dend <- distance(ps, "horn") %>% 
        hclust(method = "ward.D") %>% 
        as.dendrogram
# plot the dendogram
dend %>% plot

label_order <- dend %>% labels

# rearrange the sample order in the otutable26C
otutable26C$Sample <- factor(otutable26C$Sample, levels = label_order)

#extract the metadata to order the metadata sample accordingly to the cluster order 
metadata <- microbiome::meta(ps)
metadata$Idf..Nummer <- factor(metadata$Idf..Nummer, levels = label_order)
#merge with the clusternumber (df cluster was produced above)
metadata <- dplyr::left_join(metadata, cluster , by = c("Idf..Nummer"= "X"))

barplot <- ggplot(otutable26C,
                  aes(x = Sample, y = Abundance, fill = Genus, reorder(arrange(otutable26C, Hclust)$Hclust)))
barplot + 
  geom_bar(stat = "identity")+
  my_theme_ppt()+
  scale_fill_manual(values= colorsPalette, name = names(colorsPalette))+
 scale_x_discrete(breaks = metadata$Idf..Nummer, 
                   labels = metadata$Hclust)


