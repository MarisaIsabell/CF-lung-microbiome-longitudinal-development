library(phyloseq)
library(microbiome)
library(ggplot2)
library(dplyr)

ps <- readRDS(file= "~/data/pspaper.rds")

# Prevotella mean ####
genus <- "Prevotella"
ps_R <- microbiome::transform(ps, transform = "compositional")
ps_R <- tax_glom(ps_R, taxrank = "Genus")

microbe <- subset_taxa(ps_R, Genus == genus)

otu <- as.data.frame(t(otu_table(microbe)))
otu$sample <- rownames(otu)
colnames(otu) <- c("genus", "sample")

metad <- meta(microbe)
metad <- left_join(metad, otu, by=c("Idf..Nummer"="sample"))

metadata <- metad %>%
  group_by(clinic_id) %>% 
  mutate(microbemean = mean(genus))

## boxplots
fillvector <- c("decreasing"= "orangered", 
                "stable"= "yellowgreen")

plot <- ggplot(metadata, aes(x = FEV1year_group, y = genus, fill = FEV1year_group))
plot +
  geom_boxplot()+
  my_theme_ppt()+
  ylab("rel. Abundance Prevotella")+
  scale_x_discrete(limits=c("stable", "decreasing"))+
  scale_y_log10()+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = TRUE,
              test = "wilcox.test",
              textsize = 4,
              y_position = 3,
              vjust = 0.5)+
  scale_fill_manual(values = fillvector)
  
  
  
# Veillonella mean ####
ps_R <- microbiome::transform(ps, transform = "compositional")
ps_R <- tax_glom(ps_R, taxrank = "Genus")

veillonella <- subset_taxa(ps_R, Genus == "Veillonella")

otu_veill <- as.data.frame(t(otu_table(veillonella)))
otu_veill$sample <- rownames(otu_veill)

meta_veill <- meta(veillonella)
meta_veill <- left_join(meta_veill, otu_veill, by=c("Idf..Nummer"="sample"))

meta_veill <- meta_veill %>%
  group_by(clinic_id) %>% 
  mutate(Veillonella_mean = mean(RSV3038))


## boxplots
fillvector <- c("decreasing"= "orangered", 
                "stable"= "yellowgreen")

plot <- ggplot(meta_veill, aes(x = FEV1year_group, y = RSV3038, fill = FEV1year_group))
plot +
  geom_boxplot()+
  my_theme_ppt()+
  ylab("rel. Abundance Veillonella")+
  scale_x_discrete(limits=c("stable", "decreasing"))+
  scale_y_log10()+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = TRUE,
              test = "wilcox.test",
              textsize = 4,
              y_position = 3,
              vjust = 0.5)+
  scale_fill_manual(values = fillvector)
  
  
# Pseudomonas mean  ####
genus <- "Pseudomonas"
ps_R <- microbiome::transform(ps, transform = "compositional")
ps_R <- tax_glom(ps_R, taxrank = "Genus")

microbe <- subset_taxa(ps_R, Genus == genus)

otu <- as.data.frame(t(otu_table(microbe)))
otu$sample <- rownames(otu)
colnames(otu) <- c("genus", "sample")

metad <- meta(microbe)
metad <- left_join(metad, otu, by=c("Idf..Nummer"="sample"))

metadata <- metad %>%
  group_by(clinic_id) %>% 
  mutate(microbemean = mean(genus))


## boxplots
fillvector <- c("decreasing"= "orangered", 
                "stable"= "yellowgreen")

plot <- ggplot(metadata, aes(x = FEV1year_group, y = genus, fill = FEV1year_group))
plot +
  geom_boxplot()+
  my_theme_ppt()+
  ylab("rel. Abundance Pseudomonas")+
  scale_x_discrete(limits=c("stable", "decreasing"))+
  scale_y_log10()+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = TRUE,
              test = "wilcox.test",
              textsize = 4,
              y_position = 3,
              vjust = 0.5)+
  scale_fill_manual(values = fillvector)
