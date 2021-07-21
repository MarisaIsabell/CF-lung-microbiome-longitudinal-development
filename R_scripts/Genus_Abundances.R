library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggsignif)
library(dplyr)
library(here)

ps <- readRDS(file=paste0(here(), "/data/pspaper.rds"))

  
# Prevotella mean ####
genus <- "Prevotella"
ps_R <- microbiome::transform(ps, transform = "compositional")
ps_R <- tax_glom(ps_R, taxrank = "Genus")

microbe <- subset_taxa(ps_R, Genus == genus)

otu <- as.data.frame(t(otu_table(microbe)))
otu$sample <- rownames(otu)
colnames(otu) <- c("genus", "sample")

metad <- meta(microbe)
metad <- dplyr::left_join(metad, otu, by=c("Idf..Nummer"="sample"))

metadata <- metad %>%
  group_by(clinic_id) %>% 
  mutate(microbemean = mean(genus))

md_plot <- metadata %>% 
  select(clinic_id, FEV1year_group, microbemean) %>% distinct()

## boxplots
fillvector <- c("decreasing"= "#E21214", 
                "stable"= "#377DB8")

plot <- ggplot(md_plot, aes(x = FEV1year_group, y = microbemean, fill = FEV1year_group))
plot +
  geom_boxplot()+
  theme_classic()+
  ylab("rel. Abundance Prevotella")+
  scale_x_discrete(limits=c("stable", "decreasing"))+
 # scale_y_log10()+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = F,
              test = "wilcox.test",
              textsize = 4,
              y_position = ,
              vjust = 0.5)+
  scale_fill_manual(values = fillvector)
  
ggsave(paste0(here(), "/outputs/PrevotellaBoxplot_patientMean.pdf"),
       device = "pdf", 
       width = 4, 
       height = 8,
       units = "cm")
  
  
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

md_veill <- meta_veill %>% 
  select(clinic_id, FEV1year_group, Veillonella_mean) %>% distinct()

## boxplots
fillvector <- c("decreasing"= "#E21214", 
                "stable"= "#377DB8")

plot <- ggplot(md_veill, aes(x = FEV1year_group, y = Veillonella_mean, fill = FEV1year_group))
plot +
  geom_boxplot()+
  theme_classic()+
  ylab("rel. Abundance Veillonella")+
  scale_x_discrete(limits=c("stable", "decreasing"))+
  #scale_y_log10()+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = F,
              test = "wilcox.test",
              textsize = 3,
              y_position = 0.15,
              vjust = 0.1)+
  scale_fill_manual(values = fillvector)
 
ggsave(paste0(here(), "/outputs/VeillonellaBoxplot_patientMean.pdf"),
       device = "pdf", 
       width = 4, 
       height = 8,
       units = "cm")


# #test: as subject wise boxplot 
# plot <- ggplot(meta_veill, aes(x = FEV1year_group, y = meta_veill$RSV3038, group = clinic_id, fill = FEV1year_group))
# plot +
#   geom_boxplot()+
#   theme_classic()+
#   ylab("rel. Abundance Veillonella")+
#   scale_x_discrete(limits=c("stable", "decreasing"))+
#   #scale_y_log10()+
#   theme(legend.position = "None")+
#   xlab("")+
#   geom_signif(comparison = list(c("decreasing", "stable")), 
#               map_signif_level = F,
#               test = "wilcox.test",
#               textsize = 3,
#               y_position = 0.5,
#               vjust = 0.1)+
#   scale_fill_manual(values = fillvector)
  
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
  mutate(microbemean = mean(genus)) %>% 
  select(clinic_id, FEV1year_group, microbemean) %>% distinct()

## boxplots
fillvector <- c("decreasing"= "#E21214", 
                "stable"= "#377DB8")

plot <- ggplot(metadata, aes(x = FEV1year_group, y = microbemean, fill = FEV1year_group))
plot +
  geom_boxplot()+
  theme_classic()+
  ylab("rel. Abundance Pseudomonas")+
  scale_x_discrete(limits=c("stable", "decreasing"))+
 # scale_y_log10()+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = T,
              test = "wilcox.test",
              textsize = 3,
              y_position = 0.8,
              vjust = 0.1)+
  scale_fill_manual(values = fillvector)

ggsave(paste0(here(), "/outputs/PseudomonasBoxplot_patientMean.pdf"),
       device = "pdf", 
       width = 4, 
       height = 8,
       units = "cm")
#


# Staphylococcus mean  ####
genus <- "Staphylococcus"
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
fillvector <- c("decreasing"= "#E21214", 
                "stable"= "#377DB8")

plot <- ggplot(metadata, aes(x = FEV1year_group, y = genus, fill = FEV1year_group))
plot +
  geom_boxplot()+
  theme_classic()+
  ylab("rel. Abundance Staphylococcus")+
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

ggsave(paste0(here(), "/outputs/StaphylococcusBoxplot.pdf"),
       device = "pdf", 
       width = 4, 
       height = 8,
       units = "cm")

