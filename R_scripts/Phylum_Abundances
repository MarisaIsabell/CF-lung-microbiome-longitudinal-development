library(plyr)
library(dplyr)
library(phyloseq)
library(microbiome)

ps <- readRDS(file= "~/data/pspaper.rds")
# load functions 
source("~/R_scripts/Analysis_functions.R")

# glom to phylum level

ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_phylum <- microbiome::transform(ps_phylum, "compositional")

phyl <- "Actinobacteria" # change this Phylum accordingly to the Phylum you want to plot
phylum_sub <- subset_taxa(ps_phylum, Phylum == phyl)

otu_ph <- as.data.frame(t(otu_table(phylum_sub)))
otu_ph$sample <- rownames(otu_ph)
colnames(otu_ph) <- c("Phylum", "sample" )

meta_ph <- meta(phylum_sub)
meta_ph <- left_join(meta_ph, otu_ph, by=c("Idf..Nummer"="sample"))


fillvector <- c("decreasing"= "orangered", 
                "stable"= "yellowgreen")

library(ggsignif)

plot <- ggplot(meta_ph, aes(x = FEV1year_group, y = Phylum, fill = FEV1year_group))
plot +
  geom_boxplot()+
  ylab(paste0("rel. Abundance | ", phyl))+
  xlab("FEV1 Slope/year group")+
  scale_x_discrete(limits=c("stable", "decreasing"))+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = TRUE,
              test = "wilcox.test",
              textsize = 4,
              y_position = 1.05,
              vjust = 0.5)+
  scale_fill_manual(values = fillvector)
