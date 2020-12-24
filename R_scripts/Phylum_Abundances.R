library(plyr)
library(dplyr)
library(phyloseq)
library(microbiome)
library(here)

ps <- readRDS(file=paste0(here(), "/data/pspaper.rds"))
# load functions 
source("~/R_scripts/Analysis_functions.R")

# glom to phylum level

ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_phylum <- microbiome::transform(ps_phylum, "compositional")

phyl <- "Proteobacteria" # change this Phylum accordingly to the Phylum you want to plot
phylum_sub <- subset_taxa(ps_phylum, Phylum == phyl)

otu_ph <- as.data.frame(t(otu_table(phylum_sub)))
otu_ph$sample <- rownames(otu_ph)
colnames(otu_ph) <- c("Phylum", "sample" )

meta_ph <- meta(phylum_sub)
meta_ph <- left_join(meta_ph, otu_ph, by=c("Idf..Nummer"="sample"))

# removing samples with chronically infected Pseudomonas
#meta_ph <- meta_ph %>% 
#  filter(!clinic_id %in% c("RBB01b0197", "RBB01a0025", "RBB01a0013", "RBB01a0084"))


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
