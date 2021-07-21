library(plyr)
library(dplyr)
library(phyloseq)
library(microbiome)
library(here)
library(ggsignif)
library(ggplot2)
library(patchwork)

ps <- readRDS(file=paste0(here(), "/data/pspaper.rds"))
# test: prune samples in case you want to remove the chronic patients 
#ps <- prune_samples(!sample_data(ps)$clinic_id %in% c("patient 2", "patient 4", "patient 9", "patient 12"), ps)
# load functions 
source(paste0(here(), "/R_scripts/Analysis_functions.R"))

# glom to phylum level

ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_phylum <- microbiome::transform(ps_phylum, "compositional")

# select for the specific Phylum 
md.list <- list()
phyl <- c("Proteobacteria", "Firmicutes", "Fusobacteria", "Actinobacteria", "Bacteroidetes")  # change this Phylum accordingly to the Phylum you want to plot
for(p in phyl){
phylum_sub <- subset_taxa(ps_phylum, Phylum == p)

otu_ph <- as.data.frame(t(otu_table(phylum_sub)))
otu_ph$sample <- rownames(otu_ph)
colnames(otu_ph) <- c("Phylum", "sample" )

meta_ph <- meta(phylum_sub)
meta_ph <- left_join(meta_ph, otu_ph, by=c("Idf..Nummer"="sample"))

# calculate the mean relative abundance for this phylum for each patient 
tmp <- meta_ph %>% 
  group_by(clinic_id) %>% 
  mutate(phyl_mean = mean(Phylum)) %>% 
  select(clinic_id, phyl_mean, FEV1year_group) %>% 
  mutate(Phylum = p) %>% 
  distinct()
md.list[[p]] <- tmp
}

md.total <- do.call("rbind", md.list) %>% 
  mutate(Phylum = factor(Phylum, levels = phyl))
# Plot 
  
fillvector <- c("decreasing"= "#E21214", 
                "stable"= "#377DB8")

plot <- ggplot(md.total, aes(x = FEV1year_group, y = phyl_mean, fill = FEV1year_group))+
  geom_boxplot()+
  theme_classic()+
  ylab(paste0("rel. Abundance [mean / patient]"))+
  xlab("FEV1 Slope/year group")+
  scale_x_discrete(limits=c("stable", "decreasing"))+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = F,
              test = "wilcox.test",
              textsize = 2.5,
              y_position = 0.82,
              vjust = 0.1)+
  scale_fill_manual(values = fillvector)+
  facet_wrap(~ Phylum, nrow=1, strip.position = "bottom")
plot




ggsave(paste0(here(), "/outputs/Phylum_allPhyla_patientmean_sigvalues.pdf"),
       device = "pdf", 
       width = 14, 
       height = 8,
       units = "cm")
