library(phyloseq)
library(microbiome)
library(vegan)
library(dplyr)
library(here)

# load data
ps <- readRDS(file=paste0(here(),"/data/pspaper.rds"))
# load functions 
source(paste0(here(),"/R_scripts/Analysis_functions.R"))


# PERMANOVA 
set.seed(1000)
# calculate bray curtis distance matrix 
bray <- phyloseq::distance(ps, method = "bray")
bray <- as.matrix(bray)

metadata <- microbiome::meta(ps)

# load antibiotic information 
# anti <- read.csv(file =paste0(here(),"/data/Antibiotika_table.csv"))
# 
# anti <- anti %>% 
#   mutate(reason = If.any..reason.for.iv.antibiotic.therapy) %>% 
#   dplyr::select(Pseudonym, reason) 
#metadata <- dplyr::left_join(metadata, anti, by =c("clinic_id" = "Pseudonym"))


# Adonis test: Permutational Multivatiate Analysis of Variance Using Distance Matrices
# based on ANderson 2001 
vegan::adonis(bray ~  clinic_id + mutation + Geschlecht + Age, data = metadata, by="terms")
vegan::adonis2(bray ~ clinic_id + mutation + Geschlecht + Age, data = metadata, by="terms")
vegan::adonis2(bray ~ clinic_id + reason,  data =metadata, by = "terms")

# for exacerbations 
# merge the routine visits to one level for the factor, so it is treated as identical in the PERMANOVA.
metadata <- metadata %>% 
  dplyr::mutate(Reason.for.visit = case_when(.$Reason.for.visit == "routine with yearly check-up" ~ "routine", 
                                      .$Reason.for.visit == "routine without yearly check-up" ~ "routine",
                                      .$Reason.for.visit == "unclear" ~ as.character(NA),
                                      TRUE ~ Reason.for.visit)) %>% 
  filter(!is.na(Reason.for.visit))

bray <- bray[rownames(metadata), rownames(metadata)] 
tmp <- vegan::adonis2(bray ~ clinic_id + mutation + Geschlecht + Age + Reason.for.visit, data= metadata, by = "margin")

knitr::kable(tmp)
