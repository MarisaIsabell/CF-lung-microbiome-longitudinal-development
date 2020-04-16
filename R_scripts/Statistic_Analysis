library(phyloseq)
library(microbiome)
library(vegan)
library(dplyr)

# load data
ps <- readRDS(file= "~/data/pspaper.rds")
# load functions 
source("~/R_scripts/Analysis_functions.R")


# PERMANOVA 
set.seed(1000)
# calculate bray curtis distance matrix 
bray <- phyloseq::distance(ps, method = "bray")

metadata <- microbiome::meta(ps)

# load antibiotic information 
anti <- read.csv(file ="C:/Users/Marisa/Google Drive/Masterarbeit/Clinical_Data/Antibiotic_time_period.csv")

anti <- dplyr::select(anti, Pseudonym, percentage)
metadata <- dplyr::left_join(metadata, anti, by =c("clinic_id" = "Pseudonym"))


# Adonis test: Permutational Multivatiate Analysis of Variance Using Distance Matrices
# based on ANderson 2001 
vegan::adonis(bray ~  clinic_id + mutation + Geschlecht + Age, data = metadata, by="terms")
vegan::adonis2(bray ~ clinic_id + mutation + Geschlecht + Age, data = metadata, by="terms")

