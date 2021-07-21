# Alpha Diversity, DOminance und Eveness Analysis

library(plyr)
library(dplyr)
library(phyloseq)
library(microbiome)
library(here)
library(ggsignif)
library(ggridges)

ps <- readRDS(file=paste0(here::here(), "/data/pspaper.rds"))

# test: prune samples in case you want to remove the chronic patients 
#ps <- prune_samples(!sample_data(ps)$clinic_id %in% c("patient 2", "patient 4", "patient 9", "patient 12"), ps)
md <- meta(ps)

md$FEV1year_group <- factor(md$FEV1year_group, levels = c("stable", "decreasing"))

fillvector <- c("decreasing"= "#E21214", 
                "stable"= "#377DB8")

# ALPHA DIVERSITY 
# calculate the alpha diversity mean for each patient
alpha_md <- md %>% 
  group_by(clinic_id) %>% 
  mutate(alpha_mean = mean(Shannon)) %>% 
  select(clinic_id, alpha_mean, FEV1year_group) %>% 
  distinct()

plot <- ggplot(alpha_md, aes(x = FEV1year_group, y = alpha_mean, fill = FEV1year_group))
plot +
  geom_boxplot()+
  theme_classic()+
  ylab(paste0("Alpha Diversity | Shannon"))+
 # scale_x_discrete(limits=c("stable", "decreasing"))+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = F,
              test = "wilcox.test",
              textsize = 3,
              y_position = 3.3,
              vjust = 0.1)+
  scale_fill_manual(values = fillvector)

ggsave(paste0(here(), "/outputs/AlphaDIversity_mean.pdf"),
       device = "pdf", 
       width = 4, 
       height = 8,
       units = "cm")


# EVENESS 
# calculate the evenesss  mean for each patient
even_md <- md %>% 
  group_by(clinic_id) %>% 
  mutate(even_mean = mean(simpson.Eveness)) %>% 
  select(clinic_id, even_mean, FEV1year_group) %>% 
  distinct()

plot <- ggplot(even_md, aes(x = FEV1year_group, y = even_mean, fill = FEV1year_group))
plot +
  geom_boxplot()+
  theme_classic()+
  ylab(paste0("Eveness | Simpson"))+
  # scale_x_discrete(limits=c("stable", "decreasing"))+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = F,
              test = "wilcox.test",
              textsize = 3,
              y_position = 0.17,
              vjust = 0.1)+
  scale_fill_manual(values = fillvector)

ggsave(paste0(here(), "/outputs/Eveness_mean.pdf"),
       device = "pdf", 
       width = 4, 
       height = 8,
       units = "cm")

# DOMINANCE 
# calculate the dominance  mean for each patient
dom_md <- md %>% 
  group_by(clinic_id) %>% 
  mutate(dom_mean = mean(relative.Dominance)) %>% 
  select(clinic_id, dom_mean, FEV1year_group) %>% 
  distinct()

plot <- ggplot(dom_md, aes(x = FEV1year_group, y = dom_mean, fill = FEV1year_group))
plot +
  geom_boxplot()+
  theme_classic()+
  ylab(paste0("Dominance | relative"))+
  # scale_x_discrete(limits=c("stable", "decreasing"))+
  theme(legend.position = "None")+
  xlab("")+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = F,
              test = "wilcox.test",
              textsize = 3,
              y_position = 0.8,
              vjust = 0.1)+
  scale_fill_manual(values = fillvector)

ggsave(paste0(here(), "/outputs/Dominance_mean.pdf"),
       device = "pdf", 
       width = 4, 
       height = 8,
       units = "cm")


# Bacterial Burdon 
# as boxplot 
burdon_md <- md %>% 
  filter(!is.na(copies.µl)) %>% 
  group_by(clinic_id) %>% 
  mutate(burdon_mean = mean(copies.µl)) %>% 
  select(clinic_id, burdon_mean, FEV1year_group) %>% 
  distinct()


plot <- ggplot(burdon_md, aes(x = FEV1year_group, y = burdon_mean, fill = FEV1year_group))
plot +
  geom_boxplot()+
  theme_classic()+
  ylab(paste0("Burdon [mean / patient] | copies / µl"))+
  # scale_x_discrete(limits=c("stable", "decreasing"))+
  theme(legend.position = "None")+
  xlab("")+
  scale_y_log10()+
  geom_signif(comparison = list(c("decreasing", "stable")), 
              map_signif_level = TRUE,
              test = "wilcox.test",
              textsize = 3,
              y_position = 8.4,
              vjust = 0.1)+
  scale_fill_manual(values = fillvector)

ggsave(paste0(here(), "/outputs/Burdon_mean.svg"),
       device = "svg", 
       width = 4, 
       height = 8,
       units = "cm")

# as burdon over time 
tmp <- md %>% rename("burdon"  = copies.µl) %>% 
  filter(!is.na(burdon)) %>% 
  mutate(Shannon = factor(Shannon))

ggplot(tmp, aes(x=Time_index, y=burdon,  group = clinic_id, color = FEV1year_group))+
  geom_line()+
  theme_classic()+
  theme(legend.position = "bottom")+
  scale_color_manual(values = fillvector)+
  scale_y_log10()

ggsave(paste0(here(), "/outputs/Burdon_overtime.svg"),
       device = "svg", 
       width = 16, 
       height = 8,
       units = "cm")
