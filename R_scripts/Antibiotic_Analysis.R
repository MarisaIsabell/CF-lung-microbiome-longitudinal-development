library(phyloseq)
library(dplyr)
library(here)


Clean_range <- read.csv(file=paste0(here(), "/data/Antibiotika_table.csv"))
FEV <- read.csv(file=paste0(here(), "/data/Regression_FEV1_year.csv"))
ps <- readRDS(file=paste0(here(), "/data/pspaper.rds"))
source(paste0(here(), "/R_scripts/Analysis_Functions.R"))
# Number of Antibiotic Therapy 
clean_nb <- select(Clean_range, Pseudonym, begin_index, If.any..reason.for.iv.antibiotic.therapy)
clean_nb <- distinct(clean_nb)

nb_iv_antibiotic <- as.data.frame(table(clean_nb$Pseudonym))
colnames(nb_iv_antibiotic) <- c("Pseudonym", "number_iv_antibiotic")

# number of exacerbations/elevtive antibiotic usage 
nbanti_FEV <- left_join(nb_iv_antibiotic, select(FEV, Pseudonym, Slope, Group), by = "Pseudonym")

# calculate regression 
reg <- calculate_regression(df = nbanti_FEV, y = nbanti_FEV$Slope,
                            x = nbanti_FEV$number_iv_antibiotic)
nbanti_FEV_wooutlier <- filter(nbanti_FEV, !Pseudonym == "RBB01a0064" )
reg_out <- lm(nbanti_FEV_wooutlier$Slope ~ nbanti_FEV_wooutlier$number_iv_antibiotic)
reg_out_f <- calculate_regression(df = nbanti_FEV_wooutlier, 
                                  y = nbanti_FEV_wooutlier$Slope, 
                                  x = nbanti_FEV_wooutlier$number_iv_antibiotic)

plot <- ggplot(nbanti_FEV, aes(x = number_iv_antibiotic, y = Slope, color = Pseudonym))
plot+
  scale_color_manual(values = colors_clinic_id)+
  geom_smooth(method = "lm", color = "black", size = 1, se=FALSE)+
  geom_abline(slope = reg_out$coefficients[[2]], intercept = reg_out$coefficients[[1]], color = "gray", linetype = "dashed")+
   geom_point(size =4 )+
  xlab("# iv Antibiotic therapy [Quantity]")+
  ylab("Slope FEV1pred Regression")+
  theme(legend.position = "none")+
  annotate("text", x = 10.8, y = -0.0072,
           label = paste0("Slope= ", round(reg$Slope, digits=6), 
               ", R²= ",round(reg$R.squared, digits=2), 
               ", p= ", round(reg$P.value, digits =4)),
           size = 3)+
  annotate("text", x = 10.5, y = -0.0078, 
           label = paste0("Slope= ", round(reg_out_f$Slope, digits = 6), 
                          ", R²= ", round(reg_out_f$R.squared, digits = 2), 
                          ", p= ", round(reg_out_f$P.value, digits =4)), 
            size = 3, color = "gray")

    # boxplot between the two patient group
plot <- ggplot(nbanti_FEV, aes(x = Group, y = number_iv_antibiotic, fill = Group))
plot+
  geom_boxplot()+
  geom_jitter(alpha = 0.5, width = 0.1)+
  xlab("")+
  ylab("# iv Antibiotic therapy [Quantity]")+
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_blank())+
  ggsignif::geom_signif(comparisons = list(c("decliner", "stable")),
                        test = "t.test", 
                        map_signif_level= TRUE)

#boxplot with number of antibiotics per year
#calculate the duration for each patient (total number of day on iv antobiotica)
Clean_range$duration <- Clean_range$end_index -  Clean_range$begin_index
Duration <- select(Clean_range, Pseudonym, begin_index, end_index, duration) %>%
  #since the time index are in months, i have to multiply the month with the mean number of days per month (29.53)
  mutate(duration = duration * 29.53) %>% 
  distinct()
# for some start days we have multiple different end days, we need the longest timeperiod, the other row can be discarded
remove_sec_duration <- function(df){
  df <- df[order(df$begin_index, -abs(df$duration)),] # sort first
  df <- df[!duplicated(df$begin_index), ] # keep highest
  return(df)
}
longer <- Duration %>%
  group_by(Pseudonym) %>%
  do(., remove_sec_duration(.))


# sup up the duration for each patient
antidays <- aggregate(longer$duration, by = list(Pseudonym = longer$Pseudonym), FUN = sum)
colnames(antidays) <- c("Pseudonym", "anti.days")


# Total number of days for each patient in the study

metadata <- microbiome::meta(ps)
totaldays <- select(metadata, clinic_id, Eingangsdatum)
totaldays <- totaldays %>% group_by(clinic_id) %>% mutate(period = (difftime(max(Eingangsdatum), min(Eingangsdatum))))
period <- select(totaldays, clinic_id, period)
period$period <- as.integer(period$period)

# merge information 
percent <- left_join(antidays, period, by = c("Pseudonym"  = "clinic_id")) %>% 
  select(Pseudonym, period, anti.days) %>% 
    distinct()

# calculate the percentage
percent$percentage <- (percent$anti.days/ percent$period)*100

# how many days per year was the patient in the mean on in antibiotics?
perYear = percent %>%
  mutate(NbYear = period/365,
         anti.year = anti.days/NbYear)


antiPerYear <- left_join(perYear, select(FEV, Pseudonym, Slope, Group), by = "Pseudonym")# %>% 

antiPerYear$Group <- factor(antiPerYear$Group, levels = c("stable", "decliner"))
plot <- ggplot(antiPerYear, aes(x = Group, y = anti.year, fill = Group))
plot+
  geom_boxplot()+
  theme_classic()+
  geom_jitter(alpha = 0.5, width = 0.1)+
  xlab("")+
  ylab("iv Antibiotic therapy [days/year]")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_blank())+
  ggsignif::geom_signif(comparisons = list(c("decliner", "stable")),
                        test = "wilcox.test",
                        map_signif_level= FALSE)


ggsave("Antibiotic_Boxplot_antibioticdaysPerYear.pdf",
       device = "pdf",
       path = here(),
       scale = 1,
       widt = 9,
       height =9,
       units = c("cm"),
       dpi = 300)
