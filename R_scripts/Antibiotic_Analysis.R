library(phyloseq)
library(dplyr)
library(here)


Clean_range <- read.csv(file=paste0(here(), "/data/Antibiotika_table.csv"))
FEV <- read.csv(file=paste0(here(), "/data/Regression_FEV1_year.csv"))
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
        axis.text.x = element_text(angle = 0, hjust = 0.5), axis.ticks.x = element_blank())
  ggsignif::geom_signif(comparisons = list(c("decliner", "stable")),
                        test = "t.test", 
                        map_signif_level= TRUE)
  
  
