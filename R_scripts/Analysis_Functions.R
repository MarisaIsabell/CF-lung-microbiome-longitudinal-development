# Marisa Metzger 2019  
# Master thesis 
# functions for lung microbiome Analysi for usage in other R scripts 
# based on R packages phyloseq and microbiome 


colors_clinic_id <- c("patient 1" = "#ff0e70", 
                      "patient 2" = "#9b0016",
                      "patient 3" = "#f45300",
                      "patient 4" = "#eeef02",
                      "patient 5" = "#77d600",
                      "patient 6" = "#00c89a",
                      "patient 7" = "#5ddcd7",
                      "patient 8" = "#029fc7",
                      "patient 9" = "#666cff",
                      "patient 10" = "#500069",
                      "patient 11" = "#ffaafa",
                      "patient 12" = "#ff939a")



# Alpha Diversity 
# creates in the sample_data from the phyloseq object an column with alpha-diversity index, default-method: "Shannon"
# but you can also use other indexes (see ?estimate_richness)

alpha_column_to_ps <- function(ps,measures ="Shannon",...){
  richness <- phyloseq::estimate_richness(ps, split=TRUE, measures = measures)
  richness$`Idf..Nummer`<- rownames(richness)
  metadata <- meta(ps)
  metadata <- left_join(metadata, richness, by = "Idf..Nummer")
  rownames(metadata) <- metadata$`Idf..Nummer`
  ps <- merge_phyloseq(ps, sample_data(metadata))
  return(ps)
}


# make linear regression for alpha diversity. 

# this function works with the previously calculated Alpha-Diversity and Time_index 
# it creates an extra dataframe with the Slope, Intercept, R squared value, and p value in the columns for each patient (rows) 

regression_for_each_patient <- function(phyloseq.object = ps, alpha.measure = "Shannon",Time.index = "Time_index",...){
  metadata <- meta(phyloseq.object)
 colnames(metadata)[which(colnames(metadata) == alpha.measure)] <- "alpha"
  colnames(metadata)[which(colnames(metadata) == Time.index)] <- "index"
  
  regression_for_one_patient <- function(df){
    library(stats)
    #calculate the regression
    reg_tmp <- lm(alpha ~ index, data = df)
    #extract swlope and intercept value from the regression
    intercept_tmp <- coef(reg_tmp)[[1]]
    slope_tmp <- coef(reg_tmp)[[2]]
    # extract R_squared value from the regression
    r_squared <- summary(reg_tmp)$r.squared
    # make an anova test for the p value 
    p_value <- anova(reg_tmp)$"Pr(>F)"[1]
    # make a dataframe with all the information
    regression_df <- data.frame(Intercept = intercept_tmp, Slope = slope_tmp, R.squared = r_squared, P.value = p_value)
    return(regression_df)
  }
  Regression <- metadata %>% group_by(clinic_id) %>% do(.,regression_for_one_patient(.))
  #assign the patients to a group
  Regression$Group <- NA
  Regression <- within(Regression, Group[P.value < 0.05 && Slope > 0.0000000] <- "increasing")
  Regression <- within(Regression, Group[P.value < 0.05 && Slope < 0.0] <- "decreasing")
  Regression <- within(Regression, Group[P.value > 0.05] <- "consistent")

  return(Regression)
}

# calculate the regression
calculate_regression<- function(df,y,x,...){
  library(stats)
  #calculate the regression
  reg_tmp <- lm(y ~ x, data = df)
  #extract swlope and intercept value from the regression
  intercept_tmp <- coef(reg_tmp)[[1]]
  slope_tmp <- coef(reg_tmp)[[2]]
  # extract R_squared value from the regression
  r_squared <- summary(reg_tmp)$r.squared
  # make an anova test for the p value 
  p_value <- anova(reg_tmp)$"Pr(>F)"[1]
  
  
  # make a dataframe with all the information
  regression_df <- data.frame(Intercept = intercept_tmp, 
                              Slope = slope_tmp, 
                              R.squared = r_squared, 
                              P.value = p_value)
  return(regression_df)
}


# Most Abundant species for each patient
# this function creates for each patient a plot with the genus with at least 5 % abundance in at least one sample
# and stores it in a seperate path 

plot_5percent_RSV_for_each_patient <- function(phyloseq.object = ps, 
                                               path_for_storage = "C:/Users/Marisa/Google Drive/Masterarbeit"){
  ps <- phyloseq.object
  metadata <- meta(ps)
  patients <- unique(metadata$clinic_id)
  PMY <- vector("list", length(unique(sample_data(ps)$clinic_id)))
  names(PMY) <- unique(sample_data(ps)$clinic_id)
  
  for (i in patients){
    # prune the ps object for p 
    ps_tmp <- prune_samples(sample_data(ps)$clinic_id == i, ps)
    # make relative abundance 
    ps_tmpR <- microbiome::transform(ps_tmp, "compositional")
    # filter for the taxa wich are above the threshold
    ps_above <- filter_taxa(ps_tmpR, function(x) max(x) > 0.05, TRUE)
    top <- taxa_names(ps_above)
    # filter  the ps for the top taxa in the taxtable 
    taxtable <- tax_table(ps_tmpR) %>%
      data.frame(stringsAsFactors = FALSE)
    taxtable$seq <- rownames(taxtable)
    taxtable_top <- dplyr::filter(taxtable, rownames(taxtable) %in% top)
    # filter for all the "others"
    taxtable_others <- dplyr::setdiff(taxtable, taxtable_top)
    othersNames <- taxtable_others$seq
    # merge all the "others" in the phyloseq object, results in a phyloseq object with only 26 taxas
    ps_top <- merge_taxa(ps_tmpR, othersNames)
    
    # make the right order of the RSVs in the ps
    # make a vector with the order of the genus
    right_order <-c("Actinomyces", "Rothia", "Alloprevotella", "Bergeyella", "Capnocytophaga", "F0058", "Porphyromonas", "Prevotella", "Campylobacter", "Abiotrophia", "Gemella", "Granulicatella", "Megasphaera", "Oribacterium", "Parvimonas", "Staphylococcus", "Streptococcus", "Veillonella", "Fusobacterium", "Leptotrichia", "Candidatus_Saccharimonas", "Achromobacter", "Acinetobacter", "Actinobacillus", "Aggregatibacter", "Bordetella", "Eikenella", "Escherichia/Shigella", "Haemophilus", "Kingella", "Lautropia", "Moraxella", "Neisseria", "Pseudomonas", "Stenotrophomonas", "Treponema", "Mycoplasma", "Others")
    # extract taxatable as a dataframe   
    taxtabletop <- tax_table(ps_top) %>% 
      data.frame(stringsAsFactors = FALSE)
    # rename the NA Genus to "Others"
    taxtabletop$Genus[is.na(taxtabletop$Genus)] <- "Others"
    # change the factororder of the Genus column according to the rsv_level vector 
    taxtabletop$Genus <- gdata::reorder.factor(taxtabletop$Genus, new.order=right_order)
    # arrange the table accoring to the factors of the Genus levels 
    taxtabletop$names <- rownames(taxtabletop) # because the rownames get lost by arranging
    taxtabletop <- arrange(taxtabletop, Genus)
    rownames(taxtabletop) <- taxtabletop$names
    taxtabletop$names <- NULL
    taxtabletop <- as.matrix(taxtabletop)
    # remerge the net taxtable to the phyloseq object  
    tax_table(ps_top) <- taxtabletop
    
    # create a color pal accorrding to the Phylum membership 
    library(RColorBrewer)
    nb.color <- length(unique(tax_table(ps_top)[,"Genus"]))
    myColor <- colorRampPalette(brewer.pal(nb.color, "Set1"))
    genusList <- unique(tax_table(ps_top)[,"Genus"])
    genusPalette <- myColor(length(genusList))
    names(genusPalette) <- genusList
    genusPalette["Achromobacter"] <- "coral1"
    genusPalette["Abiotrophia"] <- "lightgoldenrod1"
    genusPalette["Acinetobacter"] <- "mediumpurple1"
    genusPalette["Actinobacillus"] <- "palevioletred"
    genusPalette["Actinomyces"] <- "turquoise4"
    genusPalette["Aggregatibacter"] <- "brown4"
    genusPalette["Alloprevotella"] <- "green3"
    genusPalette["Atopobium"] <- "steelblue"
    genusPalette["Bergeyella"] <- "olivedrab"
    genusPalette["Bordetella"] <- "hotpink"
    genusPalette["Campylobacter"] <- "lavender"
    genusPalette["Candidatus_Saccharimonas"] <- "lavenderblush"
    genusPalette["Capnocytophaga"] <-"forestgreen"
    genusPalette["Eikenella"] <- "indianred1"
    genusPalette["Escherichia/Shigella"] <- "lightpink"
    genusPalette["F0058"] <- "greenyellow"
    genusPalette["Fusobacterium"] <- "turquoise3"
    genusPalette["Gemella"] <- "chocolate1"
    genusPalette["Granulicatella"] <- "darkgoldenrod1"
    genusPalette["Haemophilus"] <- "mediumorchid3"
    genusPalette["Kingella"] <- "maroon3"
    genusPalette["Lautropia"] <- "salmon3"
    genusPalette["Leptotrichia"] <- "deepskyblue"
    genusPalette["Megasphaera"] <- "lightsalmon"
    genusPalette["Moraxella"] <- "plum3"
    genusPalette["Mycoplasma"] <- "slategray1"
    genusPalette["Neisseria"] <- "rosybrown3"
    genusPalette["Oribacterium"] <- "burlywood1"
    genusPalette["Parvimonas"] <- "gold2"
    genusPalette["Porphyromonas"] <-"mediumspringgreen"
    genusPalette["Prevotella"] <- "mediumaquamarine"
    genusPalette["Pseudomonas"] <- "red2"
    genusPalette["Rothia"] <- "mediumblue"
    genusPalette["Staphylococcus"] <- "yellow"
    genusPalette["Stenotrophomonas"] <- "violetred2"
    genusPalette["Streptococcus"] <- "orange2"
    genusPalette["Treponema"] <- "lightcyan"
    genusPalette["Veillonella"] <- "yellow3"
    genusPalette["Others"] <- "gray"
    
    
    
    # don't care about previous plots 
    p <- NULL 
    # create a plot and store it in a temporary variable 
    p <- plot_bar(ps_top, x = "Time_index", fill = "Genus")
    # changes the order of genera in the plot 
    df <- p$data
    # change the level order of the genus
    df$Genus <- gdata::reorder.factor(df$Genus, new.order=right_order)
    # replace the data in the plot object 
    p$data <- df
    
    p+
      geom_bar(aes(color= Genus, fill = Genus),
               stat = "identity", 
               position = "stack",
               width = 0.5)+
      xlim(0,60)+
      #scale_y_continuous(limits=c(0,1.0), breaks=(c(0,0.2,0.4,0.6,0.8,1.0)))+
      my_theme_ppt()+
      xlab("")+
      scale_color_manual(values=genusPalette)+
      scale_fill_manual(values= genusPalette)+
      theme(legend.position = "bottom")
    #geom_point(aes(y="FEV1pred", fill = "white", color ="black" ))
    
    #text(aes(y = -0.5, y = Time_index , label = Idf..Nummer))
    
    # store the plot
    PMY[[i]] <- p 
    # save the plot to the directory 
    ggsave(paste0(i, ".pdf"), 
           device = "pdf", 
           path = path_for_storage, 
           scale = 1, 
           width = 18, 
           height = 8,  
           units = c("cm"), 
           dpi = 300)
  }
}

# make tax table with most abunbdant genus(each patient, like in the plot) 
table_most_abundant_genus <- function(phyloseq.object = ps, 
                                              number_of_genus = 10, 
                                              path_for_storage = "C:/Users/Marisa/Google Drive/Masterarbeit"){
  ps <- phyloseq.object
  metadata <- meta(ps)
  patients <- unique(metadata$clinic_id)
  tax_df <- data.frame(matrix(ncol = 7, nrow=0))
  colnames(tax_df) <- colnames(tax_table(ps))
  
    for (i in patients){
    # prune the ps object for p 
    ps_tmp <- prune_samples(sample_data(ps)$clinic_id == i, ps)
    # make relative abundance 
    ps_tmpR <- microbiome::transform(ps_tmp, "compositional")
    #ps_tmpR <- tax_glom(ps_tmpR, taxrank = "Genus")
    # filter for the most abundant genus
   # top <- microbiome::top_taxa(ps_tmpR, n = number_of_genus)
  #  ps_top <- prune_taxa(top, ps_tmpR)
    # filter for genus > 1% 
    ps_top <- filter_taxa(ps_tmpR, function(x) max(x) > 0.01, TRUE)
    tax_tmp <- tax_table(ps_top) %>%
      data.frame(stringsAsFactors = FALSE)
    
    #tax_tmp <- as.data.frame(tax_table(ps_top))
    tax_df <- bind_rows(tax_df, tax_tmp)
    }
  # save the df
  tax_df <- distinct(tax_df)
  write.csv(tax_df, file = paste0(path_for_storage, "/Tax_table_mostabundant.csv"))
}

## TOP 25 Genus, merge others to "Others ## 
# ps should have absolute abundances, for searching for the most abundant taxa, the data will transformed to relative abundances
# but the output phyloseq object will habe absolute abundances again. 

topTaxa_and_others <- function(phyloseq.object = ps, 
                              number.of.top.taxa=10,
                              taxa.level = "Genus"){
  # relative abundance (to filter for the top25 species)
  ps_R <- microbiome::transform(phyloseq.object, transform = "compositional")
  ps_R <- tax_glom(ps_R, taxrank = taxa.level)
  # filter the most abundant 25 genus (later can be species)
  top <- microbiome::top_taxa(ps_R, n = number.of.top.taxa)

  # filter  the ps for the top taxa in the taxtable (with absolute abundance data!!)
  taxtable <- tax_table(ps_R) %>%
      data.frame(stringsAsFactors = FALSE)
  taxtable$seq <- rownames(taxtable)
  taxtable_top <- dplyr::filter(taxtable, rownames(taxtable) %in% top)
  
  # filter for all the "others"
  taxtable_others <- dplyr::setdiff(taxtable, taxtable_top)
  othersNames <- taxtable_others$seq
  # merge all the "others" in the phyloseq object, results in a phyloseq object with only 26 taxas
  ps_top <- merge_taxa(ps_R, othersNames)
  return(ps_top)
}

Taxa_abovepercent_and_others <- function(phyloseq.object = ps, 
                               percentage = 0.05){
  # relative abundance
  ps_R <- microbiome::transform(phyloseq.object, transform = "compositional")
  #ps_glom <- tax_glom(ps_R, taxrank = taxa.level)
  # filter for the taxa wich are above the threshold
  ps_above <- filter_taxa(ps_R, function(x) max(x) > 0.05, TRUE)
  top <- taxa_names(ps_above)
  
  # filter  the ps for the top taxa in the taxtable 
  taxtable <- tax_table(ps_R) %>%
    data.frame(stringsAsFactors = FALSE)
  taxtable$seq <- rownames(taxtable)
  taxtable_top <- dplyr::filter(taxtable, rownames(taxtable) %in% top)
  
  # filter for all the "others"
  taxtable_others <- dplyr::setdiff(taxtable, taxtable_top)
  othersNames <- taxtable_others$seq
  # merge all the "others" in the phyloseq object, results in a phyloseq object with only 26 taxas
  ps_top <- merge_taxa(ps_R, othersNames)
  return(ps_top)
}


# BETA Diversity 

# this function is used inside the function "beta_for_each_patient_over_time"

beta_individual_over_time <- function(phyloseq.object = ps, clinic_id, method = "horn"  ){
  ps <- phyloseq.object
  ps_ind <- prune_samples(sample_data(ps)$clinic_id == clinic_id , ps)
  rightorder <- sample_data(ps_ind)[order(sample_data(ps_ind)$Eingangsdatum),1]
  otu_table(ps_ind) <- otu_table(ps_ind)[, rightorder$Idf..Nummer]
  
  # calculate distance matrix
  dist <- phyloseq::distance(ps_ind, method = method, type = "samples") 
  # extract the diagonal 
  distmat <- as.matrix(dist)
  beta <- distmat[row(distmat) == col(distmat)-1]
  # assign to each value the compared samples 
  beta <- cbind(beta,rownames(distmat)[2:ncol(distmat)], rownames(distmat)[1:(ncol(distmat)-1)])
  colnames(beta) <- c("beta", "sample.x", "sample.y")
  beta <- as.data.frame(beta)
  # assign the timepoint to the beta diversity table 
  time <- dplyr::select(meta(ps_ind), Time_index)
  time$sample <- rownames(time)
  library(dplyr)
  beta <- left_join(beta, time, by = c("sample.x" = "sample"))
  beta <- left_join(beta, time, by = c("sample.y" = "sample"))
  beta$index <- clinic_id
  beta$beta <- as.numeric(levels(beta$beta))[beta$beta]
  return(beta)
}



beta_for_each_patient_over_time <- function(phyloseq.obejct = ps, method = "horn"){
  ps <- phyloseq.obejct
  patients <- unique(sample_data(ps)$clinic_id)
  beta_tmp <- data.frame()
  for (p in patients){
    beta_ind <- beta_individual_over_time(ps, p, method = method)
    beta_tmp <- plyr::rbind.fill(beta_tmp, beta_ind)
  }
  return(beta_tmp)
}

## beta diversity withing and between patients 

beta_diversity_within_and_between <- function(phyloseq.object = ps, method = "horn", type = "samples", Col_with_sample_index = "Idf..Nummer"){
  ps <- phyloseq.object
  # make the distance matrix
  dist <- phyloseq::distance(ps, method = method, type = type)
  # melt the matrix to a df
  library(reshape2)
  library(dplyr)
  dist_df <- melt(as.matrix(dist), varnames = c("sample1", "sample2"))
  # assign the patient_id to the samples 
  metadata <- meta(ps)
  names(metadata)[names(metadata) == Col_with_sample_index] <- "Idf.Nummer"
  patients <- metadata[, c( "Idf.Nummer" , "clinic_id")]
  dist_df <- left_join(dist_df, patients, by = c("sample1" = "Idf.Nummer"))
  colnames(dist_df)[colnames(dist_df) == "clinic_id"] <- "patient1"
  dist_df <- left_join(dist_df, patients, by = c("sample2" = "Idf.Nummer"))
  colnames(dist_df)[colnames(dist_df) == "clinic_id"] <- "patient2"
  # assign the within and between index, based on whether the patients 1 and 2 are the same or not 
  dist_df$group <- apply(dist_df, 1, function(x) 
    {if (x["patient1"] == x["patient2"]){dist_df$group <- "within"} else{dist_df$group <- "between"}})
  return(dist_df)
}
  

## % Percentage of Stability, added to the sample_data
# adds two columns to the phyloseq object: the Change column which shows you if the cluster will change until the next visit (1) or if it will stay the same(0)
# and the Stability column which tells you the percentage of stability 

stability_percentage <- function(phyloseq.object = ps, Cluster.column = "Hclust", Time.index = "Time_index"){
  metadata <- meta(phyloseq.object)
  metadata$Cluster <- metadata[Cluster.column][[1]]
  colnames(metadata)[colnames(metadata) == Time.index] <- "Time_index"
  metadata$Change <- NA

  stability_for_individuals <- function(df){
    metadata <- df
    metadata <- dplyr::arrange(metadata, Time_index)
    metadata$Cluster_next <- c(metadata$Cluster[-1], tail(metadata$Cluster, n =1))
    for (i in 1:nrow(metadata)){
      if (metadata$Cluster[i] == metadata$Cluster_next[i]){
        metadata$Change[i] <- 0 
      }else{
        metadata$Change[i] <- 1
      }}
    metadata$Stability <- 100-(sum(metadata$Change)/length(metadata$Change))*100
    return(metadata)
  }
  
  new_metadata <- metadata %>% group_by(clinic_id) %>% do(.,stability_for_individuals(.))
  new_metadata <- as.data.frame(new_metadata)
  rownames(new_metadata) <- new_metadata$`Idf..Nummer`
  
  ps<- merge_phyloseq(phyloseq.object, sample_data(new_metadata))
  return(ps)
}





