#August 2023
#Arista Fourie
#Analysis of beta diversity of metagenomic data, based on kaiju classification of reads on genus rank.
#This includes individual beta-diversity analyses per dataset
#Read abundance of genera were converted to a phyloseq object to analyse data in a similar manner as metabarcoding data

library(dplyr)
library(tidyr)
library(stringr)
library (phyloseq)
library(metagMisc)
library("vegan")
library(ggpubr)
library(pairwiseAdonis)

path <- "Selected directory"
Kaiju_files <- sort(list.files(path, pattern="summary_genus_rn.tsv",full.names = TRUE))

AllKaiju <- data.frame()
disper_data_sum <- data.frame()
anova_sum <- data.frame()
pairwise_sum <- data.frame()
for (i in Kaiju_files){
  Kaiju_data_plant <- read.table(i, sep="\t", header = TRUE) %>% 
    dplyr::select(-Kingdom)
  plant <- sapply(strsplit(basename(i), "_"), `[`, 1)
  AllKaiju <- mutate(Kaiju_data_plant,Genus=ifelse((Family %in% "Propionibacteriaceae" & Genus %in% "Ponticoccus"),"Enemella",Genus))

  AllKaiju_OTUs <- AllKaiju %>% 
    dplyr::rename(sample=file) %>% 
    dplyr::select(sample,reads,Genus)

  Kaiju_data_pivot <- pivot_wider(AllKaiju_OTUs, names_from = sample, values_from = reads) 
  Kaiju_data_pivot[is.na(Kaiju_data_pivot)] <- 0
  Kaiju_data_pivot <- mutate(Kaiju_data_pivot,Genus=(paste0(Genus,"_gen"))) %>% 
    dplyr::rename(Genus_ID=Genus)
  
  Taxonomy_data <- AllKaiju %>%
    dplyr::select(-file,-reads)
  Taxonomy_data <- unique(Taxonomy_data)
  Taxonomy_data <- Taxonomy_data %>% 
    mutate(Genus_ID=(paste0(Genus,"_gen")))

  metadata_p <- read.table ("Metadata_metagenomes.tsv", header=TRUE, row.names = 1, sep = "\t") %>% 
    mutate(host_source=paste0(origin,"_",source))
  otu_table_p <- data.frame(Kaiju_data_pivot, row.names = 1)
  taxonomy_p <-  data.frame(Taxonomy_data, row.names = 6)
  otu_table_p = as.matrix(otu_table_p)
  taxonomy_p = as.matrix(taxonomy_p)
  FinalData_sampledata = sample_data(metadata_p)
  Final_otu = otu_table(otu_table_p, taxa_are_rows = TRUE)
  Final_tax = tax_table(taxonomy_p)
  ps <- phyloseq(Final_otu, Final_tax, FinalData_sampledata)

  #Prevalence should be selected based on how many samples are considered - Should be in at least three samples
  ps_filt <- phyloseq_filter_prevalence(ps,prev.trh = 0.04, abund.trh = 20)

  #Relative abundance values for PCoA
  ps_filt_r = transform_sample_counts(ps_filt, function(x) x/sum(x))

  #Do addiontal filter for low RA genera removed
  AllKaiju_rFilt = filter_taxa(ps_filt_r, function(x) mean(x) > 1e-5, TRUE)
  bray.dist = phyloseq::distance(AllKaiju_rFilt, method="bray")
  ord.pcoa.bray <- ordinate(AllKaiju_rFilt, method="PCoA", distance=bray.dist)
  p2 <-plot_ordination(AllKaiju_rFilt, ord.pcoa.bray, color="source", shape = "source", title="Bray-Curtis PCoA Rel. abundance") +
    theme(panel.grid = element_blank(),
          panel.background = element_rect(fill = "white", colour = "white"),
          axis.line = element_line(colour = "black"))
  p2 + geom_point(size=5) + scale_colour_manual(values = c("Soil"="#945C23" , "Root"="#397D09"))
  ggsave(filename = paste0(plant,"_Bray-Curtis PCoA Kaiju RA 1e-5.jpeg"), device="jpeg", width=20, height =16, units="cm")

  #PERMANOVA
  #Test homogeneity of the data, if dispersion around the centroid is equal the PERMANOVA is more reliable
  disper_data <- anova(betadisper(bray.dist, sample_data(AllKaiju_rFilt)$source))
  disper_data$host <- rep(plant,times=length(disper_data$Df))
  disper_data_sum <- rbind(disper_data_sum, disper_data)
  #Then do permutational ANOVA test
  anova_res <- adonis(bray.dist ~ sample_data(AllKaiju_rFilt)$source)
  anova_res <- anova_res$aov.tab
  anova_res$host <- rep(plant,times=length(anova_res$Df))
  anova_sum <- rbind(anova_sum,anova_res)
  pair.mod <- pairwise.adonis(bray.dist,sample_data(AllKaiju_rFilt)$source,p.adjust.m ="BH")
  pair.mod$host <- rep(plant,times=length(pair.mod$Df))
  pairwise_sum <- rbind(pairwise_sum,pair.mod)
}
write.table(disper_data_sum, file="Beta div dispersion.tsv", sep="\t", quote = FALSE)
write.table(anova_sum, file="Beta div anova.tsv", sep="\t", quote = FALSE)
write.table(pairwise_sum, file="Beta div pairwise comparison.tsv", sep="\t", quote = FALSE)
