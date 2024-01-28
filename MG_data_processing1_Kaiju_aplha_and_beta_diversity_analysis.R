#August 2023
#Arista Fourie
#Analysis of alpha and beta diversity of metagenomic data, based on kaiju classification of reads on genus rank.
#This includes beta-diversity analysis of the combined data from all plant datasets, individual beta-diversity analyses were performed per dataset in antoher script (MG data processing2_Kaiju beta diversity indiv analysis.R)
#Read abundance of genera were converted to a phyloseq object to analyse data in a similar manner as metabarcoding data
#The identified genera were also used to construct an UpSet plot to summarise the core genera in the rhizosphere

library(dplyr)
library(tidyr)
library(stringr)
library (phyloseq)
library(metagMisc)
library("vegan")
library(ggpubr)
library(pairwiseAdonis)
library(UpSetR)
Kaiju_files <- sort(list.files(path, pattern="summary_genus_rn.tsv",full.names = TRUE))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#7729C4","#5CD1B6")

#Combine kaiju classification of reads for all plant datasets
AllKaiju <- data.frame()
for (i in Kaiju_files){
  Kaiju_data_plant <- read.table(i, sep="\t", header = TRUE) %>% 
    dplyr::select(-Kingdom)
  AllKaiju <- rbind(AllKaiju,Kaiju_data_plant)
}
AllKaiju <- mutate(AllKaiju,Genus=ifelse((Family %in% "Propionibacteriaceae" & Genus %in% "Ponticoccus"),"Enemella",Genus)) #Renaming changes in taxonomy

#Convert data to correct table formats to create a phyloseq object
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

#Alpha diversity
AlphaDiv <- estimate_richness(ps_filt)

#Add the origin info as extra column
AlphaDiv$host_source <- sample_data(ps_filt)$host_source

write.table(AlphaDiv, file="Alpha diversity.tsv", sep="\t", quote = FALSE)

Host.sources <- list(c("Athal_Stringlis_Root","Athal_Stringlis_Soil"),c("Wheat_Zhou_Soil","Chickpea_Zhou_Root"),c("Wheat_Zhou_Soil","Wheat_Zhou_Root"),
                     c("Wheat_Ofek_Soil","Wheat_Ofek_Root"),c("Cucumber_Ofek_Soil","Cucumber_Ofek_Root"),
                     c("Citrus_Italy_Xu_Soil","Citrus_Italy_Xu_Root"),c("Citrus_Brazil_Xu_Soil","Citrus_Brazil_Xu_Root"),
                     c("Athal_Sanchez_Soil","Athal_Sanchez_Root"),c("Citrus_Spain_Xu_Soil","Citrus_Spain_Xu_Root"),
                     c("Citrus_China_Xu_Soil","Citrus_China_Xu_Root"))
Fig_labels <-  metadata_p$host_source %>% 
  unique()
Fig_labels <-  c("Athal_Sanchez_Root","Athal_Sanchez_Soil","Athal_Stringlis_Root","Athal_Stringlis_Soil",
                 "Citrus_Brazil_Xu_Root","Citrus_Brazil_Xu_Soil","Citrus_China_Xu_Root","Citrus_China_Xu_Soil",
                 "Citrus_Italy_Xu_Root","Citrus_Italy_Xu_Soil","Citrus_Spain_Xu_Root","Citrus_Spain_Xu_Soil",
                 "Cucumber_Ofek_Root","Cucumber_Ofek_Soil","Wheat_Ofek_Root","Wheat_Ofek_Soil","Chickpea_Zhou_Root",
                 "Wheat_Zhou_Soil","Wheat_Zhou_Root")
ps_filt@sam_data$host_source <- factor(ps_filt@sam_data$host_source, levels=Fig_labels)
plot_richness(ps_filt, x="host_source", measures=c("Shannon"), color = "origin") + 
  geom_point(aes(fill=origin), size=3) + 
  theme(panel.grid = element_blank(),
                         panel.background = element_rect(fill = "white", colour = "white"),
                         axis.line = element_line(colour = "black"),
                         axis.text.x = element_text(angle = 90, hjust=1,size = 12),
                         axis.text.y = element_text(size =12))+
  scale_colour_manual(values=cbPalette) +
  stat_compare_means(comparisons = Host.sources, method = "t.test", 
                     p.adjust.methods="BH",label = "p.signif",label.y =6.3, size=5)
ggsave(filename = "Shannon div. index all plants.jpeg", device="jpeg", width=30, height =16, units="cm")

#pairwise t-test test, is there significant difference between populations
kruskal.test(Shannon ~ host_source, data = AlphaDiv)
#NB - this setting is similar to ggplot t-test. Normal parameters is pooled=TRUE (assumes equal variance)
pairwise.t.test(AlphaDiv$Shannon, AlphaDiv$host_source, p.adjust.method = "BH",pool.sd=FALSE)

#Do addiontal filter to remove low RA genera and analyse beta-diversity in the combined dataset
AllKaiju_rFilt = filter_taxa(ps_filt_r, function(x) mean(x) > 1e-5, TRUE)
bray.dist = phyloseq::distance(AllKaiju_rFilt, method="bray")
ord.pcoa.bray <- ordinate(AllKaiju_rFilt, method="PCoA", distance=bray.dist)

p2 <-plot_ordination(AllKaiju_rFilt, ord.pcoa.bray, color="origin", shape = "source", title="Bray-Curtis PCoA Rel. abundance") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black"))
p2 + geom_point(size=5) + scale_colour_manual(values=cbPalette)
ggsave(filename = "Bray-Curtis PCoA Rel. abundance all plants.jpeg", device="jpeg", width=20, height =16, units="cm")

#PERMANOVA
#Test homogeneity of the data, if dispersion around the centroid is equal the PERMANOVA is more reliable
anova(betadisper(bray.dist, sample_data(AllKaiju_rFilt)$host_source))
#Then do permutational ANOVA test
adonis(bray.dist ~ metadata_p$host_source + metadata_p$method)
pair.mod <- pairwise.adonis(bray.dist,sample_data(AllKaiju_rFilt)$host_source,p.adjust.m ="BH")
pairwise.adonis(bray.dist,sample_data(AllKaiju_rFilt)$host_source,p.adjust.m ="BH")

#UpSet plot - visualise core species
AllSamples_merge <- AllKaiju_rFilt %>%
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.005) %>% 
  filter(source == "Root") %>% 
  arrange(Genus)
Input_table <- AllSamples_merge %>% 
  dplyr::select(Sample,Abundance,host_source,Genus) 
Input_table <- Input_table %>% 
  group_by(host_source,Genus) %>% 
  summarise(mean=mean(Abundance))
Input_table_ed <- pivot_wider(Input_table, names_from = host_source, values_from = mean)
Input_table_ed[is.na(Input_table_ed)] <- 0 
Input_table_ed[,2:11][Input_table_ed[,2:11] > 0] <- 1 #Create True/False matrix and set true to 1
Input_table_ed <- as.data.frame(Input_table_ed)
jpeg(file="UpSet Plot_only rhizosphere.jpeg", width=25, height =18, units="cm", res=300)
upsetfig <- upset(Input_table_ed, nsets=10, nintersects=NA,
                  text.scale = 1.5,
                  order.by = "freq",mainbar.y.label="Number of shared genera")
upsetfig
dev.off()

#Select core genera, in 80% of datasets
Input_table_ed$Total <- rowSums(Input_table_ed[,2:11])
All_samples_genera <- Input_table_ed %>% filter(Total == 10)
Core_genera <- c(All_samples_genera$Genus)
Core_gen_data <- Input_table %>% filter(Genus %in% Core_genera)
All_genera_80 <- Input_table_ed %>% filter(Total >= 8)
All_genera_80 <- c(All_genera_80$Genus)
All_genera_80_data <- Input_table %>% filter(Genus %in% All_genera_80)

ggplot(All_genera_80_data, aes(x = host_source, y = mean, fill = Genus)) + 
  geom_bar(stat = "identity", colour = "black") +
  ylab("Relative Abundance of Genera") +
  ggtitle("Core Genera in >80% of the samples") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size=12),
        axis.text.y = element_text(size=12),
        legend.text = element_text(face = "italic", size=12), panel.grid = element_blank(),
        panel.background = element_rect(fill = "white", colour = "white"),
        axis.line = element_line(colour = "black")) +
  scale_y_continuous(expand = c(0,0))
ggsave(filename = "Rel Abund of core genera.jpeg", device="jpeg", width=20, height =15, units="cm")
