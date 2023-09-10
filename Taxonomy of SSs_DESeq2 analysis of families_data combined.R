#August 2023
#Arista Fourie
#Combine DESeq2 results from all datasets to identify families with SS's enriched in multiple studies.
#Create a heatmap to illustrate the log fold enrichment of a SS of a specific family in rhizosphere or soil
library(dplyr)
library(ggplot2)
library(ggpubr)
path <- "Summarise_all_plants/DESeq/"
SS_family_abundance<- sort(list.files(path, pattern="AllSS_res_summary_filter",full.names = TRUE))

SS_summaries_All <- data.frame()
for (i in SS_family_abundance){
  SS_per_plant <- read.table(i,header = TRUE, sep = "\t")
  plant <- sapply(strsplit(basename(i), "_AllSS"), `[`, 1)
  SS_per_plant$plant <- rep(plant, times=length(SS_per_plant$SS))
  SS_summaries_All <- rbind(SS_summaries_All,SS_per_plant)
}

SS_list <- unique(SS_summaries_All$SS)
SS_list <- SS_list[SS_list != "T4SS.I"]

#All combined  
SS_summaries_select <- SS_summaries_All %>% 
  dplyr::filter(SS=="T2SS" | SS=="T3SS" |SS=="T4SS.T" |SS=="T5aSS" |SS=="T5bSS" |SS=="T5cSS" |SS=="T6SSi" | SS=="T6SSiii")
#Set to min and max 5 in order to see lower values clearer
SS_summaries_select <- mutate(SS_summaries_select,log2FoldChange=ifelse(log2FoldChange > 5,5,log2FoldChange))
SS_summaries_select <- mutate(SS_summaries_select,log2FoldChange=ifelse(log2FoldChange < -5,-5,log2FoldChange))

Fig_labels <- SS_summaries_select %>% 
  group_by(SS,Family) %>% 
  count(Family)
Fig_labels <-  Fig_labels[order(Fig_labels$n,Fig_labels$Family, decreasing = TRUE),]
SS_summaries_select$Family <- factor(SS_summaries_select$Family, levels=unique(Fig_labels$Family))
ggplot(SS_summaries_select, aes(plant, Family, fill=log2FoldChange)) + 
    geom_tile() +
    scale_fill_gradient2(low='#945C23',high="#397D09",mid="white",midpoint = 0) +
    facet_wrap(vars(SS), ncol=11) +
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = "#ABB0B8", colour = "#ABB0B8"),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 12), 
          axis.text.y = element_text(hjust=1, vjust = 0.5, size = 10))
    
ggsave(filename = "Heatmap select SS families all plants_maxLog5.jpeg", device="jpeg", width=38, height =46, units="cm")
write.table(SS_summaries_select, "Select SS's logfold values in all plants.tsv",
            sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)  
