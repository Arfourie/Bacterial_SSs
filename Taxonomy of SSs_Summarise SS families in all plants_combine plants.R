#August 2023
#Arista Fourie
#Combine results from all datasets to identify trends where the same families with the same SS's were identified in the majority of datasets.
#The abundance of a SS per family is summarised over all datasets to be visualised in a heatmap. Only families observed in 8-10 of the studies are shown. Output file was used to construct heatmap in iTOL.
#A horizontal barplot is also created to summarise the number of families in which a SS was oberved.
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
library(UpSetR)
path <- "SS_families"
SSFam_files <- sort(list.files(path, pattern="abundance_combined",full.names = TRUE))

Allplants <- data.frame()
for (i in SSFam_files){
  SSFam_data_plant <- read.table(i, sep="\t", header = TRUE)
  plant <- sapply(strsplit(basename(i), "_All"), `[`, 1)
  SSFam_data_plant$plant <- rep(plant, times=length(SSFam_data_plant$source))  
  Allplants <- rbind(Allplants,SSFam_data_plant)
}

Sum_AllPlants <- Allplants %>% 
  filter(str_detect(source, "Root")) %>% 
  group_by(plant,SS,family) %>% 
  summarise(mean=mean(SS_avg))
  
#Calculate total nr of families per SS
SS_count <- ungroup(Sum_AllPlants) %>% 
  select(SS,family) 
SS_count <-  unique(SS_count) %>% 
  count(SS)  
  
SS_list <- unique(Sum_AllPlants$SS)
SS_list <- SS_list[-4]
Core_fam_data_combined <- data.frame()
for (g in SS_list){
  SS_name=g
  SS_select <- filter(Sum_AllPlants, SS==g)
  Input_table_ed <- pivot_wider(SS_select, names_from = plant, values_from = mean)
  Input_table_ed[is.na(Input_table_ed)] <- 0 
  Input_table_ed[,3:12][Input_table_ed[,3:12] > 0] <- 1 #Create True/False matrix and set true to 1
  Input_table_ed <- as.data.frame(Input_table_ed)
  #Look at core families
  Input_table_ed$Total <- rowSums(Input_table_ed[,3:12])
  Core80perc <- Input_table_ed %>% filter(Total > 7)
  Core_families <- c(Core80perc$family)
  Core_fam_data <- SS_select %>% filter(family %in% Core_families)
  Core_fam_data_combined <- rbind(Core_fam_data_combined,Core_fam_data)
}

Core_families_all <- unique(Core_fam_data_combined$family)
Core80_inAll <- Sum_AllPlants %>% 
  filter(family %in% Core_families_all)
Core80_inAll_summarise <- Core80_inAll 
Core80_inAll_summarise$plant <- 1
Core80_inAll_summarise <- Core80_inAll_summarise %>% 
  group_by(SS,family) %>% 
  summarise(Nr_hosts = sum(plant))
Core80_inAll_summarise$Nr_hosts[Core80_inAll_summarise$Nr_hosts<8] <- 0
write.table(Core80_inAll_summarise, "Core80 families_heatmap.txt",sep="\t", quote = FALSE, row.names = FALSE)

Family_count <- Core80_inAll_summarise
Family_count <- subset(Family_count, Nr_hosts > 0) %>% 
  count(SS)

SS_Family_count <- ggplot(Family_count, aes(n,SS)) +
  geom_bar(stat = "identity") +
  theme(panel.grid = element_blank(),panel.background = element_rect(fill = "white"), axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20)) + 
  xlab("Nr. of families")
ggsave(filename = "Barplot SS family counts.jpeg", device="jpeg", width=20, height =15, units="cm")
