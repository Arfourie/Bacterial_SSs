#August 2023
#Arista Fourie
#Perform t-test comparisons in each plant study, adjusting for multiple testing among all SS's within a study
#Results of SS abundance are also concatenated to produce a final summary figure of which SS's are significantly different between rhizosphere and soil

library(dplyr)
library(ggplot2)
library("ggpubr")
library(tidyr)
library(rstatix)
library(plotrix)
library(OTUtable)
path <- "Summarise_all_plants"
SS_summaries <- sort(list.files(path, pattern = "_abundance_GTDB_MGs",full.names = TRUE))
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7","#7729C4","#5CD1B6")

SS_summaries_All <- data.frame()
for (i in SS_summaries) {
  SS_per_plant <- read.table(i,header = TRUE, sep = "\t")
  plant <- sapply(strsplit(basename(i), "_SS"), `[`, 1)
  SS_per_plant$plant <- rep(plant, times=length(SS_per_plant$med_TPM))
  SS_summaries_All <- rbind(SS_summaries_All,SS_per_plant) 
}

root <- c("SRR908273","SRR908275","SRR908276","SRR908208","SRR908211","SRR908272","SRR5195114","SRR5195112","SRR5195135","SRR5195133","SRR5195131","SRR5195137","SRR5195139","SRR5195141","SRR5195119","SRR5195121","SRR5195123",
          "C1A","C2A","C3A","D1A","D2A","D3A","SRR6797242","SRR6797249","SRR6797250","S-1-20","S-2-24",
          "S-3-19","S-4-15","S-4-19","S-1-2","S-1-3","S-2-6","S-2-9","S-3-10","S-3-13")
soil <-c("SRR908290","SRR908291","SRR908279","SRR908281","SRR5195136","SRR5195134","SRR5195132","SRR5195115","SRR5195113","SRR5195138","SRR5195140","SRR5195142","SRR5195120","SRR5195124","SRR5195122",
         "A1A","A2A","A3A","A1A_Ch","A2A_Ch","A3A_Ch","SRR6797243","SRR6797244","SRR6797246","S-1-9","S-2-11","S-2-21","S-3-24","S-4-23")

Summary_ed <- mutate(SS_summaries_All, source=ifelse((sample %in% root),"Root",sample))
Summary_ed <- mutate(Summary_ed, source=ifelse((sample %in% soil),"Soil",source))
Summary_ed <- Summary_ed %>% 
  select(-med_TPM) %>% 
  filter(SS!= "T4SS.I")

plant_types <- unique(Summary_ed$plant)

#Perform t-test per study, adjusting for multiple testing among all SS's
Plants.ttest.summary <- data.frame()
for (i in plant_types){
  Plant.indiv_ed <- filter(Summary_ed,plant == i)
  t_test_out <- Plant.indiv_ed %>% 
    group_by(SS) %>% 
    pairwise_t_test(Adjusted_abund ~ source, p.adjust.method = "none")
  Output <- t_test_out %>%
    dplyr::select(SS,p)  
  Output = Output[order(Output$p),]
  Output$BH = p.adjust(Output$p, method = "BH")
  Output$plant <- rep(i, times=length(Output$p))
  summary_statsR <- Plant.indiv_ed %>% 
    filter(source == "Root") %>% 
    group_by(SS) %>% 
    summarise(Root_meanTPM = mean(Adjusted_abund))
  summary_statsS <- Plant.indiv_ed %>% 
    filter(source == "Soil") %>% 
    group_by(SS) %>% 
    summarise(Soil_meanTPM = mean(Adjusted_abund))
  Output <- left_join(Output,summary_statsR, by="SS")
  Output <- left_join(Output,summary_statsS, by="SS")
  Plant.ttest.summary <- rbind(Plant.ttest.summary,Output)
  }
write.table(Plant.ttest.summary, "All plants t-test BH plants tested seperately.tsv", 
            sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Summarise the log fold diffference in SS abundance between rhizosphere and soil for visualisation purposes
Combined_logFold <- data.frame()
for (i in plant_types){
  Logfold_dat <- filter(Summary_ed,plant == i)
  Logfold_dat <- Logfold_dat %>% 
    group_by(SS,source) %>% 
    summarise(avg_TPM=mean(Adjusted_abund))
  Logfold_dat_ed <- pivot_wider(Logfold_dat, names_from = source, values_from = avg_TPM)
  Logfold_dat_ed$ratio <- Logfold_dat_ed$Root / Logfold_dat_ed$Soil
  Logfold_dat_ed$log2 <- log2(Logfold_dat_ed$ratio)
  Logfold_dat_ed$plant <- rep(i, times=length(Logfold_dat_ed$log2))
  Combined_logFold <- rbind(Combined_logFold,Logfold_dat_ed)
}

#Combine info on log fold differences and p-values for final figure construction
Pval_info <- Plants.ttest.summary %>% 
  dplyr::select(SS,BH,plant)
Pval_info$Plant_SS <- paste0(Pval_info$plant,"_",Pval_info$SS)
Pval_info <- mutate(Pval_info,signif=ifelse(BH <= 0.051, "Yes", "No"))
Pval_info <- Pval_info %>% 
  dplyr::select(-SS,-plant)
Combined_logFold$Plant_SS <- paste0(Combined_logFold$plant,"_",Combined_logFold$SS)
Combined_logFold <- left_join(Combined_logFold,Pval_info, by="Plant_SS")

#Use the abundance of the SS in the rhizosphere for relative size of dots in the dotplot, where larger dots will indicate higher abundance
TPM_select <- Combined_logFold
TPM_select$TPM_abund <- ifelse((TPM_select$Root > TPM_select$Soil),TPM_select$Root,TPM_select$Soil)

#Create final figure, including the logfold data, t-test p-values and abundance of SS in the rhizosphere
Fig_labels <-  TPM_select$plant %>% 
    unique()
  Fig_labels <-  sort(Fig_labels, decreasing = FALSE)
  TPM_select$plant <- factor(TPM_select$plant, levels=Fig_labels)
Logfold_plot <- ggplot(TPM_select, aes(x = log2, y = plant, fill=signif, size=TPM_abund)) +
  geom_point(alpha=0.9, shape=21, color="black") +   
  scale_fill_manual(values = c("No"="white" , "Yes"="#009900")) +
  facet_wrap(vars(SS), scales = "free_x", ncol=11) +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "#D3D3D3", colour = "#D3D3D3"),
        axis.line = element_line(colour = "black"))+
  geom_vline(xintercept = 0)
ggsave(filename = "SS abundance Logfold rhizosphere vs soil all plants.png", device="png", width=35, height =8, units="cm", dpi=400)