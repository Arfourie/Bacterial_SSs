#August 2023
#Arista Fourie
#This script produces a barplot to summarise the abundance of each of the genes of the respective SSs in the metatranscriptomic data. 
#Data processed and summarised in the previous script is used as input (Transcriptome Summarise read total in TPM)

library(dplyr)
library(stringr)
library(ggplot2)
library(plotrix)
path <- "nrCDS_HtSeq"
countfiles <- sort(list.files(path, pattern="TPMfromAllgenes",full.names = TRUE))

Comb_counts_table <- data.frame()
for (i in countfiles){
  sample_name <- sapply(strsplit(basename(i), "_TPM"), `[`, 1)
  Count_table <- read.table(i,header = TRUE, sep = "\t")
  colnames(Count_table) <- c("DB_gene","TPM_count")
  Count_table$SS <- sapply(strsplit(basename(Count_table$DB_gene), "_"), `[`, 1)
  Count_table$sample <- rep(sample_name, times=length(Count_table$SS))
  Comb_counts_table <- rbind(Comb_counts_table, Count_table)
}
write.table(Comb_counts_table,"AllSamples_SScounts_summary.tsv", sep="\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

root <- c("S-1-2","S-1-3","S-2-6","S-2-9","S-3-10","S-3-13","S-1-20","S-2-24", "S-3-19","S-4-15","S-4-19")
soil <-c("S-1-9","S-2-11","S-2-21","S-3-24","S-4-23")
Comb_counts_table_rn <- mutate(Comb_counts_table, source=ifelse((sample %in% root),"Rhizosphere",sample))
Comb_counts_table_rn <- mutate(Comb_counts_table_rn, source=ifelse((sample %in% soil),"Soil",source))
Comb_counts_table_rn$SS <- sub("K02113", "atpD", Comb_counts_table_rn$SS) 
Comb_counts_table_rn$SS <- sub("K03040", "rpoA", Comb_counts_table_rn$SS) 
Comb_counts_table_rn$SS <- sub("K03531", "ftsZ", Comb_counts_table_rn$SS) 
Comb_counts_table_rn$SS <- sub("K03628", "rho", Comb_counts_table_rn$SS)

Comb_counts_m_se <- Comb_counts_table_rn %>% 
  select(SS,DB_gene,TPM_count,source) %>% 
  group_by(SS,DB_gene,source) %>% 
  summarise(mean_abund = mean(TPM_count),
            se_abund = std.error(TPM_count))

diamond_counts <- read.table("diamond_SS_gene_avg_summary_readcount.tsv", sep="\t", header = TRUE)
All_combined <- rbind(Comb_counts_m_se,diamond_counts)

m_se_plot <- ggplot(data = Comb_counts_m_se,
                    aes(x = DB_gene,
                        y = mean_abund, fill= source)) +
  geom_bar(stat="identity", width = 0.5, position = position_dodge(width=0.5)) +
  scale_fill_manual(values = c("#397D09","#945C23")) +
  geom_errorbar(aes(ymin = mean_abund - se_abund,ymax = mean_abund + se_abund),width = 0.2,position = position_dodge(0.5)) +
  facet_wrap(vars(SS), scales = "free", ncol = 3) +
  ylab("Mean abundance (TPM)")+
  theme_bw() +
  theme(legend.position = "right", axis.text.x = element_text(angle = 60, hjust = 0.95, vjust = 0.95))
ggsave(filename = "SS gene expression_TPM counts.jpeg", device="jpeg", width=20, height =30, units="cm")
