#August 2023
#Arista Fourie
#Summarise the taxonomic classification of the reads for each sample dataset, as produced by Kaiju. 
#This script is based on phylum rank summaries and also filters for the specific phyla of interest to obtain the total number of reads assigned to each phylum 
#The output was used in the next script to incorporate for calculation of SS abundance (Identification of SSs2_Summarise_diamond_results_convert_to_TPM.R)


library(dplyr)
kaiju_phyla <- read.table("All_MGs_kaiju_phylum.tsv", header = TRUE, sep = "\t", dec = ".") %>%   #Output from kaiju analysis table output on phylum level
  filter(Kingdom=="Bacteria") %>% 
  dplyr::rename(sample = file)
MG_samples <- unique(kaiju_phyla$sample)

#Run script for all sample files
All_MGs_phyla <- data.frame()
for (i in MG_samples){
  name=i
  Sample_phyla <- filter(kaiju_phyla,sample==i)
  Total_reads <- sum(Sample_phyla$reads)
  Sample_phyla$RA <- Sample_phyla$reads / Total_reads
  Phyla_info <- filter(Sample_phyla, Phylum=="Proteobacteria" | Phylum=="Bacteroidetes" | 
                            Phylum=="Chlamydiae")
  All_MGs_phyla <- rbind(All_MGs_phyla,Phyla_info)
}
All_MGs_phyla_select <- All_MGs_phyla %>% 
  select(sample,RA,Phylum)
write.table(All_MGs_phyla_select,"Kaiju_phyla_summary.tsv", sep = "\t",
            col.names = TRUE, row.names = FALSE, quote = FALSE)

