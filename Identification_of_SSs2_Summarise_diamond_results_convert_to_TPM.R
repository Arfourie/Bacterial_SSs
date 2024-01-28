#August 2023
#Arista Fourie
#Process the results produced by DIAMOND blastx to identify genes with >60% read coverage
#Convert read counts to Transcripts per kilobase million (TPM) and adjust the counts based on the abundance of speficic phyla that encode for the specific SS

library(dplyr)
library(ggplot2)
library(ggpubr)
library(stringr)
#Provide location and folder name for specific plant species
setwd("plant")
path <- "plant/diamond"
path2 <- "plant/diamond/summaries"

#Information on phylum adundance in eahc sample dataset, as produced by the script Identification_of_SSs1_Summarise_kaiju_table_for_SS_Summaries.R
Phyla_abund <- read.table("plant/Kaiju_phyla_summary.tsv", header = TRUE, sep = "\t", dec = ".") 

#Provide the SS gene information of the genes identified in the metagenomic contigs
#Use as input the SS info file produced from the script 2_summarise_data.sh
Contig_data_prep <- read.table("All_MGs_SSgenes_info_renamed.tsv", header = FALSE, sep = "\t")  
colnames(Contig_data_prep) <- c("DB_gene","SS_gene","SS") 

#Provide the SS gene information of the genes from the GTDB database
SS_type_GTDB_genes <- read.table("gtdb_All_SS_ID.tsv", header = TRUE, sep = "\t")
colnames(SS_type_GTDB_genes) <- c("DB_gene","SS","SS_gene")
SS_type_GTDB_genes <- select(SS_type_GTDB_genes,DB_gene,SS_gene,SS) %>%
	distinct()
All_genes_DB_SSinfo <- rbind(Contig_data_prep,SS_type_GTDB_genes)

#Read in the diamond result file for each sample, as processed by the script 6_prep_diamond_final_files_R.sh
diamond_files<- sort(list.files(path, pattern="ID40_GTDB_MGs_Rfile.tsv",full.names = TRUE))

#For each metagenome file, Calculate the fraction of the gene covered by reads,
#summarise this for all genes in one sample, then do this for all MG samples
for (f in diamond_files) {
  sample.name <- sapply(strsplit(basename(f), "_ID40"), `[`, 1)
  diamond_result <- read.table(f, header = FALSE, sep = " ")
  colnames(diamond_result) <- c("qseqid","sseqid","slen","sstart","send")
  Summary_gene_cov = data.frame()
  GeneList <- unique(diamond_result$sseqid)
  for (g in GeneList) {
    SS_gene=filter(diamond_result,sseqid==g)
    GeneLen=SS_gene[1,3]
    Gene_matrix=matrix(data=NA,nrow=nrow(SS_gene), ncol=GeneLen)
    for (i in 1:nrow(SS_gene)){
      Gene=(1:GeneLen)
      read_match=(SS_gene[i,4]:SS_gene[i,5])
      Gene_matrix[i,]=match(Gene,read_match,nomatch = 0)
    }
    Gene_matrix[Gene_matrix>0]=1
    Gene_covg=colSums(Gene_matrix)
    Gene_covg_frac=length(Gene_covg[Gene_covg>0])/length(Gene_covg)
    Read_mapping=nrow(SS_gene)
    row_info <- c(g,GeneLen,Read_mapping,Gene_covg_frac)
    Summary_gene_cov=rbind(Summary_gene_cov,row_info)
  }
  colnames(Summary_gene_cov)=c("DB_gene","Gene_length","Nr_reads_mapped","Gene_coverage")
  Summary_gene_60cov <- filter(Summary_gene_cov,Gene_coverage>=0.6)
  file.out <- file.path(path2,paste0(sample.name,"_Summary_genes_GTDB_MGs_60cov.tsv"))
  write.table(Summary_gene_60cov, file.out, sep="\t", col.names=TRUE, row.names=FALSE, quote = FALSE)
}

#Add SS info to genes and calculate their TPM from the read mapping data
diamond_Summ_files<- sort(list.files(path2, pattern="_Summary_genes_GTDB_MGs_60cov.tsv",full.names = TRUE))
for (h in diamond_Summ_files) {
  sample.name <- sapply(strsplit(basename(h), "_Summary"), `[`, 1)
  input <- read.table(h, header = TRUE, sep = "\t")
  input$RPK <- input$Nr_reads_mapped / (input$Gene_length/1000)
  PerMil <- sum(input$RPK) / 1000000
  input$TPM <- input$RPK / PerMil
  Final_TMP_result <- merge(input, All_genes_DB_SSinfo, by.x = "DB_gene", all.x=TRUE, incomparables=NULL, no.dups=FALSE)
  Final_TMP_result <- Final_TMP_result %>% 
    select(SS,SS_gene,DB_gene,TPM)
  file.out <- file.path(path2,paste0(sample.name,"_TPM_SS_genes_AllReads.tsv"))
  write.table(Final_TMP_result, file.out, sep="\t", col.names=TRUE, row.names=FALSE, quote = FALSE)
  SS_abund_summary <- Final_TMP_result %>%
    group_by(SS_gene) %>%
    summarise(sum_TPM=sum(TPM))
  SS_abund_summary$SS <- sapply(strsplit(basename(SS_abund_summary$SS_gene), "_"), `[`, 1)
  file2.out <- file.path(path2,paste0(sample.name,"_TPM_SS_abundance_summary.tsv"))
  write.table(SS_abund_summary, file2.out, sep="\t", col.names=TRUE, row.names=FALSE, quote = FALSE)
}
#Determine abundance of specifically revelant phyla in all samples
All_Proteobacteria <- filter(Phyla_abund, phylum=="Proteobacteria") %>% 
  group_by(sample) %>% 
  summarise(Proteobacteria_RA=sum(RA))
Proteo_Chlam <- filter(Phyla_abund, phylum=="Proteobacteria" | phylum=="Chlamydiae") %>% 
  group_by(sample) %>% 
  summarise(Proteo_Chlam_RA=sum(RA))
Bacteroidetes <- filter(Phyla_abund, phylum=="Bacteroidetes") %>% 
  group_by(sample) %>% 
  summarise(Bacteroidetes_RA=sum(RA))
  
#Remove genes not part of the calculation for SS abundance
#Summarise for each SS the abundance and also adjust abundance based on abudance of the phyla the encode for the SS
SS_abund_files<- sort(list.files(path2, pattern="_TPM_SS_abundance_summary",full.names = TRUE))
SS_abund_files_allSamples = data.frame()
for (i in SS_abund_files) {
  MGsample <- read.table(i,header = TRUE, sep="\t")
  MGsample <- filter(MGsample, !str_detect(SS_gene,"T4P") & !str_detect(SS_gene,"Tad") & !str_detect(SS_gene,"Flg") & SS_gene != "T2SS_gspC" & SS_gene !="T2SS_gspD"
                     & SS_gene !="T2SS_gspE" & SS_gene !="T2SS_gspF" & SS_gene !="T2SS_gspG" & SS_gene !="T2SS_gspL"
                     & SS_gene !="T2SS_gspM" & SS_gene !="T2SS_gspN" & SS_gene !="T2SS_gspO" & !str_detect(SS_gene,"t4cp1") 
                     & !str_detect(SS_gene,"virb4") & SS_gene != "T4SS.I_traI" & SS_gene != "T4SS.I_traK" 
                     & SS_gene != "T4SS.I_traL" & SS_gene != "T4SS.I_trbB" & SS_gene != "T4SS.I_traW" & SS_gene != "T4SS.I_traY"
                     & SS_gene != "T9SS_gldJ_TIGR03524" & SS_gene != "T9SS_porQ" & SS_gene !="T9SS_porU"                     
                     & SS_gene !="T3SS_sctC" & SS_gene !="T3SS_sctN" & SS_gene !="T1SS_abc" & SS_gene !="T1SS_omf"
                     & SS_gene != "T6SSi_evpJ" & SS_gene != "T6SSi_tssA" & SS_gene != "T6SSi_tssB" & SS_gene != "T6SSi_tssC"
                     & SS_gene != "T6SSi_tssD" & SS_gene != "T6SSi_tssE" & SS_gene != "T6SSi_tssF" & SS_gene != "T6SSi_tssG"
                     & SS_gene != "T6SSi_tssH" & SS_gene != "T6SSi_tssI" & SS_gene != "T6SSi_tssK")
  MGsample.name <- sapply(strsplit(basename(i), "_"), `[`, 1) 
  MGsample_Proteo <- filter(All_Proteobacteria, sample==MGsample.name)
  MGsample_Proteo_chlam <- filter(Proteo_Chlam, sample==MGsample.name) 
  MGsample_Bacteroid <- filter(Bacteroidetes, sample==MGsample.name)
  SS_TPM <- MGsample %>% 
    group_by(SS) %>% 
    summarise(med_TPM=median(sum_TPM))
  T2and6SS <- filter(SS_TPM, SS=="T2SS" | SS=="T6SSi")
  T2and6SS$Adjusted_abund <- T2and6SS$med_TPM * (MGsample_Proteo$Proteobacteria_RA)
  T2and6SS$sample <- rep(MGsample.name, times= length(T2and6SS$Adjusted_abund))
  T1_4andT5SS <- filter(SS_TPM, SS=="T1SS" | SS=="T4SS.T" | SS=="T4SS.I" | SS=="T5aSS" | SS=="T5bSS" | SS=="T5cSS") 
  T1_4andT5SS$Adjusted_abund <- T1_4andT5SS$med_TPM / 1
  T1_4andT5SS$sample <- rep(MGsample.name, times= length(T1_4andT5SS$Adjusted_abund))
  T3SS <- filter(SS_TPM, SS=="T3SS")
  T3SS$Adjusted_abund <- T3SS$med_TPM * (MGsample_Proteo_chlam$Proteo_Chlam_RA)
  T3SS$sample <- rep(MGsample.name, times= length(T3SS$Adjusted_abund))
  T6SSiii_T9 <- filter(SS_TPM, SS=="T6SSiii" | SS=="T9SS")
  T6SSiii_T9$Adjusted_abund <- T6SSiii_T9$med_TPM * (MGsample_Bacteroid$Bacteroidetes_RA)
  T6SSiii_T9$sample <- rep(MGsample.name, times= length(T6SSiii_T9$Adjusted_abund))
  SS_abund_files_allSamples <- rbind(SS_abund_files_allSamples,T2and6SS,T1_4andT5SS,T3SS,T6SSiii_T9)
}
write.table(SS_abund_files_allSamples, "plant_SS_abund_Relative_to_Phyla_abundance_GTDB_MGs.tsv", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)
