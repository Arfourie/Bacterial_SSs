#August 2023
#Arista Fourie
#Identify families with >40% of SS protein components present to determine the presence and abundance of the SS's of different families in each dataset
#Both the GTDB genome taxonomy and CAT classification of metagenome contig SS genes are included here to get overall taxonomic contents
library(dplyr)
library(ggplot2)
library("ggpubr")
library(stringr)
library(tidyr)
path <- "Plant"
SS_genes_MRPM<- sort(list.files(path, pattern="MRPM_SS_genes_AllReads",full.names = TRUE))
GTDB_taxon <- read.table("GTDBbac120_ncbitaxonomy.tsv",
                          header = TRUE, sep = "\t", colClasses = "character", na.strings = c("","\t")) %>% 
  dplyr::rename(ID = ident)
MGcontigs_taxon <- read.table("plant_SS_genes_CAT_lineage_info.tsv", 
                              header = TRUE, sep = "\t") %>%
  dplyr::rename(DB_gene = SS_gene)

Summary_SSs_family <- data.frame()
for (i in SS_genes_mapping) {
  MG.name <- sapply(strsplit(basename(i), "_"), `[`, 1)
  SS_genes_MRPM <- read.table(i, header = TRUE, sep = "\t")
  SS_genes_MRPM <- filter(SS_genes_MRPM, !str_detect(SS_gene,"T4P") & !str_detect(SS_gene,"Tad") & !str_detect(SS_gene,"Flg") & SS_gene != "T2SS_gspC" & SS_gene !="T2SS_gspD"
                                & SS_gene !="T2SS_gspE" & SS_gene !="T2SS_gspF" & SS_gene !="T2SS_gspG" & SS_gene !="T2SS_gspL"
                                & SS_gene !="T2SS_gspM" & SS_gene !="T2SS_gspN" & SS_gene !="T2SS_gspO" & !str_detect(SS_gene,"t4cp1") 
                                & !str_detect(SS_gene,"virb4") & SS_gene != "T4SS.I_traI" & SS_gene != "T4SS.I_traK" 
                                & SS_gene != "T4SS.I_traL" & SS_gene != "T4SS.I_trbB" & SS_gene != "T4SS.I_traW" & SS_gene != "T4SS.I_traY"
                                & SS_gene != "T9SS_gldJ_TIGR03524" & SS_gene != "T9SS_porQ" & SS_gene !="T9SS_porU"		     
                                & SS_gene !="T3SS_sctC" & SS_gene !="T3SS_sctN" & SS_gene !="T1SS_abc" & SS_gene !="T1SS_omf"
                                & SS_gene != "T6SSi_evpJ" & SS_gene != "T6SSi_tssA" & SS_gene != "T6SSi_tssB" & SS_gene != "T6SSi_tssC"
                                & SS_gene != "T6SSi_tssD" & SS_gene != "T6SSi_tssE" & SS_gene != "T6SSi_tssF" & SS_gene != "T6SSi_tssG"
                                & SS_gene != "T6SSi_tssH" & SS_gene != "T6SSi_tssI" & SS_gene != "T6SSi_tssK")
  SS_genes_MRPM$ID <- sub("_([^_]*)$", "", SS_genes_MRPM$DB_gene)
  #Split genes from GTDB and those from metagenomes
  GTDB_match <- SS_genes_MRPM %>%
    filter(!str_detect(ID,"^SRR")) 
  MGs_DB_match <- SS_genes_MRPM %>% 
   filter(str_detect(ID,"^SRR"))
  
  #Merge taxonomy info from each database with diamond results
  GTDB_match_taxonID <- merge(GTDB_match, GTDB_taxon, by.x = "ID", all.x=TRUE, incomparables=NULL, no.dups=FALSE) %>% 
    select(SS,SS_gene,DB_gene,MRPM,phylum,class,order,family,genus)
  MGs_DB_match_taxonID <- merge(MGs_DB_match, MGcontigs_taxon, by.x = "DB_gene", all.x=TRUE, incomparables=NULL, no.dups=FALSE) %>% 
    select(SS,SS_gene,DB_gene,MRPM,phylum,class,order,family,genus)
  Combined_SS <- rbind(GTDB_match_taxonID,MGs_DB_match_taxonID)
  Combined_SS[is.na(Combined_SS)] <- "unclassified"
  Combined_SS[Combined_SS == "no support"] <- "unclassified"
  Combined_SS[Combined_SS == ""] <- "unclassified"
  filename=paste0(MG.name,"_SS_genes_family_taxon_ID.tsv")
  write.table(Combined_SS,filename, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  #Determine abundance of family for each SS
  MG_Group <- Combined_SS %>% 
    group_by(SS, SS_gene, family) %>% 
    summarise(Avg_SS_gene_fam=sum(MRPM))
  MG_Group$name <- rep(MG.name, times= length(MG_Group$family))
  Summary_SSs_family <- rbind(Summary_SSs_family,MG_Group)
}  
write.table(Summary_SSs_family,"Plant_Summary_SS_family_abundance.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)

#Remove families with <40% of the genes for a SS, gapfill NA's and get average
SS_list <- unique(Summary_SSs_family$SS)  
Combined_SS_avgs <- data.frame()
for (g in SS_list){
  SS_name=g
  SS_select <- filter(Summary_SSs_family, SS==g)
  SS_gene_cols <- unique(SS_select$SS_gene)
  Nr_genes <- (length(SS_gene_cols) * 0.4)
  Nr_genes <- sapply(Nr_genes, ceiling)
  Last_col <- length(SS_gene_cols) + 3
  SS_pivot <- pivot_wider(SS_select, names_from=SS_gene, values_from=Avg_SS_gene_fam)
  SS_pivot[is.na(SS_pivot)] <- 0
  SS_pivot$Zeros <- rowSums(SS_pivot==0)
  SS_pivot <- SS_pivot %>% filter(!Zeros > Nr_genes) %>% 
    select(-Zeros)
  SS_pivot[SS_pivot == 0] <- NA
  #The 1 indicates to perform calculation on rows. Only apply to numeric columns and select lowest
  SS_pivot$minval <- apply(SS_pivot[,4:Last_col], 1, function(x) sort(x,na.last=NA)[1])
  SS_pivot <- apply(SS_pivot, 2, function(x) ifelse(is.na(x),SS_pivot$minval,x))
  SS_pivot <- as.data.frame(SS_pivot)
  numeric_cols <- colnames(SS_pivot[,4:ncol(SS_pivot)])
  SS_pivot <- SS_pivot %>% mutate_at(vars(numeric_cols), as.numeric)
  SS_pivot$SS_avg <- rowMeans(subset(SS_pivot,select=SS_gene_cols), na.rm = TRUE)
  SS_summary <- SS_pivot %>% select(SS, family, name, SS_avg)
  file.out <- paste0(g,"_abundance_Plant.tsv")
  write.table(SS_summary,file.out,sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  Combined_SS_avgs <- rbind(Combined_SS_avgs,SS_summary)
}   

root <- c("SRR908208","SRR908211","SRR908272")
soil <-c("SRR908279","SRR908281")

Combined_SS_avgs_rn <- mutate(Combined_SS_avgs, source=ifelse((name %in% root),"Root",name))
Combined_SS_avgs_rn <- mutate(Combined_SS_avgs_rn, source=ifelse((name %in% soil),"Soil",source)) %>% 
  filter(family != "unclassified")
write.table(Combined_SS_avgs_rn,"Plant_AllSS_abundance_combined.tsv", sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)  

#Create modified output for DESeq2 analysis. Also rounded values of read counts.
for (g in SS_list){
  SS_name=g
  SS_summary <- filter(Combined_SS_avgs_rn, SS==g) %>%
   select(family,SS_avg,name)
  SS_summary <- pivot_wider(SS_summary, names_from = name, values_from = SS_avg)
  is.num <- sapply(SS_summary, is.numeric)
  SS_summary[is.num] <- lapply(SS_summary[is.num], round, 0)
  SS_summary[is.na(SS_summary)] <- 0
  fileout=paste0(path,"/DESeq/",SS_name,"_MRPMGeneAvg_avg per family_rounded.tsv")
  write.table(SS_summary, fileout, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
}
