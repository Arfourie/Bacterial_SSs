#August 2023
#Arista Fourie
#The read mappings to CDSs were counted and converted to TPM by HiSeq. 
#Here results are combined for all samples and modified to create barplots that can summarise transcript abundance of the SS genes.

library(dplyr)
library(stringr)
library(ggplot2)
library(ggpubr)
All_genes_DB_SSinfo <- read.table("All_MGs_SSgenes_info.tsv", sep="\t", header = TRUE)
colnames(All_genes_DB_SSinfo) <- c("DB_gene","SS_gene","SS")
gtf.file <- read.table("trinity_comb.transdecoder_cdsShift.gtf", sep="\t") %>% 
  dplyr::rename(DB_gene = V1)
Select_list <- c("K03040","K03628","K03531","K02113","K02134")
SS_genes <-unique(All_genes_DB_SSinfo$SS_gene)
SS_genes <- SS_genes[SS_genes !="T4P_pilM" & SS_genes !="Tad_tadZ" & SS_genes != "Flg_flgB" & SS_genes != "Flg_flgC" & SS_genes !="Flg_fliE"]
Select_list <- c(Select_list, SS_genes)
gtf.file <- gtf.file %>% mutate(gene_length = V5 - V4)
gtf_length.file <- gtf.file %>% select(DB_gene,gene_length)
KO_terms <- read.table("comb_KOterms.tsv"), header = FALSE)
colnames(KO_terms) <- c("DB_gene","KO_term")

Total_mapped <- sort(list.files(pattern="Comb_AllnrCDS.count", full.names=TRUE))
for (i in Total_mapped){
  sample.name <- sapply(strsplit(basename(i), "_Comb_AllnrCDS.count"), '[',1)
  map.file <- read.table(i, header = FALSE, sep = "\t")
  colnames(map.file) <- c("DB_gene","read_count")
  map.file <- head(map.file, -5)
  map.file <- map.file %>% 
    filter(read_count>0)
  input <- merge(map.file,gtf_length.file, by="DB_gene",)
  #Add KO and SS info
  input_functions <- left_join(input,KO_terms, by="DB_gene")
  input_functions <- left_join(input_functions, All_genes_DB_SSinfo, by="DB_gene")
  input_functions <- mutate(input_functions,annot=ifelse(!is.na(SS_gene),SS_gene, KO_term))
  input_select <- input_functions %>% 
    select(DB_gene,read_count,annot) %>% 
    filter(!is.na(annot))
  input_select <- input_select %>% 
    group_by(annot) %>% 
    summarise(read_total=sum(read_count))
  input_select$sample <- rep(sample.name,times=length(input_select$annot))
  file2.out <- file.path(paste0(sample.name,"_ReadTotalAllgenes_functionsAnnot.tsv"))
  write.table(input_select, file2.out, sep="\t", col.names=TRUE, row.names=FALSE, quote = FALSE)
  SS_abund_summary <- input_functions %>% 
    filter(annot %in% Select_list) %>% 
    group_by(annot) %>% 
    summarise(read_total=sum(read_count))
  file3.out <- file.path(paste0(sample.name,"_ReadTotal_SSgenesandHK.tsv"))
  write.table(SS_abund_summary, file3.out, sep="\t", col.names=TRUE, row.names=FALSE, quote = FALSE)
}
