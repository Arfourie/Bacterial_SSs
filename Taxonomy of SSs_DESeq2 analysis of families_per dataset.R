#August 2032
#Arista Fourie
#Identification of families with a SS that were enriched in the rhizosphere or soil, using DESeq2 analysis
library(dplyr)
library(DESeq2)
library(ggpubr)
path <- "DeSeq/"
SS_family_abundance<- sort(list.files(path, pattern="TPMGeneAvg_avg per family_rounded.tsv",full.names = TRUE))

#Added manual size factor of 1, since values are already TPM values, thus it should not be further adjusted to a relative abundance
AllSS_res_summary <- data.frame()
for (i in SS_family_abundance){
  SS.name <- sapply(strsplit(basename(i), "_"), `[`, 1)
  countData <- read.table(i, sep="\t",header = TRUE, row.names = 1) 
  countData <- countData[rowSums(countData == 0) <= 2, ] 
  countData <- as.matrix(countData)
  colData <- read.table("Col Info.tsv", sep="\t",header = TRUE)
  dds <- DESeqDataSetFromMatrix(countData=countData,
                                colData=colData,
                                design=~condition)
  dds$condition <- relevel(dds$condition, ref="soil")
  size_val <- rep(1,times=ncol(countData))
  sizeFactors(dds) <- size_val
  dds <- DESeq(dds,fitType='mean')
  res <- results(dds)
  file.out <- paste0(path,SS.name,"_DeSeq2_result.tsv")
  write.table(res, file.out, sep="\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  res_summary <- as.data.frame(res) %>% 
    select(baseMean,log2FoldChange,padj)
  res_summary$SS <- rep(SS.name, times= length(res_summary$log2FoldChange))
  res_summary <- tibble::rownames_to_column(res_summary, "Family")
  file.out2 <- paste0(path,SS.name,"_DeSeq2_result_summary.tsv")
  write.table(res_summary, file.out2, sep="\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  AllSS_res_summary <- rbind(AllSS_res_summary,res_summary)
}
write.table(AllSS_res_summary, "AllSS_results_family_summary.tsv",col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")

AllSS_res_summary_filt <- AllSS_res_summary %>% 
  filter(padj < 0.05) 
write.table(AllSS_res_summary_filt, "Plant_AllSS_res_summary_filter.tsv",col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")
