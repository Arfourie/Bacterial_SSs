#August 2032
#Arista Fourie
#Identification of families with a SS that were enriched in the rhizosphere or soil, using DESeq2 analysis
path <- "DeSeq/"
SS_family_abundance<- sort(list.files(path, pattern="abundance_Plant_AvgReads.tsv",full.names = TRUE))

AllSS_res_summary <- data.frame()
for (i in SS_family_abundance){
  countData <- read.table(i, sep="\t",header = TRUE) %>% 
    dplyr::select(family,SS_avg,name)
  countData <- pivot_wider(countData, names_from = name, values_from = SS_avg)
  countData[is.na(countData)] <- 0
  countData <- countData[rowSums(countData == 0) <= 4, ]
  countData <- filter(countData, family!="unclassified")
  is.num <- sapply(countData, is.numeric)
  countData[is.num] <- lapply(countData[is.num], round, 0)
  countData <- countData %>% remove_rownames %>% column_to_rownames(var='family')
  countData <- as.matrix(countData)
  colData <- read.table("Col Info_order.tsv", sep="\t",header = TRUE)
  dds <- DESeqDataSetFromMatrix(countData=countData,
                                colData=colData,
                                design=~condition)
  dds$condition <- relevel(dds$condition, ref="soil")
  dds <- DESeq(dds,fitType='mean')
  res <- results(dds)
  res_summary <- as.data.frame(res) %>% 
    select(baseMean,log2FoldChange,padj)
  res_summary$SS <- rep(SS.name, times= length(res_summary$log2FoldChange))
  res_summary <- tibble::rownames_to_column(res_summary, "Family")
  file.out2 <- paste0(path2,SS.name,"_DeSeq2_result_rawReadsAvg.tsv")
  write.table(res_summary, file.out2, sep="\t", col.names = TRUE, row.names = TRUE, quote = FALSE)
  AllSS_res_summary <- rbind(AllSS_res_summary,res_summary)
}
write.table(AllSS_res_summary, "AllSS_results_family_summary.tsv",col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")

AllSS_res_summary_filt <- AllSS_res_summary %>% 
  filter(padj < 0.05) 
write.table(AllSS_res_summary_filt, "Plant_AllSS_res_summary_filter.tsv",col.names = TRUE, row.names = TRUE, quote = FALSE, sep="\t")
