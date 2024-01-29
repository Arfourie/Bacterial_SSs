#!/bin/bash

DATAROOT="genomes_SS_ID.tsv"
OUT="/home/Effectors/T5SS"

#1. Collect the gene ID's that were identified as T5bSS genes, in each genome of the At genome database
for name in $(cat /home/Effectors/Genomes.txt)
do
  grep ${name}_ $DATAROOT | \
  grep 'T5bSS_translocator' | awk -F"\t" '{print $1}' > "${OUT}/T5b/Data_files/${name}_T5b.txt"
done
cat "${OUT}/T5b/Data_files/*_T5b.txt" | sort > "${OUT}/T5b/Data_files/Select_genomes_T5bGenes_sort.txt"

find ${OUT}/T5b/Data_files/ -type f -empty -print -delete

#2. Extract the genomic location of T5bSS genes from GFF files and prepare input file for bedtools, perform this per genome
DATAFILES="/home/Effectors/T5SS/T5b/Data_files"
GFF="/home/Effectors/Data_files/gff"

for file in ${DATAFILES}/*_T5b.txt
do 
  name=$(basename ${file} _T5b.txt)
  for gene in $(cat ${file})
  do
    grep "${gene};" ${GFF}/${name}.gff | awk -F'[\t;]' '{print $1 "\t" $4-1 "\t" $5-1 "\t" $9 "\t" $7}' >> ${DATAFILES}/${name}_T5b.bed
  done
  
  sort -k1,1 -k2,2n ${DATAFILES}/${name}_T5b.bed > ${DATAFILES}/${name}_T5b_sorted.bed
  rm ${DATAFILES}/${name}_T5b.bed
done

#3. Use bedtools to identify the genes flanking the T5bSS genes in each genome
OUT="/home/Effectors/T5SS/T5b/Data_files/Flanking_genes"

for file in ${DATA}/*_T5b_sorted.bed
do
  name=$(basename ${file} _T5b_sorted.bed)
  bedtools closest -k 2 -iu -a ${file} -b ${GFF}/${name}_gff.bed -D ref > ${OUT}/${name}_T5b1geneDown.bed
  bedtools closest -k 2 -id -a ${file} -b ${GFF}/${name}_gff.bed -D ref > ${OUT}/${name}_T5b1geneUp.bed
  cat ${OUT}/${name}_T5b1geneDown.bed ${OUT}/${name}_T5b1geneUp.bed > ${OUT}/${name}_T5b1geneAll.bed
  rm ${OUT}/${name}_T5b1geneDown.bed ${OUT}/${name}_T5b1geneUp.bed
  
  awk -F"\t" '{print $9}' ${OUT}/${name}_T5b1geneAll.bed | sed 's/ID=//g' | sort | uniq > ${OUT}/${name}_1flankGenes.txt
  comm -13 ${DATAFILES}/Select_genomes_T5bGenes_sort.txt ${OUT}/${name}_1flankGenes.txt > ${OUT}/${name}_T5bFlankGenes.txt
done

#4. Extract the identified flanking proteins from the .faa files
#Need to create combined file with all proteins from all databases

PROTFILE="/All_At_proteins.faa"
FLANK_GENES="/home/Effectors/T5SS/T5b/Data_files/Flanking_genes"


for file in ${FLANK_GENES}/*_T5bFlankGenes.txt 
do
  name=$(basename $file _T5bFlankGenes.txt)
  seqtk subseq -l 60 "${PROTFILE}" $file > "${FLANK_GENES}/${name}_T5bFlankGenes.faa"
done

#5. Run iprscan on proteins to predict function

for file in $FLANK_GENES/*_T5bFlankGenes.faa
do
	name=$(basename ${file} _T5bFlankGenes.faa)

	/home/Software/interproscan-5.54-87.0/interproscan.sh \
	-appl Pfam \
	-f TSV, GFF3 \
	-i ${file} \
	-T /linuxhome/tmp \
	-b "/home/Effectors/T5SS/pfam/${name}"

done