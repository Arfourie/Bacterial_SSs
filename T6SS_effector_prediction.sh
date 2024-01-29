#!/bin/bash

DATAROOT="/home/Effectors"

#1. Extract list of predicted SS proteins only for the genomes of interest
for name in $(cat ${DATAROOT}/Genomes.txt)
do
  grep ${name}_ ${DATAROOT}/All_genomes_genes_SS_ID.tsv >> ${DATAROOT}/Data_files/Select_genomes_SS_ID.tsv
done

#2. Extract the names of the vgrG and hcp proteins from the list of predicted SS proteins
SS_GENES="/home/Effectors/Data_files/Select_genomes_SS_ID.tsv"
DATAFILES="/home/Effectors/Data_files"

for name in $(cat ${DATAROOT}/Genomes.txt)
do
  grep ${name}_ $SS_GENES | \
  grep 'T6SSi_tssI\|T6SSi_tssD' | awk -F"\t" '{print $1}' > "${DATAFILES}/vgrG.hcp/${name}_vgrG.hcp.txt"
done

find ${DATAFILES}/vgrG.hcp/ -type f -empty -print -delete

#3. Prepare bedtools format files by extrating the genomic location of vgrG and hcp genes from gff files

for file in $DATAFILES/vgrG.hcp/*_vgrG.hcp.txt
do 
  name=$(basename ${file} _vgrG.hcp.txt)
  for gene in $(cat ${file})
  do
    grep ${gene} ${DATAFILES}/gff/${name}.gff | awk -F'[\t;]' '{print $1 "\t" $4-1 "\t" $5-1 "\t" $9 "\t" $7}' >> ${DATAFILES}/vgrG.hcp/${name}_vgrG.hcp.bed
  done
  
  sort -k1,1 -k2,2n ${DATAFILES}/vgrG.hcp/${name}_vgrG.hcp.bed > ${DATAFILES}/vgrG.hcp/${name}_vgrG.hcp_sorted.bed
  rm ${DATAFILES}/vgrG.hcp/${name}_vgrG.hcp.bed
done

for file in $DATAFILES/gff/*.gff
do
  name=$(basename ${file} .gff)
  grep "CDS" $file | awk -F'[\t;]' '{print $1 "\t" $4-1 "\t" $5-1 "\t" $9 "\t" $7}' | sort -k1,1 -k2,2n > ${DATAFILES}/gff_bed/${name}_gff.bed
done

#4. Identify the 10 genes flanking each vgrG and hcp gene in the gff files and save gene names to obtain list of putative effectors
OUT="/home/arista/Effectors/Data_files/Flanking_genes"

grep "T6SSi" ${DATAFILES}/Select_genomes_SS_ID.tsv | sort > ${DATAFILES}/Select_genomes_SS_structureGenes_sort.txt

for file in ${DATAFILES}/vgrG.hcp/*_vgrG.hcp_sorted.bed
do
  name=$(basename ${file} _vgrG.hcp_sorted.bed)
  bedtools closest -k 10 -iu -a ${file} -b ${DATAFILES}/gff_bed/${name}_gff.bed -D ref > ${OUT}/${name}_vgrG.hcp10genesDown.bed
  bedtools closest -k 10 -id -a ${file} -b ${DATAFILES}/gff_bed/${name}_gff.bed -D ref > ${OUT}/${name}_vgrG.hcp10genesUp.bed
  cat ${OUT}/${name}_vgrG.hcp10genesDown.bed ${OUT}/${name}_vgrG.hcp10genesUp.bed > ${OUT}/${name}_vgrG.hcp10genesAll.bed
  rm ${OUT}/${name}_vgrG.hcp10genesDown.bed ${OUT}/${name}_vgrG.hcp10genesUp.bed
  
  awk -F"\t" '{print $9}' ${OUT}/${name}_vgrG.hcp10genesAll.bed | sed 's/ID=//g' | sort | uniq > ${OUT}/${name}_10flankGenes.txt
  comm -13 ${DATAFILES}/Select_genomes_SS_structureGenes_sort.txt ${OUT}/${name}_10flankGenes.txt > ${OUT}/${name}_putativeEffectors.txt
done

#5. Use gene name list to extract putative effector proteins
#Need to create combined file with all proteins from all databases
PROTFILE="All_At_proteins.faa"

for file in $DATAFILES/Flanking_genes/*_putativeEffectors.txt
do
  name=$(basename $file _putativeEffectors.txt)
  seqtk subseq -l 60 "${PROTFILE}" $file | seqtk seq -l 60 -L 50 > "${DATAFILES}/Effectors/${name}_putativeEffectors.faa"
done

#Also extract the vgrG and hcp proteins to identify functional domains on these proteins, cargo proteins might be attached to vgrG or hcp protein
for file in $DATAFILES/vgrG.hcp/*_vgrG.hcp.txt
do
  name=$(basename $file _vgrG.hcp.txt)
  seqtk subseq -l 60 "${PROTFILE}" $file | seqtk seq -l 60 > "${DATAFILES}/Effectors/${name}_vgrG.hcp.faa"
done

#6. Predict functional domains of putative effectors with iprscan

for file in $DATAFILES/Effectors/All*.faa
do
	name=$(basename ${file} .faa)

	/home/Software/interproscan-5.54-87.0/interproscan.sh \
	-appl Pfam \
	-f TSV, GFF3 \
	-i ${file} \
	-T /linuxhome/tmp \
	-b "/home/Effectors/pfam/${name}"

done
