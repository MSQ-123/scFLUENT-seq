#!/bin/bash
#SBATCH -J scbed
#SBATCH -p shenxhlab
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
source /WORK/Samples/bio.sh

barcode_file="CTRL_barcode.txt"

# Read barcodes from the file into an array
mapfile -t barcode_list < "$barcode_file"

#samtools view -H $bam_file > _header.sam

for barcode in "${barcode_list[@]}"; do
    #output_file="Extracted_reads_${barcode}.bam"
    #samtools view -h $bam_file | awk '$1 ~ /^@/ {next} {print $1}' > read_names.txt
    bedtools bamtobed -i splitBamCTRL_new/${barcode}_split.sort.bam > splitBedCTRL_new/${barcode}_split.sort.bed
    bedtools merge -i splitBedCTRL_new/${barcode}_split.sort.bed -s -c 6 -o distinct > splitBedCTRL_new/${barcode}_split.sort.merged.bed
    sort splitBedCTRL_new/${barcode}_split.sort.merged.bed | uniq > splitBedCTRL_new/${barcode}_split.sort.merged.uniq.bed
    awk -F $'\t' '{$(NF+1)=++i;}1' OFS=$'\t' splitBedCTRL_new/${barcode}_split.sort.merged.uniq.bed |awk -F $'\t' '{$(NF+1)=++i;}1' OFS=$'\t' | awk -F $'\t' '{ print $1, $2, $3, $5,$6,$4}' OFS=$'\t' > splitBedCTRL_new/${barcode}_union_expand.bed
# single cell also using -d 1000
    sort -k 1,1 -k2,2n splitBedCTRL_new/${barcode}_union_expand.bed > splitBedCTRL_new/${barcode}_union_sort_expand.bed
    bedtools merge -i splitBedCTRL_new/${barcode}_union_sort_expand.bed -s -d 1000 -c 6 -o distinct > splitBedCTRL_new/${barcode}_union_merged.bed
    sort -k 1,1 -k2,2n splitBedCTRL_new/${barcode}_union_merged.bed | uniq > splitBedCTRL_new/${barcode}_union_merged_uniq.bed
    rm splitBedCTRL_new/${barcode}_split.sort.bed
    rm splitBedCTRL_new/${barcode}_split.sort.merged.bed
    rm splitBedCTRL_new/${barcode}_split.sort.merged.uniq.bed
    rm splitBedCTRL_new/${barcode}_union_sort_expand.bed
    rm splitBedCTRL_new/${barcode}_union_expand.bed
    rm splitBedCTRL_new/${barcode}_union_merged.bed
done


#---------------
