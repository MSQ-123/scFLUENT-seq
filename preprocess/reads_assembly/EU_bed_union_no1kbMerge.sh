#!/bin/bash
#SBATCH -J scbed
#SBATCH -p shenxhlab
#SBATCH --ntasks-per-node=1
#SBATCH -N 1
source /WORK/Samples/bio.sh

barcode_file="False_IAA.txt"

# Read barcodes from the file into an array
mapfile -t barcode_list < "$barcode_file"

#samtools view -H $bam_file > _header.sam

for barcode in "${barcode_list[@]}"; do
    #output_file="Extracted_reads_${barcode}.bam"
    #samtools view -h $bam_file | awk '$1 ~ /^@/ {next} {print $1}' > read_names.txt
    bedtools bamtobed -i splitBamIAA_new/${barcode}_split.sort.bam > splitBedIAA_new/${barcode}_split.sort.bed
    bedtools merge -i splitBedIAA_new/${barcode}_split.sort.bed -s -c 6 -o distinct > splitBedIAA_new/${barcode}_split.sort.merged.bed
    sort splitBedIAA_new/${barcode}_split.sort.merged.bed | uniq > splitBedIAA_new/${barcode}_split.sort.merged.uniq.bed
    awk -F $'\t' '{$(NF+1)=++i;}1' OFS=$'\t' splitBedIAA_new/${barcode}_split.sort.merged.uniq.bed |awk -F $'\t' '{$(NF+1)=++i;}1' OFS=$'\t' | awk -F $'\t' '{ print $1, $2, $3, $5,$6,$4}' OFS=$'\t' > splitBedIAA_new/${barcode}_union_expand.bed
# single cell also using -d 1000
    sort -k 1,1 -k2,2n splitBedIAA_new/${barcode}_union_expand.bed > splitBedIAA_new/${barcode}_union_sort_expand.bed
    bedtools merge -i splitBedIAA_new/${barcode}_union_sort_expand.bed -s -c 6 -o distinct > splitBedIAA_new/${barcode}_union_merged.bed
    sort -k 1,1 -k2,2n splitBedIAA_new/${barcode}_union_merged.bed | uniq > splitBedIAAFalse_new_noMerge/${barcode}_union_merged_uniq.bed
    rm splitBedIAA_new/${barcode}_split.sort.bed
    rm splitBedIAA_new/${barcode}_split.sort.merged.bed
    rm splitBedIAA_new/${barcode}_split.sort.merged.uniq.bed
    rm splitBedIAA_new/${barcode}_union_sort_expand.bed
    rm splitBedIAA_new/${barcode}_union_expand.bed
    rm splitBedIAA_new/${barcode}_union_merged.bed
done


#---------------
#EU start
#bedtools bamtobed -i ../R2_bam/IAA_10min.assigned_sorted.bam > IAA_10min.assigned_sorted.bed

#bedtools merge -i IAA_10min.assigned_sorted.bed -s -c 6 -o distinct > IAA_10min.assigned_sorted_merged.bed

#EU finished
#-------------------

#sort IAA_10min.assigned_sorted_merged.bed | uniq > IAA_10min_merged_noDuplicates.bed


#cat IAA_10min_merged_noDuplicates.bed CTRL_10min_merged_noDuplicates.bed IAA_10min_mature_merged_noDuplicates.bed CTRL_10min_mature_merged_noDuplicates.bed | sort -k 1,1 -k2,2n > IAA_CTRL_bamtobed_union_sort.bed

#----------------------

#EU start


#EU finished

#-------------------



#-----------------
#EU start 

#awk -F $'\t' '{$(NF+1)=++i;}1' OFS=$'\t' IAA_10min_merged_noDuplicates.bed |awk -F $'\t' '{$(NF+1)=++i;}1' OFS=$'\t' | awk -F $'\t' '{ print $1, $2, $3, $5,$6,$4}' OFS=$'\t' > IAA_bamtobed_union_sort_expand.bed
# single cell also using -d 1000
#bedtools merge -i IAA_bamtobed_union_sort_expand.bed -s -d 1000 -c 6 -o distinct > IAA_bamtobed_union_merged.bed
#sort -k 1,1 -k2,2n IAA_bamtobed_union_merged.bed | uniq > IAA_bamtobed_union_merged_noDuplicates.bed


#EU finished
#----------------------

#---------------------------



