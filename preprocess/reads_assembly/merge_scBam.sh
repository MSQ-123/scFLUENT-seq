#!/bin/bash
#SBATCH -J mergeBam
#SBATCH -p shenxhlab
#SBATCH --ntasks-per-node=4
#SBATCH -N 1
source /WORK/Samples/bio.sh

# Assuming your barcode list is in a file named "barcode_list.txt"
# Each line contains a barcode to merge

# Read each barcode from the file and merge corresponding BAM files
#while IFS= read -r barcode; do
    # Collect BAM files corresponding to the current barcode
#    bam_files=$(ls ../split_R2_barcode/splitBamIAA_new/${barcode}_split.sort.bam)
    
    # Check if there are BAM files for the current barcode
#    if [ -n "$bam_files" ]; then
        # Extract header from the first BAM file for the current barcode
#        samtools view -H ${bam_files%% *} > header_${barcode}.sam
        
        # Merge BAM files for the current barcode along with the header
#        samtools merge -h header_${barcode}.sam merged_${barcode}.bam $bam_files
        
        # Remove temporary header file
#        rm header_${barcode}.sam
#    else
#        echo "No BAM files found for barcode $barcode"
#    fi
#done < False_IAA.txt

# Assuming your barcode list is in a file named "barcode_list.txt"
# Each line contains a barcode to match against BAM filenames

# Collect all matching BAM files
matching_bams=()
while IFS= read -r barcode; do
   # matching_bams+=( ../split_R2_barcode/splitBamIAA_new/"${barcode}"_split.sort.bam )
   #matching_bams+=( splitBamIAA_mature_False/"${barcode}"_split.sort.bam )
     matching_bams+=( raw_reads_split/splitBamFalseR2_mature/"${barcode}"_split.sort.bam )
done < False_IAA.txt

# Merge all matching BAM files
if [ ${#matching_bams[@]} -gt 0 ]; then
    #samtools merge True_IAA_merged_all_new.bam "${matching_bams[@]}"
    #samtools sort True_IAA_merged_all_new.bam > True_IAA_merged_all_new.sort.bam
    #samtools index True_IAA_merged_all_new.sort.bam
    #rm True_IAA_merged_all_new.bam

    #samtools sort False_IAA_merged_all_new.bam >False_IAA_merged_all_new.sort.bam
    #samtools index False_IAA_merged_all_new.sort.bam
    #rm False_IAA_merged_all_new.bam
	#samtools merge False_IAA_merged_all.bam "${matching_bams[@]}"
    #samtools merge False_IAA_merged_cage_mature.bam "${matching_bams[@]}"
    #samtools sort False_IAA_merged_cage_mature.bam > False_IAA_merged_cage_mature.sort.bam
    #samtools index False_IAA_merged_cage_mature.sort.bam

    samtools merge raw_reads_split/False_IAA_merged_mature_R2.bam "${matching_bams[@]}"
    samtools sort raw_reads_split/False_IAA_merged_mature_R2.bam > raw_reads_split/False_IAA_merged_mature_R2.sort.bam
    samtools index raw_reads_split/False_IAA_merged_mature_R2.sort.bam
    #rm False_IAA_merged_cage_mature.bam
    rm raw_reads_split/False_IAA_merged_mature_R2.bam
    echo "Merging completed."
else
    echo "No matching BAM files found."
fi

