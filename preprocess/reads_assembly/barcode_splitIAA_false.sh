#!/bin/bash
#SBATCH -J splitbam
#SBATCH -p shenxhlab
#SBATCH --ntasks=4
#SBATCH -N 2
source /WORK/Samples/bio.sh
source /WORK/Samples/bio2.sh
#bam_file="../pair_bam/IAA_10min_pair_Aligned.sortedByCoord.out.bam"
#barcode_file="barcode_list.txt"  # File containing one barcode per line


#mkdir splitBamIAA_mature_False_R2
barcode_file="False_IAA.txt"

#umi_tools dedup -I $bam_file --paired --per-cell -S splitBam/IAA_10min_pair_split.sort.deduplicated.bam
#samtools index splitBam/IAA_10min_pair_split.sort.deduplicated.bam

#bam_file="../R2_bam/IAA_10min_R2_Aligned.sortedByCoord.out.bam"

umi_tools dedup -I ../R2_bam/IAA_10min_mature_R2_Aligned.sortedByCoord.out.bam --per-cell -S raw_reads_split/IAA_10min_mature_R2.sort.deduplicated.bam

bam_file="raw_reads_split/IAA_10min_mature_R2.sort.deduplicated.bam"
samtools index $bam_file

#bam_file="switched_IAA_10min_mature/IAA_10min_mature.R1switched.R2.CAGE.sort.bam"
#rm switched_IAA_10min/IAA_10min.R1switched.R2.CAGE.pseudoLong.bam
# Check if the barcode file exists
if [ ! -f "$barcode_file" ]; then
    echo "Barcode file not found: $barcode_file"
    exit 1
fi

# Read barcodes from the file into an array
mapfile -t barcode_list < "$barcode_file"

samtools view -H $bam_file > false_header.sam

for barcode in "${barcode_list[@]}"; do
    #output_file="Extracted_reads_${barcode}.bam"
    #samtools view -h $bam_file | awk '$1 ~ /^@/ {next} {print $1}' > read_names.txt

    samtools view $bam_file | grep "_${barcode}_"  | cat false_header.sam - | samtools view -b > raw_reads_split/splitBamFalseR2_mature/${barcode}_split.bam
    # deduplicate using umi-tools

    samtools sort raw_reads_split/splitBamFalseR2_mature/${barcode}_split.bam > raw_reads_split/splitBamFalseR2_mature/${barcode}_split.sort.bam
    samtools index raw_reads_split/splitBamFalseR2_mature/${barcode}_split.sort.bam
    #umi_tools dedup -I splitBam/${barcode}_split.sort.bam --paired --per-cell -S splitBam/${barcode}_split.sort.deduplicated.bam

    rm raw_reads_split/splitBamFalseR2_mature/${barcode}_split.bam
    #rm splitBam/${barcode}_split.sort.bam
    #samtools index splitBam/${barcode}_split.sort.deduplicated.bam

    #samtools view -h -bS "$bam_file" | grep "_${barcode}_" > "$output_file"
    #samtools view -h -bS "$bam_file" | awk -v barcode="${barcode}" '$0 ~ "_" barcode "_"' > "$output_file"

    # If you want to convert the SAM file to BAM format, uncomment the following line:
    # samtools view -bS "$output_file" > "extracted_reads_${barcode}.bam"
done

