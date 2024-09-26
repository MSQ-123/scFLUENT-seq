#!/bin/bash
#SBATCH -J switchbed
#SBATCH -p shenxhlab
#SBATCH --ntasks=4
#SBATCH -N 2
source /WORK/Samples/bio.sh

#samtools sort -n ../CTRL_10min_pair_nojunction.bam -o CTRL_10min_pair_sortedByname.bam
#samtools index CTRL_10min_pair_sortedByname.bam
#pairToBed -abam CTRL_10min_pair_sortedByname.bam -b TSS_mouse.bed -ubam -s > CTRL_10min_pair_TSS_mouse.bed.bam 

# Input and output file names
# input_bed="CTRL_10min_pair_nojunction.bed12"
# output_bed="CTRL_10min_pair_nojunction.R1switched_test.bed12"

# Switch strand information for entries with '/2' suffix in column 4

bam_file="../pair_bam/CTRL_10min_pair_Aligned.sortedByCoord.out.bam"
barcode_file="CTRL_barcode.txt"

umi_tools dedup -I $bam_file --paired --per-cell -S splitBamCTRL_new/CTRL_10min_pair_split.sort.deduplicated.bam
samtools index splitBamCTRL_new/CTRL_10min_pair_split.sort.deduplicated.bam

bedtools bamtobed -bed12 -i splitBamCTRL_new/CTRL_10min_pair_split.sort.deduplicated.bam > CTRL_10min_pair.bed12

awk 'BEGIN{OFS="\t"} $4 ~ "/1$"{if($6 == "+") $6 = "-"; else if($6 == "-") $6 = "+"; print; next} {print}' CTRL_10min_pair.bed12 > switched_CTRL_new/CTRL_10min_pair.R1switched.bed12


bedtools bedtobam -bed12 -i switched_CTRL_new/CTRL_10min_pair.R1switched.bed12 -g mm10.chrom.sizes > switched_CTRL_new/CTRL_10min_pair.R1switched.bam


#samtools view -H switched/IAA_10min_pair.R1switched.bam > _header.sam
#samtools view -h switched/IAA_10min_pair.R1switched.bam | awk '$1 ~ "/1$"' | cat _header.sam - | samtools view -b > switched/IAA_10min.R1switched.R1.bam

#samtools view -h switched/IAA_10min_pair.R1switched.bam | awk '$1 ~ "/2$"' | cat _header.sam - | samtools view -b > switched/IAA_10min.R1switched.R2.bam
#bedtools intersect -abam IAA_10min.R1switched.R2.bam -b mm10_liftover+new_CAGE_peaks_phase1and2.bed -wa -s > IAA_10min.R1switched.R2.CAGE.bam


# Input BAM file names
#input_bam="CTRL_10min.R1switched.R2.CAGE.bam"
#output_read1_bam="CTRL_10min.R1switched.R1.CAGE.bam"

# Step 1: Extract read names with "/2" suffix
#echo "samtools view $id.R1switched.R2.bam | cut -f1 | sed 's/\\/2$//' > $id.read2_names.txt">>$id.switch.sh
# Step 2: Extract read 1 using the read names
#echo "samtools view $id.R1switched.R1.bam | grep -f $id.read2_names.txt | cat ${id}_header.sam - | samtools view -b > ${id}.R1switched.R1.CAGE.bam">>$id.switch.sh

#echo "samtools cat ${id}.R1switched.R2.CAGE.bam ${id}.R1switched.R1.CAGE.bam -o ${id}.R1switched.cat.CAGE.bam">>$id.switch.sh
samtools sort switched_CTRL_new/CTRL_10min_pair.R1switched.bam > switched_CTRL_new/CTRL_10min_pair.R1switched.sort.bam

samtools index switched_CTRL_new/CTRL_10min_pair.R1switched.sort.bam
bedtools bamtobed -bed12 -i switched_CTRL_new/CTRL_10min_pair.R1switched.sort.bam > switched_CTRL_new/CTRL_10min_pair.R1switched.sort.bed12

Rscript pseudolong.R switched_CTRL_new/CTRL_10min_pair.R1switched.sort.bed12 switched_CTRL_new/CTRL_10min_pair.R1switched.R2.pseudoLong.bed12

bedtools bedtobam -bed12 -i switched_CTRL_new/CTRL_10min_pair.R1switched.R2.pseudoLong.bed12 -g mm10.chrom.sizes > switched_CTRL_new/CTRL_10min_pair.R1switched.R2.pseudoLong.bam

samtools sort switched_CTRL_new/CTRL_10min_pair.R1switched.R2.pseudoLong.bam > switched_CTRL_new/CTRL_10min_pair.R1switched.R2.pseudoLong.sort.bam
samtools index switched_CTRL_new/CTRL_10min_pair.R1switched.R2.pseudoLong.sort.bam
#pairToBed -abam CTRL_10min_pair_nojunction.R1switched.sortByname.bam -b mm10_liftover+new_CAGE_peaks_phase1and2.bed -ubam -s > CTRL_10min_pair_nojunction.R1switched_CAGE_peaks.bam

#sbatch $id.switch.sh                              
#done<id

