while read id
do
echo "#!/bin/bash
#SBATCH -J umi-tools
#SBATCH -p shenxhlab
#SBATCH --ntasks-per-node=4
#SBATCH -N 1">>$id.sh



echo "source /WORK/Samples/bio2.sh">>$id.sh
echo "source /WORK/Samples/bio.sh">>$id.sh
echo "umi_tools whitelist --stdin raw/Rawdata/$id.R1.fastq.gz \
                    --bc-pattern=CCCCCCCCCCCCCCCCCNNNNNNNNNNNN \
                   --set-cell-number=4000 \
		     --plot-prefix=$id.barcode_qc \
		     --log2stderr > $id.txt">>$id.sh
echo "umi_tools extract --bc-pattern=CCCCCCCCCCCCCCCCCNNNNNNNNNNNN \
                   --stdin raw70G/${id}_R1_001.fastq.gz \
                   --stdout processed_data/extracted/${id}_R1_extracted.fastq.gz \
                   --read2-in raw70G/${id}_R2_001.fastq.gz \
                   --read2-out=processed_data/extracted/${id}_R2_extracted.fastq.gz \
		   --whitelist=processed_data/$id.txt">>$id.sh
echo "trim_galore -q 5 --phred33 --stringency 12 --length 15 -e 0.1 --no_report_file --trim-n -a CCCCATGTCTCTGCGTTGATACCACTGCTT --gzip --clip_R2 30 --clip_R1 36 --paired processed_data/extracted/${id}_R1_extracted.fastq.gz processed_data/extracted/${id}_R2_extracted.fastq.gz -o processed_data/trimed">>$id.sh


sbatch $id.sh
done<id



#perform STARsolo alignment 

#featureCounts

while read id
do
echo "#!/bin/bash
#SBATCH -J umi-tools
#SBATCH -p shenxhlab
#SBATCH --ntasks-per-node=4
#SBATCH -N 1">>$id.sh
echo "source /WORK/Samples/bio.sh">>$id.sh
echo "source /WORK/Samples/bio2.sh" >> $id.sh


echo "featureCounts -a processed_data/counts/new_bed/bed_union_new_anno_no_nearest_refseq_stringent_df_refadd_gadd_first.gtf -s 1 -o processed_data/counts/new_bed/$id-custom_assigned processed_data/R2_bam/${id}_R2_Aligned.sortedByCoord.out.bam -R BAM -T 4 --ignoreDup -M --primary -g gene_id --extraAttributes Set_id,Subset_id">>$id.sh
echo "samtools sort processed_data/counts/new_bed/${id}_R2_Aligned.sortedByCoord.out.bam.featureCounts.bam -o processed_data/counts/new_bed/$id.assigned_sorted.bam">>$id.sh
echo "samtools index processed_data/counts/new_bed/$id.assigned_sorted.bam">>$id.sh
echo "samtools index processed_data/counts/Antisense_exon/$id.assigned_sorted.bam">>$id.Excount.sh

echo "umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell --wide-format-cell-counts -I processed_data/counts/new_bed/$id.assigned_sorted.bam -S processed_data/counts/new_bed/$id.custom_counts.tsv.gz" >> $id.sh

sbatch $id.sh
done<id
