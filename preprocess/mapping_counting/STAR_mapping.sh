while read id
do
echo "#!/bin/bash
#SBATCH -J STAR
#SBATCH -p shenxhlab
#SBATCH --ntasks-per-node=4
#SBATCH -N 1">>$id.R2.sh

# for pair end mapping ###
# for R2 mapping (assembly) ###

echo "/NFSdata01/msq/10x_zUMI_pipe/assembly/STAR-2.7.3a/source/STAR --runThreadN 5 \
       --genomeDir /NFSdata01/msq/10x_zUMI_pipe/AID_data/STARsolo/mm10_STARidx_GTF \
       --readFilesIn processed_data/trimed/${id}_R2_extracted_val_2.fq.gz \
       --readFilesCommand zcat \
       --outBAMsortingThreadN 2 --outFilterMultimapNmax 1\
       --winAnchorMultimapNmax 1\
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix processed_data/R2_bam/${id}_R2_">>$id.R2.sh

echo "samtools index processed_data/R2_bam/${id}_R2_Aligned.sortedByCoord.out.bam">>$id.R2.sh

sbatch $id.R2.sh
done<id




