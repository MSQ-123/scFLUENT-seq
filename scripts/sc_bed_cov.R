library(data.table)
library(ggplot2)
require("GenomicFeatures")
require("rtracklayer")
require("data.table")
require("ggplot2")
require("ggpubr")
require("cowplot")
require("ggthemes")
require("ggsci")
require("ggforce")
require("ggExtra")
require("ggrepel")
require("scales")
require("DT")
require("circlize")
#require("BSgenome.Hsapiens.UCSC.hg38")
require(BSgenome.Mmusculus.UCSC.mm10)
require("ggbio")
#require("phastCons100way.UCSC.hg38")
#require(phast)
require("GenomicAlignments")
require("genomation")
require("VennDiagram")
require("viridis")
library(tidyverse)


x = list.files('./Paper/rds/bed_union_sc/strand/')
x = x[grep(x = x, pattern = '*bed')]
bed_list<- lapply(x, function(i)fread(paste('Paper/rds/bed_union_sc/strand/', sep = '',i)))

# Define the set of chromosome to use
chromosomes = genomeStyles()[["Mus_musculus"]]$UCSC
chromosomes = setdiff(chromosomes, c("chrY", "chrM"))
chromosomes.lengths = seqlengths(Mmusculus)
chromosomes.lengths = chromosomes.lengths[names(chromosomes.lengths)%in%chromosomes]

# Load the pre-processed GRanges and print the list
load(file.path('../IntergenicTranscription-master/IAA_gffcompare/GencodeReference/', "mm10_GencodeM22_genes_annotations.RData"))

genome.length = data.table(as.data.frame(seqinfo(genes)), keep.rownames="seqnames")[seqnames%in%c(chromosomes,'chrY'), sum(seqlengths)]

gr_list<- lapply(bed_list, function(i)GRanges(seqnames = i$V1, IRanges(i$V2, i$V3),i$V6))
# bulk
gr_list2 <- c()
gr_list2<- c(gr_list2,gr_list[[1]])[[1]]
for(i in 2:length(gr_list)){
  gr_list2 <- c(gr_list2,gr_list[[i]])
}

# for bulk level coverage , use reduce ####
or<- findOverlaps(gr_CTRL, new.anno.genes)
idx <- as.integer(names(which(table(or@from)==1)))
or<- or[or@from %in% idx]
gr_tmp <- gr_CTRL[or@from]
gr_tmp$Set <- new.anno.genes[or@to]$Set
gr_tmp$Subset <- new.anno.genes[or@to]$Subset
gr_tmp <-gr_tmp[!is.na(gr_tmp$Subset),]
nSet <- names(table(gr_tmp$Subset))
widthSet <- list()
for (i in 1:length(nSet)){
  widthSet[[nSet[i]]] <- sum(width(GenomicRanges::reduce(gr_tmp[gr_tmp$Subset == nSet[i]], ignore.strand=FALSE)))
}

Total_cov <- sum(unlist(widthSet))

unlist(widthSet)/Total_cov

gr_tmp%>%group_by(Subset)%>%summarise(cov = sum(width))%>%mutate(perc = cov/Total_cov)



# merged gr for IAA #

bed.plus<- lapply(bed_list, function(i)i[i$V6 == '+',])
bed.plus<- do.call(rbind, bed.plus)
bed.plus.gr<- reduce(GRanges(seqnames = bed.plus$V1, IRanges(bed.plus$V2, bed.plus$V3),bed.plus$V6))
sum(width(bed.plus.gr))/genome.length # 69.3%

bed.minus<- lapply(bed_list, function(i)i[i$V6 == '-',])
bed.minus<- do.call(rbind, bed.minus)
bed.minus.gr<- reduce(GRanges(seqnames = bed.minus$V1, IRanges(bed.minus$V2, bed.minus$V3),bed.minus$V6))
sum(width(bed.minus.gr))/genome.length # 68.65%

bed<- do.call(rbind, bed_list)
bed.gr<- reduce(GRanges(seqnames = bed$V1, IRanges(bed$V2, bed$V3),bed$V6), ignore.strand=T)
sum(width(bed.gr))/genome.length # 83.56% for IAA, 86.56% for external CTRL

# noMerged for False CTRL (550 cells): up to 65%, merged for 65%
# noMerged for True IAA (974 cells): up to 74.5%, merged for 74.5%
# combined for False and True IAA: 81.8%




# individual percent ####
Perc_list <- unlist(lapply(gr_list, function(i)sum(width(i[strand(i) == '+',]))/genome.length))
names(Perc_list) <- unlist(lapply(strsplit(x, split = '_'),function(i)i[[1]]))

Perc_list_df <- data.frame(Cells = 'IAA 1hr', Coverage = Perc_list, row.names = names(Perc_list))

#
Perc_list_total <- unlist(lapply(gr_list, function(i)sum(width(reduce(i, ignore.strand=T)))/genome.length))
names(Perc_list_total) <- unlist(lapply(strsplit(x, split = '_'),function(i)i[[1]]))
Perc_list_total_df <- data.frame(Cells = 'IAA 1hr', Coverage = Perc_list_total, row.names = names(Perc_list_total))
Perc_list_total_df <- Perc_list_total_df[IAA_barcode,]
Perc_list_total_df[rownames(Perc_list_total_df)%in% True_IAA,]$Cells <- 'IAA 1hr'
Perc_list_total_df[rownames(Perc_list_total_df)%in% False_IAA,]$Cells <- 'CTRL'
Perc_list_total_df$Cells <- factor(Perc_list_total_df$Cells, levels = c('CTRL','IAA 1hr'))


p1 <- ggplot(Perc_list_total_df, aes(x=Cells, y=Coverage, color = Cells)) + 
  geom_violin(trim=FALSE)+ #scale_y_log10() + 
  #ylim(c(10,200))+
  scale_color_manual(values = c('CTRL' = "#009999",
                               'IAA 1hr' = "#D6604D"))+
  geom_boxplot(width=0.1, position = position_dodge(width = 0.9), fill = 'white')+
  theme_classic()+
  stat_compare_means(comparisons = list(c('CTRL','IAA 1hr')), label = 'p.signif')+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"), legend.position = 'none')

library(eoffice)
plot.width = 8
plot.height = 8
font.size = 12
topptx(p2, filename = './Paper/Figures/Fig1/sc_cov_total.pptx')

out.plot <- './Paper/Figures/Fig1/batch_corr.pdf'
pdf(out.plot, width=plot.width, height=plot.height, pointsize=font.size)

out.dev <- dev.off()

saveRDS(Perc_list_df_rev, file = './Paper/rds/scPerc_list_IAA_nascent_rev.rds')
saveRDS(Perc_list_total_df, file = './Paper/rds/scPerc_list_IAA_nascent_total.rds')

# calculate genome coverage versus cell number ####

cell_n <- seq(from = 10, to = 1500)
names(bed_list) <- unlist(lapply(strsplit(x, split = '_'),function(i)i[[1]]))
seq_cov<- c()
ncells<- length(bed_list)

IAA_list <- bed_list[names(bed_list)%in%True_IAA]
CTRL_list <- bed_list[names(bed_list)%in%False_IAA]

for (n in cell_n) {
  idx <- sample(ncells, size = n)
  bed_tmp <- bed_list[idx]
  bed<- do.call(rbind, bed_tmp)
  bed.gr<- reduce(GRanges(seqnames = bed$V1, IRanges(bed$V2, bed$V3),bed$V6), ignore.strand=T)
  seq_cov[n-9]<- sum(width(bed.gr))/genome.length # 69.3%
}

seq_cov_df<- data.frame(cell_n = cell_n, genome_cov = seq_cov)

seq_cov_df$Sample <- 'IAA 1hr'
seq_cov_df_CTRL$Sample <- 'CTRL'
tab <- rbind(seq_cov_df_CTRL, seq_cov_df)

p2 <-  ggplot(seq_cov_df, aes(x=cell_n, y=genome_cov, color = Sample)) + 
  geom_jitter(alpha = 0.5)+
   geom_smooth(color = '#B2182B', lwd = 1)+
   xlab('Number of cells')+
   ylab('Genome coverage')+
  scale_color_manual(values = c('IAA mix' = "black"))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

p4 <-  ggplot(tab, aes(x=cell_n, y=genome_cov, color = Sample)) + 
  geom_jitter(alpha = 0.5)+
  geom_smooth(color = 'darkgrey', lwd = 1.5)+
  xlab('Number of cells')+
  ylab('Genome coverage')+
  scale_color_manual(values = c('CTRL' = "#009999",
                                'IAA' = "#D6604D"))+
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))


out.plot <- './Paper/Figures/Fig1/ncell_cov.pdf'
pdf(out.plot, width=plot.width, height=plot.height, pointsize=font.size)

out.dev <- dev.off()

# compute correlation between UMI and gene length ####

CTRL_new_df <- as.matrix(Genic_seu[, Genic_seu$orig.ident == 'CTRL_10min_beads']@assays$RNA@counts)
CTRL_mat_df <- as.matrix(Genic_seu[, Genic_seu$orig.ident == 'CTRL_10min_mature']@assays$RNA@counts)

rownames(CTRL_new_df)
CTRL_new_meta<- read.table('./Paper/counts/new_custom_counts/CTRL_10min-custom_assigned', 
           header = T, sep = '\t', row.names = 1, skip = 1, stringsAsFactors = F)
CTRL_mat_meta<- read.table('./Paper/counts/new_custom_counts/CTRL_10min_mature-custom_assigned', 
                           header = T, sep = '\t', row.names = 1, skip = 1, stringsAsFactors = F)

CTRL_new_meta <- CTRL_new_meta[grep(x = rownames(CTRL_new_meta), pattern = 'Distal|Proximal', invert = T),]
CTRL_mat_meta <- CTRL_mat_meta[grep(x = rownames(CTRL_mat_meta), pattern = 'Distal|Proximal', invert = T),]
# transform id of Genic ####

g2s <- fread('../IntergenicTranscription-master/IAA_gffcompare/GencodeReference/g2s_vm22_gencode.txt',header = F,data.table = F) #disable data.table mode
colnames(g2s) <- c("geneid","symbol")

ids <- data.frame(geneid=rownames(CTRL_new_meta), symbol = 'NULL')
ids <- ids[ids$geneid %in% g2s$geneid,] #
ids$symbol <- g2s[match(ids$geneid,g2s$geneid),2] #
CTRL_new_meta$gene <- rownames(CTRL_new_meta)
CTRL_new_meta$gene[CTRL_new_meta$gene %in% ids[!duplicated(ids$symbol),]$geneid] <- ids[!duplicated(ids$symbol),]$symbol
CTRL_new_meta<- CTRL_new_meta%>%distinct(gene, .keep_all = TRUE)
rownames(CTRL_new_meta) <- CTRL_new_meta$gene
CTRL_new_meta <- CTRL_new_meta%>%select(-gene)

# filter out no expressed features ###

CTRL_new_meta<- CTRL_new_meta[CTRL_new_meta$processed_data.R2_bam.CTRL_10min_R2_Aligned.sortedByCoord.out.bam > 0,]

CTRL_10min <- as.matrix(Genic_seu[, Genic_seu$orig.ident == 'CTRL_10min_beads']@assays$RNA@counts)
#CTRL_10min_mat <- as.matrix(Genic_seu[, Genic_seu$orig.ident == 'CTRL_10min_mature']@assays$RNA@counts)

gene_median_UMI<- apply(CTRL_10min, 1, median)
gene_median_UMI<- gene_median_UMI[gene_median_UMI > 0]

tab<- data.frame(UMI = gene_median_UMI, 
                 Length = CTRL_new_meta[match(names(gene_median_UMI),rownames(CTRL_new_meta)),]$Length,
                 Bulk_sum = rowSums(CTRL_10min[names(gene_median_UMI),]))

gene_median_UMI2<- apply(CTRL_10min_mat, 1, median)
gene_median_UMI2 <- gene_median_UMI2[gene_median_UMI2 > 0]
tab2<- data.frame(UMI = gene_median_UMI2, 
                 Length = CTRL_mat_meta[match(names(gene_median_UMI2),rownames(CTRL_mat_meta)),]$Length,
                 Bulk_sum = rowSums(CTRL_10min_mat[names(gene_median_UMI2),]))

p1 <- ggplot(tab, aes(x=Length, y=UMI)) + # N = 1418
  geom_jitter(color = 'black', alpha = 0.5)+
  geom_smooth(color = 'red', lwd = 1.5)+
  xlab('Feature length')+
  ylab('Median UMI')+scale_x_log10()+scale_y_log10()+
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

p2 <- ggplot(tab, aes(x=Length, y=Bulk_sum)) + 
  geom_jitter(color = 'black', alpha = 0.5)+
  geom_smooth(color = 'red', lwd = 1.5)+
  xlab('Feature length')+
  ylab('Sum of UMI')+scale_x_log10()+scale_y_log10()+
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

p3 <- ggplot(tab2, aes(x=Length, y=UMI)) + # CTRL mature N = 2464
  geom_jitter(color = 'black', alpha = 0.5)+
  geom_smooth(color = 'red', lwd = 1.5)+
  xlab('Feature length')+
  ylab('Median UMI')+scale_x_log10()+scale_y_log10()+
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

p4 <- ggplot(tab2, aes(x=Length, y=Bulk_sum)) + 
  geom_jitter(color = 'black', alpha = 0.5)+
  geom_smooth(color = 'red', lwd = 1.5)+
  xlab('Feature length')+
  ylab('Sum of UMI')+scale_x_log10()+scale_y_log10()+
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

p <- plot_grid(p1,p2, nrow = 2)
topptx(p, filename = './Paper/Figures/Fig1/newRNA_UMI_vs_featureLength.pptx')
p <- plot_grid(p3,p4, nrow = 2)
topptx(p, filename = './Paper/Figures/Fig1/oldRNA_UMI_vs_featureLength.pptx')

# plot exon/intron, intergenic percentage (bulk)

df<- data.frame(intergenic  = c(0.16, 0.085, 0.178, 0.105),
           intron = c(0.524,0.565,0.36,0.394),
           exon = c(0.177, 0.204,0.304,0.347),
           Sample = c('IAA_10min_beads','CTRL_10min_beads','IAA_10min_mature','CTRL_10min_mature'))

df2<- gather(df, key = 'Regions', value = "Percentage", -Sample)
df2$Regions <- factor(df2$Regions, levels = c('exon','intron','intergenic'))

p <- ggplot(data=df2, aes(x=Sample, y=Percentage, fill=Regions)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

topptx(p, filename = './Paper/Figures/Fig1/Exon_intron_intergenic_percent.pptx')

#[1] GRCh38 0.01612265 in beads 0.6979718 in mature

df2 <- data.frame(GRCh38 = c(0.016,0.698), mm10 = c(1-0.016,1-0.698), Group = c('beads', 'supernatant'))

df2<- gather(df2, key = 'Genome', value = 'Percent',-c('Group'))

p <- ggplot(data=df2, aes(x=Group, y=Percent, fill=Genome)) +
  geom_bar(stat="identity") + 
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

topptx(p, filename = './Paper/Figures/Fig1/contamination_percent.pptx')

# single cell level exon and intron counts ####

IAA_10min<- fread('./Paper/counts/constitutiveExon/IAA_10min.constExon_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
IAA_10min<- as.data.frame(IAA_10min)
rownames(IAA_10min) <- IAA_10min$gene
IAA_10min<- IAA_10min%>%dplyr::select(-gene)
IAA_10min <- colSums(IAA_10min)
IAA_10min <- IAA_10min[IAA_barcode]

CTRL_10min<- fread('./Paper/counts/constitutiveExon/CTRL_10min.constExon_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
CTRL_10min<- as.data.frame(CTRL_10min)
rownames(CTRL_10min) <- CTRL_10min$gene
CTRL_10min<- CTRL_10min%>%dplyr::select(-gene)
CTRL_10min <- colSums(CTRL_10min)
CTRL_10min <- CTRL_10min[CTRL_barcode]

IAA_10min_m<- fread('./Paper/counts/constitutiveExon/IAA_10min_mature.constExon_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
IAA_10min_m<- as.data.frame(IAA_10min_m)
rownames(IAA_10min_m) <- IAA_10min_m$gene
IAA_10min_m<- IAA_10min_m%>%dplyr::select(-gene)
IAA_10min_m <- colSums(IAA_10min_m)
IAA_10min_m <- IAA_10min_m[IAA_barcode]

CTRL_10min_m<- fread('./Paper/counts/constitutiveExon/CTRL_10min_mature.constExon_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
CTRL_10min_m<- as.data.frame(CTRL_10min_m)
rownames(CTRL_10min_m) <- CTRL_10min_m$gene
CTRL_10min_m<- CTRL_10min_m%>%dplyr::select(-gene)
CTRL_10min_m <- colSums(CTRL_10min_m)
CTRL_10min_m <- CTRL_10min_m[CTRL_barcode]

# for consti intron
if(FALSE){
  IAA_10min_I<- fread('./Paper/counts/constitutiveIntron/IAA_10min.constIntron_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  IAA_10min_I<- as.data.frame(IAA_10min_I)
  rownames(IAA_10min_I) <- IAA_10min_I$gene
  IAA_10min_I<- IAA_10min_I%>%dplyr::select(-gene)
  IAA_10min_I <- colSums(IAA_10min_I)
  IAA_10min_I <- IAA_10min_I[IAA_barcode]
  
  CTRL_10min_I<- fread('./Paper/counts/constitutiveIntron/CTRL_10min.constIntron_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  CTRL_10min_I<- as.data.frame(CTRL_10min_I)
  rownames(CTRL_10min_I) <- CTRL_10min_I$gene
  CTRL_10min_I<- CTRL_10min_I%>%dplyr::select(-gene)
  CTRL_10min_I <- colSums(CTRL_10min_I)
  CTRL_10min_I <- CTRL_10min_I[CTRL_barcode]
  
  IAA_10min_Im<- fread('./Paper/counts/constitutiveIntron/IAA_10min_mature.constIntron_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  IAA_10min_Im<- as.data.frame(IAA_10min_Im)
  rownames(IAA_10min_Im) <- IAA_10min_Im$gene
  IAA_10min_Im<- IAA_10min_Im%>%dplyr::select(-gene)
  IAA_10min_Im <- colSums(IAA_10min_Im)
  IAA_10min_Im <- IAA_10min_Im[IAA_barcode]
  
  CTRL_10min_Im<- fread('./Paper/counts/constitutiveIntron/CTRL_10min_mature.constIntron_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  CTRL_10min_Im<- as.data.frame(CTRL_10min_Im)
  rownames(CTRL_10min_Im) <- CTRL_10min_Im$gene
  CTRL_10min_Im<- CTRL_10min_Im%>%dplyr::select(-gene)
  CTRL_10min_Im <- colSums(CTRL_10min_Im)
  CTRL_10min_Im <- CTRL_10min_Im[CTRL_barcode]
}

# for intergenic ###

if(FALSE){
  IAA_10min_<- fread('./Paper/counts/new_custom_counts/IAA_10min.custom_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  IAA_10min_<- as.data.frame(IAA_10min_)
  rownames(IAA_10min_) <- IAA_10min_$gene
  idx <- grep(x = IAA_10min_$gene, pattern = 'Distal|Proximal')
  IAA_10min_<- IAA_10min_[idx,]%>%dplyr::select(-gene)
  IAA_10min_ <- colSums(IAA_10min_)
  IAA_10min_ <- IAA_10min_[IAA_barcode]
  
  CTRL_10min_<- fread('./Paper/counts/new_custom_counts/CTRL_10min.custom_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  CTRL_10min_<- as.data.frame(CTRL_10min_)
  rownames(CTRL_10min_) <- CTRL_10min_$gene
  idx <- grep(x = CTRL_10min_$gene, pattern = 'Distal|Proximal')
  CTRL_10min_<- CTRL_10min_[idx,]%>%dplyr::select(-gene)
  CTRL_10min_ <- colSums(CTRL_10min_)
  CTRL_10min_ <- CTRL_10min_[CTRL_barcode]
  
  IAA_10min_m_<- fread('./Paper/counts/new_custom_counts/IAA_10min_mature.custom_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  IAA_10min_m_<- as.data.frame(IAA_10min_m_)
  rownames(IAA_10min_m_) <- IAA_10min_m_$gene
  idx <- grep(x = IAA_10min_m_$gene, pattern = 'Distal|Proximal')
  IAA_10min_m_<- IAA_10min_m_[idx,]%>%dplyr::select(-gene)
  IAA_10min_m_ <- colSums(IAA_10min_m_)
  IAA_10min_m_ <- IAA_10min_m_[IAA_barcode]
  
  CTRL_10min_m_<- fread('./Paper/counts/new_custom_counts/CTRL_10min_mature.custom_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  CTRL_10min_m_<- as.data.frame(CTRL_10min_m_)
  rownames(CTRL_10min_m_) <- CTRL_10min_m_$gene
  idx <- grep(x = CTRL_10min_m_$gene, pattern = 'Distal|Proximal')
  CTRL_10min_m_<- CTRL_10min_m_[idx,]%>%dplyr::select(-gene)
  CTRL_10min_m_ <- colSums(CTRL_10min_m_)
  CTRL_10min_m_ <- CTRL_10min_m_[CTRL_barcode]
}

sc_count<- data.frame(Exon = c(CTRL_10min, CTRL_10min_m, IAA_10min, IAA_10min_m),
           Intron = c(CTRL_10min_I, CTRL_10min_Im, IAA_10min_I, IAA_10min_Im),
           Intergenic = c(CTRL_10min_, CTRL_10min_m_, IAA_10min_, IAA_10min_m_),
           Sample = c(rep('CTRL new', length(CTRL_barcode)), 
                      rep('CTRL mature', length(CTRL_barcode)), 
                      rep('IAA new', length(IAA_barcode)), 
                      rep('IAA mature', length(IAA_barcode))),
           barcode = c(rep(CTRL_barcode,2), rep(IAA_barcode,2))
           )

sc_count<- sc_count%>%group_by(Sample, barcode)%>%
  summarise(Exon_r = Exon/sum(Exon,Intron, Intergenic),
            Intron_r = Intron/sum(Exon,Intron, Intergenic),
            Intergenic_r = Intergenic/sum(Exon,Intron, Intergenic))


sc_count<- sc_count%>%gather(key = 'Regions', value = 'ratio', -c(Sample, barcode))
sc_count$Sample <- factor(sc_count$Sample, levels = c('CTRL new','CTRL mature', 'IAA new', 'IAA mature'))
sc_count$Regions <- factor(sc_count$Regions, levels = c('Exon_r','Intron_r', 'Intergenic_r'))

p1 <- ggplot(sc_count, aes(x=Regions, y=ratio, color = Sample)) + 
  geom_violin(trim=FALSE)+ #scale_y_log10() + 
  #geom_hline(yintercept = median(CV2_ctrl), color = 'grey', linetype=2)+
  scale_color_manual(values = c('CTRL new' = '#999999',
                                'CTRL mature' = "#009999",
                                'IAA new' = "#B2182B",
                                'IAA mature' = "#D6604D"))+
  geom_boxplot(width=0.3, position = position_dodge(width = 0.9), fill = 'white')+
  theme_classic()+
  labs(x = 'Regions',y = "Counts ratio")+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))# + 
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif')

  # single cell level annotated, distal, UpG,DoG,Linker counts ####

  IAA_10min<- fread('./Paper/counts/new_custom_counts/IAA_10min.custom_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  IAA_10min<- as.data.frame(IAA_10min)
  rownames(IAA_10min) <- IAA_10min$gene
  IAA_10min<- IAA_10min%>%dplyr::select(-gene)
  Distals <- str_replace(Distals, pattern = '-',replacement = '_')
  IAA_10min_dT <- colSums(IAA_10min[rownames(IAA_10min)%in% Distals,True_IAA])
  IAA_10min_dF <- colSums(IAA_10min[rownames(IAA_10min)%in% Distals,False_IAA])
  
  IAA_10min_pT <- colSums(IAA_10min[rownames(IAA_10min)%in% proximal_refseq_stringent_refined$gene_id,True_IAA])
  IAA_10min_pF <- colSums(IAA_10min[rownames(IAA_10min)%in% proximal_refseq_stringent_refined$gene_id,False_IAA])
  
  new.anno.genes_refadd[new.anno.genes_refadd$Set == 'refseq']$Set <- 'Annotated'
  IAA_10min_aT <- colSums(IAA_10min[rownames(IAA_10min)%in% new.anno.genes_refadd[new.anno.genes_refadd$Set == 'Annotated',]$gene_id,True_IAA])
  IAA_10min_aF <- colSums(IAA_10min[rownames(IAA_10min)%in% new.anno.genes_refadd[new.anno.genes_refadd$Set == 'Annotated',]$gene_id,False_IAA])
  
  IAA_10min_m<- fread('./Paper/counts/new_custom_counts/IAA_10min_mature.custom_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
  IAA_10min_m<- as.data.frame(IAA_10min_m)
  rownames(IAA_10min_m) <- IAA_10min_m$gene
  IAA_10min_m<- IAA_10min_m%>%dplyr::select(-gene)

  IAA_10min_mdT <- colSums(IAA_10min_m[rownames(IAA_10min_m)%in% Distals,True_IAA])
  IAA_10min_mdF <- colSums(IAA_10min_m[rownames(IAA_10min_m)%in% Distals,False_IAA])
  
  IAA_10min_mpT <- colSums(IAA_10min_m[rownames(IAA_10min_m)%in% proximal_refseq_stringent_refined$gene_id,True_IAA])
  IAA_10min_mpF <- colSums(IAA_10min_m[rownames(IAA_10min_m)%in% proximal_refseq_stringent_refined$gene_id,False_IAA])

  IAA_10min_maT <- colSums(IAA_10min_m[rownames(IAA_10min_m)%in% new.anno.genes_refadd[new.anno.genes_refadd$Set == 'Annotated',]$gene_id,True_IAA])
  IAA_10min_maF <- colSums(IAA_10min_m[rownames(IAA_10min_m)%in% new.anno.genes_refadd[new.anno.genes_refadd$Set == 'Annotated',]$gene_id,False_IAA])
  
  
  
  sc_count<- data.frame(Annotated = c(IAA_10min_aF, IAA_10min_maF, IAA_10min_aT, IAA_10min_maT),
                        Proximal = c(IAA_10min_pF, IAA_10min_mpF, IAA_10min_pT, IAA_10min_mpT),
                        Intergenic = c(IAA_10min_dF, IAA_10min_mdF, IAA_10min_dT, IAA_10min_mdT),
                        Sample = c(rep('CTRL new', length(False_IAA)), 
                                   rep('CTRL mature', length(False_IAA)), 
                                   rep('IAA new', length(True_IAA)), 
                                   rep('IAA mature', length(True_IAA))),
                        barcode = c(rep(False_IAA,2), rep(True_IAA,2))
  )
  
  sc_count<- sc_count%>%group_by(Sample, barcode)%>%
    summarise(Annotated_r = Annotated/sum(Annotated,Proximal, Intergenic),
              Proximal_r = Proximal/sum(Annotated,Proximal, Intergenic),
              Intergenic_r = Intergenic/sum(Annotated,Proximal, Intergenic))
  
  
  sc_count<- sc_count%>%gather(key = 'Regions', value = 'ratio', -c(Sample, barcode))
  sc_count$Sample <- factor(sc_count$Sample, levels = c('CTRL new','CTRL mature', 'IAA new', 'IAA mature'))
  sc_count$Regions <- factor(sc_count$Regions, levels = c('Annotated_r','Proximal_r', 'Intergenic_r'))

  
  p1 <- ggplot(sc_count[sc_count$Regions == 'Intergenic_r',], aes(x=Regions, y=ratio, color = Sample)) + 
    geom_violin(trim=FALSE)+ scale_y_log10() + 
    #geom_hline(yintercept = median(CV2_ctrl), color = 'grey', linetype=2)+
    scale_color_manual(values = c('CTRL new' = '#999999',
                                  'CTRL mature' = "#009999",
                                  'IAA new' = "#B2182B",
                                  'IAA mature' = "#D6604D"))+
    geom_boxplot(width=0.1, position = position_dodge(width = 0.9), fill = 'white')+
    theme_classic()+
    labs(x = 'Regions',y = "Counts ratio")+
   # facet_zoom(y = ratio <= 0.02, zoom.size = 1, split = T) +
    #facet_zoom(xy = Regions == 'Intergenic_r',split = TRUE)+
    theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
          axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
          legend.text = element_text(size = 18, lineheight = 2),
          plot.title = element_text(color="black", face="bold"))# + 
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif')
  














  
# single cell coverage #####
idx <- grep(x = new.anno.genes_refadd$gene_id, pattern = 'Distal')
Distal_gr<- new.anno.genes_refadd[idx]
Distal_gr <- Distal_gr[!Distal_gr$gene_id %in% Distals]

new.anno.genes <- new.anno.genes_refadd[!new.anno.genes_refadd$gene_id %in% Distal_gr$gene_id]
new.anno.genes <- new.anno.genes[!is.na(new.anno.genes$Subset)]
names(gr_list) <- unlist(lapply(strsplit(x, split = '_'),function(i)i[[1]]))
#new.anno.genes[new.anno.genes$Set == 'refseq']$Set <- 'Annotated'

IAA_list <- gr_list[True_IAA]
CTRL_list <- gr_list[False_IAA]
seq_cov <- list()
for (n in False_IAA) {
  gr_tmp <- CTRL_list[[n]]
  or<- findOverlaps(gr_tmp, new.anno.genes)
  idx <- as.integer(names(which(table(or@from)==1)))
  or<- or[or@from %in% idx]
  gr_tmp <- gr_tmp[or@from]
  gr_tmp$Set <- new.anno.genes[or@to]$Set
  gr_tmp$Subset <- new.anno.genes[or@to]$Subset
  gr_tmp <- as.data.frame(gr_tmp)
  Total_cov <- sum(gr_tmp$width)
  seq_cov[[n]]<- gr_tmp%>%group_by(Subset)%>%summarise(cov = sum(width))%>%mutate(perc = cov/Total_cov)
}
  
dt_CTRL<- do.call(rbind, seq_cov)
dt_CTRL$cell_barcode <- unlist(lapply(str_split(rownames(dt_CTRL),pattern = '\\.'), FUN = function(x)x[[1]]))
dt_CTRL$FeatureType <- ifelse(dt_CTRL$Subset %in% c("protein_coding", "long_ncRNA", "ncRNA", "other annotated","refseq"),'Annotated','Discovered')

factBiotype2 = c("UpstreamOfGene", "LinkerOfGene", "DownstreamOfGene", "< 10Kb", "> 10Kb",
                 "protein_coding", "long_ncRNA", "ncRNA", "other annotated")
factBiotype3 <- c(factBiotype2,'refseq')
dt_CTRL$Subset <- factor(dt_CTRL$Subset, levels=factBiotype3)
dt$Sample <- 'IAA 1hr'

tab<- rbind(dt_CTRL,dt)

ggplot(tab, aes(x=Subset, y=perc, color=Sample)) +
  ggtitle("Genome coverage of annotated and unannotated transcribed regions") +
  #geom_violin(trim=FALSE, position=position_dodge(width = 0.9))+ #scale_y_log10() + 
  #ylim(c(10,200))+
  scale_color_manual(values = c('CTRL' = "#009999",
                                'IAA 1hr' = "#D6604D"))+
  # scale_color_manual(values = c('CTRL' = "#999999",
  #                               'IAA 1hr' = "#B2182B"))+
  geom_boxplot(width=0.1, position = position_dodge(width = 0.9))+
  theme_classic()+
  scale_x_discrete("") +
  scale_y_continuous("Percentage of genome covered") +
  # stat_compare_means(comparisons = list(c('CTRL','IAA 1hr')), label = 'p.signif', data = tab[tab$Subset == '> 10Kb'])+
  #scale_fill_manual(values=setNames(c("#009999", "#D6604D"), c("Annotated", "Discovered"))) +
  # scale_fill_tableau(palette = "Classic 10 Medium", guide=FALSE) +
  facet_zoom(y = perc <= 0.1) +
  theme(legend.position=c(0.9, 0.85), axis.text=element_text(size=12), axis.title=element_text(size=14), axis.text.x=element_text(angle=45, hjust=1, vjust=1))

# correlation between annotated protein coding and intergenic Txn ####

reIdx<- unlist(lapply(seq_cov_IAA, function(x)nrow(x)))
dt$Cell <- rep(x = True_IAA, times = as.integer(reIdx))
reIdx<- unlist(lapply(seq_cov, function(x)nrow(x)))
dt_CTRL$Cell <- rep(x = False_IAA, times = as.integer(reIdx))

IAA_nc<- unlist(lapply(seq_cov_IAA, function(x)c(sum(x[x$Subset %in% c('< 10Kb','> 10Kb','DownstreamOfGene','LinkerOfGene','long_ncRNA',
                                               'UpstreamOfGene','ncRNA'),]$cov))))
IAA_nc<- unlist(lapply(seq_cov_IAA, function(x)c(sum(x[x$Subset %in% c('> 10Kb'),]$cov))))
IAA_nc<- unlist(lapply(seq_cov_IAA, function(x)c(sum(x[x$Subset %in% c('< 10Kb'),]$cov))))
IAA_nc<- unlist(lapply(seq_cov_IAA, function(x)c(sum(x[x$Subset %in% c('LinkerOfGene'),]$cov))))
IAA_nc<- unlist(lapply(seq_cov_IAA, function(x)c(sum(x[x$Subset %in% c('DownstreamOfGene'),]$cov))))
IAA_nc<- unlist(lapply(seq_cov_IAA, function(x)c(sum(x[x$Subset %in% c('long_ncRNA'),]$cov))))
IAA_c<- unlist(lapply(seq_cov_IAA, function(x)c(sum(x[x$Subset %in% c('protein_coding'),]$cov))))


nc_IAA <- data.frame(coding = IAA_c, noncoding = IAA_nc, cell = names(IAA_nc))

CTRL_nc<- unlist(lapply(seq_cov, function(x)c(sum(x[x$Subset %in% c('< 10Kb','> 10Kb','DownstreamOfGene','LinkerOfGene','long_ncRNA',
                                                                       'UpstreamOfGene','ncRNA'),]$cov))))
CTRL_nc<- unlist(lapply(seq_cov, function(x)c(sum(x[x$Subset %in% c('> 10Kb'),]$cov))))
CTRL_nc<- unlist(lapply(seq_cov, function(x)c(sum(x[x$Subset %in% c('< 10Kb'),]$cov))))
CTRL_nc<- unlist(lapply(seq_cov, function(x)c(sum(x[x$Subset %in% c('LinkerOfGene'),]$cov))))
CTRL_nc<- unlist(lapply(seq_cov, function(x)c(sum(x[x$Subset %in% c('DownstreamOfGene'),]$cov))))
CTRL_nc<- unlist(lapply(seq_cov, function(x)c(sum(x[x$Subset %in% c('long_ncRNA'),]$cov))))
CTRL_c<- unlist(lapply(seq_cov, function(x)c(sum(x[x$Subset %in% c('protein_coding'),]$cov))))


nc_CTRL <- data.frame(coding = CTRL_c, noncoding = CTRL_nc, cell = names(CTRL_nc))

nc_CTRL$noncoding_perc <- unlist(lapply(seq_cov, 
                                       function(x)c(sum(x[x$Subset %in% c('> 10Kb'),]$perc))))
nc_CTRL$coding_perc <- unlist(lapply(seq_cov, 
                                    function(x)c(sum(x[x$Subset %in% c('protein_coding'),]$perc))))

ggplot(nc_IAA, aes(x=coding, y=noncoding, fill = noncoding_perc)) + 
  #ylim(c(0,100))+
  xlab('protein coding sequence coverage (NCell = 972)') +
  ylab('> 10Kb sequence coverage')+geom_hline(yintercept = 25000, lty = 3)+
  #scale_color_grey() +
  scale_fill_gradient2(high = '#B2182B',low = 'white',limits = c(0,0.2), na.value ='white')+#labs(fill='coding_perc')+
  scale_x_log10() +
  scale_y_log10()+
  stat_cor(method = 'pearson', label.sep = '\n', size = 5)+
  geom_point(shape=21, alpha = 0.5, color = 'lightgrey', size = 4) + theme_classic() +
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

nc_IAA$states <- ifelse(nc_IAA$noncoding>25000, 'upper',
                        ifelse(nc_IAA$noncoding>0,'lower','none'))
nc_IAA$noncoding_perc <- unlist(lapply(seq_cov_IAA, 
                                       function(x)c(sum(x[x$Subset %in% c('> 10Kb'),]$perc))))
nc_IAA$coding_perc <- unlist(lapply(seq_cov_IAA, 
                                       function(x)c(sum(x[x$Subset %in% c('protein_coding'),]$perc))))


ggplot(nc_CTRL, aes(x=coding, y=noncoding, fill = noncoding_perc)) + 
  #ylim(c(0,100))+
  xlab('protein coding sequence coverage (NCell = 545)') +
  ylab('> 10Kb sequence coverage')+
  #scale_color_grey() +
  scale_fill_gradient2(high = '#D6604D')+labs(fill='coding_perc')+     # "#2166AC" for coding
  scale_x_log10() +geom_hline(yintercept = 25000, lty = 3)+
  scale_y_log10()+
  stat_cor(method = 'pearson', label.sep = '\n', size = 5)+
  geom_point(alpha = 0.5,color = 'black') + theme_classic() +
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

nc_CTRL$states <- ifelse(nc_CTRL$noncoding>25000, 'upper',
                        ifelse(nc_CTRL$noncoding>0,'lower','none'))
nc_CTRL$noncoding_perc <- unlist(lapply(seq_cov, 
                                       function(x)c(sum(x[x$Subset %in% c('> 10Kb'),]$perc))))
nc_CTRL$coding_perc <- unlist(lapply(seq_cov, 
                                    function(x)c(sum(x[x$Subset %in% c('protein_coding'),]$perc))))

save(seq_cov_IAA,seq_cov, IAA_list, CTRL_list, file = 'Paper/Figures/Fig1/seq_cov_TU.rda')

# coding potential analysis ####
old_TrueIAA$states <- ifelse(old_TrueIAA$barcode %in% nc_IAA[nc_IAA$states == 'upper',]$cell, 'upper',ifelse(old_TrueIAA$barcode %in% nc_IAA[nc_IAA$states == 'lower',]$cell, 'lower', 'none'))
table(old_TrueIAA$states)


qcparams <- c("nFeature_RNA", "nCount_RNA")
VlnPlot(object = new_FalseIAA[,!new_FalseIAA$states == 'none'], features = qcparams[1], group.by = "states", pt.size = 0, log = T, 
        cols = c("#009999","#D6604D")) + NoLegend()

new_falseMETA <- new_FalseIAA@meta.data
new_TrueMETA <- new_TrueIAA@meta.data
old_falseMETA <- old_FalseIAA@meta.data
old_TrueMETA <- old_TrueIAA@meta.data

ggplot(new_falseMETA[new_falseMETA$states != 'none',], aes(x=states, y=nCount_RNA, color = states)) + 
  geom_violin(trim=FALSE)+ #scale_y_log10() + 
  #ylim(c(10,200))+
  scale_color_manual(values = c('lower' = "#009999",
                                'upper' = "#D6604D"))+
  geom_boxplot(width=0.1, position = position_dodge(width = 0.9), fill = 'white')+
  theme_classic()+
  stat_compare_means(comparisons = list(c('lower','upper')), label = 'p.signif')+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"), legend.position = 'none')

# expression table, cytotrace ####
TrueIAA_results <- CytoTRACE(TrueIAA_tab)
pheno<- new_TrueIAA[,match(True_IAA, new_TrueIAA$barcode)]$states
names(pheno) <- unlist(lapply(str_split(names(pheno),pattern = '_'), FUN = function(x)x[[4]]))
plotCytoTRACE(TrueIAA_results, phenotype = pheno)
truecyto<- TrueIAA_results$CytoTRACE
nc_IAA$cytoTrace<- truecyto[match(nc_IAA$cell,names(truecyto))]

ggplot(nc_IAA[nc_IAA$states != 'none',], aes(x=states, y=cytoTrace, color = states)) + 
  geom_violin(trim=FALSE)+ #scale_y_log10() + 
  #ylim(c(10,200))+
  scale_color_manual(values = c('lower' = "#009999",
                                'upper' = "#D6604D"))+
  geom_boxplot(width=0.1, position = position_dodge(width = 0.9), fill = 'white')+
  theme_classic()+
  stat_compare_means(comparisons = list(c('lower','upper')), label = 'p.signif')+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"), legend.position = 'none')

# bulk pro-seq and TT-seq coverage ####

TTseq<- read.table('Paper/rds/PRO_TTseq_assembly/TTseq_bamtobed_union_merged_noDuplicates.bed')
TTseq <- TTseq[grep(x = TTseq$V1,pattern = 'chr'),]
TTseq.plus<- TTseq[TTseq$V4 =='+',]

TTseq.plus.gr<- reduce(GRanges(seqnames = TTseq.plus$V1, IRanges(TTseq.plus$V2, TTseq.plus$V3),TTseq.plus$V4))
sum(width(TTseq.plus.gr))/genome.length # 56.1%

TTseq.minus<- TTseq[TTseq$V4 =='-',]
TTseq.minus.gr<- reduce(GRanges(seqnames = TTseq.minus$V1, IRanges(TTseq.minus$V2, TTseq.minus$V3),TTseq.minus$V4))
sum(width(TTseq.minus.gr))/genome.length # 55.7%

TTseq.gr<- reduce(GRanges(seqnames = TTseq$V1, IRanges(TTseq$V2, TTseq$V3),TTseq$V4), ignore.strand=T)
sum(width(TTseq.gr))/genome.length # 74.9%

PROseq<- read.table('Paper/rds/PRO_TTseq_assembly/pro-seq_bamtobed_union_merged_noDuplicates.bed')
PROseq <- PROseq[grep(x = PROseq$V1,pattern = 'chr'),]
PROseq.plus<- PROseq[PROseq$V4 =='+',]

PROseq.plus.gr<- reduce(GRanges(seqnames = PROseq.plus$V1, IRanges(PROseq.plus$V2, PROseq.plus$V3),PROseq.plus$V4))
sum(width(PROseq.plus.gr))/genome.length # 12.6%

PROseq.minus<- PROseq[PROseq$V4 =='-',]
PROseq.minus.gr<- reduce(GRanges(seqnames = PROseq.minus$V1, IRanges(PROseq.minus$V2, PROseq.minus$V3),PROseq.minus$V4))
sum(width(PROseq.minus.gr))/genome.length # 13.5%

PROseq.gr<- reduce(GRanges(seqnames = PROseq$V1, IRanges(PROseq$V2, PROseq$V3),PROseq$V4), ignore.strand=T)
sum(width(PROseq.gr))/genome.length # 21.5%


cov_count<- data.frame(sc_cov = Perc_list_total, UMI_count = meta[match(names(Perc_list_total), meta$barcode),]$nCount_RNA)

ggplot(cov_count, aes(x=UMI_count, y=sc_cov)) + 
  #ylim(c(0,100))+
  xlab('UMI count (NCell = 550)') +
  ylab('sequence coverage')+
  #scale_color_grey() +
  #scale_fill_gradient2(high = '#B2182B',low = 'white',limits = c(0,0.2), na.value ='white')+#labs(fill='coding_perc')+
  #scale_x_log10() +
  #scale_y_log10()+
  stat_cor(method = 'pearson', label.sep = '\n', size = 5)+
  geom_point(alpha = 0.5, size = 2) + theme_classic() +
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))  

nc_CTRL$states2 <- ifelse(nc_CTRL$noncoding_perc>0.005, 'upper',
                         ifelse(nc_CTRL$noncoding>0,'lower','none'))


nc_CTRL$cc <- factor(nc_CTRL$cc, levels = c('g1s_1','g1s_2','g1s_3','s_1','s_2','s_3','g2m_1','g2m_2','g2m_3'))
ggplot(data=nc_CTRL[nc_CTRL$states !='none' & !is.na(nc_CTRL$cc),], aes(x=cc, fill=states2)) +
  geom_bar(position = 'fill')+
  scale_fill_manual(values=setNames(c("#009999", "#D6604D"), 
                                    c("lower", "upper"))) +
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))




