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
require(BSgenome.Mmusculus.UCSC.mm10)
require("ggbio")
require("GenomicAlignments")
require("genomation")
require("VennDiagram")
require("viridis")

species = "Mus musculus"
setwd("./R_Gencode_Reference-master/")
TxDb = makeTxDbFromGFF(file = './Mus_musculus/mm10/mm10.ncbiRefSeq.gtf.gz',
                       format = "gtf",
                       dataSource = "UCSC",
                       organism = species)

# Genes, transcripts and exons
g <- genes(TxDb)
g <- keepStandardChromosomes(g, pruning.mode="coarse")

gtf<- fread('../IAA_gffcompare/bed_union_new_anno_no_nearest_refseq_stringent_df_final_total.gtf', col.names=c("seqnames", "source","type", "start", "end", "score", 
                                                                                                               "strand", "phase", "attributes"))
#load("../IAA_gffcompare/GencodeReference/intergenic_stringent/proximal_intergenic_antisense_overlap_refseq_stringent.Rda")

table(proximal.dt$Subset)


#refine antisense overlap
library(plyranges)
tmp <- GRanges(proximal.dt)
#if the tmp overlap multiple times, it will occur more than 1 time in the output
tmp<- tmp %>% join_overlap_left(c(geneset, g))
table(is.na(tmp$gene_id.y))
#the proximal that overlaps geneid:
proximal_antisense<- unique(tmp[!is.na(tmp$gene_id.y),]$gene_id.x)
proximal.dt$antisense_overlap <- rep('FALSE', nrow(proximal.dt))
#antisense overlaps' nearest gene is not trustable:
proximal.dt[proximal.dt$gene_id %in% proximal_antisense, ]$antisense_overlap <- TRUE

#warning: the subset of proximal is inaccurate in refseq_bamtobed classification ()

proximal.dt[, `:=`(PreviousGene=as.numeric(NA), NextGene=as.numeric(NA),
                   PreviousDistance=as.integer(NA), NextDistance=as.integer(NA))]

# Add the indexes and distances of neighbouring genes

proximal = with(proximal.dt, GRanges(seqnames, IRanges(start, end), strand, gene_id, Set, Subset, antisense_overlap))
proximal.dt[, `:=`(PreviousGene=as.numeric(NA), NextGene=as.numeric(NA),
                   PreviousDistance=as.integer(NA), NextDistance=as.integer(NA))]


#saveRDS(geneset, file = './ensembl_geneset.rds')
#saveRDS(g, file = './refseq_geneset.rds')



proximal.dt = data.table(as.data.frame(proximal), index=seq_along(proximal))
proximal.dt[, `:=`(PreviousGene = follow(proximal, c(geneset, g)),
                   NextGene = precede(proximal, c(geneset,g)))]
#intergenic.dt = data.table(as.data.frame(intergenic), index=seq_along(intergenic))

#distance(proximal[seq_along(proximal[!is.na(proximal.dt$PreviousGene)]),], c(geneset, g)[proximal.dt[!is.na(proximal.dt$PreviousGene),]$PreviousGene])
#distance(proximal[seq_along(proximal[!is.na(proximal.dt$NextGene)]),], c(geneset, g)[proximal.dt[!is.na(proximal.dt$NextGene),]$NextGene])

proximal.dt[!is.na(PreviousGene), PreviousDistance := distance(proximal[seq_along(proximal[!is.na(proximal.dt$PreviousGene)]),], c(geneset, g)[proximal.dt[!is.na(proximal.dt$PreviousGene),]$PreviousGene])]

proximal.dt[!is.na(NextGene), NextDistance := distance(proximal[seq_along(proximal[!is.na(proximal.dt$NextGene)]),], c(geneset, g)[proximal.dt[!is.na(proximal.dt$NextGene),]$NextGene])]

#great!


# Assign additional columns to the proximal set

proximal.dt[!is.na(PreviousDistance) & PreviousDistance <= max.gapwidth | is.na(NextDistance), Subset := "DownstreamOfGene"]
proximal.dt[!is.na(NextDistance) & NextDistance <= max.gapwidth | is.na(PreviousDistance), Subset := "UpstreamOfGene"]
proximal.dt[PreviousDistance <= max.gapwidth & NextDistance <= max.gapwidth, Subset := "LinkerOfGene"]

# Close the gaps with the nearest gene(s)  the gap has been closed before


# # Remove features that are below the minimum threshold
proximal.dt = proximal.dt[width>=min.feature.width,]

head(proximal.dt)

# reassign paired gene id for DoG and UoG
proximal.dt[, `:=`(PairedGene=as.character(NA))]
proximal.dt[Subset == "DownstreamOfGene", PairedGene :=c(geneset, g)[proximal.dt[proximal.dt$Subset == "DownstreamOfGene",]$PreviousGene]$gene_id]

proximal.dt[proximal.dt$Subset == "UpstreamOfGene" & !is.na(proximal.dt$NextGene), PairedGene :=c(geneset, g)[proximal.dt[proximal.dt$Subset == "UpstreamOfGene" & !is.na(proximal.dt$NextGene),]$NextGene]$gene_id]

# for linker:

proximal.dt[proximal.dt$Subset == "LinkerOfGene", PairedGene :=paste(c(geneset, g)[proximal.dt[proximal.dt$Subset == "LinkerOfGene",]$PreviousGene]$gene_id,c(geneset, g)[proximal.dt[proximal.dt$Subset == "LinkerOfGene",]$NextGene]$gene_id, sep = ',')]
proximal.dt$antisense_overlapped_gene<- tmp[match(proximal.dt$gene_id, tmp$gene_id.x),]$gene_id.y


saveRDS(proximal.dt, file = '../../scEU-seq/Paper/rds/bed_union_no1kb_threshold/proximal_refseq_stringent_refined.rds')


tmp2 <- GRanges(intergenic.dt)
tmp2<- tmp2 %>% join_overlap_left(c(geneset, g))
table(is.na(tmp2$gene_id.y))
#the intergenics that overlaps geneid:
intergenic_antisense<- unique(tmp2[!is.na(tmp2$gene_id.y),]$gene_id.x) #79978 intergenics with geneid overlap
intergenic.dt$antisense_overlap <- rep('FALSE', nrow(intergenic.dt))
intergenic.dt[intergenic.dt$gene_id %in% intergenic_antisense, ]$antisense_overlap <- TRUE
intergenic.dt$antisense_overlapped_gene<- tmp2[match(intergenic.dt$gene_id, tmp2$gene_id.x),]$gene_id.y
#refine intergenic subset

tmp2 <- GRanges(intergenic.dt)
tmp<- nearest(tmp2, c(geneset,g))
intergenic.dt[tmp[!is.na(tmp)], `:=`(Subset = ifelse(distance(tmp2[!is.na(tmp)], c(geneset,g)[tmp[!is.na(tmp)]])<min.gapwidth, "< 10Kb", "> 10Kb"))]

# only keep > 1000 bp distal
intergenic.dt<- intergenic.dt[width>=min.feature.width,]
idx<- overlapsAny(GRanges(intergenic.dt), g, maxgap=max.gapwidth)
intergenic.dt = intergenic.dt[!idx]

saveRDS(intergenic.dt, '../../scEU-seq/Paper/rds/bed_union_no1kb_threshold/intergenic_refseq_stringent_refined.rds')




