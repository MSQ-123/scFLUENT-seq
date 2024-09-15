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


#### setParameters ####
# Whether or not to regenerate all the objects
REBUILD = FALSE

# Define the species
#species = "Homo_sapiens"
species = "Mus_musculus"

# Define the distance to consider two regions and separated
#max.gapwidth = 1e3 # was 2e2
max.gapwidth = 2e2

# Define the distance to consider a region distal rather than proximal
min.gapwidth = 1e4

# Define the minimum width for a novel region
#min.frag.width = 1e3
min.frag.width = 0
#min.feature.width = 1e3
min.feature.width = 0

# Define the minimum TPM expression threshold
min.expr = 1 

# Define the quantile threshold for i) Removing lowly expressed fragments (frags) and ii) Adding unassigned regions (gaps) 
qnt_th = 0.1

# Define the factor levels for the subset
factSubset = c("GencodeGene", "UpstreamOfGene", "LinkerOfGene", "DownstreamOfGene", "< 10Kb", "> 10Kb")

# Define the biotype levels
factSubset2 = c("UpstreamOfGene", "LinkerOfGene", "DownstreamOfGene", "< 10Kb", "> 10Kb",
                "Protein_coding (expressed)", "LincRNA", "Protein_coding (not expressed)")

# Define the biotype levels
factBiotype = c("UpstreamOfGene", "LinkerOfGene", "DownstreamOfGene", "< 10Kb", "> 10Kb",
                "Protein_coding", "Long_non-coding", "Non-coding", "Other annotated")

# Define the biotype levels
factBiotype2 = c("UpstreamOfGene", "LinkerOfGene", "DownstreamOfGene", "< 10Kb", "> 10Kb",
                 "protein_coding", "long_ncRNA", "ncRNA", "other annotated")

# Define the colour set
colSet = setNames(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#EDC948"),
                  c("Protein coding", "Long ncRNA", "Small ncRNA", "Other annotated", "Gene-associated TU", "Independent TU"))

# Define the colour set
colSet2 = setNames(c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2", "#59A14F", "#59A14F", "#59A14F", "#EDC948"),
                   c("Protein coding", "Long ncRNA", "Small ncRNA", "Other annotated",
                     "UpstreamOfGene", "LinkerOfGene", "DownstreamOfGene", "Independent TU"))


### load anno ####
# Define the set of chromosome to use
chromosomes = genomeStyles()[["Mus_musculus"]]$UCSC
chromosomes = setdiff(chromosomes, c("chrY", "chrM"))
chromosomes.lengths = seqlengths(Mmusculus)
chromosomes.lengths = chromosomes.lengths[names(chromosomes.lengths)%in%chromosomes]

# Load the pre-processed GRanges and print the list
working.folder    = getwd()
setwd('../IntergenicTranscription-master//IAA_gffcompare/')
annotation.folder = file.path(working.folder, "GencodeReference")

load(file.path(annotation.folder, "mm10_GencodeM22_genes_annotations.RData"))
tx.regions <- readRDS(file.path(annotation.folder, "mm10_GencodeM22_annotations.all.genes.transcript.regions.rds"))

rm(exonsByTxs, txs, txs.metadata, txs.longest.pc.rna.granges, txs.longest.nc.rna.granges)

exons.pc.granges = updateObject(exons.pc.granges, verbose=TRUE)
exons.nc.granges = updateObject(exons.nc.granges, verbose=TRUE)

# Get the introns
introns = genes
names(introns) = NULL
introns = split(introns, introns$gene_id)
exons = c(exons.pc.granges, exons.nc.granges)
introns = GenomicRanges::setdiff(introns[match(names(exons), names(introns))],
                                 c(exons.pc.granges, exons.nc.granges))  

lnc.gff = import.gff3(file.path(annotation.folder, "gencode.vM22.long_noncoding_RNAs.gff3.gz"))
lnc.gff = keepSeqlevels(lnc.gff, chromosomes, pruning.mode="coarse")

# Remove the duplicated genes on the chromosome Y
x.genes = unique(lnc.gff[seqnames(lnc.gff)%in%"chrX"]$gene_id)
y.genes = unique(lnc.gff[seqnames(lnc.gff)%in%"chrY"]$gene_id)
y.genes = y.genes[y.genes%in%x.genes]
lnc.gff = lnc.gff[!(seqnames(lnc.gff)%in%"chrY" & lnc.gff$gene_id%in%y.genes)]

# Group gene biotypes into classes
protein.coding = gene.metadata[gene_type%in%c("protein_coding", paste(rep(c("IG","TR"),each=4), c("C","D","J","V"), "gene", sep="_"), "IG_LV_gene"), gene_id]
long.ncRNA = gene.metadata[gene_type%in%unique(names(sort(-table(lnc.gff$gene_type)))), gene_id]
ncRNA = c(setdiff(gene.metadata[grepl("RNA", gene_type), gene_id], long.ncRNA), "ribozyme")
other = setdiff(gene.metadata[, gene_id], c(ncRNA, long.ncRNA, protein.coding))

rm(lnc.gff)

tab.biotype = data.table(Gene_id = c(protein.coding, long.ncRNA, ncRNA, other),
                         Biotype = rep(c("protein_coding", "long_ncRNA", "ncRNA", "other annotated"),
                                       c(length(protein.coding), length(long.ncRNA), length(ncRNA), length(other))),
                         key="Gene_id")

entrez.genes = fread(cmd=paste("gunzip -c", file.path(annotation.folder, "gencode.vM22.metadata.EntrezGene.gz")))

geneset = keepSeqlevels(genes, chromosomes, pruning.mode="coarse")

TxDb = makeTxDbFromGFF(file = '../R_Gencode_Reference-master/Mus_musculus/mm10/mm10.ncbiRefSeq.gtf.gz',
                      format = "gtf",
                      dataSource = "UCSC",
                      organism = "Mus musculus")

# Genes, transcripts and exons
g <- genes(TxDb)
g <- keepStandardChromosomes(g, pruning.mode="coarse")


### load gtf ####
# Import the StringTie (or rather gffcompare) GFF output
#here I load the bed file

#gff<- fread('./total_bambed/downsampled/IAA-24h-wo-36h-3_bamtobed_union_merged_noDuplicates.bed', col.names=c("seqname", "start", "end", "strand"))
gff<- fread('../../scEU-seq/Paper/stringtie/IAA_CTRL_bamtobed_union_merged_noDuplicates.bed', col.names=c("seqname", "start", "end", "strand"))

gff[, feature := paste('read', seq_along(1:nrow(gff)), sep = '_')]
gff[, c(1,5, 2:4)]

# Create the final data.table and GRanges
gff.genes = gff[, list(seqname=unique(seqname), start=min(start), end=max(end), strand=unique(strand)), by="feature"]
gff.genes.gr = sort.GenomicRanges(with(gff.genes, GRanges(seqname, IRanges(start, end), strand, feature)), ignore.strand=TRUE)

gff.genes.gr = sort.GenomicRanges(with(gff.genes, GRanges(seqname, IRanges(start, end), strand, feature)))

rm(gff)


# Remove the Mitochondrial chromosome
gff.genes.gr = reduce(gff.genes.gr[!seqnames(gff.genes.gr)%in%"chrM",], min.gapwidth=max.gapwidth)

annotated.gr = geneset

mcols(annotated.gr) = data.frame(Set = factor("Annotated", levels=c("Annotated", "New")))
mcols(gff.genes.gr) = data.frame(Set = factor("New", levels=c("Annotated", "New")))
#refseq 
#mcols(g) = data.frame(Set = factor("Annotated", levels=c("Annotated", "New")))
gff.genes.gr<- keepStandardChromosomes(gff.genes.gr, pruning.mode="coarse")
### prepare refseq and ensembl genesets####

# Remove the intragenic fragments and separate proximal from distal
#I added refseq gtf to be more stringent



#if only consider transcribed protein coding regions 
if(FALSE){
  idx<- findOverlaps(intragenic, c(GRanges(annotated.dt), GRanges(g_add)))
  transcribed_pc.gr<-c(GRanges(annotated.dt), GRanges(g_add))[idx@to]
  transcribed_pc.gr<- unique(transcribed_pc.gr)
}


anno_g_add<- c(GRanges(annotated.dt), GRanges(g_add))
transcribed_pc_range<- intersect_ranges_directed(intragenic, anno_g_add)

idx<- findOverlaps(transcribed_pc_range, anno_g_add)
transcribed_pc_range <- as.data.frame(transcribed_pc_range[idx@from])
transcribed_pc_range$gene_id <- anno_g_add[idx@to,]$gene_id
transcribed_pc_range$Set <- anno_g_add[idx@to,]$Set
transcribed_pc_range$Subset <- anno_g_add[idx@to,]$Subset

new.anno.genes.cov = rbindlist(list(transcribed_pc_range[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")], 
                                proximal.dt[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")], 
                                intergenic.dt[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")]))

#add refseq only genes to make the reference more complete
add_refseq_gene<- keepSeqlevels(g[!overlapsAny(g, annotated.gr)], chromosomes, pruning.mode="coarse")

# Create data.table objects for each set
mcols(annotated.gr) = NULL
#mcols(g) = NULL

# transcribed_pc.dt = data.table(as.data.frame(transcribed_pc.gr), index=seq_along(transcribed_pc.gr), gene_id=names(transcribed_pc.gr))
# transcribed_pc.dt[, `:=`(Set = "Annotated", Subset = tab.biotype[match(gene_id, tab.biotype[, Gene_id]), Biotype])]

annotated.dt = data.table(as.data.frame(annotated.gr), index=seq_along(annotated.gr), gene_id=names(annotated.gr))
annotated.dt[, `:=`(Set = "Annotated", Subset = tab.biotype[match(gene_id, tab.biotype[, Gene_id]), Biotype])]

add_refseq_gene.dt = data.table(as.data.frame(add_refseq_gene), index=seq_along(add_refseq_gene), gene_id=names(add_refseq_gene))
add_refseq_gene.dt[, `:=`(Set = "Annotated", Subset = 'refseq_added')]

idx<- overlapsAny(g,GRanges(annotated.dt), maxgap=max.gapwidth)
g_add = g[!idx]

g_add$Set <- 'refseq'
g_add$Subset <- 'refseq'
intragenic = gff.genes.gr[overlapsAny(gff.genes.gr, c(annotated.gr,GRanges(g_add)), maxgap=max.gapwidth)]

proximal = setdiff(intragenic, c(annotated.gr,GRanges(g_add)))
intergenic = reduce(gff.genes.gr[!overlapsAny(gff.genes.gr, intragenic, maxgap=max.gapwidth)])


proximal.dt = data.table(as.data.frame(proximal), index=seq_along(proximal))
intergenic.dt = data.table(as.data.frame(intergenic), index=seq_along(intergenic))

# Find the genes preceding and/or following the region, then following genic_nearbyProximals_combine.R script ! ####


  # NB: precede() returns the index of the GRange that is preceded by (i.e., comes after) the query;
  # follow() returns the index of the GRange that is followed by (comes before) the query
  proximal.dt[, `:=`(PreviousGene=as.numeric(NA), NextGene=as.numeric(NA),
                     PreviousDistance=as.integer(NA), NextDistance=as.integer(NA))]
  
  # Add the indexes and distances of neighbouring genes
  # proximal.dt[, `:=`(PreviousGene = follow(proximal, geneset),
  #                    NextGene = precede(proximal, geneset))]
  proximal.dt[, `:=`(PreviousGene = follow(proximal, c(geneset)),
                   NextGene = precede(proximal, c(geneset)))]  

  proximal.dt[!is.na(PreviousGene), PreviousDistance := distance(proximal[index,], geneset[PreviousGene,])]

  proximal.dt[!is.na(NextGene), NextDistance := distance(proximal[index,], geneset[NextGene,])]

  
  # Assign additional columns to the proximal set
  proximal.dt[, `:=`(gene_id = paste("Proximal", seq_along(seqnames), sep="_"), Set = "Proximal", Subset = as.character(NA))]
  proximal.dt[!is.na(PreviousDistance) & PreviousDistance <= max.gapwidth | is.na(NextDistance), Subset := "DownstreamOfGene"]
  proximal.dt[!is.na(NextDistance) & NextDistance <= max.gapwidth | is.na(PreviousDistance), Subset := "UpstreamOfGene"]
  proximal.dt[PreviousDistance <= max.gapwidth & NextDistance <= max.gapwidth, Subset := "LinkerOfGene"]

  # Close the gaps with the nearest gene(s)
  proximal.dt<- na.omit(proximal.dt)
  proximal.dt[Subset%in%"DownstreamOfGene" & strand=="+", start := as.integer(end(geneset[PreviousGene,])+1)]
  proximal.dt[Subset%in%"DownstreamOfGene" & strand=="-", end := as.integer(start(geneset[PreviousGene,])-1)]
  proximal.dt[Subset%in%"UpstreamOfGene" & strand=="+", end := as.integer(start(geneset[NextGene,])-1)]
  proximal.dt[Subset%in%"UpstreamOfGene" & strand=="-", start := as.integer(end(geneset[NextGene,])+1)]
  proximal.dt[Subset%in%"LinkerOfGene" & strand=="+", `:=`(start = as.integer(end(geneset[PreviousGene,])+1), end = as.integer(start(geneset[NextGene,])-1))]
  proximal.dt[Subset%in%"LinkerOfGene" & strand=="-", `:=`(start = as.integer(end(geneset[NextGene,])+1), end = as.integer(start(geneset[PreviousGene,])-1))]
  
  # # Remove features that are below the minimum threshold
  proximal.dt = proximal.dt[width>=min.feature.width,]
  
  # Check that the linker of genes are actually connecting pairs of expressed genes
  linkers.dt = proximal.dt[Subset%in%"LinkerOfGene",]
  linkers.dt[, `:=`(NextID=names(geneset[NextGene]), PreviousID=names(geneset)[PreviousGene])]
  
  # proximal.dt[Subset%in%"LinkerOfGene", Subset := linkers.dt[, Subset]]
  
  # Remove unneeded columns
  proximal.dt[, `:=`(PreviousGene=names(geneset[PreviousGene]), NextGene=names(geneset[NextGene]))]
  saveRDS(proximal.dt, file = '../../scEU-seq/Paper/rds/bed_union_no1kb_threshold/proximal.dt.rds')
  #proximal.dt[, `:=`(PreviousGene=NULL, NextGene=NULL, PreviousDistance=NULL, NextDistance=NULL)]
  
  # Assign additional columns to the intergenic set
  intergenic.dt[, `:=`(gene_id = paste("Distal", seq_along(seqnames), sep="_"), Set = "Intergenic",
                       Subset = "NA")]
  
  idx<- distanceToNearest(intergenic, geneset)
  intergenic.dt[idx@from]$Subset <- ifelse(mcols(distanceToNearest(intergenic, geneset))$distance<min.gapwidth, "< 10Kb", "> 10Kb")
  
  intergenic.dt[Subset == 'NA',]$Subset <-  NA

  
  # # Remove features that are below the minimum threshold
  # intergenic.dt = intergenic.dt[width>=min.feature.width,]
  
  # begin genic_nearbyProximals_combine.R pipeline ####
  
  # Merge all annotation sets
  
  # proximal.dt<- proximal.dt[,c(1:5,8, 6,7,9)]
  # intergenic.dt<- intergenic.dt[,c(1:5,8, 6,7,9)]
  
  
  

  # new.anno.genes = rbindlist(list(annotated.dt[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")], 
  #                                 proximal.dt[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")], 
  #                                 intergenic.dt[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")]))
  
  # or 
  #g_add<- as.data.frame(g_add)[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")]
  new.anno.genes = rbindlist(list(annotated.dt[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")],
                                  as.data.frame(g_add)[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")],
                                  proximal.dt[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")], 
                                  intergenic.dt[,c("seqnames","start","end","width","strand", "gene_id","Set","Subset")]))
  #check
  table(new.anno.genes$Subset)
  
  # Create and sort the GRanges
  new.anno.genes = with(new.anno.genes, GRanges(seqnames, IRanges(start, end), strand, gene_id, Set, Subset))
  new.anno.genes = sort.GenomicRanges(new.anno.genes, ignore.strand=TRUE)
  #saveRDS(new.anno.genes, file = '../../scEU-seq/Paper/rds/new.anno.genes.rds')
  saveRDS(new.anno.genes, file = '../../scEU-seq/Paper/rds/bed_union_no1kb_threshold/new.anno.genes_refadd.rds')

# Merge all annotation sets ####
#export.bed(new.anno.genes, "annotatedRegions/bed_union_new.anno.genes_no_nearest_refseq_stringent.bed")

  
# Save the object
### previous write gtf file####
new.anno.df<- as.data.frame(new.anno.genes)
new.anno.df$source <- 'custom'
new.anno.df$type <- 'exon'
new.anno.df$score <- '.'
new.anno.df$phase <- '.'
new.anno.df$attributes <- paste('gene_id ', '"', new.anno.df$gene_id, '"; ',  
                                'transcript_id ',  '"', new.anno.df$gene_id, '"; ',
                                'Set_id ', '"', new.anno.df$Set, '"; ',
                                'Subset_id ', '"', new.anno.df$Subset, '"; ', sep = '')

new.anno.df <- new.anno.df[,c("seqnames", "source","type", "start", "end", "score", 
                              "strand", "phase", "attributes")]

# write.table(new.anno.df, file = '../../scEU-seq/Paper/rds/bed_union_no1kb_threshold/bed_union_new_anno_no_nearest_refseq_stringent_df.gtf', sep = '\t', quote = F,row.names = F, 
#             col.names = F)
# write.table(new.anno.df, file = '../../scEU-seq/Paper/rds/bed_union_no1kb_threshold/bed_union_new_anno_no_nearest_refseq_stringent_df_refadd.gtf', sep = '\t', quote = F,row.names = F, 
#             col.names = F)
write.table(new.anno.df, file = '../../scEU-seq/Paper/rds/bed_union_no1kb_threshold/bed_union_new_anno_no_nearest_refseq_stringent_df_refadd_gadd_first.gtf', sep = '\t', quote = F,row.names = F, 
            col.names = F)

#save(proximal.dt, intergenic.dt, new.anno.genes, file = 'GencodeReference/proximal_intergenic_new_refseq_stringent.Rda')

### plot #####
















