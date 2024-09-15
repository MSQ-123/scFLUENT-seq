# sc
seq_covIAA <- list()
seq_countIAA <- list()
for (n in False_IAA) {
  gr_tmp <- CTRL_list[[n]]
  or<- findOverlaps(gr_tmp, new.comp)
  idx <- as.integer(names(which(table(or@from)==1)))
  or<- or[or@from %in% idx]
  gr_tmp <- gr_tmp[or@from]
  gr_tmp$Comp <- new.comp[or@to]$comp
  #gr_tmp<- GenomicRanges::reduce(gr_tmp, ignore.strand=T)
  gr_tmp <- as.data.frame(gr_tmp)
  #Total_cov <- sum(gr_tmp$width)
  Total_count <- nrow(gr_tmp)
  #seq_covIAA[[n]]<- gr_tmp%>%group_by(Comp)%>%summarise(cov = sum(width))%>%mutate(perc = cov/Total_cov)
  seq_count[[n]]<- gr_tmp%>%group_by(Comp)%>%summarise(count = length(strand))%>%mutate(perc = count/Total_count)
}

dt_CTRL<- do.call(rbind, seq_cov)
dt_CTRL$cell_barcode <- unlist(lapply(str_split(rownames(dt_CTRL),pattern = '\\.'), FUN = function(x)x[[1]]))

# bulk

bed<- do.call(rbind, bed_list)
# for count contrib do not reduce
bed.gr<- GenomicRanges::reduce(GRanges(seqnames = bed$V1, IRanges(bed$V2, bed$V3),bed$V6), ignore.strand=T)
sum(width(bed.gr))/genome.length

or<- findOverlaps(bed.gr, new.comp)
idx <- as.integer(names(which(table(or@from)==1)))
or<- or[or@from %in% idx]
bed.gr <- bed.gr[or@from]
bed.gr$Comp <- new.comp[or@to]$comp
bed.gr <- as.data.frame(bed.gr)
Total_cov <- sum(bed.gr$width)
bulk<- bed.gr%>%group_by(Comp)%>%summarise(cov = sum(width))%>%mutate(perc = cov/Total_cov)

# compare

ggplot(bulk_tab, aes(x=Comp, y=perc, fill=Sample)) +
  ggtitle("Annotated and unannotated transcribed regions") +
  geom_bar(stat="identity", position="dodge", col="black") +
  #geom_text(aes(x = Type, y = N+2e3, label = paste0(N, "\n(", Perc*100,"%)")), col="black") +
  scale_x_discrete("") + #scale_y_continuous("Count", limits=c(0, 24e3)) +
  scale_fill_manual(values=setNames(c("#009999", "#D6604D"), c("CTRL", "IAA"))) +
  # scale_fill_tableau(palette = "Classic 10 Medium", guide=FALSE) +
  theme_classic() + theme(legend.position=c(0.9, 0.9), axis.text=element_text(size=12), axis.title=element_text(size=14), axis.text.x=element_text(angle=45, hjust=1, vjust=1))

ggplot(dt_tab, aes(x=Comp, y=perc, color = Sample)) + 
  geom_violin(trim=FALSE)+ #scale_y_log10() + 
  #ylim(c(10,200))+
  scale_color_manual(values = c('CTRL' = "#009999",
                                'IAA' = "#D6604D"))+
  geom_boxplot(width=0.1, position = position_dodge(width = 0.9), fill = 'white')+
  theme_classic()+
  #stat_compare_means(comparisons = list(c('A','B')), label = 'p.signif')+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

# filter low features and calculate coverage

bed<- do.call(rbind, bed_list)
cov_fil2 <- c()
#ctrl5 <- ctrl[rowSums(ctrl[,] >0) > 5,]
for( i in 1:40){
  bed.gr<- GRanges(seqnames = bed$V1, IRanges(bed$V2, bed$V3),bed$V6)
  rm(tmp)
  gc()
  tmp <- ctrl[rowSums(ctrl[,] >0) > i,]
  # for count contrib do not reduce
  or<- findOverlaps(bed.gr, new.anno.genes_refadd_gadd_first[new.anno.genes_refadd_gadd_first$gene_id %in% rownames(tmp)])
  # idx <- as.integer(names(which(table(or@from)==1)))
  # or<- or[or@from %in% idx]
  bed.gr <- bed.gr[or@from]
  cov_fil2[i] <-  sum(width(GenomicRanges::reduce(bed.gr, ignore.strand=T)))/genome.length
}


# 85.4 UMI >3,  85.2 UMI > 4, 85.1 UMI > 5

ggplot(tab, aes(x=count, y=cov)) + 
  geom_jitter(alpha = 0.5)+
  geom_smooth(color = '#B2182B', lwd = 1)+
  xlab('UMI threshold')+
  ylab('Pseudobulk genome coverage')+
  scale_color_manual(values = c('IAA mix' = "black"))+
  theme_classic()+
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

# for protein coding #####

pc.gr <- new.anno.genes_refadd_gadd_first[new.anno.genes_refadd_gadd_first$Subset == 'protein_coding',]

or<- findOverlaps(pc.gr, new.comp)
idx <- as.integer(names(which(table(or@from)==1)))
or<- or[or@from %in% idx]
pc.gr <- pc.gr[or@from]
pc.gr$comp <- new.comp[or@to]$comp

pc <- intersect(pc.gr$gene_id, names(k_IAA_false_syn))
k_IAA_false_syn <- k_IAA_false_syn[names(k_IAA_false_syn) %in% pc]

comp <- pc.gr[match(names(k_IAA_false_syn),pc.gr$gene_id),]$comp
length = width(pc.gr[match(names(k_IAA_false_syn),pc.gr$gene_id),])

k_IAA_true_syn <- k_IAA_true_syn[match(names(k_IAA_false_syn), names(k_IAA_true_syn))]

dt <- data.frame(comp = comp, length = length, k_syn = k_IAA_false_syn)


for (i in unique(dt$comp)) {
  dens <- density(log10(dt[dt$comp == i,]$length))
  df <- data.frame(x=dens$x, y=dens$y)
  probs <- c(0, 0.25, 0.5, 0.75, 1)
  quantiles <- quantile(log10(dt[dt$comp == i,]$length), prob=probs)
  df$quant <- factor(findInterval(df$x,quantiles))
  df$comp <- i
  df_list[[i]] <- df
}
df_list<- do.call(rbind, df_list)
quantiles <- quantile(df_list$x)
p03 <- ggplot(df_list, aes(x,y, color = comp)) +
  ggtitle("Summarized compartment-level pc length") + 
  geom_line(lwd = 1) +
  geom_vline(xintercept = median(log10(dt[dt$comp == 'A',]$length)), lty = 3, lwd = 0.5, color = '#009999')+
  geom_vline(xintercept = median(log10(dt[dt$comp == 'B',]$length)), lty = 3, lwd = 0.5, color = '#D6604D')+
  #geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
  scale_x_continuous("Distance (bp)", breaks=quantiles, labels=round(10^quantiles, 0)) +
  scale_y_continuous("Density") +
  scale_color_manual(values=setNames(c("#009999", "#D6604D"), c("A", "B"))) +
  # geom_text(data=as.data.table(quantiles, keep.rownames="label"), 
  #           aes(x=round(quantiles, 1), y=-0.02, label=label), size=3) +
  theme_classic() + 
  theme(axis.text=element_text(size=14), title = element_text(size=16),
        axis.title=element_text(size=18), 
        #strip.text=element_text(size=20), 
        strip.background=element_rect(fill = NA, colour = NA))

df <- data.frame(CTRL = k_IAA_false_syn[comp == 'B'], IAA = k_IAA_true_syn[comp == 'B'])

p0 <- ggpaired(df, cond1 = "CTRL", cond2 = "IAA",
               color = "condition", line.color = "gray", line.size = 0.2) +
  scale_y_log10()+
  stat_compare_means(paired = TRUE,label = 'p.signif')+
  scale_color_manual(values=setNames(c("#009999", "#D6604D"),c("CTRL", "IAA")))

#  ####


