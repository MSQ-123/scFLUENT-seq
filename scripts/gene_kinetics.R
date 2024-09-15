library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)

# estimate gene-level syn and deg first ####

IAA_10min<- fread('./Paper/counts/constitutiveExon/IAA_10min.constExon_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
rownames(IAA_10min) <- IAA_10min$gene
IAA_10min <- IAA_10min%>%select(-gene)

IAA_10min_mature<- fread('./Paper/counts/constitutiveExon/IAA_10min_mature.constExon_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
rownames(IAA_10min_mature) <- IAA_10min_mature$gene
IAA_10min_mature <- IAA_10min_mature%>%select(-gene)

genes<- intersect(rownames(IAA_10min), rownames(IAA_10min_mature)) # 19277 genes

IAA_10min <- as.data.frame(IAA_10min[match(genes,rownames(IAA_10min)),])
IAA_10min_mature <- as.data.frame(IAA_10min_mature[match(genes,rownames(IAA_10min_mature)),])
rownames(IAA_10min) <- genes
rownames(IAA_10min_mature) <- genes

idx <- genes[which(rowSums(IAA_10min) > 10)]
IAA_10min<- IAA_10min[idx,]
IAA_10min_mature<- IAA_10min_mature[idx,]

IAA_10min_true <- IAA_10min[,True_IAA]
IAA_10min_false <- IAA_10min[,False_IAA]
IAA_10min_true_mature <- IAA_10min_mature[,True_IAA]
IAA_10min_false_mature <- IAA_10min_mature[,False_IAA]

constitutiveExons_noOverlap_length <- as.data.frame(constitutiveExons_noOverlap_length)
rownames(constitutiveExons_noOverlap_length) <- constitutiveExons_noOverlap_length$gene_id
peak_width <- constitutiveExons_noOverlap_length[rownames(IAA_10min_true),]$len

# for gene level kinetics

IAA_10min_true<- IAA_10min_true/peak_width/0.2*0.98
IAA_10min_true_mature<- IAA_10min_true_mature/peak_width/0.2
IAA_10min_false<- IAA_10min_false/peak_width/0.2*0.98
IAA_10min_false_mature<- IAA_10min_false_mature/peak_width/0.2


IAA_true_theta <- rowSums(IAA_10min_true)/(rowSums(IAA_10min_true_mature)+rowSums(IAA_10min_true))
IAA_false_theta <- rowSums(IAA_10min_false)/(rowSums(IAA_10min_false_mature)+rowSums(IAA_10min_false))

k_IAA_true <- - log(1 - IAA_true_theta)/10
k_IAA_false <- - log(1 - IAA_false_theta)/10


k_IAA_true <- k_IAA_true[is.finite(k_IAA_true)& (k_IAA_true !=0) ]
k_IAA_false <- k_IAA_false[is.finite(k_IAA_false)& (k_IAA_false !=0)]

quantile(log(2)/k_IAA_true, na.rm = T)
# 0%        25%        50%        75%       100% 
# 6.710621  45.225413  62.280384  87.783197 654.169214
quantile(log(2)/k_IAA_false, na.rm = T)
# 0%        25%        50%        75%       100% 
# 3.180615  26.230026  34.698389  45.815963 468.502365

dat<- data.frame(false_hv = c(log(2)/k_IAA_true,log(2)/k_IAA_false), 
                 comp = c(rep('IAA', length(k_IAA_true)), rep('CTRL',length(k_IAA_false))))
sis, add length info
# constitutiveExons_noOverlap_length <- as.data.frame(constitutiveExons_noOverlap_length)
# rownames(constitutiveExons_noOverlap_length) <- constitutiveExons_noOverlap_length$gene_id
# Exon_length <- constitutiveExons_noOverlap_length[names(k_IAA_true),]
# IAA_10min_true <- IAA_10min_true[names(k_IAA_true),]
# IAA_10min_false <- IAA_10min_false[names(k_IAA_false),]
# IAA_10min_true_mature <- IAA_10min_true_mature[names(k_IAA_true),]
# IAA_10min_false_mature <- IAA_10min_false_mature[names(k_IAA_false),]

# k_IAA_true_syn <- (rowSums(IAA_10min_true)/Exon_length$len + rowSums(IAA_10min_true_mature)/Exon_length$len) * k_IAA_true
# k_IAA_false_syn <-(rowSums(IAA_10min_false)/Exon_length$len + rowSums(IAA_10min_false_mature)/Exon_length$len) * k_IAA_false

k_IAA_true_syn <- (rowSums(IAA_10min_true[names(k_IAA_true),]) + rowSums(IAA_10min_true_mature[names(k_IAA_true),])) * k_IAA_true
k_IAA_false_syn <-(rowSums(IAA_10min_false[names(k_IAA_false),]) + rowSums(IAA_10min_false_mature[names(k_IAA_false),])) * k_IAA_false

df_CTRL_syn <- data.frame(k_IAA_false_syn = k_IAA_false_syn, group = 'IAA_false')
df_IAA_syn <- data.frame(k_IAA_true_syn = k_IAA_true_syn, group = 'IAA_true')

dat<- data.frame(false_hv = c(k_IAA_true_syn,k_IAA_false_syn), 
                 comp = c(rep('IAA', length(k_IAA_true_syn)), rep('CTRL',length(k_IAA_false_syn))))

df_list <- list()
for (i in unique(dat$comp)) {
  dens <- density(log10(dat[dat$comp == i,]$false_hv))
  df <- data.frame(x=dens$x, y=dens$y)
  probs <- c(0, 0.25, 0.5, 0.75, 1)
  quantiles <- quantile(log10(dat[dat$comp == i,]$false_hv), prob=probs)
  df$quant <- factor(findInterval(df$x,quantiles))
  df$comp <- i
  df_list[[i]] <- df
}
df_list<- do.call(rbind, df_list)
quantiles <- quantile(df_list$x)
p2 <- ggplot(df_list, aes(x,y, color = comp)) +
  ggtitle("Synthesis rate of annotated features") + 
  geom_line(lwd = 1) +
  geom_vline(xintercept = median(log10(dat[dat$comp == 'CTRL',]$false_hv)), lty = 3, lwd = 0.5, color = '#009999')+
  geom_vline(xintercept = median(log10(dat[dat$comp == 'IAA',]$false_hv)), lty = 3, lwd = 0.5, color = '#D6604D')+
  #geom_ribbon(aes(ymin=0, ymax=y, fill=quant)) +
  scale_x_continuous("Synthesis rate", breaks=quantiles, labels=round(10^quantiles, 3)) +
  scale_y_continuous("Density") +
  scale_color_manual(values=setNames(c("#009999", "#D6604D"), c("CTRL", "IAA"))) +
  # geom_text(data=as.data.table(quantiles, keep.rownames="label"), 
  #           aes(x=round(quantiles, 1), y=-0.02, label=label), size=3) +
  theme_classic() + 
  theme(axis.text=element_text(size=14), title = element_text(size=16),
        axis.title=element_text(size=18), 
        #strip.text=element_text(size=20), 
        strip.background=element_rect(fill = NA, colour = NA))

my_comparisons <- list(c("k_IAA_false_syn","k_IAA_true_syn"))
df_IAA2_syn<- gather(df_IAA_syn, key = 'class', value = 'k_syn', -group)
df_CTRL2_syn<- gather(df_CTRL_syn, key = 'class', value = 'k_syn', -group)

tab2 <- rbind(df_IAA2_syn, df_CTRL2_syn)

tab2$class <- factor(tab2$class, levels = c("k_IAA_false_syn", "k_IAA_true_syn"))
gg2 <- ggplot(tab2, aes(x=class, y=k_syn, color = class)) + 
  geom_violin(trim=FALSE)+ scale_y_log10() + 
  scale_color_manual(values=setNames(c("#009999", "#D6604D"), 
                                     c("k_IAA_false_syn", "k_IAA_true_syn"))) +
  geom_boxplot(width=0.1, position = position_dodge(width = 0.9), fill = 'white')+
  theme_classic()+
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif')+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"),legend.position="none")

save(k_IAA_false, k_IAA_false_syn, k_IAA_true, k_IAA_true_syn, file = './Paper/rds/Txn_burst/gene_kinetics_deg_syn.rda')
save(IAA_10min_true, IAA_10min_true_mature, IAA_10min_false, IAA_10min_false_mature, file = 'Paper/rds/Txn_burst/gene_kinetics_counts.rda')

# also estimate sc level.  ######
genes <- rownames(IAA_10min_true)
genes <- genes[genes%in% new.anno.genes_refadd$gene_id]
pcs <- new.anno.genes_refadd[match(genes, new.anno.genes_refadd$gene_id),]$Subset
idx <- genes[pcs == 'protein_coding']


genes <- rownames(IAA_10min_false)
genes <- genes[genes%in% new.anno.genes_refadd$gene_id]
pcs <- new.anno.genes_refadd[match(genes, new.anno.genes_refadd$gene_id),]$Subset
idx <- genes[pcs == 'protein_coding']

IAA_true_theta <- colSums(IAA_10min_true[idx,],na.rm = T)/(colSums(IAA_10min_true_mature[idx,],na.rm = T)+colSums(IAA_10min_true[idx,],na.rm = T))
IAA_false_theta <- colSums(IAA_10min_false[idx,],na.rm = T)/(colSums(IAA_10min_false_mature[idx,],na.rm = T)+colSums(IAA_10min_false[idx,],na.rm = T))

k_IAA_true <- - log(1 - IAA_true_theta)/10
k_IAA_false <- - log(1 - IAA_false_theta)/10

quantile(log(2)/k_IAA_true, na.rm = T)
quantile(log(2)/k_IAA_false, na.rm = T)

k_IAA_true_syn <- (colSums(IAA_10min_true) + colSums(IAA_10min_true_mature)) * k_IAA_true
k_IAA_false_syn <-(colSums(IAA_10min_false) + colSums(IAA_10min_false_mature)) * k_IAA_false


df_CTRL <- data.frame(k_IAA_false = k_IAA_false,  group = 'IAA_false')
df_IAA <- data.frame(k_IAA_true = k_IAA_true, group = 'IAA_true')

my_comparisons <- list(c("k_IAA_false","k_IAA_true"))
df_IAA2<- gather(df_IAA, key = 'class', value = 'k_deg', -group)
df_CTRL2<- gather(df_CTRL, key = 'class', value = 'k_deg', -group)

tab <- rbind(df_IAA2, df_CTRL2)

tab$class <- factor(tab$class, levels = c("k_IAA_false", "k_IAA_true"))
gg3 <- ggplot(tab, aes(x=class, y=k_deg, color = class)) + 
  geom_violin(trim=FALSE)+ #scale_y_log10() + 
  scale_color_manual(values=setNames(c("#009999", "#D6604D"), 
                                     c("k_IAA_false", "k_IAA_true"))) +
  geom_boxplot(width=0.1, position = position_dodge(width = 0.9), fill = 'white')+
  theme_classic()+
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif')+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"),legend.position="none")

my_comparisons <- list(c("k_IAA_false_syn","k_IAA_true_syn"))

df_CTRL_syn <- data.frame(k_IAA_false_syn = k_IAA_false_syn, group = 'IAA_false')
df_IAA_syn <- data.frame(k_IAA_true_syn = k_IAA_true_syn, group = 'IAA_true')

df_IAA2_syn<- gather(df_IAA_syn, key = 'class', value = 'k_syn', -group)
df_CTRL2_syn<- gather(df_CTRL_syn, key = 'class', value = 'k_syn', -group)

tab2 <- rbind(df_IAA2_syn, df_CTRL2_syn)

tab2$class <- factor(tab2$class, levels = c("k_IAA_false_syn", "k_IAA_true_syn"))
gg4 <- ggplot(tab2, aes(x=class, y=k_syn, color = class)) + 
  geom_violin(trim=FALSE)+ scale_y_log10() + 
  scale_color_manual(values=setNames(c("#009999", "#D6604D"), 
                                     c("k_IAA_false_syn", "k_IAA_true_syn"))) +
  geom_boxplot(width=0.1, position = position_dodge(width = 0.9), fill = 'white')+
  theme_classic()+
  stat_compare_means(comparisons = my_comparisons, label = 'p.signif')+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"),legend.position="none")

save(k_IAA_false_syn, k_IAA_true_syn, k_IAA_false, k_IAA_true, file = './Paper/rds/Txn_burst/sc_kinetics_deg_syn.rda')

# estimate normalized kinetics and corr between comp and genes ####

df_IAA_ <- data.frame(gene_syn = k_IAA_true_syn, 
                      gene_deg = k_IAA_true,
                      A_comp_syn = k_IAA_true_A_syn,
                      A_comp_deg = k_IAA_true_A,
                      B_comp_syn = k_IAA_true_B_syn,
                      B_comp_deg = k_IAA_true_B,
                      group = 'IAA_true')

df_CTRL_ <- data.frame(gene_syn = k_IAA_false_syn, 
                       gene_deg = k_IAA_false,
                       A_comp_syn = k_IAA_false_A_syn,
                       A_comp_deg = k_IAA_false_A,
                       B_comp_syn = k_IAA_false_B_syn,
                       B_comp_deg = k_IAA_false_B,
                       group = 'IAA_false')


ctrl_gene_syn_CV <- sd(df_CTRL_$gene_syn)/mean(df_CTRL_$gene_syn)
ctrl_gene_deg_CV <- sd(df_CTRL_$gene_deg)/mean(df_CTRL_$gene_deg)
ctrl_compA_syn_CV <- sd(df_CTRL_$A_comp_syn)/mean(df_CTRL_$A_comp_syn)
ctrl_compA_deg_CV <- sd(df_CTRL_$A_comp_deg)/mean(df_CTRL_$A_comp_deg)
ctrl_compB_syn_CV <- sd(df_CTRL_$B_comp_syn)/mean(df_CTRL_$B_comp_syn)
ctrl_compB_deg_CV <- sd(df_CTRL_$B_comp_deg)/mean(df_CTRL_$B_comp_deg)

IAA_gene_syn_CV <- sd(df_IAA_$gene_syn)/mean(df_IAA_$gene_syn)
IAA_gene_deg_CV <- sd(df_IAA_$gene_deg)/mean(df_IAA_$gene_deg)
IAA_compA_syn_CV <- sd(df_IAA_$A_comp_syn)/mean(df_IAA_$A_comp_syn)
IAA_compA_deg_CV <- sd(df_IAA_$A_comp_deg)/mean(df_IAA_$A_comp_deg)
IAA_compB_syn_CV <- sd(df_IAA_$B_comp_syn)/mean(df_IAA_$B_comp_syn)
IAA_compB_deg_CV <- sd(df_IAA_$B_comp_deg)/mean(df_IAA_$B_comp_deg)

df <- data.frame(class = c('gene_syn', 'gene_deg', 'compA_syn','compA_deg','compB_syn', 'compB_deg'),
                 CTRL = c(ctrl_gene_syn_CV,ctrl_gene_deg_CV,ctrl_compA_syn_CV,ctrl_compA_deg_CV,
                           ctrl_compB_syn_CV, ctrl_compB_deg_CV),
                 IAA = c(IAA_gene_syn_CV, IAA_gene_deg_CV, IAA_compA_syn_CV, IAA_compA_deg_CV,
                         IAA_compB_syn_CV, IAA_compB_deg_CV))
df<- df%>%gather(key = 'Sample',value = 'CV',-class)

gg7 <- ggplot(data=df, aes(x=class, y=CV, fill=Sample)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=setNames(c("#009999", "#D6604D"), 
                                     c("CTRL", "IAA"))) +
  theme_classic()+
  theme(panel.grid = element_blank(), axis.title.x = element_text(size = 18), axis.text = element_text(size = 14),
        axis.title.y = element_text(size = 18), title = element_text(size = 18), strip.text = element_text(size = 18),
        legend.text = element_text(size = 18, lineheight = 2),
        plot.title = element_text(color="black", face="bold"))

topptx(gg7, filename = 'Paper/Figures/Fig2/syn_deg_CV.pptx')




