library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)

A_B_merged <- rbind(A_Overlapped_final[,], B_Overlapped_final[,])
saveRDS(A_B_merged, file = './Paper/rds/Distals_DNasePeaks/A_B_merged_DNase_distal_gr.rds')
A_B_merged_width <- as.data.frame(A_B_merged%>%group_by(comp.ID)%>%
  summarise(width = sum(end - start)))

#这里没有对DNase peak排除gene nearby

# single-cell kinetics ####

IAA_10min<- fread('./Paper/rds/counts/DNasePeaks/merged/IAA_10min.merged_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
IAA_10min <- as.data.frame(IAA_10min)
rownames(IAA_10min) <- IAA_10min$gene
IAA_10min <- IAA_10min%>%select(-gene)

IAA_10min_mature<- fread('./Paper/rds/counts/DNasePeaks/merged/IAA_10min_mature.merged_counts.tsv.gz', sep = '\t', header = T, stringsAsFactors = F)
IAA_10min_mature <- as.data.frame(IAA_10min_mature)
rownames(IAA_10min_mature) <- IAA_10min_mature$gene
IAA_10min_mature <- IAA_10min_mature%>%select(-gene)

peaks<- intersect(rownames(IAA_10min), rownames(IAA_10min_mature)) # 3916 peaks

IAA_10min <- as.data.frame(IAA_10min[match(peaks,rownames(IAA_10min)),])
IAA_10min_mature <- as.data.frame(IAA_10min_mature[match(peaks,rownames(IAA_10min_mature)),])
rownames(IAA_10min) <- peaks
rownames(IAA_10min_mature) <- peaks

cells <- intersect(colnames(IAA_10min), colnames(IAA_10min_mature))
cells <- intersect(cells, IAA_barcode)

idx <- match(cells, colnames(IAA_10min))
IAA_10min<- IAA_10min[,idx]
idx <- match(cells, colnames(IAA_10min_mature))
IAA_10min_mature<- IAA_10min_mature[,idx]

idx <- peaks[which(rowSums(IAA_10min) > 10)] # filtering !
IAA_10min<- IAA_10min[idx,]
IAA_10min_mature<- IAA_10min_mature[idx,]

#IAA_10min.fwd.sum <- rowSums(IAA_10min.fwd)
True_IAA <- intersect(cells, True_IAA)
False_IAA <- intersect(cells, False_IAA)
IAA_10min_true <- IAA_10min[,True_IAA]
IAA_10min_false <- IAA_10min[,False_IAA]
IAA_10min_true_mature <- IAA_10min_mature[,True_IAA]
IAA_10min_false_mature <- IAA_10min_mature[,False_IAA]

#try using colSum to calculate each cell kinetics
rownames(A_B_merged_width) <- A_B_merged_width$comp.ID
peak_width <- A_B_merged_width[rownames(IAA_10min_true),]$width
dat<- data.frame(compID = A_B_merged_width[rownames(IAA_10min_true),]$comp.ID, length = peak_width)
dat$comp <- NA
dat$comp <- ifelse(dat$compID %in% A.compartment.gr$comp.ID, 'A', 'B')


IAA_10min_true<- IAA_10min_true/peak_width/0.2*0.98
IAA_10min_true_mature<- IAA_10min_true_mature/peak_width/0.2
IAA_10min_false<- IAA_10min_false/peak_width/0.2*0.98
IAA_10min_false_mature<- IAA_10min_false_mature/peak_width/0.2

IAA_true_theta <- colSums(IAA_10min_true)/(colSums(IAA_10min_true_mature)+colSums(IAA_10min_true))
IAA_false_theta <- colSums(IAA_10min_false)/(colSums(IAA_10min_false_mature)+colSums(IAA_10min_false))


k_IAA_true <- - log(1 - IAA_true_theta)/10
k_IAA_false <- - log(1 - IAA_false_theta)/10

quantile(log(2)/k_IAA_true, na.rm = T)
quantile(log(2)/k_IAA_false, na.rm = T)




# try to separate compartment ####
idx <- rownames(IAA_10min_true)%in%A.compartment.gr$comp.ID   # idx <- rownames(IAA_10min_true)%in%A_B_merged[A_B_merged$V8 == 'A',]$comp.ID
IAA_true_theta_A <- colSums(IAA_10min_true[idx,])/(colSums(IAA_10min_true_mature[idx,])+colSums(IAA_10min_true[idx,]))
IAA_false_theta_A <- colSums(IAA_10min_false[idx,])/(colSums(IAA_10min_false_mature[idx,])+colSums(IAA_10min_false[idx,]))

# IAA_true_theta_A <- rowSums(IAA_10min_true[idx,])/(rowSums(IAA_10min_true_mature[idx,])+rowSums(IAA_10min_true[idx,]))
# IAA_false_theta_A <- rowSums(IAA_10min_false[idx,])/(rowSums(IAA_10min_false_mature[idx,])+rowSums(IAA_10min_false[idx,]))

k_IAA_true_A <- - log(1 - IAA_true_theta_A)/10
k_IAA_false_A <- - log(1 - IAA_false_theta_A)/10

# k_IAA_true_A <- k_IAA_true_A[is.finite(k_IAA_true_A)& (k_IAA_true_A !=0) ]
# k_IAA_false_A <- k_IAA_false_A[is.finite(k_IAA_false_A)& (k_IAA_false_A !=0)]

#idx <- intersect(names(k_IAA_true_A),names(k_IAA_false_A))

quantile(log(2)/k_IAA_true_A, na.rm = T) # intergenic in A has longer half life than B
quantile(log(2)/k_IAA_false_A, na.rm = T)

idx <- rownames(IAA_10min_true)%in%B.compartment.gr$comp.ID # # idx <- rownames(IAA_10min_true)%in%A_B_merged[A_B_merged$V8 == 'A',]$comp.ID
IAA_true_theta_B <- colSums(IAA_10min_true[idx,])/(colSums(IAA_10min_true_mature[idx,])+colSums(IAA_10min_true[idx,]))
IAA_false_theta_B <- colSums(IAA_10min_false[idx,])/(colSums(IAA_10min_false_mature[idx,])+colSums(IAA_10min_false[idx,]))

#IAA_true_theta_B <- rowSums(IAA_10min_true[idx,])/(rowSums(IAA_10min_true_mature[idx,])+rowSums(IAA_10min_true[idx,]))
#IAA_false_theta_B <- rowSums(IAA_10min_false[idx,])/(rowSums(IAA_10min_false_mature[idx,])+rowSums(IAA_10min_false[idx,]))

k_IAA_true_B <- - log(1 - IAA_true_theta_B)/10
k_IAA_false_B <- - log(1 - IAA_false_theta_B)/10

# k_IAA_true_B <- k_IAA_true_B[is.finite(k_IAA_true_B)& (k_IAA_true_B !=0) ]
# k_IAA_false_B <- k_IAA_false_B[is.finite(k_IAA_false_B)& (k_IAA_false_B !=0)]

# idx <- intersect(names(k_IAA_true_B),names(k_IAA_false_B))

quantile(log(2)/k_IAA_true_B, na.rm = T) # intergenic in A has longer half life than B
quantile(log(2)/k_IAA_false_B, na.rm = T)


df_CTRL <- data.frame(k_IAA_false_A = k_IAA_false_A, k_IAA_false_B=k_IAA_false_B, group = 'IAA_false')
df_IAA <- data.frame(k_IAA_true_A = k_IAA_true_A, k_IAA_true_B=k_IAA_true_B, group = 'IAA_true')

save(k_IAA_false_A, k_IAA_false_B, k_IAA_true_A, k_IAA_true_B, file = './Paper/Figures/Figure_patch/Cell/plot_rds/sc_comp_kinetics_deg_len.rda')



# also estimate synthesis rate (I ignored length here, since the unit is not of defined length) #####
idx <- rownames(IAA_10min_true)%in%A.compartment.gr$comp.ID # 
k_IAA_true_A_syn <- colSums(IAA_10min_true[idx,] + IAA_10min_true_mature[idx,]) * k_IAA_true_A
k_IAA_false_A_syn <- colSums(IAA_10min_false[idx,] + IAA_10min_false_mature[idx,]) * k_IAA_false_A

idx <- rownames(IAA_10min_true)%in%B.compartment.gr$comp.ID
k_IAA_true_B_syn <- colSums(IAA_10min_true[idx,] + IAA_10min_true_mature[idx,]) * k_IAA_true_B
k_IAA_false_B_syn <- colSums(IAA_10min_false[idx,] + IAA_10min_false_mature[idx,]) * k_IAA_false_B


# plot correlation between syn and deg
save(k_IAA_false_A_syn, k_IAA_false_B_syn, k_IAA_true_A_syn, k_IAA_true_B_syn, file = './Paper/Figures/Figure_patch/Cell/plot_rds/sc_comp_kinetics_syn_len.rda')

# also plot for feature level #####

idx <- names(k_IAA_true_A)
k_IAA_true_A_syn <- rowSums(IAA_10min_true[idx,] + IAA_10min_true_mature[idx,]) * k_IAA_true_A
idx <- names(k_IAA_false_A)
k_IAA_false_A_syn <- rowSums(IAA_10min_false[idx,] + IAA_10min_false_mature[idx,]) * k_IAA_false_A

idx <- names(k_IAA_true_B)
k_IAA_true_B_syn <- rowSums(IAA_10min_true[idx,] + IAA_10min_true_mature[idx,]) * k_IAA_true_B
idx <- names(k_IAA_false_B)
k_IAA_false_B_syn <- rowSums(IAA_10min_false[idx,] + IAA_10min_false_mature[idx,]) * k_IAA_false_B


# quick view of intergenic distribution among cells, prepare for burst estimate and CV estimate #####

# density of rowSum for A and B ####

quantile(apply(IAA_10min_false, 2, FUN = function(x)sum(x>0))) # per cell 
quantile(apply(IAA_10min_false, 1, FUN = function(x)sum(x>0))) # per feature

tmp <- data.frame(activity = c(apply(DNase_10min, 1, FUN = function(x)sum(x>0)), apply(distal_10min, 1, FUN = function(x)sum(x>0))),
                 comp = c(rep('DHS', nrow(DNase_10min)), rep('distal', nrow(distal_10min))),
                 comp.ID = c(rownames(DNase_10min), rownames(distal_10min)) )


df <- data.frame(activity = c(rowSums(CTRL_10min_A_n),rowSums(CTRL_10min_B_n)), 
                 comp = c(rep('A', nrow(CTRL_10min_A_n)), rep('B', nrow(CTRL_10min_B_n))),
                 comp.ID = c(rownames(CTRL_10min_A_n), rownames(CTRL_10min_B_n)))

df2 <- data.frame(activity = c(colSums(CTRL_10min_A_n),colSums(CTRL_10min_B_n)), 
                 comp = c(rep('A', ncol(CTRL_10min_A_n)), rep('B', ncol(CTRL_10min_B_n))))

ggplot(df2, aes(x=activity, color = comp))+
  #scale_x_log10()+
  geom_histogram(fill="white", binwidth = 0.05, alpha=0.2, position = 'identity')+
  scale_color_manual(values=c( 'A' = "#999999",'B' = "#E69F00")) +
  theme_classic() + 
  facet_grid(rows = vars(comp))+
  # facet_zoom(x >= 100) +
  theme(axis.text=element_text(size=14), title = element_text(size=16),
        axis.title=element_text(size=18), 
        strip.text=element_text(size=14), 
        strip.background=element_rect(fill = NA, colour = NA))




