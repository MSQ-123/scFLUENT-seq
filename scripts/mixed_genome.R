library(ggplot2)

# contamination estimate ####

beads_mix<- read.table('./2312_2/CTRL_30min_mix_beads_idxstats.txt')
mature_mix<- read.table('./2312_2/CTRL_30min_mix_mature_idxstats.txt')

beads_mix$species <- unlist(lapply(strsplit(beads_mix$V1, split = '_'), function(x)x[[1]]))
mature_mix$species <- unlist(lapply(strsplit(mature_mix$V1, split = '_'), function(x)x[[1]]))
# beads on human 

sum(beads_mix[beads_mix$species == 'GRCh38',]$V3)/sum(beads_mix$V3) 
#[1] 0.01612265

#(1-0.016) 

sum(mature_mix[mature_mix$species == 'GRCh38',]$V3)/sum(mature_mix$V3) 
#[1] 0.6979718
