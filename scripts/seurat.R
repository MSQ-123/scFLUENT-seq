### load libraries
library(Seurat)
library(dplyr)
library(reticulate)
library(sctransform)
library(cowplot)
library(ggplot2)
library(viridis)
library(tidyr)
library(magrittr)
library(reshape2)
library(readxl)
library(progeny)
library(readr)
library(stringr)
library(tidydr)
library(tidyverse)
library(data.table)


### preprocess ####

  nFeature_lower <- 500
  nFeature_upper <- NA
  nCount_lower <- 500
  nCount_upper <- 100000
  
  # for 2312_2 data, IAA threshold: nFeature_lower 1000 for mature, CTRL unchanged. To define True/False IAA, use 500 threshold
  
  theme_set(theme_cowplot())
  options(future.globals.maxSize = 4000 * 1024^6)
  
  ### sample list
  samples <- c('./Paper//counts/new_custom_counts/CTRL_10min.custom_counts.tsv.gz',
               './Paper//counts/new_custom_counts/CTRL_10min_mature.custom_counts.tsv.gz',
               './Paper//counts/new_custom_counts/IAA_10min.custom_counts.tsv.gz',
               './Paper//counts/new_custom_counts/IAA_10min_mature.custom_counts.tsv.gz')

  
  samplenames <- c('CTRL_10min_beads','CTRL_10min_mature',
                   'IAA_10min_beads', 'IAA_10min_mature')
  objpath <- './Paper/rds/counts/'

  g2s <- fread('../IntergenicTranscription-master/IAA_gffcompare/GencodeReference/g2s_vm22_gencode.txt',header = F,data.table = F) #disable data.table mode
  colnames(g2s) <- c("geneid","symbol")
  # ids <- data.frame(geneid=IAA$gene, symbol = 'NULL')
  # table(ids$geneid %in% g2s$geneid) #
  # ids <- ids[ids$geneid %in% g2s$geneid,] #
  # ids$symbol <- g2s[match(ids$geneid,g2s$geneid),2] #
  
  
  creatObj_save.helper<- function(filePath, objPath, g2s, samplename, saveFile = F){
    require(data.table)
    require(dplyr)
    require(Seurat)
    
    if(nargs() < 3)stop('must supply g2s table!')
    
    if(file.exists(filePath)){
      tmp<- fread(filePath, sep = '\t', header = T, stringsAsFactors = F)
    }
    else{
      stop('file does not exist!')
    }
    
    tmp <- as.data.frame(tmp)
    #rownames(tmp) <- tmp$gene
    ids <- data.frame(geneid=tmp$gene, symbol = 'NULL')
    ids <- ids[ids$geneid %in% g2s$geneid,] #
    ids$symbol <- g2s[match(ids$geneid,g2s$geneid),2] #
    tmp$gene[tmp$gene %in% ids[!duplicated(ids$symbol),]$geneid] <- ids[!duplicated(ids$symbol),]$symbol
    tmp<- tmp%>%distinct(gene, .keep_all = TRUE)
    rownames(tmp) <- tmp$gene
    tmp <- tmp%>%dplyr::select(-gene)
    
    message('filtering the dataframe by excluding features no more than 3 cells expressed')
    
    tmp <- tmp[rowSums(tmp[,] >0) > 3,]
    #saveRDS(tmp, file = objPath)
    myobj<- paste(samplename,'obj', sep = '_')
    
    
    if(!exists(myobj, where = .GlobalEnv)){
      assign(myobj, CreateSeuratObject(counts = tmp, project = samplename), envir = .GlobalEnv)
    }
    else{
      warning(paste('Seurat obj', myobj, 'is overwritten!', sep = ' '))
    }
    if(saveFile){
      if(!file.exists(file.path(objPath,paste(myobj,'rds', sep = '.')))){
        message('save your seurat obj as rds')
        saveRDS(myobj, file = file.path(objPath,paste(myobj,'rds', sep = '.')))
      }
    }
    
  }
  
  creatObj_save.helper(filePath = samples[1], objPath = NULL, g2s = g2s, samplename = samplenames[1], saveFile = F)
  
  CTRL_10min_beads_obj$barcode <- colnames(CTRL_10min_beads_obj)
  CTRL_10min_mature_obj$barcode <- colnames(CTRL_10min_mature_obj)
  IAA_10min_beads_obj$barcode <- colnames(IAA_10min_beads_obj)
  IAA_10min_mature_obj$barcode <- colnames(IAA_10min_mature_obj)

  for (i in 1:length(samples)) {
  message(paste('Processing file', i, sep = ' '))
  creatObj_save.helper(filePath = samples[i], objPath = objpath, g2s = g2s, samplename = samplenames[i], saveFile = T)
  }
  
  CTRL_10min_beads_obj <- subset(CTRL_10min_beads_obj, subset = nFeature_RNA > nFeature_lower & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper)
  CTRL_10min_mature_obj <- subset(CTRL_10min_mature_obj, subset = nFeature_RNA > nFeature_lower & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper)
  IAA_10min_beads_obj <- subset(IAA_10min_beads_obj, subset = nFeature_RNA > nFeature_lower & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper)
  IAA_10min_mature_obj <- subset(IAA_10min_mature_obj, subset = nFeature_RNA > nFeature_lower & nCount_RNA > nCount_lower & nCount_RNA < nCount_upper)
  
  IAA_10min_beads_obj$barcode <- colnames(IAA_10min_beads_obj)
  IAA_10min_mature_obj$barcode <- colnames(IAA_10min_mature_obj)
  CTRL_10min_beads_obj$barcode <- colnames(CTRL_10min_beads_obj)
  CTRL_10min_mature_obj$barcode <- colnames(CTRL_10min_mature_obj)
  
  
  IAA_barcode <- intersect(IAA_10min_beads_obj$barcode, IAA_10min_mature_obj$barcode)
  CTRL_barcode <- intersect(CTRL_10min_beads_obj$barcode, CTRL_10min_mature_obj$barcode)
  
  seu_obj <- merge(CTRL_10min_beads_obj[,CTRL_10min_beads_obj$barcode %in%CTRL_barcode],  
                   y = c(CTRL_10min_mature_obj[,CTRL_10min_mature_obj$barcode %in%CTRL_barcode],
                         IAA_10min_beads_obj[,IAA_10min_beads_obj$barcode %in% IAA_barcode], 
                         IAA_10min_mature_obj[,IAA_10min_mature_obj$barcode%in% IAA_barcode]), 
                   add.cell.ids = c('CTRL_10min_beads', 'CTRL_10min_mature', 'IAA_10min_beads','IAA_10min_mature'),  project = "scEU")
  
  IAA_10min_beads_obj <- IAA_10min_beads_obj[,IAA_10min_beads_obj$barcode %in% IAA_barcode]
  IAA_10min_mature_obj <- IAA_10min_mature_obj[,IAA_10min_mature_obj$barcode %in% IAA_barcode]

  intersect(colnames(Genic_seu[,Genic_seu$orig.ident %in%c('IAA_10min_beads','IAA_10min_mature')]), colnames(Genic_seu_IAA_10min))
  
  
  ### calculate mitochondrial, hemoglobin and ribosomal gene counts
  
  Idents(seu_obj) <- factor(Idents(seu_obj), levels = c('CTRL_10min_beads','CTRL_10min_mature', 'IAA_10min_beads', 'IAA_10min_mature'))
  qcparams <- c("nFeature_RNA", "nCount_RNA")
  p_list <- list()
  for (i in seq_along(qcparams)){
    p_list[[i]] <- VlnPlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident", pt.size = 0, log = T, 
                           cols = c('#999999',"#009999","#B2182B", "#D6604D")) + NoLegend()
  }
  for (i in seq_along(qcparams)){
    print(RidgePlot(object = seu_obj, features = qcparams[i], group.by = "orig.ident"))
  }
  my_levels <- c('CTRL_10min_beads','CTRL_10min_mature', 
                 'IAA_10min_beads', 'IAA_10min_mature')
  unfil_seu <- seu_obj
  
  unfil_seu@active.ident <- factor(x = unfil_seu@active.ident, levels = my_levels)
  unfil_seu$orig.ident <- factor(x = unfil_seu$orig.ident, levels = my_levels)
  
 
  saveRDS(seu_obj, file = 'Paper//rds/counts/matched_seu_filtered100G.rds')
  

  
  #only minor filtered ####
  Distal_seu<- seu_obj[grep(rownames(seu_obj), pattern = 'Distal', value = T),]
  Proximal_seu <- seu_obj[grep(rownames(seu_obj), pattern = 'Proximal', value = T),]
  Genic_seu <- seu_obj[grep(rownames(seu_obj), pattern = 'Proximal|Distal', value = T, invert = T),]
  
  preprocess_pca.helper <- function(seu_obj, nfeatures = 3000, ndims = 15){
    seu_obj <- NormalizeData(seu_obj)
    seu_obj <- FindVariableFeatures(seu_obj, nfeatures = nfeatures)
    seu_obj <- ScaleData(seu_obj,verbose = T)
    
    seu_obj <- RunPCA(seu_obj)
    #ElbowPlot(seu_obj, ndims = 50)
    
    seu_obj <- FindNeighbors(seu_obj, reduction = "pca", dims = 1:ndims)
    seu_obj <- FindClusters(seu_obj, resolution = 0.2)
    seu_obj <- RunUMAP(seu_obj, reduction = "pca", dims = 1:ndims)
    return(seu_obj)
  }

  
  Proximal_seu<- preprocess_pca.helper(Proximal_seu, nfeatures = 3000, ndims = 20)
  Genic_seu<- preprocess_pca.helper(Genic_seu, nfeatures = 3000,ndims = 20)


  ### ####
 library(ggsc)

  
  
p0 <-   DimPlot(Genic_seu, reduction = "umap", group.by = "orig.ident", cols = c('#999999',"#009999","#B2182B", "#E69F00")) + 
    labs(title = "UMAP plot of Genic") + 
    theme_dr(xlength = 0.2, ylength = 0.2, arrow = grid::arrow(length = unit(0.1, 'inches'), type = 'closed')) +
    #NoLegend() + 
    theme(aspect.ratio = 1, panel.grid = element_blank(), axis.title.x = element_text(size = 18, hjust = 0.03),
          axis.title.y = element_text(size = 18, hjust = 0.03), title = element_text(size = 18), legend.text = element_text(size = 18, lineheight = 2),
          plot.title = element_text(color="black", face="bold")) 

