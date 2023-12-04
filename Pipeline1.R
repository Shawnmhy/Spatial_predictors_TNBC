#-------------------------------------------#
# METABRIC TNBC processing pipeline 1-------#
# @Author: Haoyang Mi ----------------------#
# Date: April 22nd 2023------- --------------#
# REF: 


# The object is to study the difference between IMMC and IDC
# IMMC: Invasive Mixed Mucinous Carcinoma
# IDC: Invasive Ductal Carcinoåma
library(ggplot2); library(ComplexHeatmap); library(circlize); library(ggthemes); library(ggprism)
library(tidyr); library(dplyr); library(plyr); library(readr); library(data.table)
library(survival); library(survminer)
library(RANN); library(Rtsne)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd('..')
source("./Codes/Function.r")

#--------------- Useful info --------------------#
# Clinical Groups:
#                1. ER+/PR-/HER2-:
#                2. ER+/PR+/HER2+: N = 17      
# read dataset
sg_cell_expr <- fread('./SingleCells.csv')


sg_cell_expr_ft <- sg_cell_expr[, c(12:50)] 




normalize_imc_data <- function(data, cofactor = 0.8) {
  # 1. arc-hyperbolic-sine 转换，使用 0.8 作为 cofactor
  arcsinh_transformed_data <- asinh(data / cofactor)
  
  # 2. clipped at 99th centile
  clipped_data <- apply(arcsinh_transformed_data, 2, function(x) {
    p99 <- quantile(x, 0.99)
    x[x > p99] <- p99
    return(x)
  })
  
  # 3. z-score 
  # 3. scale to
  normalized_data <- apply(clipped_data, 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  #featuers
  
  return(normalized_data)
}

# normalization
sg_cell_ft_normalized <- normalize_imc_data(sg_cell_expr_ft) %>%
  as.data.frame() 


sg_cell_expr[, c(12:50)] <- sg_cell_ft_normalized

#sg_cell_expr <- sg_cell_expr1

# which patients has too little number of cells
mid_to_include <- sg_cell_expr %>%
  data.frame() %>%
  group_by(metabric_id, ImageNumber) %>%
  tally() %>%
  data.frame() %>%
  dplyr::filter(n >= 1000)


sg_cell_expr <- sg_cell_expr %>%
  dplyr::filter(metabric_id %in% mid_to_include$metabric_id) %>%
  collect() %>%
  collect()


# formatting the datasert
colnames(sg_cell_expr)[2] <- 'CellID' # rename

sg_cell_expr$Phenotype <- ifelse(sg_cell_expr$is_epithelial, 'Tumor', sg_cell_expr$cellPhenotype) # Discard the old classification criteria


sg_cell_expr[sg_cell_expr$cellPhenotype == 'CD4^{+} T cells', 'Phenotype'] <- 'CD4+ T cells'
sg_cell_expr[sg_cell_expr$cellPhenotype == 'CD4^{+} T cells & APCs', 'Phenotype'] <- 'CD4+ T cells & APCs'
sg_cell_expr[sg_cell_expr$cellPhenotype == 'CD8^{+} T cells', 'Phenotype'] <- 'CD8+ T cells'
sg_cell_expr[sg_cell_expr$cellPhenotype == 'T_{Reg} & T_{Ex}', 'Phenotype'] <- 'Tregs and Tex'
sg_cell_expr[sg_cell_expr$cellPhenotype == 'Fibroblasts FSP1^{+}', 'Phenotype'] <- 'Fibroblasts'
sg_cell_expr[sg_cell_expr$cellPhenotype == 'Myofibroblasts PDPN^{+}', 'Phenotype'] <- 'Myofibroblasts'
sg_cell_expr[sg_cell_expr$cellPhenotype == 'Ki67^{+}', 'Phenotype'] <- 'Ki67+'
sg_cell_expr[sg_cell_expr$cellPhenotype == 'CD57^{+}', 'Phenotype'] <- 'CD57+'
sg_cell_expr[sg_cell_expr$cellPhenotype == 'CD38^{+} lymphocytes', 'Phenotype'] <- 'CD38+ lymphocytes'

# read clinical data and rename columns
clinical <- read_tsv('clinical_data.tsv') 

colnames(clinical)[10] <- 'subtype'
colnames(clinical)[14] <- 'Grade'
colnames(clinical)[24] <- 'Index'
colnames(clinical)[25] <- 'Code'
colnames(clinical)[26] <- 'Survival'
colnames(clinical)[27] <- 'Status'
colnames(clinical)[38] <- 'Stage'


# Get all IDC

selectedDF <- clinical %>%
  dplyr::filter(`Sample ID` %in% sg_cell_expr$metabric_id) %>%
  dplyr::filter(`Cancer Type Detailed` == 'Breast Invasive Ductal Carcinoma') %>%
  dplyr::filter(`ER Status` == 'Negative' & `PR Status` == 'Negative' & `HER2 Status` == 'Negative')

selectedDF$Status <- sapply(strsplit(selectedDF$Status, ':'), "[[", 1) %>%
  as.numeric()

clinical$`PR Status`
unique(selectedDF$`Patient ID`)

#####################################################
#---------------------------------------------------#
#------------ ANALYSIS BEGIN HERE ------------------#
#---------------------------------------------------#
#####################################################
#write.csv(TA_area, 'TA_area.csv')


sg_cell_USETHIS <- sg_cell_expr %>%
  dplyr::filter(metabric_id %in% selectedDF$`Patient ID`) 



# Some patients contain multiple TMAs
# We randomly pick one of these for downstream analysis
# First, identify which patients contain multiple regions
multiPatients <- sg_cell_USETHIS %>%
  group_by(metabric_id, ImageNumber) %>%
  tally() %>% # Until here, the output would be the number of cells for each image from eahc patient
  group_by(metabric_id) %>%
  tally() %>%# Until here, the number of subregion is computed
  data.frame() %>%
  dplyr::filter(n > 1)


TA_area <- data.frame(matrix(nrow = 0, ncol = 0))
for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  #mid <- 'MB-3277'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid)
  
  TA <- ifelse(max(sg_metabricobj$Location_Center_X) > 1000, pi*(0.5)^2, pi*(0.3)^2)
  
  TA_area <- rbind.data.frame(TA_area, cbind.data.frame(mid, TA))
}

# Finally, for these patients, pick up the subregion with the highest number of cells
for(pt in multiPatients$metabric_id){
  
  imagenumber <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == pt) %>%
    group_by(ImageNumber) %>%
    tally() 
  
  
  selected_subregion <- imagenumber[which.max(imagenumber$n), 'ImageNumber'] %>%
    dplyr::select(ImageNumber) %>%
    as.character()
  
  # filter with the randomly selected image number
  sg_cell_USETHIS <- sg_cell_USETHIS %>%
    dplyr::filter(!(metabric_id == pt & ImageNumber != selected_subregion))
}




#########################################################
# Plot Cell Boundaries for illustrative purposes
boundaryFiles <- list.files('./SingleCell_Boundary') %>%
  data.frame() %>%
  `colnames<-` (c('filename'))




for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
 # mid <- 'MB-3277'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid)
  
  
  
  # get the serial number to match the boundary data
  serial <- paste0("MB", strsplit(mid, '-')[[1]][2], '_', unique(sg_metabricobj$ImageNumber))
  
  
  # read single cell boundary data
  sg_boundary <- read.csv(paste0('./SingleCell_Boundary/', 
                                 boundaryFiles[boundaryFiles$filename %like% serial,]))
  

  NewID_vec <- c()
  for(pl in seq_len(max(sg_boundary$CellID))){
    


    polygon_pl <- sg_boundary %>%
      dplyr::filter(CellID == pl)
    
  
    inPoly_vec <- point.in.polygon(sg_metabricobj$Location_Center_X, sg_metabricobj$Location_Center_Y,
                     polygon_pl$X, polygon_pl$Y)
    
    newIndex <- which(inPoly_vec == 1) %>%
      as.numeric()
    
    newID <- sg_metabricobj[newIndex, ] %>%
      slice(1) %>%
      dplyr::select(CellID) %>%
      as.numeric()
    
    #
    NewID <- ifelse(length(newID) == 0, 0, newID)
    
    NewID_vec <- c(NewID_vec, NewID)
    
    
    #ggplot() +
    #  geom_point(data = polygon_pl, aes(X, Y)) #+
    #  geom_point(data = newID, aes(Location_Center_X, Location_Center_Y)) 
      
  }
  
  testxx <- NewID_vec %>%
    data.frame() %>%
    mutate(falseid = seq_len(max(sg_boundary$CellID))) %>%
    replace(is.na(.), 0) %>%
    `colnames<-` (c('TrueID', 'CellID'))
    

  
  test <- merge(sg_boundary, testxx, by = 'CellID') %>%
    dplyr::select(Y, X, TrueID) %>%
    `colnames<-` (c('Y', 'X', 'CellID'))
  
  
  # assign cell type to boundary file
  sg_boundary_w_ctype <- merge(sg_metabricobj, test, by = 'CellID') %>%
    dplyr::select(X, Y, CellID, 'Phenotype')
  
  
  
  
  p <- ggplot() +
    #theme_void() +
    geom_polygon(data = sg_boundary_w_ctype, aes(X, Y, group = CellID, fill = Phenotype), size = 0.1) +
    #geom_point(data = sg_metabricobj, aes(Location_Center_X, Location_Center_Y, color = Phenotype), size = 0.1) +

    scale_fill_manual(values = c('Tumor' = '#e31a1c','Ki67+' = '#ec8a49','Tregs and Tex' = '#f2ad80', 'CD4+ T cells & APCs' = '#ae46cd','CD4+ T cells' = '#d581fa', 
                                 'CD38+ lymphocytes' = '#ddaaf7', 'CD57+' = '#fdc5f6',
                                 'B cells' = '#193ddc', 
                                 'CD8+ T cells' = '#3b75db', 
                                 'Macrophages' = '#6594e9', 
                                 'Macrophages & granulocytes' = '#95cee9', 
                                 'Granulocytes' = '#5eccd0',
                                 'Endothelial' = '#32a02d',
                                 'Fibroblasts' = '#b2df8a',
                                 'Myofibroblasts' = '#b9b8b8')) +
    theme(plot.background = element_rect(fill = "black"),
          panel.background = element_rect(fill = "black"),
          legend.position = 'NA'
    ) +
    ylim(max(sg_boundary$Y) + 5,0) +
    xlim(max(sg_boundary$X) + 5,0) +
    coord_fixed(ratio = 1)
  p
  
  
  
  
  ggsave(p, file=paste0("Figures/CellMap/", mid, ".png"), width = 8, height = 8, units = "in", dpi = 300)
  
  
}


#########################################################
#-------tSNE analysis------#
tsneData <- Rtsne(sg_cell_USETHIS[,c('Histone H3', 'panCK', 'CD38', 'HLA-DR', 'CD15', 'FSP1', 'CD163', 'ICOS', 'OX40',
                                     'CD68', 'CD3', 'CD11c', 'PD-1', 'GITR', 'CD16', 'CD45RA', 'CD45RO', 'FOXP3', 'CD20', 'CD8', 'CD57',
                                     'PDGFRB', 'Caveolin-1', 'CD4', 'CD31-vWF', 'HLA-ABC', 'c-Caspase3', 'DNA1', 'DNA2')]) 

tsneDataDF <- tsneData$Y %>%
  data.frame() %>%
  `colnames<-` (c('x', 'y')) %>%
  mutate(Phenotype = sg_cell_USETHIS$Phenotype)


#-------UMAP analysis------#
names(sg_cell_USETHIS)
umapData <- umap(sg_cell_USETHIS[,c('Histone H3', 'panCK', 'CD38', 'HLA-DR', 'CD15', 'FSP1', 'CD163', 'ICOS', 'OX40',
                                    'CD68', 'CD3', 'CD11c', 'PD-1', 'GITR', 'CD16', 'CD45RA', 'CD45RO', 'FOXP3', 'CD20', 'CD8', 'CD57',
                                    'PDGFRB', 'Caveolin-1', 'CD4', 'CD31-vWF', 'HLA-ABC', 'c-Caspase3', 'DNA1', 'DNA2')])
umapDataDF <- umapData$layout %>%
  data.frame() %>%
  `colnames<-` (c('x', 'y')) %>%
  mutate(Phenotype = sg_cell_USETHIS$Phenotype)


p <- ggplot(umapDataDF, aes(x, y, color = Phenotype)) +
  geom_point() +
  theme_void() +
  theme(legend.position = 'none') +
  scale_color_manual(values = c('Tumor' = '#e31a1c','Ki67+' = '#ec8a49','Tregs and Tex' = '#f2ad80', 'CD4+ T cells & APCs' = '#ae46cd','CD4+ T cells' = '#d581fa', 
                               'CD38+ lymphocytes' = '#ddaaf7', 'CD57+' = '#fdc5f6',
                               'B cells' = '#193ddc', 
                               'CD8+ T cells' = '#3b75db', 
                               'Macrophages' = '#6594e9', 
                               'Macrophages & granulocytes' = '#95cee9', 
                               'Granulocytes' = '#5eccd0',
                               'Endothelial' = '#32a02d',
                               'Fibroblasts' = '#b2df8a',
                               'Myofibroblasts' = '#b9b8b8')) +
  coord_fixed(ratio = 1)
p
ggsave(p, file=paste0("Figures/UMAP.png"), width = 10, height = 10, units = "in", dpi = 300)





##---------- Plot Vessel Mask ------------##
vessel <- read.csv('./Vessel_Boundary/MB3277_685_VesselMask.csv')

pos <- holePolygon(vessel)
gg3 <- ggplot(pos, aes(x=x, y=y)) + 
  theme_void() +
  geom_polygon(aes(group=id), fill = 'grey') +
  geom_path(data=pos, aes(x=x, y=y, group=draw), color = 'black', size=1) +  #draw is the original ring
  scale_fill_discrete("ID") #+
  #ylim(550, 0)
gg3
ggsave(gg3, file=paste0("Figures/MB3277_Vessel.png"), width = 10, height = 10, units = "in", dpi = 300, bg="black")



Epithelium <- read.csv('./Epithleial_Boundary/MB0340_462_EpithleialMask.csv')

pos <- holePolygon(Epithelium)
gg3 <- ggplot(pos, aes(x=x, y=y)) + 
  theme_void() +
  geom_polygon(aes(group=id), fill = 'grey') +
  geom_path(data=pos, aes(x=x, y=y, group=draw), color = 'black', size=1) +  #draw is the original ring
  scale_fill_discrete("ID") +
  ylim(550, 0)
gg3
ggsave(gg3, file=paste0("Figures/MB0340_Epitheleial.png"), width = 10, height = 10, units = "in", dpi = 300, bg="black")





#------------------------------------------------------#
#------------ Cell Proportion Stack -------------------#
#------------------------------------------------------#



# Part 1: calculate immune cell counts 
#-----------------------------------------#

selectedDF %>%
  group_by(`Cancer Type Detailed`) %>%
  tally()



# get the single cell (immune) data for the current metabric ID
sg_cell_USETHIS$Phenotype <- as.factor(sg_cell_USETHIS$Phenotype)
sg_cell_USETHIS$metabric_id <- as.factor(sg_cell_USETHIS$metabric_id)

# read tissue area
tissue_area <- read.csv('TA_area.csv')

metabricobj_im <- sg_cell_USETHIS %>%
  #dplyr::filter(metabric_id == mid) %>%
  dplyr::filter(Phenotype != 'Tumor') %>%
  dplyr::filter(Phenotype != 'Ki67+') %>%
  dplyr::filter(Phenotype != 'Myofibroblasts') %>%
  dplyr::filter(Phenotype != 'Endothelial') %>%
  dplyr::filter(Phenotype != 'Fibroblasts') %>%
  group_by(metabric_id, .drop = FALSE) %>%
  dplyr::summarise(count = n()) %>%
  merge(tissue_area, by = 'metabric_id') %>%
  mutate(density = count / TA)
  


# use what order to reorder stacked plot?
metabricobj_im$density <- as.numeric(metabricobj_im$density)




orders <- metabricobj_im %>% 
  arrange(desc(density))
 

# split into two groups: LumA versus Her2

# Part 2: LumA stacked plot
#-----------------------------------------#

#------- Immune infiltrate -----#


p <-ggplot(data=metabricobj_im, aes(y= reorder(metabric_id, -density), x=density)) +
  geom_bar(stat="identity", fill="grey")+
  theme_bw() +
  xlab('Immune density') +
  #xlim(0, 3500) +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 25),
        axis.title.x = element_text(size = 25),
        axis.ticks.x = element_blank()) +
  scale_x_continuous(breaks = c(0, 9000))

p
ggsave(p, file=paste0("Figures/Immune_sorted_low_to_high.png"), width = 2, height = 15, units = "in", dpi = 300)



# stack plot
BCdata_stacked <- sg_cell_USETHIS %>%
  group_by(Phenotype,metabric_id, .drop = FALSE) %>%
  dplyr::summarise(count = n()) 

# reorder sample according to immune density

orders <- metabricobj_im %>% arrange(desc(density))

BCdata_stacked$metabric_id <- factor(BCdata_stacked$metabric_id, levels = orders$metabric_id)

p <- ggplot(BCdata_stacked, aes(fill=Phenotype, x=count, y=as.factor(metabric_id))) + 
  scale_fill_manual(values = c('Tumor' = '#e31a1c','Ki67+' = '#ec8a49','Tregs and Tex' = '#f2ad80', 'CD4+ T cells & APCs' = '#ae46cd','CD4+ T cells' = '#d581fa', 
                               'CD38+ lymphocytes' = '#ddaaf7', 'CD57+' = '#fdc5f6',
                               'B cells' = '#193ddc', 
                               'CD8+ T cells' = '#3b75db', 
                               'Macrophages' = '#6594e9', 
                               'Macrophages & granulocytes' = '#95cee9', 
                               'Granulocytes' = '#5eccd0',
                               'Endothelial' = '#32a02d',
                               'Fibroblasts' = '#b2df8a',
                               'Myofibroblasts' = '#b9b8b8')) +
  geom_bar(position="fill", stat="identity") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.ticks = element_blank(),
        legend.position = 'NA') +
  xlab('Distribution of cell types (%)')
p

ggsave(p, file=paste0("Figures/Stacked_Barplot.png"), width = 5, height = 15, units = "in", dpi = 300)








longTable <- sg_cell_USETHIS %>%
  #dplyr::filter(Phenotype != 'Tumor') %>%
  group_by(metabric_id, Phenotype) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  select(-metabric_id) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate(across(everything(), ~ . / row_sum)) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols = `CD38+ lymphocytes`:`Tregs and Tex`, names_to = 'celltype') 


group_means <- longTable %>%
  group_by(celltype) %>%
  dplyr::summarise(median_value = median(value))

# reorder according to median
sorted_data <- longTable %>%
  inner_join(group_means, by = "celltype") %>%
  arrange(median_value) %>%
  mutate(group = factor(celltype, levels = unique(celltype)))

mean(sorted_data[sorted_data$celltype == 'Tumor', ]$value)

p <- ggplot(sorted_data, aes(y = group, x = value * 100), fill = 'black') +
  stat_boxplot( aes(value*100, group), 
                geom='errorbar', linetype=1, width=0.3)+  #whiskers
  geom_boxplot(aes(ymin = min(value*100),ymax = max(value*100)),  outlier.size = -1) +
  geom_jitter(shape = 21, size = 2, aes(fill = group), alpha = 0.4)+
  xlab("Proportion of non-tumor cells (%)" ) +
  scale_fill_manual(values = c('Tumor' = '#e31a1c','Ki67+' = '#ec8a49','Tregs and Tex' = '#f2ad80', 'CD4+ T cells & APCs' = '#ae46cd','CD4+ T cells' = '#d581fa', 
                               'CD38+ lymphocytes' = '#ddaaf7', 'CD57+' = '#fdc5f6',
                               'B cells' = '#193ddc', 
                               'CD8+ T cells' = '#3b75db', 
                               'Macrophages' = '#6594e9', 
                               'Macrophages & granulocytes' = '#95cee9', 
                               'Granulocytes' = '#5eccd0',
                               'Endothelial' = '#32a02d',
                               'Fibroblasts' = '#b2df8a',
                               'Myofibroblasts' = '#b9b8b8')) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 80, by = 20)) +
  theme(axis.text = element_text(size = 32),
        axis.title.y = element_blank(),
        axis.title = element_text(size = 32),
        legend.position = 'none') 

p
ggsave(p, file = './Figures/Proportion_immune.jpg', width = 13, height = 16, units = "in", dpi = 300)




#-------------------------------------------------------------#
# Myeloid versus lymphoid clustering -------------------------#
#-------------------------------------------------------------#



wideTable <- sg_cell_USETHIS %>%
  dplyr::filter(Phenotype != 'Tumor') %>%
  dplyr::filter(Phenotype != 'Ki67+') %>%
  dplyr::filter(Phenotype != 'Myofibroblasts') %>%
  dplyr::filter(Phenotype != 'Endothelial') %>%
  dplyr::filter(Phenotype != 'Fibroblasts') %>%
  dplyr::filter(Phenotype != 'CD57+') %>%
  group_by(metabric_id, Phenotype) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  #select(-metabric_id) %>%
  rowwise() %>%
  #mutate(row_sum = sum(c_across(-metabric_id))) %>%
  #mutate(across(-c(metabric_id, row_sum), ~ ./row_sum)) %>%
  replace(is.na(.), 0) %>%
  #select(-row_sum) %>%
  t() %>%
  data.matrix() 
  
colnames(wideTable) <- wideTable[1,]
wideTable <- wideTable[-1,]


#set.seed(2)
#km <- MiniBatchKmeans(matrix(as.numeric(wideTable), ncol=ncol(wideTable)), clusters = 3, batch_size = 5)
#clusters <- predict_MBatchKMeans(as.matrix(wideTable[,2:14]), km$centroids)



#nearest_neighbors_clus <- sapply(nearest_neighbors_clus, as.numeric)


library(ComplexHeatmap)
library(colorRamp2)

wideTable_clus <- matrix(as.numeric(wideTable), ncol = ncol(wideTable))
rownames(wideTable_clus) <- rownames(wideTable)
color_mapping <- colorRamp2(c(-1, 0, 1), c("#0047ab", "white", "#fb3640"))

png(file="Figures/Neighborhood_cluster_validation.png", width = 8.5, height = 9, units = "in", res = 300)
p <- Heatmap(wideTable_clus, #col = color_mapping,
             row_order = c(1,2,3,4,7,9,5,6,8),
             cluster_rows = FALSE,
             #cluster_columns = FALSE,
             column_dend_height = unit(4, "cm"),
             show_heatmap_legend = FALSE,
             #heatmap_legend_param = list(legend_direction = "horizontal", 
            #                             legend_width = unit(4, "cm"), 
             #                            at = seq(-1, 1, by = 1),
            #                             labels_gp = gpar(fontsize = 15)),
             row_names_gp = gpar(fontsize = 12), 
             column_names_gp = gpar(fontsize = 20)
)
p
dev.off()





#--------------- Correlation analysis -------------------------------#

sg_cell_USETHIS_Cor <- sg_cell_USETHIS %>%
  dplyr::filter(Phenotype != 'Tumor') %>%
  dplyr::filter(Phenotype != 'Ki67+') %>%
  dplyr::filter(Phenotype != 'Myofibroblasts') %>%
  dplyr::filter(Phenotype != 'Endothelial') %>%
  dplyr::filter(Phenotype != 'Fibroblasts')

sg_cell_USETHIS_Cor$Phenotype <- as.factor(sg_cell_USETHIS_Cor$Phenotype)
freq_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(mid in unique(sg_cell_USETHIS_Cor$metabric_id)){
  
  sg_metabricobj <- sg_cell_USETHIS_Cor %>%
    dplyr::filter(metabric_id == mid) %>%
    
    group_by(Phenotype, .drop = FALSE) %>%
    tally() 
  
  freq <- setNames(data.frame(t(sg_metabricobj[ , - 1])),  as.data.frame(sg_metabricobj[ , 1]) %>%
                     as.matrix() %>%
                     as.character()) %>%
    sweep(2, rowSums(.), '/')
  
  freq_all <- rbind.data.frame(freq_all, cbind.data.frame(freq, mid))
}

colnames(freq_all)[16] <- 'metabric_id' 

test <- merge(freq_all, metabricobj_im, by = 'metabric_id')


cor.test(test[, 2], test[, 12])
# plot the correlation plot betwenn Tregs percentage and immune count
p <- ggplot(test, aes(`Macrophages` * 100, density)) +
  theme_bw() +
  theme(legend.position = 'none',
        #axis.title.y = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26)) +
  geom_point(size = 4, shape =21, fill = 'black') +
  geom_smooth(method = "lm", se = FALSE, color = 'black') +
  ylab('Total immune') +
  xlab(expression(CD4^ "+" ~ "T Cells & APCs, %"))


p



Freq_all <- freq_all %>%
  pivot_longer(cols = `B cells`:`Tumor`, names_to = 'celltype') %>%
  dplyr::filter(celltype == 'Macrophages')

orders <- metabricobj_im %>% arrange(desc(density))

Freq_all$metabric_id <- factor(Freq_all$metabric_id, levels = orders$metabric_id)


p <- ggplot(Freq_all, aes(x=value*100, y=as.factor(metabric_id), fill = celltype)) + 
  scale_fill_manual(values = '#6594e9') + #1a3ddc
  geom_col(aes(y=as.factor(metabric_id))) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 24),
        axis.title.x = element_text(size = 24),
        axis.ticks = element_blank(),
        legend.position = 'none') +
  xlab('Proportion of immune cells (%)')
p
ggsave(p, file=paste0("Figures/Macrophage_Barplot.png"), width = 4, height = 15, units = "in", dpi = 300)





