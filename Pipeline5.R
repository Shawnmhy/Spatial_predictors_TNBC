#-------------------------------------------#
# METABRIC TNBC processing pipeline 4-------#
# @Author: Haoyang Mi ----------------------#
# Date: May 3rd 2023------------------------#
# REF: 


library(ggplot2); 
library(ComplexHeatmap); 
library(circlize); 
library(ggthemes);
library(ggprism)
library(tidyr); library(dplyr); library(plyr); library(readr); library(data.table)
library(survival); library(survminer); library(ade4); library(vegan)
library(RANN); library(Rtsne); library(rstatix)
library(rgeos); library(igraph);library(ClusterR);library(ggvoronoi)

library(tripack)
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
  # 1. arc-hyperbolic-sine cofactor
  arcsinh_transformed_data <- asinh(data / cofactor)
  
  # 2. clipped at 99th centile
  clipped_data <- apply(arcsinh_transformed_data, 2, function(x) {
    p99 <- quantile(x, 0.99)
    x[x > p99] <- p99
    return(x)
  })
  
  # 3. Normalization
  normalized_data <- apply(clipped_data, 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  #featuers
  
  return(normalized_data)
}

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



#####################################################
#---------------------------------------------------#
#------------ ANALYSIS BEGIN HERE ------------------#
#---------------------------------------------------#
#####################################################

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


# make Phenotype as factor
sg_cell_USETHIS$Phenotype <- as.factor(sg_cell_USETHIS$Phenotype)

#---------------------------------------------------#
#------------ ANALYSIS BEGIN HERE ------------------#
#---------------------------------------------------#
nearest_neighbors_all <- data.frame(matrix(nrow = 0, ncol = 0))
#sg_data <- sg_metabricobj

neighbor_freq <- function(i, sg_data){
  
  #i <- 1
  #sg_data <- sg_metabricobj
  nn_idx <- nn2(sg_data[, c('Location_Center_X', 'Location_Center_Y')], 
                query = sg_data[i, c('Location_Center_X', 'Location_Center_Y')], 
                k = 20, treetype = 'kd', searchtype = 'priority') %>%
    .$nn.idx %>%
    as.vector()
  
  # use the index to get the protein expression profile
  res <- sg_data[nn_idx,] %>%
    #slice(-1) %>%
    select(Phenotype) %>%
    group_by(Phenotype, .drop = FALSE) %>%
    dplyr::summarise(count = n()) %>%
    pivot_wider(names_from = Phenotype, values_from = count) %>%
    as.data.frame()
  
  return(res)
}

for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  print(mid)
  ##mid <- 'MB-3277'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid)
  
  
  # 
  
  
  nearest_neighbors <- sapply(seq_len(nrow(sg_metabricobj)),  neighbor_freq, sg_data = sg_metabricobj) %>%
    t() %>%
    data.frame()
  
  nearest_neighbors <- sapply(nearest_neighbors, as.numeric)
  #nearest_neighbors <- nearest_neighbors/10 
  
  nearest_neighbors <- as.matrix(nearest_neighbors) %>%
    data.frame()
  
  
  nearest_neighbors_all <- rbind.data.frame(nearest_neighbors_all, cbind.data.frame(nearest_neighbors, mid, sg_metabricobj$CellID))
  
}



# Mini-Batch KMean clustering
set.seed(2)
km <- MiniBatchKmeans(as.matrix(nearest_neighbors_all[,1:15]), clusters = 10, batch_size = 100)
clusters <- predict_MBatchKMeans(as.matrix(nearest_neighbors_all[,1:15]), km$centroids)


nearest_neighbors_clus <- nearest_neighbors_all

# assign clusters
nearest_neighbors_clus$cluster <- clusters


nearest_neighbors_clus <- nearest_neighbors_clus %>%
  group_by(cluster) %>%
  select(1:15) %>%
  summarise_all(mean, na.rm = TRUE)

nearest_neighbors_clus <- sapply(nearest_neighbors_clus, as.numeric)


# scale clusters
nearest_neighbors_clus[,2:16] <- scale(nearest_neighbors_clus[,2:16])
nearest_neighbors_clus[nearest_neighbors_clus < -1] <- -1; nearest_neighbors_clus[nearest_neighbors_clus >  1] <- 1

# rename columns
colnames(nearest_neighbors_clus) <- c("cluster", "B cells", "CD38+ lymphocytes", "CD4+ T cells", "CD4+ T cells & APCs",
                                      "CD57+", "CD8+ T cells", "Endothelial", "Fibroblasts", "Granulocytes", 'Ki67+', 'Macrophages',
                                      "Macrophages & granulocytes", "Myofibroblasts", "Tregs and Tex", "Tumor")

rownames(nearest_neighbors_clus) <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10')

library(ComplexHeatmap)
library(colorRamp2)

color_mapping <- colorRamp2(c(-1, 0, 1), c("#0047ab", "white", "#fb3640"))

png(file="Figures/Neighborhood_cluster.png", width = 8.5, height = 9, units = "in", res = 300)
p <- Heatmap(nearest_neighbors_clus[,2:16], col = color_mapping,
             cluster_rows = FALSE,
             #cluster_columns = FALSE,
             column_dend_height = unit(4, "cm"),
             show_heatmap_legend = FALSE,
             heatmap_legend_param = list(legend_direction = "horizontal", 
                                         legend_width = unit(4, "cm"), 
                                         at = seq(-1, 1, by = 1),
                                         labels_gp = gpar(fontsize = 15)),
             row_names_gp = gpar(fontsize = 12), 
             column_names_gp = gpar(fontsize = 20)
)
p
dev.off()


#-------------------------------------#
#----------- Donut plot --------------#
#-------------------------------------#


# merge with clinical data
selectedDF_dup <- selectedDF
colnames(selectedDF_dup)[2] <- 'metabric_id'
selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')

donutData_low <- merge(sg_cell_clus, selectedDF_dup, by = 'metabric_id') %>%
  select(clusters, Survival) %>%
  dplyr::filter(Survival == 'Survival low') %>%
  group_by(clusters) %>%
  tally()

donutData_high <- merge(sg_cell_clus, selectedDF_dup, by = 'metabric_id') %>%
  select(clusters, Survival) %>%
  dplyr::filter(Survival == 'Survival high') %>%
  group_by(clusters) %>%
  tally()


# Compute percentages
donutData_low$fraction = donutData_low$n / sum(donutData_low$n)

# Compute the cumulative percentages (top of each rectangle)
donutData_low$ymax = cumsum(donutData_low$fraction)

# Compute the bottom of each rectangle
donutData_low$ymin = c(0, head(donutData_low$ymax, n=-1))

# Make the plot
p <- ggplot(donutData_low, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=clusters)) +
  geom_rect() +
  theme_void() +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  scale_fill_manual(values = c('1' = '#0078ba','2' = '#afcae8','3' = '#ff7700', 
                               '4' = '#ffb96b', '5' = '#00a400', '6' = '#84e080', 
                               '7' = '#a858cd','8' = '#b88fcc', '9' = '#d50b24', 
                               '10' = '#fe969a')) +
  theme(legend.position = 'none')
p
ggsave(p, file=paste0("Figures/Donut_Survival_low.png"), width = 10, height = 10, units = "in", dpi = 300)



# Compute percentages
donutData_high$fraction = donutData_high$n / sum(donutData_high$n)

# Compute the cumulative percentages (top of each rectangle)
donutData_high$ymax = cumsum(donutData_high$fraction)

# Compute the bottom of each rectangle
donutData_high$ymin = c(0, head(donutData_high$ymax, n=-1))

# Make the plot
p <- ggplot(donutData_high, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=clusters)) +
  geom_rect() +
  theme_void() +
  coord_polar(theta="y") + 
  xlim(c(2, 4)) +
  scale_fill_manual(values = c('1' = '#0078ba','2' = '#afcae8','3' = '#ff7700', 
                               '4' = '#ffb96b', '5' = '#00a400', '6' = '#84e080', 
                               '7' = '#a858cd','8' = '#b88fcc', '9' = '#d50b24', 
                               '10' = '#fe969a')) +
  theme(legend.position = 'none')
p
ggsave(p, file=paste0("Figures/Donut_Survival_high.png"), width = 10, height = 10, units = "in", dpi = 300)







# Spatial location of clusters

sg_cell_clus <- cbind.data.frame(sg_cell_USETHIS, clusters)
sg_cell_clus$clusters <- as.factor(sg_cell_clus$clusters)

# Voronoi tesselation

for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  print(mid)
  #mid <- 'MB-0316'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_clus %>%
    dplyr::filter(metabric_id == mid)
  
  p <- ggplot(sg_metabricobj, aes(Location_Center_X, Location_Center_Y, fill = clusters)) +
    geom_voronoi() +
    theme_void() +
    #theme(legend.position = 'none') +ff7700   0078ba
    scale_fill_manual(values = c('1' = '#0078ba','2' = '#afcae8','3' = '#ff7700', 
                                 '4' = '#ffb96b', '5' = '#00a400', '6' = '#84e080', 
                                 '7' = '#a858cd','8' = '#b88fcc', '9' = '#d50b24', 
                                 '10' = '#fe969a')) 
  
  p 
  ggsave(p, file=paste0("Figures/Neighborhoods/", mid, ".png"), width = 8, height = 8, units = "in", dpi = 300)
  
}

# compare the percentage of cells locate in each CN between LTS and STS
selectedDF_dup <- selectedDF
colnames(selectedDF_dup)[2] <- 'metabric_id'
selectedDF_dup$`Relapse Free Status` <- ifelse(selectedDF_dup$`Relapse Free Status` == '1:Recurred', 1, 0)

for(cid in seq_len(10)){
  
  #cid <- 10
  # frequencies
  freqCN <- sg_cell_clus %>%
    group_by(metabric_id, clusters, .drop = FALSE) %>%
    tally() %>%
    dplyr::filter(clusters == cid) 
  
  # get the median
  
  medianFreq <- median(freqCN$n)
  
  # merge dfs
  freqCN_clinical <- merge(freqCN, selectedDF_dup, by = 'metabric_id')
  freqCN_clinical$CN_status <- ifelse(freqCN_clinical$n > medianFreq, 'high', 'low')
  
  
  # 
  #fit <- survfit(Surv(Survival, Status) ~ CN_status, data = freqCN_clinical)
  fit <- survfit(Surv(`Relapse Free Status (Months)`, `Relapse Free Status`) ~ CN_status, data = freqCN_clinical)
  p <- survminer::ggsurvplot(fit, data = freqCN_clinical, 
                             palette = c('#325698', '#e83032'),
                             #legend = 'none',
                             legend.title = '',
                             #legend.labs = c('', '),
                             risk.table = TRUE,
                             surv.scale = 'percent',
                             font.tickslab = c(26),
                             font.title = c(26),
                             font.x = c(28),
                             font.y = c(28),
                             font.legend = c(16),
                             fontsize = 10,
                             tables.theme = theme(axis.text = element_text(size = 16),
                                                  axis.title = element_text(size = 16),
                                                  title = element_text(size = 14)),
                             risk.table.y.text = FALSE,
                             conf.int = FALSE,
                             size = 1.5,
                             censor.size = 8,
                              pval = TRUE
  ) +
    xlab('Time, (months)') +
    ylab('Overall survival')
  
  p
  ggsave(file= paste0("Figures/RFS_Survival", "_CN", cid, '.pdf'), plot = p$plot, width = 6, height = 5, dpi = 300)
  
  
  
}


sg_cell_clus$clusters <- as.factor(sg_cell_clus$clusters)

clusterCount <- sg_cell_clus %>%
  group_by(metabric_id, clusters, .drop = FALSE) %>%
  tally() 

TotalCount <- sg_cell_clus %>%
  group_by(metabric_id, .drop = FALSE) %>%
  tally()

# merge with clinical data
selectedDF_dup <- selectedDF
colnames(selectedDF_dup)[2] <- 'metabric_id'

selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')

clusterPercent <- merge(clusterCount, TotalCount, by = c('metabric_id')) %>%
  mutate(Percent = n.x / n.y) %>%
  merge(selectedDF_dup, by = 'metabric_id') %>%
  dplyr::select(Survival, Percent, clusters) %>%
  dplyr::filter(clusters != 7) %>%
  arrange(Survival, clusters)


clusterPercent[clusterPercent$Survival == 'Survival low' & clusterPercent$clusters == '5',]$Percent
# reorder x axis
clusterPercent$clusters <- paste0("CN", clusterPercent$clusters)
clusterPercent$clusters <- factor(clusterPercent$clusters, levels=c("CN1", "CN2", "CN3", "CN4", "CN5", "CN6", "CN8", "CN9", 'CN10'))

library(rstatix)
clusterPercent %>%
  group_by(clusters) %>%
  wilcox_test(Percent~Survival)




# bar plot
p <- ggbarplot(
  clusterPercent, x = "clusters", y = "Percent", 
  add = c("mean_se", "jitter"), add.params=list(shape = 21, size = 1, width = 0.2),
  color = "Survival", palette = c("#325698", "#e83032"),
  position = position_dodge(0.8), size = 1
)

# layer 2 contains the size of jitter point
p$layers[[2]]$aes_params$size <- 2

p <- p + theme_prism() +
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 22)) +
  ylab('Percent of cells in each CN \n (average across all patients)') 

p
ggsave(p, file=paste0("Figures/CN_percentage.png"), width = 14, height = 5, units = "in", dpi = 300)



#----------------------------------------------------#
#-------- Cox Regression Analysis -------------------#
#----------------------------------------------------#


# frequencies
phenoCount <- sg_cell_clus %>%
  group_by(metabric_id, clusters, Phenotype) %>%
  tally() 

clusterCount <- sg_cell_clus %>%
  group_by(metabric_id, clusters) %>%
  tally()


# get the median
selectedDF_dup <- selectedDF
colnames(selectedDF_dup)[2] <- 'metabric_id'

#selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')

for(cid in c(1,2,3,4,5,6,8,9,10)){
  
  print(cid)  
  cid <- 9
  try({
    for(ctype in unique(sg_cell_clus$Phenotype)){
      
      ctype <- 'Myofibroblasts'
      ###### Inside CN
      InsideDF <- merge(phenoCount, clusterCount, by = c('metabric_id', 'clusters')) %>%
        mutate(Percent = n.x / n.y) %>%
        merge(selectedDF_dup, by = 'metabric_id') %>%
        dplyr::select(Survival, Status, Percent, clusters, Phenotype) %>%
        dplyr::filter(clusters == cid) %>%
        dplyr::filter(Phenotype == ctype)
      
      # reorder x axis
      res.cox <- coxph(Surv(Survival, Status) ~ Percent, data = InsideDF) %>%
        summary()
      
      coef(res.cox)
      ####### Outside CN
      spec_cell_outside <- phenoCount %>%
        dplyr::filter(Phenotype == ctype) %>%
        dplyr::filter(clusters != cid) %>%
        group_by(metabric_id) %>%
        dplyr::summarize(total = sum(n))# celltype-specific count
      
      all_cell_outside <- clusterCount %>%
        dplyr::filter(clusters != cid) %>%
        group_by(metabric_id) %>%
        dplyr::summarize(total = sum(n))
        
      # all celltypes
      outside_DF <- merge(spec_cell_outside, all_cell_outside, by = 'metabric_id') %>%
        mutate(Percent = total.x / total.y) %>%
        merge(selectedDF_dup, by = 'metabric_id') %>%
        dplyr::select(Survival, Status, Percent)
      
      res.cox_out <- coxph(Surv(Survival, Status) ~ Percent, data = outside_DF) %>%
        summary()
      
      
      
      ####### All tissue
      spec_cell_at <- sg_cell_clus %>%
        group_by(metabric_id, Phenotype) %>%
        tally() %>%
        dplyr::filter(Phenotype == ctype)
      
      all_cell_at <- sg_cell_clus %>%
        group_by(metabric_id) %>%
        tally()
      
      # all celltypes
      at_DF <- merge(spec_cell_at, all_cell_at, by = 'metabric_id') %>%
        mutate(Percent = n.x / n.y) %>%
        merge(selectedDF_dup, by = 'metabric_id') %>%
        dplyr::select(Survival, Status, Percent)
      
      res.cox_at <- coxph(Surv(Survival, Status) ~ Percent, data = at_DF) %>%
        summary()
      
      
      
      if(res.cox$coefficients[,5] < 0.05){
        
        print(ctype)
        print(res.cox$coefficients[,5])
        print(res.cox$conf.int[, 3:4])
        print(res.cox_out$coefficients[,5])
        print(res.cox_at$coefficients[,5])
        
      }
      
      
      # cox-regression model
      
    }
    
    
  })
  
}


# frequencies
phenoCount <- sg_cell_clus %>%
  group_by(metabric_id, clusters, Phenotype) %>%
  tally() 

clusterCount <- sg_cell_clus %>%
  group_by(metabric_id, clusters) %>%
  tally()

# get the median
selectedDF_dup <- selectedDF
colnames(selectedDF_dup)[2] <- 'metabric_id'


for(ctype in unique(sg_cell_clus$Phenotype)){
  
  
  
  clusterPercent <- merge(phenoCount, clusterCount, by = c('metabric_id', 'clusters')) %>%
    mutate(Percent = n.x / n.y) %>%
    merge(selectedDF_dup, by = 'metabric_id') %>%
    dplyr::select(Survival, Status, Percent, clusters, Phenotype) %>%
    dplyr::filter(clusters == 5) %>%
    dplyr::filter(Phenotype == ctype)
  
  clusterPercent$Survival <- ifelse(clusterPercent$Survival < 66, 'Survival low', 'Survival high')
  
  high <- clusterPercent[clusterPercent$Survival == 'Survival high', ]
  low <- clusterPercent[clusterPercent$Survival == 'Survival low', ]
  
  
  
  
  tests <- wilcox.test(high$Percent, low$Percent)
  print(tests$p.value)
  
  if(tests$p.value < 0.05){
    print(ctype)
    print(tests$p.value)
    
  }
  
}


library(ggvoronoi)
print(mid)
mid <- 'MB-0269'
# get the single cell data for the current metabric object
sg_metabricobj <- sg_cell_clus %>%
  dplyr::filter(metabric_id == mid)


CD8T <- sg_metabricobj %>%
  dplyr::filter(Phenotype == 'CD8+ T cells')

Granulocytes <- sg_metabricobj %>%
  dplyr::filter(Phenotype == 'Granulocytes')

p <- ggplot(sg_metabricobj, aes(Location_Center_X, Location_Center_Y, fill = clusters, color = 'grey')) +
  geom_voronoi(alpha = 0.4) +
  geom_point(data = Granulocytes, aes(Location_Center_X, Location_Center_Y, color = Phenotype), size = 4) +
  theme_void() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c('1' = '#0078ba','2' = '#afcae8','3' = '#ff7700', 
                               '4' = '#ffb96b', '5' = '#00a400', '6' = '#84e080', 
                               '7' = '#a858cd','8' = '#b88fcc', '9' = '#d50b24', 
                               '10' = '#fe969a')) +
  scale_color_manual(values = c('grey' = 'grey', 'CD8+ T cells' = '#3b75db', 'Granulocytes' = '#3b75db')) 

p 
ggsave(p, file="Figures/Granulocytes_distribution_Voronoi_cox.png", width = 8, height = 8, units = "in", dpi = 300)



#---------------------------------------------------------#
#------------------ Spatial Analysis ---------------------#

# Examine the functional marker expression

# OX40, GITR, ICOS, PD-1

funcMarkers <- c('PD-1', 'OX40', 'GITR', 'ICOS')
marker <- 'OX40'

thresh <- sg_cell_clus %>%
  dplyr::filter(Phenotype == 'Granulocytes') %>%
  select(marker) 


# assign a label 'in cluster' versus 'out cluster'
sg_cell_clus_binary <- sg_cell_clus
sg_cell_clus_binary$inClusterCN5 <- ifelse(sg_cell_clus$clusters == '5', 'TRUE', 'FALSE')
sg_cell_clus_binary$inClusterCN10 <- ifelse(sg_cell_clus$clusters == '10', 'TRUE', 'FALSE')


CD8_clusdist <- sg_cell_clus_binary %>%
  #dplyr::filter(clusters == 3 | clusters == 9) %>%
  dplyr::filter(Phenotype == 'Granulocytes') %>%
  select(metabric_id, Phenotype, eval(funcMarkers), clusters, inClusterCN10) %>%
  pivot_longer(cols = funcMarkers, names_to = c("Markers"), values_to = "Expr") %>%
  group_by(inClusterCN10, metabric_id, Markers) %>%
  dplyr::summarize(proportion = mean(Expr > median(thresh[,1]))) %>%
  merge(selectedDF_dup, by = 'metabric_id')




# compare between incluster proportion versus not incluster propertion

test <- CD8_clusdist %>%
  dplyr::filter(inClusterCN10 == TRUE) %>%
  dplyr::filter(Survival == 'Survival low') %>%
  dplyr::filter(Markers == marker)

test1 <- CD8_clusdist %>%
  dplyr::filter(inClusterCN10 == TRUE) %>%
  dplyr::filter(Survival == 'Survival high') %>%
  dplyr::filter(Markers == marker)


test1$proportion
test$proportion


stat.test <- surv_high_inClusterCN10 %>%
  dplyr::filter(Markers == marker) %>%
  group_by(Markers) %>%
  t_test(proportion~inClusterCN10) %>%
  add_significance() %>%
  add_xy_position()

stat.test

p <- ggplot(surv_high_inClusterCN10, aes(x = inClusterCN10, y = proportion, fill = inClusterCN10), color = 'black') +
  stat_boxplot(geom='errorbar', linetype=1, width=0.3, size = 1.5)+
  geom_boxplot(size = 1.5, outlier.colour = NA) +
  #geom_jitter(aes(fill = inClusterCN3), width = 0, size = 4, shape = 21, stroke = 2) +
  theme_prism() +
  scale_color_manual(values = c('FALSE' = '#7f8082', 'TRUE' = '#7ba6da')) +
  scale_fill_manual(values = c('FALSE' = '#7f8082', 'TRUE' = '#7ba6da')) +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size = 35),
        axis.ticks = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_discrete(labels=c("FALSE" = "Inside", "TRUE" = "Outside")) +
  ylim(0, 1.2)
p
ggsave(p, file="Figures/InCluster_ICOS_CN3_SurvLow.png", width = 5, height = 10, units = "in", dpi = 300)






p <- ggplot(surv_high_inClusterCN9, aes(x = inClusterCN9, y = proportion, fill = inClusterCN9), color = 'black') +
  stat_boxplot(geom='errorbar', linetype=1, width=0.3, size = 1.5)+
  geom_boxplot(size = 1.5, outlier.colour = NA) +
  #geom_jitter(aes(fill = inClusterCN3), width = 0, size = 4, shape = 21, stroke = 2) +
  theme_prism() +
  scale_color_manual(values = c('FALSE' = '#7f8082', 'TRUE' = '#7ba6da')) +
  scale_fill_manual(values = c('FALSE' = '#7f8082', 'TRUE' = '#7ba6da')) +
  theme(legend.position = 'none',
        axis.title = element_blank(),
        axis.text = element_text(size = 35),
        axis.ticks = element_blank(),
        axis.text.x = element_blank()) +
  scale_x_discrete(labels=c("FALSE" = "Inside", "TRUE" = "Outside")) +
  ylim(0, 1.2)
p
ggsave(p, file="Figures/InCluster_ICOS_CN3_SurvLow.png", width = 5, height = 10, units = "in", dpi = 300)





# Deprecated
#-----------------------------------------------------------------------#
#-------- CD8+ T, Granulocytes, Endothelial cells interactions ---------#
#-----------------------------------------------------------------------#
ratios_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  print(mid)
  #mid <- 'MB-0316'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_clus %>%
    dplyr::filter(metabric_id == mid)
  
  
  # get CD8 T cells
  Tumor <- sg_cell_clus %>%
    dplyr::filter(metabric_id == mid) %>%
    dplyr::filter(Phenotype == 'Tumor') %>%
    #dplyr::filter(`PD-1` > thresh1) %>%
    #dplyr::filter(ICOS > thresh2) %>%
    select(Location_Center_X, Location_Center_Y)
  
  
  # Endothelial cells
  Treg <- sg_cell_clus %>%
    dplyr::filter(metabric_id == mid) %>%
    dplyr::filter(Phenotype == "Tregs and Tex") %>%
    select(Location_Center_X, Location_Center_Y)
  
  # Endothelial cells
  Myo <- sg_cell_clus %>%
    dplyr::filter(metabric_id == mid) %>%
    dplyr::filter(Phenotype == "Myofibroblasts") %>%
    select(Location_Center_X, Location_Center_Y)
  
  
  if(nrow(Myo) > 5 & nrow(Treg) > 5 & nrow(Tumor) > 5){
    Tumor_Myo <- nn2(Myo, Treg, k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.dists %>%
      as.vector() 
    
    Tumor_Treg <- nn2(Tumor, Treg, k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.dists %>%
      as.vector()
    
    ratio <- Tumor_Treg / (Tumor_Treg + Tumor_Myo)  
    
    ratios_all <- rbind.data.frame(ratios_all, cbind.data.frame(ratio, mid))  
  }
  
  
}

colnames(ratios_all)[2] <- 'metabric_id'

test <- merge(ratios_all, selectedDF_dup, by ='metabric_id')
stattest_all <- data.frame(matrix(nrow = 0, ncol = 0))
#test$Survival <- ifelse(test$Survival < 66, 'Survival low', 'Survival high')
for(to_exclude in unique(test$metabric_id)){
  
  test_excluded <- test %>%
    dplyr::filter(metabric_id != to_exclude)
  
  
  low <- test_excluded[test_excluded$Survival == 'Survival low', ]
  high <- test_excluded[test_excluded$Survival == 'Survival high', ]
  
  
  stattest <- wilcox.test(low$ratio, high$ratio)
  print(stattest$p.value)
  print(length(unique(low$metabric_id)))
  print(length(unique(high$metabric_id)))
  
  
  stattest_all <- rbind.data.frame(stattest_all, cbind.data.frame(to_exclude, stattest$p.value))
  
}


ggplot(stattest_all, aes(to_exclude, log(`stattest$p.value`))) +
  theme_bw() +
  geom_hline(yintercept = log(0.05), linetype = 'dashed', color = 'red', size = 1) +
  geom_point(size = 5) +
  theme(axis.text = element_text(size = 25),
        axis.title = element_text(size = 25),
        axis.text.x = element_text(angle = 90)) +
  xlab('') +
  ylab('log(p value)')





mycomparison <- compare_means(ratio ~ Survival,  data = test)
mycomparison

my_comparisons <- list( c("Survival low", "Survival high"))

library(Rmisc)
tgc <- summarySE(test, measurevar= "ratio", groupvars=c("Survival"))
#tgc$class <- factor(tgc$class, level = c('short', 'long'))

pd <- position_dodge(0.1) # move them .05 to the left and right
#wilcox.test(eff.scores[eff.scores$class == 'short', 'eff.score'], eff.scores[eff.scores$class == 'long', 'eff.score'])
p <- ggplot(tgc, aes(x= Survival, y = ratio, colour=Survival, group=Survival)) + 
  theme_bw() +
  geom_errorbar(aes(ymin=ratio - 10*se, ymax=ratio + 10*se, color = Survival), width=.05, size = 1) +
  geom_point(aes(fill = Survival, color = 'black'), position=pd, size=5, shape=21) + # 21 is filled circle
  
  scale_color_manual(values = c('Survival high' = '#9dc9e6', 'Survival low' = "#787878")) +
  scale_fill_manual(values = c('Survival high' = '#9dc9e6', 'Survival low' = "#787878")) +
  expand_limits(y=0) +                        # Expand y range
  ylim(0, 1) +
  geom_signif(comparisons = list(c("Survival low", "Survival high")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 0.9, annotations = '****'
  ) +
  theme_prism() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26),
        legend.position = 'none')  +
  ylab('STAIN')
p

ggsave(p, file="Figures/RISC.png", width = 5, height = 7, units = "in", dpi = 300)


#--------------------------------------#
#-------- Pan-Immune Networks ---------#
#--------------------------------------#

ratios_all <- data.frame(matrix(nrow = 0, ncol = 0))
CN_interactions <- data.frame(matrix(nrow = 0, ncol = 0))

CN_interactions_from <- data.frame(matrix(nrow = 0, ncol = 0))
CN_interactions_to <- data.frame(matrix(nrow = 0, ncol = 0))
allpairs <- data.frame(matrix(nrow = 0, ncol = 0))
for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  #print(mid)
  #mid <- 'MB-2834'
  # get the single cell data for the current metabric object
  corePts <- sg_cell_clus %>%
    dplyr::filter(metabric_id == mid) %>%
    #dplyr::filter(clusters == 5 | clusters == 9) %>%
    select(Location_Center_X, Location_Center_Y, clusters, Phenotype)
  
  
  
  
  r <- tri.mesh(corePts$Location_Center_X, corePts$Location_Center_Y)
  delaunayTri <- tripack::triangles(r)
  
  Num_tri <- nrow(delaunayTri)
  # coord length of triangle
  a <- delaunayTri[,1] # node 1 index
  b <- delaunayTri[,2] # node 2 index
  c <- delaunayTri[,3] # node 3 index
  tri_sidelength <- cbind(sqrt( ( r$x[a] - r$x[b] ) ^ 2 + ( r$y[a] - r$y[b] ) ^ 2 ), sqrt( ( r$x[a] - r$x[c] ) ^ 2 + ( r$y[a] - r$y[c] ) ^ 2 ), sqrt( ( r$x[b] - r$x[c] ) ^ 2 + ( r$y[b] - r$y[c] ) ^ 2 ))
  
  
  # area of triangles
  
  # clustering
  k <- seq_len( r$tlnew - 1 )
  i <- r$tlist[k]
  j <- r$tlist[r$tlptr[k]]
  keep <- i > 0
  
  
  i <- abs(i[keep])
  j <- abs(j[keep])
  
  
  distances <- sqrt( ( r$x[i] - r$x[j] ) ^ 2 + ( r$y[i] - r$y[j] ) ^ 2 )
  
  i <- i[distances <= 30]
  j <- j[distances <= 30]
  #--------------------------------------------#
  
  Dual_NodeList <- cbind(seq(1, nrow(corePts)),  corePts)
  colnames(Dual_NodeList) <- c('nodes', 'x', 'y', 'clusters', 'Phenotype')
  
  
  
  Dual_EdgeList <- data.frame(cbind(i, j))
  colnames(Dual_EdgeList) <- c('col1', 'col2')
  
  Dual_EdgeList <- Dual_EdgeList %>%
    mutate(from = pmin(col1, col2),
           to = pmax(col1, col2)) %>%
    distinct(from, to)
  
  
  
  # remove isolated nodes
  Dual_NodeList <- Dual_NodeList %>%
    dplyr::filter(nodes %in% Dual_EdgeList$from | nodes %in% Dual_EdgeList$to)
  
  # identify boundaries between valid pair of clusters
  vClusters <- corePts %>%
    group_by(clusters) %>%
    tally() %>%
    dplyr::filter(n > 10) %>%
    select(clusters) %>%
    as.matrix() 
  
  level <- 1 # number larger than 2, not recomend larger than 3
  
  if(nrow(vClusters) > 1){
    
    vClusters <- vClusters%>%
      combn(2) %>%
      t()
    
    for(pair in 1:nrow(vClusters)){
      
      pheno1 <- vClusters[pair,1]
      pheno2 <- vClusters[pair,2]
      
      allpairs <- rbind.data.frame(allpairs, cbind.data.frame(pheno1, pheno2))
      #pheno1 <- 1
      #pheno2 <- 9
      Dual_EdgeList_Pheno <- cbind.data.frame(Dual_EdgeList, Dual_NodeList[match(Dual_EdgeList$from, Dual_NodeList$nodes), 4],
                                              Dual_NodeList[match(Dual_EdgeList$to, Dual_NodeList$nodes), 4]) %>%
        `colnames<-` (c('from', 'to', 'from_pheno', 'to_pheno')) %>%
        dplyr::filter((from_pheno == pheno1 | from_pheno == pheno2) & (to_pheno == pheno1 | to_pheno == pheno2)) %>%
        mutate(status = ifelse(from_pheno == to_pheno, 'intra', 'inter')) %>%
        select(from, to, status) %>%
        pivot_longer(cols = from:to, values_to = 'nodes')
      
      
      interList <- Dual_EdgeList_Pheno %>%
        dplyr::filter(status == 'inter')
      
      intraList <- Dual_EdgeList_Pheno %>%
        dplyr::filter(status == 'intra')
      
      # level 1 (immediate neighbor)
      
      if(nrow(interList) >= 10){
        
        bordercells <- Dual_EdgeList_Pheno[which(ifelse(Dual_EdgeList_Pheno$nodes %in% interList$nodes, TRUE, FALSE)), 'nodes'] %>%
          unique()
        
        
        # extend to higher level neighbors
        
        
        colnames(Dual_EdgeList) <- c('from', 'to')
        
        # assign cell type
        
        border_NodeList <- Dual_NodeList[Dual_NodeList$nodes %in% bordercells$nodes, ] %>%
          mutate(level = '1')
        
        
        lvl <- 1
        while(lvl < level) {
          
          
          # find the nodes in the edgelist, but not in border_NodeList
          
          nl_Dual_EdgeList <- Dual_EdgeList %>%
            dplyr::filter(from %in% border_NodeList$nodes & !(to %in% border_NodeList$nodes) |
                            to %in% border_NodeList$nodes & !(from %in% border_NodeList$nodes))
          
          
          # get next level 
          nl_nodes <- unique(c(nl_Dual_EdgeList$from, nl_Dual_EdgeList$to))
          
          # ensure all next level are the correct phenotypes
          nl_Dual_NodeList <- Dual_NodeList %>%
            dplyr::filter(nodes %in% nl_nodes) %>%
            dplyr::filter(!(nodes %in% border_NodeList$nodes)) %>%
            dplyr::filter(clusters == pheno1 | clusters == pheno2) %>%
            mutate(level = as.character(lvl))
          
          border_NodeList <- rbind.data.frame(border_NodeList, nl_Dual_NodeList)
          
          lvl <- lvl + 1 
          
          
        }
        
        
        ## Compute number of cross-CN interactions
        
        # get border_edgelist
        border_EdgeList <- Dual_EdgeList %>%
          dplyr::filter(from %in% border_NodeList$nodes) %>%
          dplyr::filter(to %in% border_NodeList$nodes) 
        
        # replace the node id with the cluster id so that we can compute number of interactions  
        border_EdgeList <- border_EdgeList %>%
          mutate(from_cn = as.numeric(border_NodeList[match(border_EdgeList$from, border_NodeList$nodes), 'clusters'])) %>%
          mutate(to_cn   = as.numeric(border_NodeList[match(border_EdgeList$to, border_NodeList$nodes), 'clusters'])) %>%
          dplyr::filter(from_cn != to_cn)
        
        
        
        from_celltypeFreq <- table(border_NodeList[border_NodeList$nodes %in% border_EdgeList$from, 'Phenotype']) %>%
          data.frame() %>%
          spread(key = Var1, value = Freq, fill = 0)
        
        to_celltypeFreq <- table(border_NodeList[border_NodeList$nodes %in% border_EdgeList$to, 'Phenotype']) %>%
          data.frame() %>%
          spread(key = Var1, value = Freq, fill = 0)
        
        celltypeFreq <- table(border_NodeList[, 'Phenotype']) %>%
          data.frame() %>%
          spread(key = Var1, value = Freq, fill = 0)
        
        CN_interactions_from <- rbind.data.frame(CN_interactions_from, cbind.data.frame(mid, pheno1, pheno2, nrow(border_EdgeList), from_celltypeFreq / nrow(corePts[corePts$clusters == pheno1,])))
        CN_interactions_to <- rbind.data.frame(CN_interactions_to, cbind.data.frame(mid, pheno1, pheno2, nrow(border_EdgeList), to_celltypeFreq / nrow(corePts[corePts$clusters == pheno2,])))
        CN_interactions <- rbind.data.frame(CN_interactions, cbind.data.frame(mid, pheno1, pheno2, nrow(border_EdgeList), celltypeFreq / (nrow(border_NodeList))))
        
        
        
        #---------------- Plot --------------------------#
        # Plot interactions between CNs
        
        p <- ggplot() +
          #theme_bw() +
          theme_void() +
          geom_segment(aes(x = Dual_NodeList[match(Dual_EdgeList$from, Dual_NodeList$nodes), 2],
                           xend = Dual_NodeList[match(Dual_EdgeList$to, Dual_NodeList$nodes), 2],
                           y = Dual_NodeList[match(Dual_EdgeList$from, Dual_NodeList$nodes), 3],
                           yend = Dual_NodeList[match(Dual_EdgeList$to, Dual_NodeList$nodes), 3]), size = 0.3) +
          #geom_point(data = Dual_NodeList, aes(x, y), size = 3, shape = 21, fill = 'grey', stroke = 0.5) +
          #geom_point(data = Dual_NodeList, aes(x, y, fill = clusters), size = 4, shape = 21) +
          geom_point(data = border_NodeList, aes(x, y, fill = clusters), size = 3.5, shape = 21, stroke = 0.8) +
          #geom_point(data = nl_Dual_NodeList, aes(x, y, fill = 'black'), size = 3, shape = 21, stroke = 1.5) +
          
          theme(legend.position = 'none') +
          coord_fixed(ratio = 1) +
          scale_fill_manual(values = c('1' = '#0078ba','2' = '#afcae8','3' = '#ff7700', 
                                       '4' = '#ffb96b', '5' = '#00a400', '6' = '#84e080', 
                                       '7' = '#a858cd','8' = '#b88fcc', '9' = '#d50b24', 
                                       '10' = '#fe969a')) 
        
        
        p
        
        
        
        #ggsave(p, file= paste0("Figures/Cluster_boundary/MB-0158-Clusters.png"), width = 15, height = 15, units = "in", dpi = 300)
        
        #ggsave(p, file= paste0("Figures/Cluster_boundary/", mid, '/', 'CN', pheno1, '_', 'CN', pheno2,  ".png"), width = 15, height = 15, units = "in", dpi = 300)
        
        
      }
    }
    
  }
}


allpairs_unique <- unique(allpairs)
colnames(CN_interactions_to)[1] <- 'metabric_id'
colnames(CN_interactions_from)[1] <- 'metabric_id'
colnames(CN_interactions)[1] <- 'metabric_id'

# long_term survivors
ls_survivors <- selectedDF_dup[selectedDF_dup$Survival == 'Survival high', 'metabric_id']
ss_survivors <- selectedDF_dup[selectedDF_dup$Survival == 'Survival low', 'metabric_id']
pvals <- data.frame(matrix(nrow = 0, ncol =0))


for(rowid in 1:nrow(allpairs_unique)){
  #rowid <- 1
  phe1 <- allpairs_unique[rowid, 1]
  phe2 <- allpairs_unique[rowid, 2]
  
  phe1 <- 1
  phe2 <- 9
  
  from_pheno <- CN_interactions_from[CN_interactions_from$pheno1 == phe1 & CN_interactions_from$pheno2 == phe2,]
  to_pheno <- CN_interactions_to[CN_interactions_to$pheno1 == phe1 & CN_interactions_to$pheno2 == phe2,]
  all_pheno <- CN_interactions[CN_interactions$pheno1 == phe1 & CN_interactions$pheno2 == phe2,]
  
  from_surv <- merge(from_pheno, selectedDF_dup, by = 'metabric_id')
  to_surv <- merge(to_pheno, selectedDF_dup, by = 'metabric_id')
  all_surv <- merge(all_pheno, selectedDF_dup, by = 'metabric_id')
  
  
  

  
  print('-----------')
  for(id in 5:19){
    try({
      id <- 10
      #from_surv$Group <- ifelse(from_surv[,id] > median(from_surv[,id]), 'High', "Low")
      #from_surv
      test1 <- all_surv[all_surv$Survival == 'Survival high', id]
      test2 <- all_surv[all_surv$Survival == 'Survival low', id]
      
      
      #fit <- coxph(Surv(Survival, Status) ~ Group, data = from_surv)
      #pval <- summary(fit)$coefficient[,5]
      stattest <- wilcox.test(test1, test2)
      
      pvals <- rbind.data.frame(pvals, cbind.data.frame(phe1, phe2, stattest$p.value))
      
      
      if(stattest$p.value < 0.05){
        print(phe1)
        print(phe2)
        print(colnames(to_surv)[id])
        mean(test1)
        mean(test2)
        
        test1DF <- cbind.data.frame(test1, 'Survival high') %>%
          `colnames<-` (c('Percentage', 'Group'))
        test2DF <- cbind.data.frame(test2, 'Survival low') %>%
          `colnames<-` (c('Percentage', 'Group'))
        
        testDF <- rbind.data.frame(test1DF, test2DF)
        
        stat.test <- testDF %>%
          #group_by(Group) %>%
          wilcox_test(Percentage ~ Group) %>%
          add_significance() %>%
          add_xy_position(x = "Group") 
        
        stat.test$y.position <- stat.test$y.position * 100 
        #stat.test$y.position <- log(stat.test$y.position) - 2.5 # - closeness
        pd = position_dodge(width = 0.75)
        
        p <- ggplot(testDF, aes(x= Group, y= Percentage * 100, color=as.factor(Group))) + 
          
          stat_boxplot(geom = "errorbar", width = 0.25, position=pd, size = 1) + 
          scale_color_manual(values = c('Survival high' = '#9dc9e6', 'Survival low' = "#787878")) +
          stat_pvalue_manual(stat.test, size = 10) +
          rremove("x.title") + rremove('y.title') +
          
          geom_boxplot(outlier.color = NA, size = 1) +
          geom_jitter(size = 3, position = position_jitter(height = 0, width = 0)) +
          
          #ylim(-9, -4) + # degree centrality
          #ylim(2, 20) + # betweenness centrality
          #ylim(-2, 1) + # transitivity centrality
          #ylim(-5, -2) + # transitivity centrality
          theme_classic() +
          xlab('') +
          ylab('') +
          theme(axis.text = element_text(size = 26),
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title = element_text(size = 20),)+
          ylim(0 , 5) +
          #font("xy.text", size = 25) +
        #  font("xy.text", size = 25) +
          rremove('legend') 
        p
        
        ggsave(p, file="Figures/CN1CN9_CD8_comparison_infiltration.png", width = 4, height = 7, units = "in", dpi = 300, bg="transparent")
        
        

        
        print(stattest$p.value)
        print(paste('Number of test 1:', length(test1)))
        print(paste('Number of test 2:', length(test2)))

      }
    }, silent = TRUE)  
  }
  
}


#-----------------------------------------------------------------------#
# Construct two heatmaps for long and short-term survivors -------------#
#-----------------------------------------------------------------------#

# create a all_commbination grid to make sure each pair will be covered
all_combinations <- expand.grid(pheno1 = 1:10, pheno2 = 1:10)
library(purrr)

# Function to add ID and combine data frames
add_id_and_combine <- function(df, ID) {
  list_of_dfs <- map(ID, ~df %>%
                       mutate(metabric_id = .x))
  
  combined_df <- bind_rows(list_of_dfs)
  return(combined_df)
}

# Apply the function
ID <- unique(CN_interactions$metabric_id)
all_combinations <- add_id_and_combine(all_combinations, ID)




# merge with survival data
#CN_interactions_surv <- merge(CN_interactions, selectedDF_dup, by = 'metabric_id')
# convert to integers
CN_interactions$pheno1 <- as.integer(CN_interactions$pheno1)
CN_interactions$pheno2 <- as.integer(CN_interactions$pheno2)
colnames(CN_interactions)[4] <- 'total'
  
# join the real interaction matrix with the expanded, NA values will be added and replaced by 0


df_complete <- all_combinations %>%
  left_join(CN_interactions, by = c("pheno1", "pheno2", 'metabric_id')) %>%
  merge(selectedDF_dup, by = 'metabric_id') %>%
  replace_na(list(count = 0))

df_complete <- df_complete %>%
  mutate_all(~ifelse(is.na(.), 0, .))


# duplicate so that have the other way
df_complete1 <- df_complete
df_complete1$pheno1 <- df_complete$pheno2
df_complete1$pheno2 <- df_complete$pheno1

df_complete <- rbind.data.frame(df_complete, df_complete1)
df_complete$pheno1 <- as.factor(df_complete$pheno1)
df_complete$pheno2 <- as.factor(df_complete$pheno2)


# get the cluster count
clusterCount <- sg_cell_clus %>%
  group_by(metabric_id, clusters) %>%
  tally() %>%
  pivot_wider(names_from = clusters, values_from = n)
  
clusterCount[is.na(clusterCount)] <- 0

# run Survival high, then Survival low
#surv_heatmap_CN_interactions <- df_complete %>%
#  group_by(Survival, metabric_id, pheno1, pheno2, .drop = FALSE) %>%
#  dplyr::summarise(total_all = sum(total)) %>%
#  merge(clusterCount, by = 'metabric_id') %>%
#  rowwise() %>%
#  mutate(norm = total_all / (get(as.character(pheno1)) + get(as.character(pheno2)))) #%>%


#surv_high_heatmap_CN_interactions <- surv_heatmap_CN_interactions %>%
#  dplyr::filter(Survival == 'Survival high') %>%
#  group_by(pheno1, pheno2, .drop = FALSE) %>%
#  dplyr::summarise(mean_norm = mean(norm, na.rm = TRUE)) %>%
#  spread(key = pheno1, value = mean_norm) %>%
#  data_frame() %>%
#  select(-pheno2)
  

#surv_low_heatmap_CN_interactions <- surv_heatmap_CN_interactions %>%
#  dplyr::filter(Survival == 'Survival low') %>%
#  group_by(pheno1, pheno2, .drop = FALSE) %>%
#  dplyr::summarise(mean_norm = mean(norm, na.rm = TRUE)) %>%
#  spread(key = pheno1, value = mean_norm) %>%
#  data_frame() %>%
#  select(-pheno2)

#dplyr::filter(Survival == 'Survival high') %>%
#  spread(key = pheno1, value = total_all) %>%
#  data.frame() %>%
#  select(-c('Survival', 'pheno2'))

surv_low_heatmap_CN_interactions <- df_complete %>%
  group_by(Survival, pheno1, pheno2, .drop = FALSE) %>%
  dplyr::summarise(total_all = sum(total)) %>%
  dplyr::filter(Survival == 'Survival low') %>%
  spread(key = pheno1, value = total_all) %>%
  data.frame() %>%
  select(-c('Survival', 'pheno2'))

surv_high_heatmap_CN_interactions <- df_complete %>%
  group_by(Survival, pheno1, pheno2, .drop = FALSE) %>%
  dplyr::summarise(total_all = sum(total)) %>%
  dplyr::filter(Survival == 'Survival high') %>%
  spread(key = pheno1, value = total_all) %>%
  data.frame() %>%
  select(-c('Survival', 'pheno2'))




surv_heatmap_CN_interactions <- surv_high_heatmap_CN_interactions - surv_low_heatmap_CN_interactions
colnames(surv_heatmap_CN_interactions) <- seq(10)

Heatmap_DF <- melt(setDT(surv_heatmap_CN_interactions, keep.rownames = TRUE), "rn") %>%
  mutate(segment = 1) %>%
  data.frame() %>%
  `colnames<-` (c('col', 'row', 'n', 'segment'))

max_val <- max(abs(Heatmap_DF$n))
Heatmap_DF$col1_scaled <- Heatmap_DF$n / max_val

Heatmap_DF$col <- as.numeric(Heatmap_DF$col)
Heatmap_DF$row <- as.numeric(Heatmap_DF$row)
p <- ggplot() +
  geom_tile(data = Heatmap_DF, aes(x = 1, y = segment, fill = col1_scaled)) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)) +
  #geom_text(data = pval_df, aes(x = 1, y = segment, label = n), vjust = 1) + 
  facet_grid(row~col, switch = "both") +
  #scale_y_continuous(expand = c(0,0)) +
  #scale_x_continuous(expand = c(0,0)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0, "line"),
        panel.border = element_rect(fill = NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none'
        )
p
ggsave(p, file= paste0("Figures/Communication_Heatmap.png"), width = 15, height = 15, units = "in", dpi = 300)



#-------------------------------------------------------#
#--------- Draw Boxplots 
unique_combs <- unique(df_complete[, c('pheno1', 'pheno2')])

TA_areas <- read.csv('TA_area.csv')
for(id in 1:nrow(unique_combs)){
  phe1 <- unique_combs[id, 1]
  phe2 <- unique_combs[id, 2]
  
  phe1 <- 8
  phe2 <- 9
  
  # how many cells are of these clusters? (excluding at boundaries)
  # this is to do normalization
  total_num_clusters_of_interest <- sg_cell_clus %>%
    #dplyr::filter(clusters == phe1 | clusters == phe2) %>%
    group_by(metabric_id) %>%
    tally()
    
    
  try({
    data_for_bp <- df_complete %>%
      group_by(metabric_id, Survival, pheno1, pheno2, .drop = FALSE) %>%
      dplyr::summarise(total_all = sum(total)) %>%
      data.frame() %>%
      dplyr::filter(pheno1 == phe1) %>%
      dplyr::filter(pheno2 == phe2) %>%
      merge(total_num_clusters_of_interest, by = 'metabric_id') %>%
      mutate(norm = total_all / n) 
    
    
    
    test1 <- data_for_bp[data_for_bp$Survival == 'Survival low', ]
    test2 <- data_for_bp[data_for_bp$Survival == 'Survival high', ]
    
    stats <- wilcox.test(test1$norm, test2$norm)
    #print(stats$p.value)
    if(stats$p.value < 0.05){
      print(stats$p.value)    
      print(phe1)
      print(phe2)
    }
    
  })
  
}

data_for_bp <- df_complete %>%
  group_by(metabric_id, Survival, pheno1, pheno2, .drop = FALSE) %>%
  dplyr::summarise(total_all = sum(total)) %>%
  data.frame() %>%
  dplyr::filter(pheno1 == '1') %>%
  dplyr::filter(pheno2 == '8') #%>%
  #dplyr::filter(metabric_id != 0)


mean(test1$total_all)
mean(test2$total_all)
test1 <- data_for_bp[data_for_bp$Survival == 'Survival low', ]
test2 <- data_for_bp[data_for_bp$Survival == 'Survival high', ]

wilcox.test(test1$total_all, test2$total_all)
# add p values
mycomparison <- compare_means(total_all ~ Survival,  data = data_for_bp)
mycomparison

my_comparisons <- list( c("Survival low", "Survival high"))

p <- ggplot(data_for_bp, aes(x = Survival, y = total_all), fill = 'black') +
  stat_boxplot( aes(Survival, total_all, color = Survival), 
                geom='errorbar', linetype=1, width=0.3, size = 2)+  #whiskers
  geom_boxplot(aes(ymin = min(total_all),ymax = max(total_all), color = Survival),  outlier.size = -1, size = 2) +
  geom_jitter(shape = 21, size = 4, aes(fill = Survival), alpha = 1, width = 0)+
  #xlab("Proportion of non-tumor cells (%)" ) +
  scale_fill_manual(values = c('Survival high' = '#8c96c6','Survival low' = '#fdc5c0')) +
  scale_color_manual(values = c('Survival high' = '#8c96c6','Survival low' = '#fdc5c0')) +
  theme_classic() +
  #scale_x_continuous(breaks = seq(0, 80, by = 20)) +
  theme(axis.text = element_text(size = 32),
        axis.text.x = element_blank(),
        axis.title = element_blank(),
        #axis.title = element_text(size = 32),
        legend.position = 'none') 

p

 ggsave(p, file= paste0("Figures/CN1_CN9_Interaction_compare.png"), width = 5, height = 9, units = "in", dpi = 300)

#-----------------------------------------#
# Try validation using Keren et al's data #

validation_data <- readRDS('../../clusterEdge/tnbc.rds')

cellDat <- validation_data$cellDat %>%
  group_by(Key) %>%
  mutate(CellID = row_number()) %>%
  ungroup()




ratios_all <- data.frame(matrix(nrow = 0, ncol = 0))

for(iid in unique(cellDat$Key)){
  
  #print(iid)
  #iid <- 849
  #mid <- 'MB-0316'
  #unique(cellDat$celltype)
  # get CD8 T cells
  CD8T <- cellDat %>%
    dplyr::filter(Key == iid) %>%
    dplyr::filter(celltype == 'CD8 T') %>%
    select(x.um, y.um)
  
  
  # Endothelial cells
  Treg <- cellDat %>%
    dplyr::filter(Key == iid) %>%
    dplyr::filter(celltype == "Treg") %>%
    select(x.um, y.um)
  
  # Endothelial cells
  Tumor <- cellDat %>%
    dplyr::filter(Key == iid) %>%
    dplyr::filter(celltype == "Tumor") %>%
    select(x.um, y.um)
  
  
  if(nrow(CD8T) > 5 & nrow(Treg) > 5 & nrow(Tumor) > 5){
    CD8T_FoxP3 <- nn2(Treg, CD8T, k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.dists %>%
      as.vector() 
    
    CD8T_Tumor <- nn2(Tumor, CD8T, k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.dists %>%
      as.vector()
    
    ratio <- CD8T_Tumor / (CD8T_FoxP3 + CD8T_Tumor)  
    
    ratios_all <- rbind.data.frame(ratios_all, cbind.data.frame(ratio, iid))  
  }
  
  
}

colnames(ratios_all)[2] <- 'Key'
clinical_data <- validation_data$clinicalDat
clinical_data$Survival <- ifelse(clinical_data$`Survival_days_capped*` >= median(clinical_data$`Survival_days_capped*`), 'Survival high', 'Survival low')

ratios_all <- merge(ratios_all, clinical_data, by = 'Key')

library(Rmisc)
library(ggpubr)
library(ggsignif)
tgc <- summarySE(ratios_all, measurevar= "ratio", groupvars=c("Survival"))
mycomparison <- compare_means(ratio ~ Survival,  data = ratios_all)
mycomparison

pd <- position_dodge(0.1) # move them .05 to the left and right

p <- ggplot(tgc, aes(x= Survival, y = ratio, colour=Survival, group=Survival)) + 
  theme_bw() +
  geom_errorbar(aes(ymin=ratio - 10*se, ymax=ratio + 10*se, color = Survival), width=.05, size = 1) +
  geom_point(aes(fill = Survival, color = 'black'), position=pd, size=7, shape=21) + # 21 is filled circle
  
  scale_color_manual(values = c('Survival high' = '#9dc9e6', 'Survival low' = "#787878")) +
  scale_fill_manual(values = c('Survival high' = '#9dc9e6', 'Survival low' = "#787878")) +
  #scale_color_manual(values = c('pCR' = '#1e19b3', 'RD' = "black")) +
  #scale_fill_manual(values = c('pCR' = '#1e19b3', 'RD' = "black")) +
  expand_limits(y=0) +                        # Expand y range
  ylim(0, 1) +
  geom_signif(comparisons = list(c("Survival high", "Survival low")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 0.9, annotations = '****'
  ) +
  theme_prism() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26),
        legend.position = 'none')  +
  ylab('RISC')
p

ggsave(p, file="Figures/RISC_Keren.png", width = 5, height = 7, units = "in", dpi = 300)













cellDat$celltype <- as.factor(cellDat$celltype)

nearest_neighbors_all <- data.frame(matrix(nrow = 0, ncol = 0))
#sg_data <- sg_metabricobj

neighbor_freq <- function(i, sg_data){
  
  #i <- 1
  sg_data <- sg_metabricobj
  nn_idx <- RANN::nn2(sg_data[, c('x.um', 'y.um')], 
                query = sg_data[i, c('x.um', 'y.um')], 
                k = 20, treetype = 'kd', searchtype = 'priority') %>%
    .$nn.idx %>%
    as.vector()
  
  # use the index to get the protein expression profile
  res <- sg_data[nn_idx,] %>%
    #slice(-1) %>%
    select(celltype) %>%
    group_by(celltype, .drop = FALSE) %>%
    dplyr::summarise(count = n()) %>%
    pivot_wider(names_from = celltype, values_from = count) %>%
    as.data.frame()
  
  return(res)
}

for(key in unique(cellDat$Key)){
  
  print(key)
  ##mid <- 'MB-3277'
  # get the single cell data for the current metabric object
  sg_metabricobj <- cellDat %>%
    dplyr::filter(Key == key)
  
  
  # 
  
  
  nearest_neighbors <- sapply(seq_len(nrow(sg_metabricobj)),  neighbor_freq, sg_data = sg_metabricobj) %>%
    t() %>%
    data.frame()
  
  nearest_neighbors <- sapply(nearest_neighbors, as.numeric)
  #nearest_neighbors <- nearest_neighbors/10 
  
  nearest_neighbors <- as.matrix(nearest_neighbors) %>%
    data.frame()
  
  
  nearest_neighbors_all <- rbind.data.frame(nearest_neighbors_all, cbind.data.frame(nearest_neighbors, key, sg_metabricobj$CellID))
  
}



# Mini-Batch KMean clustering
set.seed(2)
km <- MiniBatchKmeans(as.matrix(nearest_neighbors_all[,1:16]), clusters = 10, batch_size = 100)
clusters <- predict_MBatchKMeans(as.matrix(nearest_neighbors_all[,1:16]), km$centroids)


nearest_neighbors_clus <- nearest_neighbors_all

# assign clusters
nearest_neighbors_clus$cluster <- clusters


nearest_neighbors_clus <- nearest_neighbors_clus %>%
  group_by(cluster) %>%
  select(1:16) %>%
  summarise_all(mean, na.rm = TRUE)

nearest_neighbors_clus <- sapply(nearest_neighbors_clus, as.numeric)


# scale clusters
nearest_neighbors_clus[,2:17] <- scale(nearest_neighbors_clus[,2:17])
nearest_neighbors_clus[nearest_neighbors_clus < -1] <- -1; nearest_neighbors_clus[nearest_neighbors_clus >  1] <- 1

# rename columns
colnames(nearest_neighbors_clus)
colnames(nearest_neighbors_clus) <- c("cluster", "B cells", "CD3+ T cells", "CD4+ T cells", "CD8+ T cells",
                                      "DC", "DC/Monocytes", "Endothelial", "Keratin+ Tumor", "Macrophages", 'Mesenchymal-like', 'Mono/Neu',
                                      "Neutrophils", "NK cells", "Other immune", "Treg", "Tumor")

rownames(nearest_neighbors_clus) <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10')

library(ComplexHeatmap)
library(colorRamp2)

color_mapping <- colorRamp2(c(-1, 0, 1), c("#0047ab", "white", "#fb3640"))

png(file="Figures/Neighborhood_cluster_validation.png", width = 8.5, height = 9, units = "in", res = 300)
p <- Heatmap(nearest_neighbors_clus[,2:17], col = color_mapping,
             cluster_rows = FALSE,
             #cluster_columns = FALSE,
             column_dend_height = unit(4, "cm"),
             show_heatmap_legend = FALSE,
             heatmap_legend_param = list(legend_direction = "horizontal", 
                                         legend_width = unit(4, "cm"), 
                                         at = seq(-1, 1, by = 1),
                                         labels_gp = gpar(fontsize = 15)),
             row_names_gp = gpar(fontsize = 12), 
             column_names_gp = gpar(fontsize = 20)
)
p
dev.off()




#------------------------------------------------------------------#
#---------------------- Construct networks ------------------------#

#cellDat <- cbind.data.frame(cellDat, clusters)
cellDat$clusters <- as.character(cellDat$clusters)

ratios_all <- data.frame(matrix(nrow = 0, ncol = 0))
CN_interactions <- data.frame(matrix(nrow = 0, ncol = 0))

CN_interactions_from <- data.frame(matrix(nrow = 0, ncol = 0))
CN_interactions_to <- data.frame(matrix(nrow = 0, ncol = 0))
allpairs <- data.frame(matrix(nrow = 0, ncol = 0))
for(key in unique(cellDat$Key)){
  
  print(key)
  #mid <- 'MB-0211'
  # get the single cell data for the current metabric object
  corePts <- cellDat %>%
    dplyr::filter(Key == key) %>%
    #dplyr::filter(clusters == 5 | clusters == 9) %>%
    select(x.um, y.um, clusters, celltype)
  
  
  
  
  r <- tri.mesh(corePts$x.um, corePts$y.um)
  delaunayTri <- tripack::triangles(r)
  
  Num_tri <- nrow(delaunayTri)
  # coord length of triangle
  a <- delaunayTri[,1] # node 1 index
  b <- delaunayTri[,2] # node 2 index
  c <- delaunayTri[,3] # node 3 index
  tri_sidelength <- cbind(sqrt( ( r$x[a] - r$x[b] ) ^ 2 + ( r$y[a] - r$y[b] ) ^ 2 ), sqrt( ( r$x[a] - r$x[c] ) ^ 2 + ( r$y[a] - r$y[c] ) ^ 2 ), sqrt( ( r$x[b] - r$x[c] ) ^ 2 + ( r$y[b] - r$y[c] ) ^ 2 ))
  
  
  # area of triangles
  
  # clustering
  k <- seq_len( r$tlnew - 1 )
  i <- r$tlist[k]
  j <- r$tlist[r$tlptr[k]]
  keep <- i > 0
  
  
  i <- abs(i[keep])
  j <- abs(j[keep])
  
  
  distances <- sqrt( ( r$x[i] - r$x[j] ) ^ 2 + ( r$y[i] - r$y[j] ) ^ 2 )
  
  i <- i[distances <= 30]
  j <- j[distances <= 30]
  #--------------------------------------------#
  
  Dual_NodeList <- cbind(seq(1, nrow(corePts)),  corePts)
  colnames(Dual_NodeList) <- c('nodes', 'x', 'y', 'clusters', 'Phenotype')
  
  
  
  Dual_EdgeList <- data.frame(cbind(i, j))
  colnames(Dual_EdgeList) <- c('col1', 'col2')
  
  Dual_EdgeList <- Dual_EdgeList %>%
    mutate(from = pmin(col1, col2),
           to = pmax(col1, col2)) %>%
    distinct(from, to)
  
  
  
  # remove isolated nodes
  Dual_NodeList <- Dual_NodeList %>%
    dplyr::filter(nodes %in% Dual_EdgeList$from | nodes %in% Dual_EdgeList$to)
  
  # identify boundaries between valid pair of clusters
  vClusters <- corePts %>%
    group_by(clusters) %>%
    tally() %>%
    dplyr::filter(n > 10) %>%
    select(clusters) %>%
    as.matrix() 
  
  level <- 1 # number larger than 2, not recomend larger than 3
  
  if(nrow(vClusters) > 1){
    
    vClusters <- vClusters%>%
      combn(2) %>%
      t()
    
    for(pair in 1:nrow(vClusters)){
      
      pheno1 <- vClusters[pair,1]
      pheno2 <- vClusters[pair,2]
      
      allpairs <- rbind.data.frame(allpairs, cbind.data.frame(pheno1, pheno2))
      #pheno1 <- 1
      #pheno2 <- 9
      Dual_EdgeList_Pheno <- cbind.data.frame(Dual_EdgeList, Dual_NodeList[match(Dual_EdgeList$from, Dual_NodeList$nodes), 4],
                                              Dual_NodeList[match(Dual_EdgeList$to, Dual_NodeList$nodes), 4]) %>%
        `colnames<-` (c('from', 'to', 'from_pheno', 'to_pheno')) %>%
        dplyr::filter((from_pheno == pheno1 | from_pheno == pheno2) & (to_pheno == pheno1 | to_pheno == pheno2)) %>%
        mutate(status = ifelse(from_pheno == to_pheno, 'intra', 'inter')) %>%
        select(from, to, status) %>%
        pivot_longer(cols = from:to, values_to = 'nodes')
      
      
      interList <- Dual_EdgeList_Pheno %>%
        dplyr::filter(status == 'inter')
      
      intraList <- Dual_EdgeList_Pheno %>%
        dplyr::filter(status == 'intra')
      
      # level 1 (immediate neighbor)
      
      if(nrow(interList) >= 10){
        
        bordercells <- Dual_EdgeList_Pheno[which(ifelse(Dual_EdgeList_Pheno$nodes %in% interList$nodes, TRUE, FALSE)), 'nodes'] %>%
          unique()
        
        
        # extend to higher level neighbors
        
        
        colnames(Dual_EdgeList) <- c('from', 'to')
        
        # assign cell type
        
        border_NodeList <- Dual_NodeList[Dual_NodeList$nodes %in% bordercells$nodes, ] %>%
          mutate(level = '1')
        
        
        lvl <- 1
        while(lvl < level) {
          
          
          # find the nodes in the edgelist, but not in border_NodeList
          
          nl_Dual_EdgeList <- Dual_EdgeList %>%
            dplyr::filter(from %in% border_NodeList$nodes & !(to %in% border_NodeList$nodes) |
                            to %in% border_NodeList$nodes & !(from %in% border_NodeList$nodes))
          
          
          # get next level 
          nl_nodes <- unique(c(nl_Dual_EdgeList$from, nl_Dual_EdgeList$to))
          
          # ensure all next level are the correct phenotypes
          nl_Dual_NodeList <- Dual_NodeList %>%
            dplyr::filter(nodes %in% nl_nodes) %>%
            dplyr::filter(!(nodes %in% border_NodeList$nodes)) %>%
            dplyr::filter(clusters == pheno1 | clusters == pheno2) %>%
            mutate(level = as.character(lvl))
          
          border_NodeList <- rbind.data.frame(border_NodeList, nl_Dual_NodeList)
          
          lvl <- lvl + 1 
          
          
        }
        
        
        ## Compute number of cross-CN interactions
        
        # get border_edgelist
        border_EdgeList <- Dual_EdgeList %>%
          dplyr::filter(from %in% border_NodeList$nodes) %>%
          dplyr::filter(to %in% border_NodeList$nodes) 
        
        # replace the node id with the cluster id so that we can compute number of interactions  
        border_EdgeList <- border_EdgeList %>%
          mutate(from_cn = as.numeric(border_NodeList[match(border_EdgeList$from, border_NodeList$nodes), 'clusters'])) %>%
          mutate(to_cn   = as.numeric(border_NodeList[match(border_EdgeList$to, border_NodeList$nodes), 'clusters'])) %>%
          dplyr::filter(from_cn != to_cn)
        
        
        
        from_celltypeFreq <- table(border_NodeList[border_NodeList$nodes %in% border_EdgeList$from, 'Phenotype']) %>%
          data.frame() %>%
          spread(key = Var1, value = Freq, fill = 0)
        
        to_celltypeFreq <- table(border_NodeList[border_NodeList$nodes %in% border_EdgeList$to, 'Phenotype']) %>%
          data.frame() %>%
          spread(key = Var1, value = Freq, fill = 0)
        
        celltypeFreq <- table(border_NodeList[, 'Phenotype']) %>%
          data.frame() %>%
          spread(key = Var1, value = Freq, fill = 0)
        
        CN_interactions_from <- rbind.data.frame(CN_interactions_from, cbind.data.frame(key, pheno1, pheno2, nrow(border_EdgeList), from_celltypeFreq / nrow(corePts[corePts$clusters == pheno1,])))
        CN_interactions_to <- rbind.data.frame(CN_interactions_to, cbind.data.frame(key, pheno1, pheno2, nrow(border_EdgeList), to_celltypeFreq / nrow(corePts[corePts$clusters == pheno2,])))
        CN_interactions <- rbind.data.frame(CN_interactions, cbind.data.frame(key, pheno1, pheno2, nrow(border_EdgeList), celltypeFreq / (nrow(border_NodeList))))
        
        
        
        #---------------- Plot --------------------------#
        # Plot interactions between CNs
        
        p <- ggplot() +
          #theme_bw() +
          theme_void() +
          geom_segment(aes(x = Dual_NodeList[match(Dual_EdgeList$from, Dual_NodeList$nodes), 2],
                           xend = Dual_NodeList[match(Dual_EdgeList$to, Dual_NodeList$nodes), 2],
                           y = Dual_NodeList[match(Dual_EdgeList$from, Dual_NodeList$nodes), 3],
                           yend = Dual_NodeList[match(Dual_EdgeList$to, Dual_NodeList$nodes), 3]), size = 0.3) +
          #geom_point(data = Dual_NodeList, aes(x, y), size = 3, shape = 21, fill = 'grey', stroke = 0.5) +
          #geom_point(data = Dual_NodeList, aes(x, y, fill = clusters), size = 4, shape = 21) +
          geom_point(data = border_NodeList, aes(x, y, fill = clusters), size = 3.5, shape = 21, stroke = 0.8) +
          #geom_point(data = nl_Dual_NodeList, aes(x, y, fill = 'black'), size = 3, shape = 21, stroke = 1.5) +
          
          theme(legend.position = 'none') +
          coord_fixed(ratio = 1) +
          scale_fill_manual(values = c('1' = '#0078ba','2' = '#afcae8','3' = '#ff7700', 
                                       '4' = '#ffb96b', '5' = '#00a400', '6' = '#84e080', 
                                       '7' = '#a858cd','8' = '#b88fcc', '9' = '#d50b24', 
                                       '10' = '#fe969a')) 
        
        
        p
        
        
        
        #ggsave(p, file= paste0("Figures/Cluster_boundary/MB-0158-Clusters.png"), width = 15, height = 15, units = "in", dpi = 300)
        
        #ggsave(p, file= paste0("Figures/Cluster_boundary/", mid, '/', 'CN', pheno1, '_', 'CN', pheno2,  ".png"), width = 15, height = 15, units = "in", dpi = 300)
        
        
      }
    }
    
  }
}




clinical <- validation_data$clinicalDat
clinical$Survival <- ifelse(clinical$`Survival_days_capped*` > median(clinical$`Survival_days_capped*`), 'Survival high', 'Survival low')

allpairs_unique <- unique(allpairs)



colnames(CN_interactions_to)[1] <- 'Key'
colnames(CN_interactions_from)[1] <- 'Key'
colnames(CN_interactions)[1] <- 'Key'


for(rowid in 1:nrow(allpairs_unique)){
  #rowid <- 1
  phe1 <- allpairs_unique[rowid, 1]
  phe2 <- allpairs_unique[rowid, 2]
  
  phe1 <- 3
  phe2 <- 9
  
  from_pheno <- CN_interactions_from[CN_interactions_from$pheno1 == phe1 & CN_interactions_from$pheno2 == phe2,]
  to_pheno <- CN_interactions_to[CN_interactions_to$pheno1 == phe1 & CN_interactions_to$pheno2 == phe2,]
  all_pheno <- CN_interactions[CN_interactions$pheno1 == phe1 & CN_interactions$pheno2 == phe2,]
  
  from_surv <- merge(from_pheno, clinical, by = 'Key')
  to_surv <- merge(to_pheno, clinical, by = 'Key')
  all_surv <- merge(all_pheno, clinical, by = 'Key')
  
  
  
  
  
  print('-----------')
  for(id in 5:19){
    try({
      id <- 7
      #from_surv$Group <- ifelse(from_surv[,id] > median(from_surv[,id]), 'High', "Low")
      #from_surv
      all_surv$Group <- ifelse(all_surv[,id] > median(all_surv[,id]), 'High', 'Low')
      fit <- survfit(Surv(`Survival_days_capped*`, Censored) ~ Group, data = all_surv)
      p <- ggsurvplot(fit, data = all_surv, 
                      palette = c('#9dc9e6', '#787878'),
                      legend = 'none',
                      legend.title = '',
                      #legend.labs = c('', '),
                      risk.table = TRUE,
                      surv.scale = 'percent',
                      font.tickslab = c(26),
                      font.title = c(26),
                      font.x = c(28),
                      font.y = c(28),
                      font.legend = c(16),
                      fontsize = 10,
                      tables.theme = theme(axis.text = element_text(size = 16),
                                           axis.title = element_text(size = 16),
                                           title = element_text(size = 14)),
                      risk.table.y.text = FALSE,
                      conf.int = FALSE,
                      size = 1.5,
                      censor.size = 8,
                      pval = TRUE
      ) +
        xlab('Time, (months)')
      p
      
      
      
      
      test1 <- all_surv[all_surv$Survival == 'Survival high', id]
      test2 <- all_surv[all_surv$Survival == 'Survival low', id]
      
      
      #fit <- coxph(Surv(Survival, Status) ~ Group, data = from_surv)
      #pval <- summary(fit)$coefficient[,5]
      
      stattest <- t.test(test1, test2)
      if(stattest$p.value < 0.05){
        print(phe1)
        print(phe2)
        print(colnames(to_surv)[id])
        print(stattest$p.value)
        print(paste('Number of test 1:', length(test1)))
        print(paste('Number of test 2:', length(test2)))
        
      }
    }, silent = TRUE)  
  }
  
}


#---------------------------------------------------------------------#
#-------------- Plot Networks for Treg - Tumor - Myo axis ------------#
#---------------------------------------------------------------------#



for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  print(mid)
  #mid <- 'MB-0211'

  Tumor <- sg_cell_clus %>%
    dplyr::filter(metabric_id == mid) %>%
    dplyr::filter(Phenotype == 'Tumor') %>%
    select(Location_Center_X, Location_Center_Y)
  
  
  # Endothelial cells
  Treg <- sg_cell_clus %>%
    dplyr::filter(metabric_id == mid) %>%
    dplyr::filter(Phenotype == "Tregs and Tex") %>%
    select(Location_Center_X, Location_Center_Y)

  # Endothelial cells
  Myo <- sg_cell_clus %>%
    dplyr::filter(metabric_id == mid) %>%
    dplyr::filter(Phenotype == "Myofibroblasts") %>%
    select(Location_Center_X, Location_Center_Y)
  
  if(nrow(Myo) > 5 & nrow(Treg) > 5 & nrow(Tumor) > 5){
    
    # assign node to data frame
    Treg$node <- seq(nrow(Treg))
    Myo$node <- seq(nrow(Myo))
    Tumor$node <- seq(nrow(Tumor))
    
    
    
    Tumor_Myo <- nn2(Myo[,1:2], Treg[,1:2], k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.idx
    
    Tumor_Treg <- nn2(Tumor[,1:2], Treg[,1:2], k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.idx
    
    
    Tumor_Myo_EdgeList <- cbind.data.frame(seq(nrow(Treg)), Tumor_Myo)
    colnames(Tumor_Myo_EdgeList) <- c('from', 'to')
    Tumor_Treg_EdgeList <- cbind.data.frame(seq(nrow(Treg)), Tumor_Treg)
    colnames(Tumor_Treg_EdgeList) <- c('from', 'to')
    
    
    # remove not nearest cells
    # "from" is alwats tumor cells 
    Tumor <- Tumor %>%
      dplyr::filter(node %in% Tumor_Treg_EdgeList$to)
    Myo <- Myo %>%
      dplyr::filter(node %in% Tumor_Myo_EdgeList$to)
    
    
    p <- ggplot() +
      #theme_bw() +
      theme_void() +
      geom_segment(aes(x = Treg$Location_Center_X,
                       xend = Myo[match(Tumor_Myo_EdgeList$to, Myo$node), 1],
                       y = Treg$Location_Center_Y,
                       yend = Myo[match(Tumor_Myo_EdgeList$to, Myo$node), 2]), size = 0.3) +
      geom_segment(aes(x = Treg$Location_Center_X,
                       xend = Tumor[match(Tumor_Treg_EdgeList$to, Tumor$node), 1],
                       y = Treg$Location_Center_Y,
                       yend = Tumor[match(Tumor_Treg_EdgeList$to, Tumor$node), 2]), size = 0.3) +
      #geom_point(data = Dual_NodeList, aes(x, y), size = 3, shape = 21, fill = 'grey', stroke = 0.5) +
      geom_point(data = Tumor, aes(Location_Center_X, Location_Center_Y, fill = 'Tumor'), size = 4, shape = 21) +
      geom_point(data = Myo, aes(Location_Center_X, Location_Center_Y, fill = 'Myofibroblasts'), size = 4, shape = 21) +
      geom_point(data = Treg, aes(Location_Center_X, Location_Center_Y, fill = 'Tregs and Tex'), size = 4, shape = 21) +
      
      #geom_point(data = border_NodeList, aes(x, y, fill = clusters), size = 3.5, shape = 21, stroke = 0.8) +
      #geom_point(data = nl_Dual_NodeList, aes(x, y, fill = 'black'), size = 3, shape = 21, stroke = 1.5) +
      
      theme(legend.position = 'none') +
      coord_fixed(ratio = 1) +
      scale_fill_manual(values = c('Tumor' = '#e31a1c','Tregs and Tex' = '#6594e9',
                                   'Myofibroblasts' = '#b9b8b8')) 
    
    
    p
    ggsave(p, file= paste0("Figures/STAIN_Network/", mid, ".png"), width = 15, height = 15, units = "in", dpi = 300) 
  }
  
  
}









