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



unique(selectedDF$`Patient ID`)

#####################################################
#---------------------------------------------------#
#------------ ANALYSIS BEGIN HERE ------------------#
#---------------------------------------------------#
#####################################################
TA_area <- data.frame(matrix(nrow = 0, ncol = 0))
#write.csv(TA_area, 'TA_area.csv')
for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  #mid <- 'MB-3277'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid)
  
  TA <- ifelse(max(sg_metabricobj$Location_Center_X) > 1000, pi*(0.5)^2, pi*(0.3)^2)
  
  TA_area <- rbind.data.frame(TA_area, cbind.data.frame(mid, TA))
}
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








#####################################################
#                 Supplementary Figure 1            #
#####################################################



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



sg_cell_USETHIS_Cor <- sg_cell_USETHIS %>%
  dplyr::filter(Phenotype != 'Tumor') %>%
  dplyr::filter(Phenotype != 'Ki67+') %>%
  dplyr::filter(Phenotype != 'Myofibroblasts') %>%
  dplyr::filter(Phenotype != 'Endothelial') %>%
  dplyr::filter(Phenotype != 'Fibroblasts') %>%
  group_by(metabric_id, Phenotype) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n) %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  mutate(across(everything(), ~ . / rowSums(across(everything())))) %>%
  merge(metabricobj_im, by = 'metabric_id') %>%
  select(-c('count', 'X', 'TA')) %>%
  pivot_longer(cols = 2:11, names_to = "variable", values_to = "value")


# Plot using ggplot
p <- ggplot(sg_cell_USETHIS_Cor, aes(x = value, y = density)) +
  geom_point(size = 4, shape =21, fill = 'gray', color = 'black') +
  geom_smooth(method = "lm", se = FALSE, color = 'black') +
  facet_wrap(~ variable, scales = "free", nrow = 3) +
  labs(y = "Immune cell density",
       x = 'Fraction to immune cells') +
  theme_prism() +
  theme(strip.text = element_text(size = 13))  

p
ggsave(p, file=paste0("Figures/Supplementary figure S1_1.png"), width = 12, height = 9, units = "in", dpi = 300)

cor_results <- sg_cell_USETHIS_Cor %>%
  group_by(variable) %>%
  dplyr::summarize(
    cor_value = cor(value, density),
    p_value = cor.test(value, density)$p.value,
    .groups = 'drop'
  )




sg_cell_USETHIS_Cor <- sg_cell_USETHIS %>%
  dplyr::filter(Phenotype != 'Tumor') %>%
  dplyr::filter(Phenotype != 'Ki67+') %>%
  dplyr::filter(Phenotype != 'Myofibroblasts') %>%
  dplyr::filter(Phenotype != 'Endothelial') %>%
  dplyr::filter(Phenotype != 'Fibroblasts') %>%
  group_by(metabric_id, Phenotype) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n) %>%
  mutate(across(everything(), ~ ifelse(is.na(.), 0, 1))) %>%
  as.data.frame()

presence_mat <- sg_cell_USETHIS_Cor[,2:11]
rownames(presence_mat) <- sg_cell_USETHIS_Cor$metabric_id

results <- list()

for(i in 1:(ncol(presence_mat)-1)) {
  for(j in (i+1):ncol(presence_mat)) {
    col1 <- presence_mat[,i]
    col2 <- presence_mat[,j]
    
    test_result <- chisq.test(table(col1, col2))
    
    pair_name <- paste0(names(presence_mat)[i], "_vs_", names(presence_mat)[j])
    results[[pair_name]] <- test_result
  }
}

# Display results
results



png(p, file=paste0("Figures/Supplementary figure S1_2.png"), width = 16, height = 8, units = "in", res = 300)
Heatmap(t(presence_mat), 
        name = "Value",
        show_row_names = TRUE, 
        show_column_names = TRUE,
        show_row_dend = FALSE,
        show_column_dend = FALSE,
        #row_title = "ID",
        col = colorRamp2(c(min(presence_mat), max(presence_mat)), c("#4026a5", "#f3f933")),
        row_names_gp = gpar(fontsize = 16),  # Adjust font size for row labels
        column_names_gp = gpar(fontsize = 16),  # Adjust font size for column labels
        row_title_gp = gpar(fontsize = 18),  # Adjust font size for row title
        column_title_gp = gpar(fontsize = 18),  # Ad
        show_heatmap_legend = FALSE  # Do not display the legend
        # Do not display the legend
        #        cluster_columns = FALSE
)
dev.off()









#####################################################
#                 Supplementary Figure 2            #
#####################################################

selectedDF_dup <- selectedDF
# CD8 / Tumor cell ratio
colnames(selectedDF_dup)[2] <- 'metabric_id'



selectedDF_dup$Stage <- ifelse(selectedDF_dup$Stage <= 2, 'Stage I-II', 'Stage III-IV')
selectedDF_dup$Age <- ifelse(selectedDF_dup$`Age at Diagnosis` <= 55, '<= 55 years of age', '>55 years of age')
selectedDF_dup$Grade <- ifelse(selectedDF_dup$Grade <= 2, 'Grade I-II', 'Grade III')
selectedDF_dup$`Relapse Free Status` <- as.numeric(ifelse(selectedDF_dup$`Relapse Free Status` == '1:Recurred', '1', '0'))
selectedDF_dup$Survival_bin <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')
selectedDF_dup$TMB <- ifelse(selectedDF_dup$`TMB (nonsynonymous)` < 5, 'TMB low', 'TMB high')
selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'High'] <- 'Cellularity high'
selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'Moderate'] <- 'Cellularity moderate'
selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'Low'] <- 'Cellularity low'
#selectedDF_dup$`Relapse Free Status (Months)`



fit <- survfit(Surv(`Relapse Free Status (Months)`, `Relapse Free Status`) ~ Cellularity, data = selectedDF_dup)
fit

pdf(file= paste("Figures/Supplementary figure S2_11" , '.pdf'), width = 9, height = 10)

p <- ggsurvplot(fit, data = selectedDF_dup, 
                palette = c('blue', '#9dc9e6', '#787878'),
                #legend = 'none',
                legend.title = '',
                #legend.labs = c('Survival high', 'Survival low'),
                legend.labs = c('Cellularity high', 'low', 'moderate'),
                risk.table = TRUE,
                risk.table.heigh = 0.2,
                surv.scale = 'percent',
                font.tickslab = c(26),
                font.title = c(26),
                font.x = c(28),
                font.y = c(28),
                font.legend = c(24),
                fontsize = 10,
                tables.theme = theme(axis.text = element_text(size = 16),
                                     axis.title = element_text(size = 16),
                                     title = element_text(size = 14)),
                risk.table.y.text = FALSE,
                risk.table.x.text = c(28),
                conf.int = FALSE,
                size = 1.5,
                censor.size = 8,
                pval = TRUE,
                pval.size = 12
) +
  xlab('Time, (months)') 
p$plot <- p$plot + 
  ylab("Relapse-free survival probability")
p$table <- p$table + theme(plot.title = element_text(size=26),
                           axis.text.x = element_text(size = 26),
                           axis.title.x = element_blank())
p
dev.off()
#$plot <- p$plot + theme(legend.text = element_text(size = 10))





#####################################################
#                 Supplementary Figure 3            #
#####################################################
selectedDF_dup <- selectedDF
# CD8 / Tumor cell ratio
colnames(selectedDF_dup)[2] <- 'metabric_id'
tissue_area <- read.csv('TA_area.csv')


selectedDF_dup$Stage <- ifelse(selectedDF_dup$Stage <= 2, 'Stage I-II', 'Stage III-IV')
selectedDF_dup$`Age at Diagnosis` <- ifelse(selectedDF_dup$`Age at Diagnosis` <= 55, '<= 55 years of age', '>55 years of age')
selectedDF_dup$Grade <- ifelse(selectedDF_dup$Grade <= 2, 'Grade I-II', 'Grade III')
selectedDF_dup$`Relapse Free Status` <- ifelse(selectedDF_dup$`Relapse Free Status` == '1:Recurred', 'Recurred', 'Not recurred')
selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')
selectedDF_dup$`TMB (nonsynonymous)` <- ifelse(selectedDF_dup$`TMB (nonsynonymous)` < 5, 'TMB low', 'TMB high')
selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'High'] <- 'Cellularity high'
selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'Moderate'] <- 'Cellularity moderate'
selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'Low'] <- 'Cellularity low'

#selectedDF_dup$`Relapse Free Status (Months)`



cell_to_evaluate <- sg_cell_USETHIS %>%
  group_by(metabric_id, Phenotype) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  merge(tissue_area, by = 'metabric_id') %>%
  select(-c(X)) %>%
  replace(is.na(.), 0) %>%
  mutate(across(2:16, ~ .x / TA)) %>%
  pivot_longer(cols = 2:16, names_to = 'CellType') %>%  
  select(-TA) 

clinical_to_evaluate <- selectedDF_dup %>%
  select(c(metabric_id, `Age at Diagnosis`, Stage, Grade, `Relapse Free Status`,
           Survival, `TMB (nonsynonymous)`)) %>%
  
  pivot_longer(cols = 2:7, names_to = 'ClinicalType') 


to_evaluate <- merge(cell_to_evaluate, clinical_to_evaluate, by = 'metabric_id') %>%
  as.data.frame() %>%
  filter(complete.cases(.)) %>%
  select(-c(ClinicalType)) %>%
  `colnames<-` (c('metabric_id', 'CellType', 'Density', 'ClinicalType'))



#colnames(cell_to_evaluate)[2] <- 'value'


#cell_to_evaluate$value <- factor(cell_to_evaluate$value, levels = c('<= 55 years of age',  '>55 years of age', 'Cellularity low', 
#                                                                  'Cellularity moderate', 'Cellularity high',
#                                                                  'Basal', 'claudin-low', 'Her2', 'Grade I-II', 'Grade III',
#                                                                  'Stage I-II', 'Stage III-IV', 'Left', 'Right', 'Not recurred', 'Recurred',
#                                                                  'TMB low', 'TMB high', 'Survival low', 'Survival high'))

p <- ggplot(to_evaluate, aes(x = as.factor(ClinicalType), y = Density)) +
  theme_bw() +
  stat_boxplot(geom='errorbar', linetype=1, width=0.3) +
  geom_boxplot(outlier.size = -1) +
  geom_jitter(shape = 21, size = 2, aes(fill = as.factor(ClinicalType)), alpha = 0.6) +
  facet_wrap(~ CellType, scales = "free", ncol = 5) +
  theme(axis.text = element_text(size = 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        panel.spacing.x = unit(0, "pt"),
        legend.position = 'none',
        axis.title = element_text(size = 35),
        strip.text = element_text(size = 25)) +
  ylab(expression('Density,mm'^-2)) +
  xlab('Clinical subpopulations')


p
ggsave(p, file=paste0("Figures/allcells_heterogeneity_subgroup.png"), width = 30, height = 20, units = "in", dpi = 300)



#####################################################
#                 Supplementary Figure 4            #
#####################################################


library(ClusterR)

Fibroblasts <- sg_cell_USETHIS %>%
  dplyr::filter(Phenotype == 'Fibroblasts') %>%
  select(metabric_id, ImageNumber, CellID, SMA, FSP1, PDGFRB, `Caveolin-1`, Location_Center_X, Location_Center_Y)

nFeature <- 5

set.seed(2)
km <- MiniBatchKmeans(as.matrix(Fibroblasts[,4:7]), clusters = 4, batch_size = 100)
clusters <- predict_MBatchKMeans(as.matrix(Fibroblasts[,4:7]), km$centroids)



Fibroblasts_clus <- Fibroblasts

# assign clusters
Fibroblasts_clus$Phenotype <- as.factor(clusters)





# Panel A



boundaryFiles <- list.files('./SingleCell_Boundary') %>%
  data.frame() %>%
  `colnames<-` (c('filename'))





# get the single cell data for the current metabric object
for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  #mid <- 'MB-0149'
  
  sg_metabricobj <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid)
  
  
  
  # get the serial number to match the boundary data
  serial <- paste0("MB", strsplit(mid, '-')[[1]][2], '_', unique(sg_metabricobj$ImageNumber))
  
  
  # read single cell boundary data
  sg_boundary <- read.csv(paste0('./SingleCell_Boundary/', 
                                 boundaryFiles[boundaryFiles$filename %like% serial,]))
  
  
  NewID_vec <- c()
  RemoveID_vec <- c()
  for(pl in seq_len(max(sg_boundary$CellID))){
    
    
    
    polygon_pl <- sg_boundary %>%
      dplyr::filter(CellID == pl)
    
    x_range <- range(polygon_pl$X)[2] - range(polygon_pl$X)[1]
    y_range <- range(polygon_pl$Y)[2] - range(polygon_pl$Y)[1]
    
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
    
    if(x_range > 100 & y_range > 100){
      RemoveID_vec <- c(RemoveID_vec, pl)
    }
    
  }
  
  
  
  testxx <- NewID_vec %>%
    data.frame() %>%
    mutate(falseid = seq_len(max(sg_boundary$CellID))) %>%
    replace(is.na(.), 0) %>%
    `colnames<-` (c('TrueID', 'CellID'))
  
  
  
  test <- merge(sg_boundary, testxx, by = 'CellID') %>%
    dplyr::filter(!(CellID %in% RemoveID_vec)) %>%
    dplyr::select(Y, X, TrueID) %>%
    `colnames<-` (c('Y', 'X', 'CellID')) 
  
  
  
  
  # assign cell type to boundary file
  sg_boundary_w_ctype <- merge(sg_metabricobj, test, by = 'CellID') %>%
    dplyr::select(X, Y, CellID, 'Phenotype', ImageNumber, metabric_id)  
  
  
  sg_boundary_w_ctype <- sg_boundary_w_ctype %>%
    left_join(Fibroblasts_clus %>%
                select(metabric_id, CellID, ImageNumber, Phenotype) %>%
                rename_at(vars(Phenotype), ~paste0("New", .)),
              by = c("metabric_id", "CellID", "ImageNumber")) %>%
    mutate(Phenotype = ifelse(is.na(NewPhenotype), Phenotype, NewPhenotype)) %>%
    select(-NewPhenotype)
  
  
  p <- ggplot() +
    theme_void() +
    geom_polygon(data = sg_boundary_w_ctype, aes(X, Y, group = CellID, fill = Phenotype), size = 0.1) +
    #geom_point(data = sg_metabricobj, aes(Location_Center_X, Location_Center_Y, color = Phenotype), size = 0.1) +
    
    scale_fill_manual(values = c('Tumor' = '#b9b8b8','Ki67+' = '#b9b8b8','Tregs and Tex' = '#b9b8b8', 'CD4+ T cells & APCs' = '#b9b8b8','CD4+ T cells' = '#b9b8b8', 
                                 'CD38+ lymphocytes' = '#b9b8b8', 'CD57+' = '#b9b8b8',
                                 'B cells' = '#b9b8b8', 
                                 'CD8+ T cells' = '#b9b8b8', 
                                 'Macrophages' = '#b9b8b8', 
                                 'Macrophages & granulocytes' = '#b9b8b8', 
                                 'Granulocytes' = '#b9b8b8',
                                 'Endothelial' = '#b9b8b8',
                                 'Fibroblasts' = '#b9b8b8',
                                 'Myofibroblasts' = '#b9b8b8',
                                 '1' = '#67a61f', '3' = '#e4191c',
                                 '2' = '#1c79b5', '4' = '#ec8a48')) +
    theme(legend.position = 'NA') +
    ylim(max(sg_boundary$Y) + 5,0) +
    xlim(max(sg_boundary$X) + 5,0) +
    coord_fixed(ratio = 1) 
  p
  
  
  # Panel B
  
  
  ggsave(p, file=paste0("Figures/FibroMap/", mid, ".png"), width = 8, height = 8, units = "in", dpi = 300)
  
}

# facet plot

sg_boundary_w_ctype_fibro <- sg_boundary_w_ctype %>%
  dplyr::filter(Phenotype == 1 |
                  Phenotype == 2 |
                  Phenotype == 3 |
                  Phenotype == 4)

sg_boundary_w_ctype_fibro[sg_boundary_w_ctype_fibro$Phenotype == 1, 'Phenotype'] <- 'Fibro-C1'
sg_boundary_w_ctype_fibro[sg_boundary_w_ctype_fibro$Phenotype == 2, 'Phenotype'] <- 'Fibro-C3'
sg_boundary_w_ctype_fibro[sg_boundary_w_ctype_fibro$Phenotype == 3, 'Phenotype'] <- 'Fibro-C2'
sg_boundary_w_ctype_fibro[sg_boundary_w_ctype_fibro$Phenotype == 4, 'Phenotype'] <- 'Fibro-C4'

p <- ggplot() +
  theme_bw() +
  geom_polygon(data = sg_boundary_w_ctype_fibro, aes(X, Y, group = CellID, fill = Phenotype), size = 0.1) +
  facet_wrap(~ Phenotype, scales = "free", ncol = 2) +
  
  #geom_point(data = sg_metabricobj, aes(Location_Center_X, Location_Center_Y, color = Phenotype), size = 0.1) +
  
  scale_fill_manual(values = c('Fibro-C1' = '#67a61f', 'Fibro-C2' = '#e4191c',
                               'Fibro-C3' = '#1c79b5', 'Fibro-C4' = '#ec8a48')) +
  theme(panel.spacing.x = unit(0, "pt"),
        legend.position = 'none',
        strip.text = element_text(size = 25),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) 
p

ggsave(p, file="Figures/Fibro_MB0149_Split.png", width = 18, height = 8, units = "in", dpi = 300)

# Panel C



Fibroblasts_clus_long <- Fibroblasts_clus %>%
  select(SMA, FSP1, PDGFRB, `Caveolin-1`, Phenotype) %>%
  pivot_longer(cols = 1:4, names_to = 'markers')

Fibroblasts_clus_long$Phenotype <- as.character(Fibroblasts_clus_long$Phenotype)
Fibroblasts_clus_long[Fibroblasts_clus_long$Phenotype == 1, 'Phenotype'] <- 'Fibro-C1'
Fibroblasts_clus_long[Fibroblasts_clus_long$Phenotype == 2, 'Phenotype'] <- 'Fibro-C3'
Fibroblasts_clus_long[Fibroblasts_clus_long$Phenotype == 3, 'Phenotype'] <- 'Fibro-C2'
Fibroblasts_clus_long[Fibroblasts_clus_long$Phenotype == 4, 'Phenotype'] <- 'Fibro-C4'


p <- ggplot(Fibroblasts_clus_long, aes(x = Phenotype, y = value), fill = as.factor(Phenotype)) +
  theme_bw() +
  stat_boxplot(geom='errorbar', linetype=1, width=0.3) +
  geom_boxplot(outlier.size = 1, aes(color = Phenotype), size = 1) +
  scale_color_manual(values = c('Fibro-C1' = '#67a61f', 'Fibro-C2' = '#e4191c',
                                'Fibro-C3' = '#1c79b5', 'Fibro-C4' = '#ec8a48')) +
  facet_wrap(~ markers, scales = "free", ncol = 4) +
  theme(axis.text = element_text(size = 25),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1),
        panel.spacing.x = unit(0, "pt"),
        legend.position = 'none',
        axis.title = element_text(size = 35),
        strip.text = element_text(size = 25)) +
  ylab('Expression') +
  xlab('')


p
ggsave(p, file="Figures/Fibro_Marker_expression.png", width = 18, height = 8, units = "in", dpi = 300)





# Panel D


Fibroblasts_clus_sub <- Fibroblasts_clus %>%
  group_by(metabric_id, ImageNumber, Phenotype) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n)

nonFibro <- sg_cell_USETHIS %>%
  dplyr::filter(Phenotype != 'Fibroblasts') %>%
  group_by(metabric_id, ImageNumber, Phenotype, .drop = FALSE) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n)


Fibro_Cor <- merge(Fibroblasts_clus_sub, nonFibro, by = 'metabric_id') %>%
  select(-c(metabric_id, ImageNumber.x, ImageNumber.y, )) %>%
  replace(is.na(.), 0) %>%
  pivot_longer(cols = 5:18, names_to = "Cell_Count", values_to = "CellCount_value") %>%
  pivot_longer(cols = 1:4, names_to = "Cluster_Count", values_to = "ClusterCount_value") %>%
  droplevels()


Fibro_Cor[Fibro_Cor$Cluster_Count == 1, 'Cluster_Count'] <- 'Fibro-C1'
Fibro_Cor[Fibro_Cor$Cluster_Count == 2, 'Cluster_Count'] <- 'Fibro-C3'
Fibro_Cor[Fibro_Cor$Cluster_Count == 3, 'Cluster_Count'] <- 'Fibro-C2'
Fibro_Cor[Fibro_Cor$Cluster_Count == 4, 'Cluster_Count'] <- 'Fibro-C4'

p <- ggplot(Fibro_Cor, aes(x = ClusterCount_value, y = CellCount_value)) +
  theme_bw() +
  geom_point(size = 2, shape =21, aes(fill = Cluster_Count), color = 'black') +
  geom_smooth(method = "lm", se = FALSE, color = 'black') +
  ggh4x::facet_grid2(Cluster_Count ~ Cell_Count, scales = "free", independent = 'all') +
  scale_color_manual(values = c('Fibro-C1' = '#67a61f', 'Fibro-C2' = '#e4191c',
                                'Fibro-C3' = '#1c79b5', 'Fibro-C4' = '#ec8a48')) +
  labs(y = "",
       x = '') +
  theme(strip.text = element_text(size = 12),
        panel.spacing.x = unit(0, "pt"),
        legend.position = 'none',
        axis.text = element_blank(),
        axis.ticks = element_blank())  

p
ggsave(p, file="Figures/Fibro_Cor.png", width = 25, height = 8, units = "in", dpi = 300)





#####################################################
#                 Supplementary Figure 5            #
#####################################################



# Plot Cell Boundaries for illustrative purposes
vesselFiles <- list.files('./Vessel_Boundary') %>%
  data.frame() %>%
  `colnames<-` (c('filename'))


boundaryFiles <- list.files('./SingleCell_Boundary') %>%
  data.frame() %>%
  `colnames<-` (c('filename'))

# read area data
tissue_area <- read.csv('TA_area.csv')

library(sp)

filter_vessel <- function(df, area_thresh) {
  # transfrom DataFrame to SpatialPolygonsDataFrame object
  #df <- pos_vessel
  df_split <- split(df, df$id)
  polygons <- lapply(df_split, function(df_group) {
    coords <- df_group[, c("x", "y")]
    p <- Polygon(coords)
    p <- Polygons(list(p), ID = df_group$id[1])
    return(p)
  })
  
  sp_polygons <- SpatialPolygons(polygons)
  areas <- gArea(sp_polygons, byid = TRUE)
  
  # filter potential artifacts
  area_filtered <- which(areas > area_thresh)
  
  #plot(df_test$x, df_test$y)
  
  return(area_filtered)
}




#------------- Vasculatre Profile -------------------------#
library(rgeos)
vasculature_area <- data.frame(matrix(nrow = 0, ncol = 0))

vasculature_density <- data.frame(matrix(nrow = 0, ncol = 0))
# factorize cell type
sg_cell_USETHIS$Phenotype <- as.factor(sg_cell_USETHIS$Phenotype)

selectedDF_dup <- selectedDF
colnames(selectedDF_dup)[2] <- 'metabric_id'
selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')

num_bins <- 200 
area_thresh <- 700 #500, 700, 800
bin_width <- (area_thresh - 0) / num_bins  

for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  print(mid)
  mid <- 'MB-3277'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid)
  
  
  # get the serial number to match the boundary data
  serial <- paste0("MB", strsplit(mid, '-')[[1]][2], '_', unique(sg_metabricobj$ImageNumber))
  
  
  # get the survival type
  Surv <- selectedDF_dup[selectedDF_dup$metabric_id == mid, 'Survival'] %>%
    as.character()
  # read vessel data
  sg_vessel <- read.csv(paste0('./Vessel_Boundary/', 
                               vesselFiles[vesselFiles$filename %like% serial,]))
  
  
  # read single cell boundary data
  sg_boundary <- read.csv(paste0('./SingleCell_Boundary/', 
                                 boundaryFiles[boundaryFiles$filename %like% serial,]))
  
  df_split <- split(sg_vessel, sg_vessel$id)
  polygons <- lapply(df_split, function(df_group) {
    coords <- df_group[, c("x", "y")]
    p <- Polygon(coords)
    p <- Polygons(list(p), ID = df_group$id[1])
    return(p)
  })
  
  sp_polygons <- SpatialPolygons(polygons)
  areas <- gArea(sp_polygons, byid = TRUE)
  
  pos_vessel <- holePolygon(sg_vessel)
  
  #---- remove vessels with very small sizes
  vesselmaskTRUE <- filter_vessel(pos_vessel, area_thresh)
  
  
  pos_vesselTRUE <- pos_vessel %>%
    dplyr::filter(id %in% vesselmaskTRUE)
  
  
  pos_vesselFALSE <- pos_vessel %>%
    dplyr::filter(!(id %in% vesselmaskTRUE))
  
  
  #----------- Vasculature density ----------------#
  #mid <- 'MB-0399'
  core_area <- tissue_area[tissue_area$metabric_id == mid, 'TA']
  
  # vasculature density
  #vasculature_density <- rbind.data.frame(vasculature_density, cbind.data.frame(mid, vasculature_area / core_area))
  vasculature_density <- rbind.data.frame(vasculature_density, cbind.data.frame(mid, length(vesselmaskTRUE) / core_area))
  #----------- Cell infiltration profile w.r.t vasculature ----------------#
  vasculature_area <- rbind.data.frame(vasculature_area, cbind.data.frame(mid, areas))
  
  
}

vasculature_area_pt <- vasculature_area %>%
  group_by(mid) %>%
  dplyr::summarise(mean_area = mean(areas, na.rm = TRUE)) %>%
  mutate(type = 'All patient')



p <- ggplot(vasculature_area_pt, aes(x = type, y = mean_area)) +
  theme_classic() +
  stat_boxplot(geom='errorbar', linetype=1, width=0.3, linewidth = 1) +
  geom_boxplot(outlier.size = -1, size = 1) +
  geom_jitter(aes(x = type), width = 0, shape = 21, size = 3 , alpha = 0.6, fill = 'black') +
  theme(axis.text = element_text(size = 25),
        panel.spacing.x = unit(0, "pt"),
        legend.position = 'none',
        axis.title = element_text(size = 35),
        strip.text = element_text(size = 25)) +
  ylab(expression('Vessel area,'~mu~'m'^2)) +
  xlab('')


p
ggsave(p, file="Figures/Vessel_area_distribution.png", width = 4, height = 12, units = "in", dpi = 300)







#--------------Vasculature density -----------------#
#vasculature_density <- vasculature_density[, -c(4,5)]
vasculature_density <- vasculature_density[complete.cases(vasculature_density),]


vasculature_density$type <- 'Sample'

colnames(vasculature_density) <- c('metabric_id', 'density', 'Sample')





selectedDF_dup <- selectedDF
# CD8 / Tumor cell ratio
colnames(selectedDF_dup)[2] <- 'metabric_id'


thresh <- median(vasculature_density$density)

#breaks <- quantile(vasculature_density$density, probs = c(0, 0.25, 0.5, 0.75, 1))
#breaks <- quantile(vasculature_density$density, probs = c(0, 0.33, 0.67, 1))
vasculature_DF <- merge(vasculature_density, selectedDF_dup, by = 'metabric_id')

vasculature_DF$vd_level <- ifelse(vasculature_DF$density > thresh, 'high', 'low')

#vasculature_DF$vd_level <- as.integer(cut(vasculature_DF$density, breaks = breaks, include.lowest = TRUE, labels = FALSE))

vasculature_DF$relapse_status <- ifelse(vasculature_DF$`Relapse Free Status` == '1:Recurred', 1, 0)


#test1 <- vasculature_DF[vasculature_DF$vd_level == 2, 'Survival']
#test2 <- vasculature_DF[vasculature_DF$vd_level == 3, 'Survival']

#mean(test1)
#mean(test2)
#wilcox.test(test1, test2)
#test1 <- vasculature_DF[vasculature_DF$`Age at Diagnosis` == '<= 55 years of age', 'density']
#test2 <- vasculature_DF[vasculature_DF$`Age at Diagnosis` == '>55 years of age', 'density']
#vasculature_DF <- vasculature_DF %>%
#  dplyr::filter(vd_level != 1)

#vasculature_DF$`Relapse Free Status (Months)`
#t.test(test1, test2)
#fit <- survfit(Surv(`Relapse Free Status (Months)`, relapse_status) ~ vd_level, data = vasculature_DF)
fit <- survfit(Surv(Survival, Status) ~ vd_level, data = vasculature_DF)
# Perform the log-rank test
result <- survdiff(Surv(`Relapse Free Status (Months)`, relapse_status) ~ vd_level, data = vasculature_DF)
#result <- survdiff(Surv(Survival, Status) ~ vd_level, data = vasculature_DF)

# Extract and print the p-value
p_value <- 1 - pchisq(result$chisq, length(result$n) - 1)
print(p_value)



p <- ggsurvplot(fit, data = vasculature_DF, 
                palette = c('#c22828', '#325698'),
                legend = 'none',
                #legend.title = '',
                #legend.labs = c('high', 'low'),
                #risk.table = TRUE,
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
                #pval = TRUE,
                #break.time.by = 2000
) +
  xlab('Time, (Months)') +
  ylab('Overall survival')
p
ggsave(print(p), file=paste0("Figures/Vasculature_Density_OS_thresh300——1.png"), width = 8, height = 6, units = "in", dpi = 300, bg="black")





p <- ggplot() +
  geom_polygon(data = pos_vesselTRUE, aes(x,y, group = id), fill = '#709fce') +
  geom_polygon(data = pos_vesselFALSE, aes(x,y, group = id), fill = 'grey') +
  theme_void()
p
ggsave(p, file=paste0("Figures/Thresh700_filteredMask.jpg"), width = 10, height = 10, units = "in", dpi = 300, bg="black")











###########################################################
#                 Supplementary Figure 6                  #
###########################################################


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



library(ggvoronoi)
print(mid)
mid <- 'MB-2857'
# get the single cell data for the current metabric object
sg_metabricobj <- sg_cell_clus %>%
  dplyr::filter(metabric_id == mid)


p <- ggplot(sg_metabricobj, aes(Location_Center_X, Location_Center_Y, fill = clusters, color = 'white')) +
  geom_voronoi(color = 'white', size = 0.2) +
  theme_void() +
  theme(legend.position = 'none') +
  scale_fill_manual(values = c('1' = '#0078ba','2' = '#afcae8','3' = '#ff7700', 
                               '4' = '#ffb96b', '5' = '#00a400', '6' = '#84e080', 
                               '7' = '#a858cd','8' = '#b88fcc', '9' = '#d50b24', 
                               '10' = '#fe969a')) 
p 
ggsave(p, file=paste0("Figures/MB", mid, "_Voronoi.jpg"), width = 10, height = 10, units = "in", dpi = 300, bg="black")



#---------------------#
#      Panel C 
#


aType <- unique(sg_cell_clus$clusters)

compHeatmap <- data.frame(matrix(nrow = 9, ncol = 10))


colorHeatmap <- data.frame(matrix(nrow = 9, ncol = 10))
colorDictionary <- data.frame(`A term` = c('#0868ac', '#ece7f2', '#ece7f2', '#a6bedb', '#9ebedb', '#8c96c6', '#4d0149', '#8c6bb2', '#31a353'),
                              `B term` = c('#a8deb5', '#3690c0', '#4db4d4', '#ffebe3', '#ffffcc', '#fdc5c0', '#7bccc4', '#fa9fb5', '#c2e699'))

sizeHeatmap <- data.frame(matrix(nrow = 9, ncol = 10))

colnames(compHeatmap) <- aType
colnames(colorHeatmap) <- aType
colnames(sizeHeatmap) <- aType


rowVariables <- c('<= 55 years of age',  '>55 years of age', 'Cellularity low', 
                  'Cellularity moderate', 'Cellularity high', 'skip', 'Grade I-II', 'Grade III',
                  'Stage I-II', 'Stage III-IV', 'Left', 'Right', 'Not recurred', 'Recurred',
                  'TMB low', 'TMB high', 'Survival low', 'Survival high')

selectedDF_dup <- selectedDF
colnames(selectedDF_dup)[2] <- 'metabric_id'

selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')
DF <- sg_cell_clus %>%
  group_by(metabric_id, clusters) %>%
  tally() %>%
  pivot_wider(names_from = clusters, values_from = n) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  dplyr::select(-metabric_id) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate(across(everything(), ~ . * 100/ row_sum)) %>%
  replace(is.na(.), 0) %>%
  mutate(metabric_id = selectedDF_dup$metabric_id) %>%
  merge(selectedDF_dup) %>%
  as.data.frame() 

t1 <- DF[DF$Survival == 'Survival low', '5'] 
  
t2 <- DF[DF$Survival == 'Survival high', '5']

wilcox.test(t1, t2)
Heatmap_rid <- 1

for(ct in aType){
  #ct <- 'Fibroblasts'
  
  #ct <- 'CD8+ T cells'
  rv <- 1
  while(rv <= length(rowVariables)-1){
    
    #rv = 9
    if(rv == 5){
      #rvid=3
      rv_1 <- rowVariables[rv-2]
      rv_2 <- rowVariables[rv]
      
      # for student t test
      
      rv_1_data <- DF %>%
        filter(apply(., 1, function(row) any(sapply(rv_1, grepl, row)))) %>%
        dplyr::select(all_of(ct)) %>%
        as.matrix()
      
      rv_1_data <- rv_1_data[rv_1_data!=0]
      
      rv_2_data <- DF %>%
        filter(apply(., 1, function(row) any(sapply(rv_2, grepl, row)))) %>%
        dplyr::select(all_of(ct)) %>%
        as.matrix()
      rv_2_data <- rv_2_data[rv_2_data!=0]
      
      Heatmap_rid <- Heatmap_rid + 1
    }
    if(rv != 5){
      rv_1 <- rowVariables[rv]
      rv_2 <- rowVariables[rv+1]
      
      # for student t test
      
      rv_1_data <- DF %>%
        filter(apply(., 1, function(row) any(grepl(all_of(rv_1), row)))) %>%
        dplyr::select(all_of(ct)) %>%
        as.matrix()
      rv_1_data <- rv_1_data[rv_1_data!=0]
      
      rv_2_data <- DF %>%
        filter(apply(., 1, function(row) any(grepl(all_of(rv_2), row))))%>%
        dplyr::select(all_of(ct)) %>%
        as.matrix()
      rv_2_data <- rv_2_data[rv_2_data!=0]
      
    }
    
    # statistical test
    if(length(rv_1_data) * length(rv_2_data) != 0){
      pval <- wilcox.test(rv_1_data, rv_2_data)$p.value
      if(pval == 'NaN'){
        pval <- 0.05
      }
    }else{
      pval <- 0.05
    }
    
    
    
    # assign color 
    if(pval >= 0.05){
      color <- 'grey'
      size <- 4
    }
    if(pval < 0.05){
      
      
      print(pval)
      print(rv)
      print(ct)
      
      if(mean(rv_1_data) > mean(rv_2_data)){
        color <- colorDictionary[(rv + 1)/2,1]
      }
      if(mean(rv_1_data) < mean(rv_2_data)){
        color <- colorDictionary[(rv + 1)/2,2]
      }
      
      size <- 8
      if(pval < 0.01){
        size <- 12
      }
      
    }
    print(rv)
    compHeatmap[(rv + 1)/2, ct] <- pval
    colorHeatmap[(rv + 1)/2, ct] <- color
    sizeHeatmap[(rv + 1)/2, ct] <- size
    
    # iterate to next identifier
    rv <- rv + 2
  }
  
}




melted_matrix <- as.data.frame(as.table(as.matrix(compHeatmap)))




melted_matrix$custom_color <- as.data.frame(as.table(as.matrix(colorHeatmap)))$Freq
melted_matrix$custom_size <- as.data.frame(as.table(as.matrix(sizeHeatmap)))$Freq
library(ggplot2)
melted_matrix$Var2 <- factor(melted_matrix$Var2, levels = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10"))

p <- ggplot(melted_matrix, aes(x = Var2, y = Var1, fill = custom_color, size = custom_size)) +
  geom_point(shape = 21) +
  scale_fill_identity() +
  scale_size_identity() +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) 
p

ggsave(p, file=paste0("Figures/CN_Freq_Comparison.png"), width = 17, height = 12, units = "in", dpi = 300)




###########################################################
#                 Supplementary Figure 7                  #
###########################################################
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

stattest_all$pvalue <- p.adjust(stattest_all$`stattest$p.value`)

p <- ggplot(stattest_all, aes(y = to_exclude, x = log(pvalue))) +
  theme_bw() +
  geom_vline(xintercept = log(0.05), linetype = 'dashed', color = 'red', size = 1) +
  geom_point(size = 5) +
  theme(axis.text.x = element_text(size = 30),
        axis.text.y = element_text(size = 25),
        axis.title = element_text(size = 32)) +
  ylab('') +
  xlab('log(p value)')
p
ggsave(p, file=paste0("Figures/Exclusion_test.png"), width = 12, height = 24, units = "in", dpi = 300)


 #------------------------------------------------------#
#------------- Validation -----------------------------#

library(ggplot2); library(ComplexHeatmap); library(circlize); library(ggthemes); library(ggprism)
library(tidyr); library(dplyr); library(plyr); library(readr); library(data.table)
library(RANN); library(usedist); library(vegan)
library(rgeos); library(igraph)


setwd(dirname(rstudioapi::getSourceEditorContext()$path))

setwd('..')
source("./Codes/Function.r")

#---------- Section 1: Count cell numbers per core per patient ----------------#
#---------- Select cores with at least 1000 cells, consistent with 
#---------- discovery cohort

cell_data <- read.csv('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/cells.csv') 


# remove patient images with less cells  
validPatients <- cell_data %>%
  #dplyr::filter(PatientID %in% clinical_data$PatientID) %>%
  group_by(PatientID, ImageNumber) %>%
  tally() %>%
  dplyr::filter(n >= 1000)



#-----------  Funciton -------------------#

filter_vessel <- function(df, area_thresh) {
  # transfrom DataFrame to SpatialPolygonsDataFrame object
  #df <- pos_vessel
  df_split <- split(df, df$id)
  polygons <- lapply(df_split, function(df_group) {
    coords <- df_group[, c("x", "y")]
    p <- Polygon(coords)
    p <- Polygons(list(p), ID = df_group$id[1])
    return(p)
  })
  
  sp_polygons <- SpatialPolygons(polygons)
  areas <- gArea(sp_polygons, byid = TRUE)
  
  # filter potential artifacts
  area_filtered <- which(areas > area_thresh)
  
  #plot(df_test$x, df_test$y)
  
  return(area_filtered)
}




#-------------- Section 2: Patient Selection ---------------#
# read clinical data
# select chemo & immunotherapy treated patients
# select patients that is processed as protocol
# select BASELINE patients for prediction
# select 20 from NR and 20 from R

set.seed(5)
clinical_data <- read.csv('../TNBC_ClinicalTrial_Validation_Cohort/clinical.csv') %>%
  dplyr::filter(#Arm == 'C&I',
    isPerProtocol == 'TRUE',
    #BiopsyPhase == 'Baseline',
    pCR %in% c("RD", "pCR")
  ) %>% # Exclude Post-treatment
  dplyr::filter(PatientID %in% validPatients$PatientID) #%>%
#group_by(pCR) %>% 
#sample_n(10) 

# select one core per patient
set.seed(1)

selecter_Cores <- validPatients %>%
  dplyr::filter(PatientID %in% clinical_data$PatientID) #%>%
#group_by(PatientID) %>% 
#sample_n(1) %>%
#select(PatientID, ImageNumber)


val_cell_USETHIS <- cell_data %>%
  inner_join(selecter_Cores, by = c("PatientID", "ImageNumber"))

val_cell_USETHIS$Phenotype <- ifelse(val_cell_USETHIS$isEpithelial, 'Tumor', val_cell_USETHIS$Label) # Discard the old classification criteria


val_cell_USETHIS_tr <- val_cell_USETHIS[, c(15:60)] 

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
  normalized_data <- apply(clipped_data, 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  #featuers
  
  return(normalized_data)
}

# 对 IMC 数据进行归一化
sg_cell_ft_normalized <- normalize_imc_data(val_cell_USETHIS_tr) %>%
  as.data.frame() 


val_cell_USETHIS[, c(15:60)] <- sg_cell_ft_normalized




ratios_all <- data.frame(matrix(nrow = 0, ncol = 0))
set.seed(97)
clinical_data <- read.csv('../TNBC_ClinicalTrial_Validation_Cohort/clinical.csv') %>%
  dplyr::filter(
    #Arm == 'C&I',
    isPerProtocol == 'TRUE',
    #BiopsyPhase == 'Post-treatment' | BiopsyPhase == "Baseline",
    BiopsyPhase == "Baseline"#,
    #pCR == 'RD'
  ) %>% # Exclude Post-treatment
  dplyr::filter(PatientID %in% validPatients$PatientID) %>%
  group_by(pCR) %>% 
  sample_n(30) 


imgnum_selected <- merge(clinical_data, val_cell_USETHIS, by = 'PatientID') %>%
  dplyr::filter(BiopsyPhase.y == 'Baseline') %>%
  #dplyr::filter(BiopsyPhase.y == 'Post-treatment') %>%
  select(ImageNumber) %>%
  unique() %>%
  as.vector()

for(iid in imgnum_selected$ImageNumber){
  
  #print(iid)
  #iid <- 849
  #mid <- 'MB-0316'
  # get the single cell data for the current metabric object
  sg_obj <- val_cell_USETHIS %>%
    dplyr::filter(ImageNumber == iid)
  
  Biopsy <- unique(sg_obj$BiopsyPhase)
  
  # get the pCR
  pCR <- merge(sg_obj, clinical_data, by = 'PatientID') %>%
    select(pCR) %>%
    unique() %>%
    as.matrix()
  # get CD8 T cells
  Tumor <- val_cell_USETHIS %>%
    dplyr::filter(ImageNumber == iid) %>%
    dplyr::filter(Phenotype == 'Tumor') %>%
    select(Location_Center_X, Location_Center_Y)
  
  
  # Endothelial cells
  Treg <- val_cell_USETHIS %>%
    dplyr::filter(ImageNumber == iid) %>%
    dplyr::filter(Phenotype == "Treg" | Phenotype == "CD8^+PD1^+T_{Ex}") %>%
    select(Location_Center_X, Location_Center_Y)
  
  # Endothelial cells
  Myo <- val_cell_USETHIS %>%
    dplyr::filter(ImageNumber == iid) %>%
    dplyr::filter(Phenotype == "Myofibroblasts") %>%
    select(Location_Center_X, Location_Center_Y)
  
  
  if(nrow(Tumor) > 5 & nrow(Treg) > 5 & nrow(Myo) > 5){
    Tumor_Treg <- nn2(Tumor, Treg, k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.dists %>%
      as.vector() 
    
    Tumor_Myo <- nn2(Myo, Treg, k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.dists %>%
      as.vector()
    
    ratio <- Tumor_Treg / (Tumor_Treg + Tumor_Myo)  
    #ratio <- sqrt(CD8T_Tumor * CD8T_FoxP3)
    ratios_all <- rbind.data.frame(ratios_all, cbind.data.frame(ratio, iid, Biopsy, pCR))  
  }
  
  
}


stattest_all <- data.frame(matrix(nrow = 0, ncol = 0))

#test$Survival <- ifelse(test$Survival < 66, 'Survival low', 'Survival high')
for(to_exclude in unique(ratios_all$iid)){
  
  #to_exclude <- 590
  test_excluded <- ratios_all %>%
    dplyr::filter(iid != to_exclude)
  
  
  low <- test_excluded[test_excluded$pCR == 'pCR', ]
  high <- test_excluded[test_excluded$pCR == 'RD', ]
  
  
  stattest <- wilcox.test(low$ratio, high$ratio)
  print(stattest$p.value)
  print(length(unique(low$metabric_id)))
  print(length(unique(high$metabric_id)))
  
  
  stattest_all <- rbind.data.frame(stattest_all, cbind.data.frame(to_exclude, stattest$p.value))
  
}

stattest_all$pvalue <- p.adjust(stattest_all$`stattest$p.value`)

p <- ggplot(stattest_all, aes(y = as.character(to_exclude), x = log(pvalue))) +
  theme_bw() +
  geom_vline(xintercept = log(0.05), linetype = 'dashed', color = 'red', size = 1) +
  geom_point(size = 5) +
  theme(axis.text.y = element_text(size = 25),
        axis.text.x = element_text(size = 30),
        axis.title = element_text(size = 32),
        #axis.text.y = element_text(angle = 90)
        ) +
  ylab('') +
  xlab('log(p value)')

p

ggsave(p, file=paste0("Figures/Exclusion_test_validation.png"), width = 12, height = 24, units = "in", dpi = 300)









t.test(c(0.635, 0.5275, 0.6950, 0.54, 0.585), c(0.6225, 0.6625, 0.6425, 0.7025, 0.605))
t.test(c(0.7325, 0.725, 0.6425, 0.7325, 0.70), c(0.6225, 0.6625, 0.6425, 0.7025, 0.605))
t.test(c(0.7325, 0.725, 0.6425, 0.7325, 0.70), c(0.635, 0.5275, 0.6950, 0.54, 0.585))

