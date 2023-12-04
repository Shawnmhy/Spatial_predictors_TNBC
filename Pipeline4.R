#-------------------------------------------#
# METABRIC TNBC processing pipeline 4-------#
# @Author: Haoyang Mi ----------------------#
# Date: May 3rd 2023------------------------#
# REF: 


library(ggplot2); library(ComplexHeatmap); library(circlize); library(ggthemes); library(ggprism)
library(tidyr); library(dplyr); library(plyr); library(readr); library(data.table)
library(survival); library(survminer); library(ade4); library(vegan)
library(RANN); library(Rtsne)
library(rgeos); library(igraph);library(ClusterR);library(ggvoronoi)


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

#----------- Calculate Difussion Distance ----------#
# For each cancer cell, compute the shortest distance to vessel boundary
# Surrogate for Oxygen supply

# Plot Cell Boundaries for illustrative purposes
vesselFiles <- list.files('./Vessel_Boundary') %>%
  data.frame() %>%
  `colnames<-` (c('filename'))

meanDV <- data.frame(matrix(nrow = 0, ncol = 0))

for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  print(mid)
  ##mid <- 'MB-3277'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid)
  
  # get the serial number to match the boundary data
  serial <- paste0("MB", strsplit(mid, '-')[[1]][2], '_', unique(sg_metabricobj$ImageNumber))
  
  # read vessel data
  sg_vessel <- read.csv(paste0('./Vessel_Boundary/', 
                               vesselFiles[vesselFiles$filename %like% serial,]))
  pos_vessel <- holePolygon(sg_vessel)
  
  
  
  #------- Difussion distances are computed for cancer cells outside vessels only
  
  # so, first filter those inside the vessel
  #---- remove vessels with very small sizes
  vesselmaskTRUE <- filter_vessel(pos_vessel, 300)
  
  pos_vesselTRUE <- pos_vessel %>%
    dplyr::filter(id %in% vesselmaskTRUE)
  
  
  pos_vesselFALSE <- pos_vessel %>%
    dplyr::filter(!(id %in% vesselmaskTRUE))
  
  
  
  try({
    
    nn_dist <- nn2(pos_vesselTRUE[, c('x','y')], query = sg_metabricobj[, c('Location_Center_X', 'Location_Center_Y')], k = nrow(pos_vesselTRUE),
                   treetype = 'kd', searchtype = 'priority') %>%
      `[[`('nn.dists') %>%
      as.data.frame()
    
    sg_metabricobj$nn_dist <- nn_dist[,1]
    

    
    # point in polygon test
    result <- apply(sg_metabricobj[, c('Location_Center_X', 'Location_Center_Y')], 1, function(x) 
      is_point_in_polygons(x, pos_vesselTRUE))
    
    # assign test result as single-cell attributes
    tm_metabricobj <- sg_metabricobj %>%
      mutate(in_vasculature = result) %>%
      dplyr::filter(in_vasculature == FALSE) %>%
      dplyr::filter(Phenotype == 'Tumor')
    
    
    #---------------------------------------#
    # Compute DD
    
    nnDistvessel <- nn2(pos_vesselTRUE[, c('x', 'y')], tm_metabricobj[, c('Location_Center_X', 'Location_Center_Y')], 
                        k = 1, treetype = 'kd', searchtype = 'priority') %>%
      .$nn.dists %>%
      median()
    
    meanDV <- rbind.data.frame(meanDV, cbind.data.frame(nnDistvessel, mid))

  })
  
}


#--------------- Survival analysis --------------#


colnames(selectedDF)[2] <- 'metabric_id'
colnames(meanDV)[2] <- 'metabric_id'
selectedDF$`Relapse Free Status` <- ifelse(selectedDF$`Relapse Free Status` == '1:Recurred', 1, 0)


test <- merge(meanDV, selectedDF, by = 'metabric_id')
test$group <- ifelse(test$nnDistvessel > median(test$nnDistvessel), 'high', 'low')
test$SurvivalGroup <- ifelse(test$Survival > 66, 'high', 'low')


high <- test[test$SurvivalGroup == 'high', 'nnDistvessel']
low <- test[test$SurvivalGroup == 'low', 'nnDistvessel']



wilcox.test(high, low)


fit <- survfit(Surv(Survival, Status) ~ group, data = test)
#fit <- survfit(Surv(`Relapse Free Status (Months)`, `Relapse Free Status`) ~ group, data = test)

p <- ggsurvplot(fit, data = test, 
                palette = c('#7688ab', '#e74d4d'),
                legend = 'none',
                legend.title = '',
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
                pval = TRUE,
                #break.time.by = 2000
) +
  xlab('Time, (Months)') +
  ylab('Overall survival')
p







