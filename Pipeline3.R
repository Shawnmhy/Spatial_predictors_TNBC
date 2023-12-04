
#-------------------------------------------#
# METABRIC TNBC processing pipeline 2-------#
# @Author: Haoyang Mi ----------------------#
# Date: May 3rd 2023------------------------#
# REF: 


library(ggplot2); library(ComplexHeatmap); library(circlize); library(ggthemes); library(ggprism)
library(tidyr); library(dplyr); library(plyr); library(readr); library(data.table)
library(survival); library(survminer); library(ade4); library(vegan)
library(RANN); library(Rtsne)
library(rgeos); library(igraph)


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
  normalized_data <- apply(clipped_data, 2, function(x) {
    (x - min(x)) / (max(x) - min(x))
  })
  #featuers
  
  return(normalized_data)
}

# 对 IMC 数据进行归一化
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



####

for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  cells <- sg_cell_USETHIS[sg_cell_USETHIS$metabric_id == mid, ]
  
  print(mid)
  print(max(cells$Location_Center_X))
}




################################################################
#--------------- Distribution related to vessel ---------------#
################################################################


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


compute_polygons_area <- function(df, area_thresh) {
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
  
  return(sum(areas[area_filtered]/1000000))
}


#------------- Vasculatre Profile -------------------------#

vasculature_density <- data.frame(matrix(nrow = 0, ncol = 0))
Tcell_profiles <- data.frame(matrix(nrow = 0, ncol = 0)) 
Bcell_profiles <- data.frame(matrix(nrow = 0, ncol = 0)) 
specific_profiles <- data.frame(matrix(nrow = 0, ncol = 0)) 
immune_profiles <- data.frame(matrix(nrow = 0, ncol = 0))
CD8Tcell_profiles <- data.frame(matrix(nrow = 0, ncol = 0))
CD4Tcell_profiles <- data.frame(matrix(nrow = 0, ncol = 0))
# parameters
num_bins <- 200 
bin_width <- (500 - 0) / num_bins  


# factorize cell type
sg_cell_USETHIS$Phenotype <- as.factor(sg_cell_USETHIS$Phenotype)

selectedDF_dup <- selectedDF
colnames(selectedDF_dup)[2] <- 'metabric_id'
selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')

for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  print(mid)
  #mid <- 'MB-3277'
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
  
  
  pos_vessel <- holePolygon(sg_vessel)
  
  #---- remove vessels with very small sizes
  vesselmaskTRUE <- filter_vessel(pos_vessel, 500)
  
  
  pos_vesselTRUE <- pos_vessel %>%
    dplyr::filter(id %in% vesselmaskTRUE)
  
  
  pos_vesselFALSE <- pos_vessel %>%
    dplyr::filter(!(id %in% vesselmaskTRUE))
  
  
  #----------- Vasculature density ----------------#
  #mid <- 'MB-0399'
  core_area <- tissue_area[tissue_area$metabric_id == mid, 'TA']
  
  vasculature_area <- compute_polygons_area(pos_vessel, 500)
  
  # vasculature density
  #vasculature_density <- rbind.data.frame(vasculature_density, cbind.data.frame(mid, vasculature_area / core_area))
  vasculature_density <- rbind.data.frame(vasculature_density, cbind.data.frame(mid, length(vesselmaskTRUE) / core_area))
  #----------- Cell infiltration profile w.r.t vasculature ----------------#
  
  try({
    
    nn_dist <- nn2(pos_vesselTRUE[, c('x','y')], query = sg_metabricobj[, c('Location_Center_X', 'Location_Center_Y')], k = nrow(pos_vesselTRUE),
                   treetype = 'kd', searchtype = 'priority') %>%
      `[[`('nn.dists') %>%
      as.data.frame()
    
    sg_metabricobj$nn_dist <- nn_dist[,1]
    
    
    
    #---------- get the filtered vessel mask --------#
    
    p <- ggplot() +
      geom_polygon(data = pos_vesselTRUE, aes(x,y, group = id), fill = '#709fce') +
      geom_polygon(data = pos_vesselFALSE, aes(x,y, group = id), fill = 'grey') +
      theme_void()
    #p  
    #ggsave(p, file=paste0("Figures/Vessel_Mask_filtered/",mid, "_filteredMask.jpg"), width = 10, height = 10, units = "in", dpi = 300, bg="black")
    
    
    
    # point in polygon test
    result <- apply(sg_metabricobj[, c('Location_Center_X', 'Location_Center_Y')], 1, function(x) 
      is_point_in_polygons(x, pos_vesselTRUE))
    
    # assign test result as single-cell attributes
    sg_metabricobj <- sg_metabricobj %>%
      mutate(in_vasculature = result) %>%
      dplyr::filter(in_vasculature == FALSE)
    
    #-------- T cell infiltration profiles ---------#
    
    # first, calculate the total cell counts at each bin
    
    
    
    sg_metabricobj$bin_indice <- cut(sg_metabricobj$nn_dist, breaks = seq(0, 500, bin_width),
                                     include.lowest = TRUE, right = FALSE) 
    
    
    bin_counts_total <- sg_metabricobj %>%
      group_by(bin_indice, .drop = FALSE) %>%
      dplyr::summarise(count = n()) %>%
      `colnames<-` (c('bin_indice', 'count_total')) 
    
    
    
    # factorize the bin
    #-------- Immune infiltration profiles ---------#
    sg_metabricobj_I <- sg_metabricobj %>%
      dplyr::filter(in_vasculature == FALSE) %>%
      dplyr::filter(Phenotype != 'Tumor' &
                      Phenotype != 'Fibroblasts' &
                      Phenotype != 'Endothelial' &
                      Phenotype != 'Myofibroblasts' &
                      Phenotype != 'Ki67+')
    
    bin_counts_immune <- sg_metabricobj_I %>%
      group_by(bin_indice, .drop = FALSE) %>%
      dplyr::summarise(count = n()) %>%
      `colnames<-` (c('bin_indice', 'count_immune')) 
    
    # compute
    bin_fractions <- left_join(bin_counts_total, bin_counts_immune, by = "bin_indice") %>%
      mutate(fraction = count_immune / count_total) %>%
      #select(bin, Phenotype, fraction) %>%
      data.frame()
    
    immune_profiles <- rbind.data.frame(immune_profiles, cbind.data.frame(bin_fractions, mid, Surv))
    
    
    #-------- T cell infiltration profiles ---------#
    
    
    sg_metabricobj_T <- sg_metabricobj %>%
      dplyr::filter(in_vasculature == FALSE) %>%
      dplyr::filter(Phenotype == 'CD4+ T cells' |
                      Phenotype == 'CD8+ T cells' |
                      Phenotype == 'CD4+ T cells & APCs' |
                      Phenotype == 'Tregs and Tex') 
    
    
    bin_counts_Tcells <- sg_metabricobj_T %>%
      group_by(bin_indice, .drop = FALSE) %>%
      dplyr::summarise(count = n()) %>%
      `colnames<-` (c('bin_indice', 'count_Tcells')) 
    
    # compute
    bin_fractions <- left_join(bin_counts_total, bin_counts_Tcells, by = "bin_indice") %>%
      mutate(fraction = count_Tcells / count_total) %>%
      #select(bin, Phenotype, fraction) %>%
      data.frame()
    
    Tcell_profiles <- rbind.data.frame(Tcell_profiles, cbind.data.frame(bin_fractions, mid, Surv))
    
    #-------- Cell type-specific infiltration profiles ---------#
    
    
    bin_counts <- sg_metabricobj %>%
      group_by(bin_indice, Phenotype, .drop = FALSE) %>%
      dplyr::summarise(count = n())
    
    
    bin_fractions <- left_join(bin_counts, bin_counts_total, by = "bin_indice") %>%
      mutate(fraction = count / count_total) %>%
      data.frame()
    
    specific_profiles <- rbind.data.frame(specific_profiles, cbind.data.frame(bin_fractions, mid, Surv))
    
    
    #-------------------- Assess exhaustion ---------------------#
    #------ PD-1 expression on B cell
    #------ ICOS on CD4+ T cells 
    #------ ICOS expression on CD8 T cell
    
    sg_metabricobj_B <- sg_metabricobj %>%
      dplyr::filter(in_vasculature == FALSE) %>%
      dplyr::filter(Phenotype == 'B cells')
    
    
    bin_counts_Bcells <- sg_metabricobj_B %>%
      group_by(bin_indice, .drop = FALSE) %>%
      dplyr::summarise(mean_ICOS = mean(`ICOS`),
                       mean_PD1 = mean(`PD-1`),
                       mean_OX40 = mean(`OX40`)) %>%
      `colnames<-` (c('bin_indice', 'mean_ICOS', 'mean_PD1', 'mean_OX40')) 

    Bcell_profiles <- rbind.data.frame(Bcell_profiles, cbind.data.frame(bin_counts_Bcells, mid, Surv))
    
    
    #--------------------CD4 T cellls ----------------------#
    sg_metabricobj_CD4T <- sg_metabricobj %>%
      dplyr::filter(in_vasculature == FALSE) %>%
      dplyr::filter(Phenotype == 'CD4+ T cells')
    
    
    bin_counts_CD4Tcells <- sg_metabricobj_CD4T %>%
      group_by(bin_indice, .drop = FALSE) %>%
      dplyr::summarise(mean_ICOS = mean(`ICOS`),
                       mean_PD1 = mean(`PD-1`),
                       mean_OX40 = mean(`OX40`)) %>%
      `colnames<-` (c('bin_indice', 'mean_ICOS', 'mean_PD1', 'mean_OX40')) 
    
    
    CD4Tcell_profiles <- rbind.data.frame(CD4Tcell_profiles, cbind.data.frame(bin_counts_CD4Tcells, mid, Surv))
    
    
    #--------------------CD4 T cellls ----------------------#
    sg_metabricobj_CD8T <- sg_metabricobj %>%
      dplyr::filter(in_vasculature == FALSE) %>%
      dplyr::filter(Phenotype == 'CD8+ T cells')
    
    
    bin_counts_CD8Tcells <- sg_metabricobj_CD8T %>%
      group_by(bin_indice, .drop = FALSE) %>%
      dplyr::summarise(mean_ICOS = mean(`ICOS`),
                       mean_PD1 = mean(`PD-1`),
                       mean_OX40 = mean(`OX40`)) %>%
      `colnames<-` (c('bin_indice', 'mean_ICOS', 'mean_PD1', 'mean_OX40')) 
    
    
    CD8Tcell_profiles <- rbind.data.frame(CD8Tcell_profiles, cbind.data.frame(bin_counts_CD8Tcells, mid, Surv))
    
  })
  
  
}







#------------- Vasculatre Profile -------------------------#

vasculature_density <- data.frame(matrix(nrow = 0, ncol = 0))
Tcell_profiles <- data.frame(matrix(nrow = 0, ncol = 0)) 
specific_profiles <- data.frame(matrix(nrow = 0, ncol = 0)) 
immune_profiles <- data.frame(matrix(nrow = 0, ncol = 0))
# parameters
num_bins <- 200 
bin_width <- (500 - 0) / num_bins  
area_thresh <- 500 #500, 700, 800

# factorize cell type
sg_cell_USETHIS$Phenotype <- as.factor(sg_cell_USETHIS$Phenotype)

selectedDF_dup <- selectedDF
colnames(selectedDF_dup)[2] <- 'metabric_id'
selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')

for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  print(mid)
  #mid <- 'MB-3277'
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
  
  
  pos_vessel <- holePolygon(sg_vessel)
  
  #---- remove vessels with very small sizes
  vesselmaskTRUE <- filter_vessel(pos_vessel, area_thresh)
  
  
  pos_vesselTRUE <- pos_vessel %>%
    dplyr::filter(id %in% vesselmaskTRUE)
  
  
  pos_vesselFALSE <- pos_vessel %>%
    dplyr::filter(!(id %in% vesselmaskTRUE))
  
  
  #----------- Vasculature density ----------------#
  #mid <- 'MB-0399'
  core_area <- tissue_area[tissue_area$mid == mid, 'TA']
  
  vasculature_area <- compute_polygons_area(pos_vessel, area_thresh)
  
  # vasculature density
  #vasculature_density <- rbind.data.frame(vasculature_density, cbind.data.frame(mid, vasculature_area / core_area))
  vasculature_density <- rbind.data.frame(vasculature_density, cbind.data.frame(mid, length(vesselmaskTRUE) / core_area))
  #----------- Cell infiltration profile w.r.t vasculature ----------------#
  
  
}


#--------------Vasculature density -----------------#
#vasculature_density <- vasculature_density[, -c(4,5)]
vasculature_density <- vasculature_density[complete.cases(vasculature_density),]


vasculature_density$type <- 'Sample'

colnames(vasculature_density) <- c('metabric_id', 'density', 'Sample')


selectedDF_dup <- selectedDF
# CD8 / Tumor cell ratio
colnames(selectedDF_dup)[2] <- 'metabric_id'


#selectedDF_dup$Stage <- ifelse(selectedDF_dup$Stage <= 2, 'Stage I-II', 'Stage III-IV')
#selectedDF_dup$`Age at Diagnosis` <- ifelse(selectedDF_dup$`Age at Diagnosis` <= 55, '<= 55 years of age', '>55 years of age')
#selectedDF_dup$Grade <- ifelse(selectedDF_dup$Grade <= 2, 'Grade I-II', 'Grade III')
#selectedDF_dup$`Relapse Free Status` <- ifelse(selectedDF_dup$`Relapse Free Status` == '1:Recurred', 'Recurred', 'Not recurred')
#selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')

#selectedDF_dup$`TMB (nonsynonymous)` <- ifelse(selectedDF_dup$`TMB (nonsynonymous)` < 5, 'TMB low', 'TMB high')
#selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'High'] <- 'Cellularity high'
#selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'Moderate'] <- 'Cellularity moderate'
#selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'Low'] <- 'Cellularity low'
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
fit <- survfit(Surv(`Relapse Free Status (Months)`, relapse_status) ~ vd_level, data = vasculature_DF)
#fit <- survfit(Surv(Survival, Status) ~ vd_level, data = vasculature_DF)
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

ggsave(print(p), file=paste0("Figures/Vasculature_Density_OS.png"), width = 7, height = 6, units = "in", dpi = 300, bg="black")


vasculature_DF$density
p <- ggplot(vasculature_DF, aes(x = Sample, y = density)) +
  geom_violin(fill = "NA", color = "black" , size = 1.5) +
  #geom_point(aes(color = `Survival`)) +
  theme_classic() +
  theme(axis.text = element_text(size = 35),
        axis.title = element_text(size = 35),
        panel.border = element_blank(),
        axis.line = element_blank(),
        axis.text.x = element_blank())+
  ylab('Vasculature density') +
  xlab('All patients') +
  geom_hline(yintercept = thresh)

p

ggsave(p, file=paste0("Figures/Vasculature_density_violin.png"), width = 4, height = 14, units = "in", dpi = 300, bg="black")

max(vasculature_DF$density)



specific_profiles[(specific_profiles == 'NaN')] <- 0
midpts <- seq(1.25, 499.5, length.out = 200)
wilcox_testDF <- data.frame(matrix(nrow = 0, ncol = 0))

#--------- Each cell type -------------#
for(ctype in unique(specific_profiles$Phenotype)){
  #ctype <- 'CD8+ T cells'
  
  result <- specific_profiles %>%
    group_by(bin_indice, Phenotype, Surv) %>%
    dplyr::summarise(
      mean = mean(fraction, na.rm = TRUE),
      n = n(),
      sd = sd(fraction, na.rm = TRUE)
    ) %>%
    mutate(sem = sd / sqrt(n)) %>%
    select(bin_indice, Phenotype, mean, sem, Surv) %>%
    dplyr::filter(Phenotype == ctype)
  
  
  result <- result[complete.cases(result),]
  
  axis_dictionary <- cbind.data.frame(unique(result$bin_indice), midpts) %>%
    `colnames<-` (c('bin_indice', 'bin'))
  
  test <- merge(axis_dictionary, result, by = 'bin_indice')
  test1 <- merge(axis_dictionary, result, by = 'bin_indice')
  
  test_high <- test[test$Surv == 'Survival high',]
  
  test_low <- test[test$Surv == 'Survival low',]
 
 
 #test_high <- merge(test[test$Surv == 'Survival high',], test1[test1$Surv == 'Survival high',], by = 'bin_indice') %>%
 #   mutate(ratio = mean.x/mean.y)
    
 #test_low <- merge(test[test$Surv == 'Survival low',], test1[test1$Surv == 'Survival low',], by = 'bin_indice') %>%
#   mutate(ratio = mean.x/mean.y)
  
  wilcox_test <- wilcox.test(test_high$mean, test_low$mean)
  #wilcox_test$p.value
  #ctype <- 'ratio'
  wilcox_testDF <- rbind.data.frame(wilcox_testDF, cbind.data.frame(wilcox_test$p.value, ctype))
  
  testx <- rbind.data.frame(test_high, test_low) %>%
    select(bin.x, ratio, Surv.x) %>%
    `colnames<-` (c('bin', 'ratio', 'Surv'))
    
  p <- ggplot(testx, aes(bin, ratio, group = Surv, color = Surv)) +
    #geom_line() +
    #geom_point() +
    theme_prism() +
    geom_smooth(aes(fill = Surv), method = "lm", formula = y ~ poly(x,20), se = TRUE, alpha = 0.3, span = 0.75, size =3) +
    geom_vline(xintercept = 0, linetype = 'dashed', size = 3) +
    scale_fill_manual(values = c('Survival high' = '#4e4cd5','Survival low' = '#e9383a')) +
    scale_color_manual(values = c('Survival high' = '#4e4cd5','Survival low' = '#e9383a')) +
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          legend.position = 'none') +
    #xlab(expression('Distance from vessel ('~ mu~'m)')) +
    #ylab('Fraction of cells') +
    xlim(0, 300)
  p
  #ctype <- 'ratio'
  ggsave(p, file=paste0("Figures/Vessel_", ctype, "_CD8TumorRatio.png"), width = 9, height = 7, units = "in", dpi = 300, bg="transparent")
  
}

wilcox_testDF$`wilcox_test$p.value` <- p.adjust(wilcox_testDF$`wilcox_test$p.value`)
#wilcox_testDF1 <- wilcox_testDF
#-------- Immune cell total ------------------#
immune_profiles[(immune_profiles == 'NaN')] <- 0
midpts <- seq(1.25, 499.5, length.out = 200)

result <- immune_profiles %>%
  group_by(bin_indice, Surv) %>%
  dplyr::summarise(
    mean = mean(fraction, na.rm = TRUE),
    n = n(),
    sd = sd(fraction, na.rm = TRUE)
  ) %>%
  mutate(sem = sd / sqrt(n)) %>%
  select(bin_indice, mean, sem, Surv)

result <- result[complete.cases(result),]

axis_dictionary <- cbind.data.frame(unique(result$bin_indice), midpts) %>%
  `colnames<-` (c('bin_indice', 'bin'))

test <- merge(axis_dictionary, result, by = 'bin_indice')


test_high <- test[test$Surv == 'Survival high',]
test_low <- test[test$Surv == 'Survival low',]

wilcox.test(test_high$mean, test_low$mean)

p <- ggplot(test, aes(bin, ratio, group = Surv, color = Surv)) +
  #geom_line() +
  #geom_point() +
  theme_classic() +
  geom_smooth(aes(fill = Surv), method = "lm", formula = y ~ poly(x,7), se = TRUE, alpha = 0.3, span = 0.75) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 1) +
  scale_fill_manual(values = c('Survival high' = '#4e4cd5','Survival low' = '#e9383a')) +
  scale_color_manual(values = c('Survival high' = '#4e4cd5','Survival low' = '#e9383a')) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 40),
        legend.position = 'none') +
  xlab(expression('Distance from vessel ('~ mu~'m)')) +
  ylab('Fraction of cells') +
  xlim(0, 300)
p
ggsave(p, file=paste0("Figures/Vessel_Immune_infiltration.png"), width = 9, height = 7, units = "in", dpi = 300, bg="transparent")




#------------------------------------------------#
# Functional molecules


#ctype <- 'CD8+ T cells'
CD4Tcell_profiles <- CD4Tcell_profiles %>%
  mutate_all(~ifelse(is.nan(.), 0, .))


result <- CD4Tcell_profiles %>%
  group_by(bin_indice, Surv) %>%
  dplyr::summarise(
    mean = mean(mean_PD1, na.rm = TRUE),
    sd = sd(mean_PD1, na.rm = TRUE)
  ) %>%
  mutate(sem = sd / sqrt(mean)) %>%
  select(bin_indice, mean, sem, Surv) %>%
  dplyr::filter(!is.na(bin_indice))

#result <- result[complete.cases(result),]

axis_dictionary <- cbind.data.frame(unique(result$bin_indice), midpts) %>%
  `colnames<-` (c('bin_indice', 'bin'))

test <- merge(axis_dictionary, result, by = 'bin_indice') 

test_high <- test[test$Surv == 'Survival high',]

test_low <- test[test$Surv == 'Survival low',]


mean(test_high$mean[1:50])
mean(test_low$mean[1:50])

wilcox_test <- wilcox.test(test_high$mean, test_low$mean, paired = TRUE)
wilcox_test$p.value
#ctype <- 'ratio'
wilcox_testDF <- rbind.data.frame(wilcox_testDF, cbind.data.frame(wilcox_test$p.value, ctype))

test <- rbind.data.frame(test_high, test_low) %>%
  select(bin, mean, Surv) %>%
  `colnames<-` (c('bin', 'mean', 'Surv'))

p <- ggplot(test, aes(bin, mean, group = Surv, color = Surv)) +
  #geom_line() +
  #geom_point() +
  theme_prism() +
  geom_smooth(aes(fill = Surv), method = "lm", formula = y ~ poly(x,20), se = TRUE, alpha = 0.3, span = 0.75, size =3) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 3) +
  scale_fill_manual(values = c('Survival high' = '#4e4cd5','Survival low' = '#e9383a')) +
  scale_color_manual(values = c('Survival high' = '#4e4cd5','Survival low' = '#e9383a')) +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        legend.position = 'none') +
  #xlab(expression('Distance from vessel ('~ mu~'m)')) +
  #ylab('Fraction of cells') +
  xlim(0, 500)
p
#ctype <- 'ratio'
ggsave(p, file=paste0("Figures/Vessel_", ctype, "_CD8TumorRatio.png"), width = 9, height = 7, units = "in", dpi = 300, bg="transparent")





#-------- T cell total ------------------#
Tcell_profiles[(Tcell_profiles == 'NaN')] <- 0
midpts <- seq(1.25, 499.5, length.out = 200)

result <- Tcell_profiles %>%
  group_by(bin_indice, Surv) %>%
  dplyr::summarise(
    mean = mean(fraction, na.rm = TRUE),
    n = n(),
    sd = sd(fraction, na.rm = TRUE)
  ) %>%
  mutate(sem = sd / sqrt(n)) %>%
  select(bin_indice, mean, sem, Surv)

result <- result[complete.cases(result),]

axis_dictionary <- cbind.data.frame(unique(result$bin_indice), midpts) %>%
  `colnames<-` (c('bin_indice', 'bin'))

test <- merge(axis_dictionary, result, by = 'bin_indice')


test_high <- test[test$Surv == 'Survival high',]
test_low <- test[test$Surv == 'Survival low',]

wilcox.test(test_high$mean, test_low$mean)

p <- ggplot(test, aes(bin, mean, group = Surv, color = Surv)) +
  theme_classic() +
  geom_smooth(aes(fill = Surv), method = "lm", formula = y ~ poly(x,7), se = TRUE, alpha = 0.3, span = 0.75) +
  geom_vline(xintercept = 0, linetype = 'dashed', size = 1) +
  scale_fill_manual(values = c('Survival high' = '#4e4cd5','Survival low' = '#e9383a')) +
  scale_color_manual(values = c('Survival high' = '#4e4cd5','Survival low' = '#e9383a')) +
  theme(axis.text = element_text(size = 40),
        axis.title = element_text(size = 40),
        legend.position = 'none') +
  xlab(expression('Distance from vessel ('~ mu~'m)')) +
  ylab('Fraction of cells') +
  xlim(0, 300)
p
ggsave(p, file=paste0("Figures/Vessel_Tcell_infiltration.png"), width = 9, height = 7, units = "in", dpi = 300, bg="transparent")


#------------- Correlations between cell types -----------------------#
specific_profiles[(specific_profiles == 'NaN')] <- 0
midpts <- seq(1.25, 499.5, length.out = 200)

result <- specific_profiles %>%
  group_by(bin_indice,Phenotype, Surv) %>%
  dplyr::summarise(
    mean = mean(fraction, na.rm = TRUE),
    n = n(),
    sd = sd(fraction, na.rm = TRUE)
  ) %>%
  mutate(sem = sd / sqrt(n)) %>%
  select(bin_indice, mean, sem, Surv)

result <- result[complete.cases(result),]

axis_dictionary <- cbind.data.frame(unique(result$bin_indice), midpts) %>%
  `colnames<-` (c('bin_indice', 'bin'))



specific_profiles_low <- result %>%
  dplyr::filter(Surv == 'Survival low') %>%
  merge(axis_dictionary, by = 'bin_indice')

specific_profiles_high <- result %>%
  dplyr::filter(Surv == 'Survival high') %>%
  merge(axis_dictionary, by = 'bin_indice')


aType <- unique(sg_cell_USETHIS$Phenotype)
Heatmap_all <- data.frame(matrix(nrow = 15, ncol = 15))

colnames(Heatmap_all) <- aType
rownames(Heatmap_all) <- aType

Heatmap_all <- replace(Heatmap_all, is.na(Heatmap_all), 0)
Heatmap_SH <- Heatmap_all
Heatmap_SL <- Heatmap_all

# switch: LumA or Her2

pval_SH <- Heatmap_all
pval_SL <- Heatmap_all

for(ctype1 in unique(sg_cell_USETHIS$Phenotype)){  
  
  ctype1 <- 'Macrophages'
  ctype1_low <-specific_profiles_low %>%
    dplyr::filter(Phenotype == ctype1) %>%
    dplyr::filter(bin < 300)
  
  ctype1_high <-specific_profiles_high %>%
    dplyr::filter(Phenotype == ctype1) %>%
    dplyr::filter(bin < 300)
  
  
  # correlation
  for(ctype2 in unique(sg_cell_USETHIS$Phenotype)){  
    
    ctype2 <- 'Tregs and Tex'
    ctype2_low <-specific_profiles_low %>%
      dplyr::filter(Phenotype == ctype2) %>%
      dplyr::filter(bin < 300)
    
    ctype2_high <-specific_profiles_high %>%
      dplyr::filter(Phenotype == ctype2) %>%
      dplyr::filter(bin < 300)
    
    
    
    cor.test(ctype2_high$mean, ctype1_high$mean)
    
    # cor and pvalue for high 
    # Low survival
    cor_test <- cor.test(ctype1_low$mean, ctype2_low$mean)
    
    
    Heatmap_SL[ctype1, ctype2] <- cor_test$estimate
    pval_SL[ctype1, ctype2] <- cor_test$p.value
    #pval_SL[ctype1, ctype2] <- symnum(cor_test$p.value, corr = FALSE, na = FALSE, 
    #                                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
    #                                  symbols = c("***", "**", "*", "", " "))
    # Low survival
    cor_test <- cor.test(ctype1_high$mean, ctype2_high$mean)
    
    
    Heatmap_SH[ctype1, ctype2] <- cor_test$estimate
    pval_SH[ctype1, ctype2] <- cor_test$p.value
    
    #pval_SH[ctype1, ctype2] <- symnum(cor_test$p.value, corr = FALSE, na = FALSE, 
    #                                  cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
    #                                  symbols = c("***", "**", "*", "", " "))
    
    
    
  }

 

   
}


Heatmap_SH_rs <- melt(setDT(Heatmap_SH, keep.rownames = TRUE), "rn") %>%
  mutate(segment = 1) %>%
  data.frame() %>%
  `colnames<-` (c('col', 'row', 'n', 'segment'))

Heatmap_SL_rs <- melt(setDT(Heatmap_SL, keep.rownames = TRUE), "rn") %>%
  mutate(segment = 2) %>%
  data.frame() %>%
  `colnames<-` (c('col', 'row', 'n', 'segment'))


#---------- multiple test correction ----------------#

p_vector <- as.vector(as.matrix(pval_SH))

p_adjusted_vector <- p.adjust(p_vector, method = "BH")
p_adjusted_vector <- symnum(p_adjusted_vector, corr = FALSE, na = FALSE, 
       cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
       symbols = c("***", "**", "*", "", " "))

SH_adjusted <- matrix(p_adjusted_vector, nrow = nrow(pval_SH), ncol = ncol(pval_SH)) %>%
  data.frame()

rownames(SH_adjusted) <- rownames(pval_SH)
colnames(SH_adjusted) <- colnames(pval_SH)







p_vector <- as.vector(as.matrix(pval_SL))

p_adjusted_vector <- p.adjust(p_vector, method = "BH")
p_adjusted_vector <- symnum(p_adjusted_vector, corr = FALSE, na = FALSE, 
                            cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1), 
                            symbols = c("***", "**", "*", "", " "))

SL_adjusted <- matrix(p_adjusted_vector, nrow = nrow(pval_SL), ncol = ncol(pval_SL)) %>%
  data.frame()

rownames(SL_adjusted) <- rownames(pval_SL)
colnames(SL_adjusted) <- colnames(pval_SL)





pval_SH_rs <- melt(setDT(SH_adjusted, keep.rownames = TRUE), "rn") %>%
  mutate(segment = 1) %>%
  data.frame() %>%
  `colnames<-` (c('col', 'row', 'n', 'segment'))

pval_SL_rs <- melt(setDT(SL_adjusted, keep.rownames = TRUE), "rn") %>%
  mutate(segment = 2) %>%
  data.frame() %>%
  `colnames<-` (c('col', 'row', 'n', 'segment'))
########

Heatmap_df <- rbind.data.frame(Heatmap_SH_rs, Heatmap_SL_rs)
pval_df <- rbind.data.frame(pval_SH_rs, pval_SL_rs)


vecnames <- c("Tumor", "Myofibroblasts", "Fibroblasts", "Macrophages", "Endothelial",
              "CD8+ T cells", "CD38+ lymphocytes", "Granulocytes", "CD57+",
              "CD4+ T cells", "CD4+ T cells & APCs", "B cells", "Ki67+", "Tregs and Tex", "Macrophages & granulocytes")
Heatmap_df$col <- factor(Heatmap_df$col,
                           levels = rev(vecnames))

Heatmap_df$row <- factor(Heatmap_df$row,
                           levels = vecnames)

pval_df$col <- factor(pval_df$col,
                      levels = rev(vecnames))
pval_df$row <- factor(pval_df$row,
                      levels = vecnames)

p <- ggplot() +
  geom_tile(data = Heatmap_df, aes(x = 1, y = segment, fill = n)) +
  scale_fill_distiller(palette = "RdBu", limits = c(-1, 1)) +
  geom_text(data = pval_df, aes(x = 1, y = segment, label = n), vjust = 1) + 
  facet_grid(row~col, switch = "both") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(panel.spacing = unit(0, "line"),
        panel.border = element_rect(fill = NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none')
p

ggsave(p, file=paste0("Figures/CorHeatmap.png"), width = 6, height = 7, units = "in", dpi = 300)
















#----------------- Vasculature heterogeneity ----------------#

library(usedist)
library(ggnetwork)
library(vegan)

pvals <- data.frame(matrix(nrow = 0, ncol = 0))
pe_clustered <- data.frame(matrix(nrow = 0, ncol = 0))
pe_random <- data.frame(matrix(nrow = 0, ncol = 0))
area_thresh <- 500

for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  print(mid)
  #mid <- 'MB-0350'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid)
  
  
  
  # get the serial number to match the boundary data
  serial <- paste0("MB", strsplit(mid, '-')[[1]][2], '_', unique(sg_metabricobj$ImageNumber))
  
  
  # read vessel data
  sg_vessel <- read.csv(paste0('./Vessel_Boundary/', 
                               vesselFiles[vesselFiles$filename %like% serial,]))
  
  
  pos_vessel <- holePolygon(sg_vessel)
  
  #---- remove vessels with very small sizes
  vesselmaskTRUE <- filter_vessel(pos_vessel, area_thresh)
  
  pos_vesselTRUE <- pos_vessel %>%
    dplyr::filter(id %in% vesselmaskTRUE)
  pos_vesselFALSE <- pos_vessel %>%
    dplyr::filter(!(id %in% vesselmaskTRUE))
  #---------- Point in polygon --------#
  # point in polygon test
  result <- apply(sg_metabricobj[, c('Location_Center_X', 'Location_Center_Y')], 1, function(x) 
    is_point_in_polygons(x, pos_vesselTRUE))
  
  sg_metabricobj <- sg_metabricobj %>%
    mutate(in_vasculature = result) %>%
    dplyr::filter(in_vasculature == FALSE)

  
  #---- First, get the centroid of vasculature area
  # Only if vesselTRUE is not emoty
  
  ids <- unique(pos_vesselTRUE$id)
  
  if(length(ids) > 2){
    
    my_polygons_sp <- lapply(ids, function(m) {
      df <- pos_vesselTRUE %>%
        dplyr::filter(id == m)
      p <- Polygon(cbind(df$x, df$y))
      Polygons(list(p), ID = as.character(m))
    }) %>% SpatialPolygons()
    
    centroids <- gCentroid(my_polygons_sp, byid = TRUE)
    
    
    centroids_df <- data.frame(id = ids,
                               x = coordinates(centroids)[, 1],
                               y = coordinates(centroids)[, 2]) 
    
    
    # physical distance
    geographic_distance <- dist(centroids_df[, 2:3])
    
    # molecular distance
    pe_mean <- sapply(ids, function(m) {
      # get df for the vessel
      #m <- 2
      df <- pos_vesselTRUE %>%
        dplyr::filter(id == m) %>%
        select(x, y)
      
      # find nearest neighbors
      
      nn_idx <- nn2(df, query = sg_metabricobj[, c('Location_Center_X', 'Location_Center_Y')], k = nrow(df),
                    treetype = 'kd', searchtype = 'radius', radius = 20) %>%
        .$nn.idx
      
      non_zero_rows <- which(apply(nn_idx, 1, function(x) all(x == 0), simplify = TRUE) == FALSE)
      
      # use the index to get the protein expression profile
      
      #----------------- For molecular phenotyping ---------------#
      pe <- sg_metabricobj[non_zero_rows,] %>%
        select(12:50) %>%
        dplyr::summarise(across(everything(), mean, na.rm = TRUE))
      
      #----------------- For cellular phenotyping ---------------#
      #pe <- sg_metabricobj[non_zero_rows,] %>%
      #  group_by(Phenotype, .drop = FALSE) %>%
      #  tally() %>%
      #  spread(key = Phenotype, value = n, fill = 0)
      
        
    }) %>%
      t() %>%
      `rownames<-` (c(ids))
    
    
    
    # first get the molecular profile 
    molecular_distance <- dist(pe_mean)
    molecular_distance <- dist_setNames(molecular_distance, ids)
    
    mantel_test <- mantel(geographic_distance, molecular_distance, permutations = 999)
    
    # plot nearest molecular distance neighbor
    edges <- melt(as.matrix(molecular_distance), varnames = c("row", "col")) %>%
      dplyr::filter(row > col) %>%
      group_by(row) %>%
      slice(which.min(value)) %>%
      data.frame()
    
    edges$x <- centroids_df[match(edges$row, centroids_df$id), 'x']
    edges$y <- centroids_df[match(edges$row, centroids_df$id), 'y']
    
    edges$xend <- centroids_df[match(edges$col, centroids_df$id), 'x']
    edges$yend <- centroids_df[match(edges$col, centroids_df$id), 'y']
    
    # the above function returns the edge, now construct networks
    #g <- graph_from_data_frame(edges[, c("row", "col")], directed = TRUE)
    
    # get the coordinates
    #colnames(centroids_df) <- c("name", "x", "y")
    
    #V(g)$x <- centroids_df$x[match(V(g)$name, centroids_df$name)]
    #V(g)$y <- centroids_df$y[match(V(g)$name, centroids_df$name)]
    
    colorcode <- ifelse(mantel_test$signif < 0.05, '#7688ab', '#e74d4d')
    label <- ifelse(mantel_test$signif < 0.05, 'signif', 'nonsignif')
    
    # Differences between clustered vs random cores in 
    # terms of proteomic/molecular domain
    
    
    if(mantel_test$signif < 0.05){
      # these are the proteomic profiles stacked for all valid vessels
      pe_clustered <- rbind.data.frame(pe_clustered, cbind.data.frame(pe_mean, mid))
    }
    if(mantel_test$signif >= 0.05){
      # these are the proteomic profiles stacked for all valid vessels
      pe_random <- rbind.data.frame(pe_random, cbind.data.frame(pe_mean, mid))
    }
    
    #colorcode <- 'red'
    
    
    #------------------- Plot ----------------------#
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
      #geom_polygon(data = pos_vesselTRUE, aes(x,y, group = id), fill = '#709fce') +
      #geom_polygon(data = pos_vesselFALSE, aes(x,y, group = id), fill = '#709fce') +
      geom_polygon(data = sg_boundary_w_ctype, aes(X, Y, group = CellID), fill = 'grey', alpha = 0.3) +
      geom_path(data = pos_vesselTRUE, aes(x,y, group = draw), color = 'black', alpha = 0.6) +
      geom_nodes(data = centroids_df, aes(x = x, y = y), color = colorcode, size = 5) +
      
      geom_edges(data = edges,arrow = arrow(length = unit(8, "pt"), type = "closed"), size = 1, 
                 aes(x = x, y = y, xend = xend, yend = yend, linetype = 'dashed'), color = colorcode) +
      theme_void() +
      theme(legend.position = 'none') 
    p  
    ggsave(p, file=paste0("Figures/Mantel_test/", mid, "_", mantel_test$signif, '_', label,"_", mantel_test$statistic, ".png"), width = 10, height = 10, units = "in", dpi = 300, bg="transparent")
    
    
    #------------ Correlation plot ---------------#
    
    
    molecular_vector <- melt(as.matrix(molecular_distance), varnames = c("row", "col")) %>%
      dplyr::filter(row > col)
    
    geographic_vector <- melt(as.matrix(geographic_distance), varnames = c("row", "col")) %>%
      dplyr::filter(row > col)
    
    Mantel_cor <- cbind.data.frame(geographic_vector$value, molecular_vector$value) %>%
      `colnames<-` (c('x', 'y'))
    
    #cor.test(Mantel_cor$x, Mantel_cor$y)

    
    p <- ggplot(Mantel_cor, aes(x, y)) + 
      geom_point(fill = colorcode, size =6, shape = 21, color = 'black', stroke = 2) +
      geom_smooth(method='lm', se = FALSE, color = 'black', size = 2) +
      theme_bw() +
      theme(axis.text = element_text(size = 26),
            axis.title = element_text(size = 26)) +
      labs(x = 'Geographic distance, \U03BCm', y = 'Proteomic distance') 
    p
    #ggsave(p, file=paste0("Figures/Random_correlation.png"), width = 6.5, height = 6, units = "in", dpi = 300)
    
    
    
    pvals <- rbind.data.frame(pvals, cbind.data.frame(mid, mantel_test$signif, length(ids)))
  }
  
  
}

mids_signif <- pvals[pvals$`mantel_test$signif` < 0.05, 'mid']



selectedDF_dup <- selectedDF
# CD8 / Tumor cell ratio

selectedDF_dup$GD <- ifelse(selectedDF_dup$`Patient ID` %in% mids_signif, 'Clustered', 'Random')
selectedDF_dup$GD <- ifelse(selectedDF_dup$`Patient ID` %in% pvals$mid, selectedDF_dup$GD, 'NA')

selectedDF_dup$relapse_status <- ifelse(selectedDF_dup$`Relapse Free Status` == '1:Recurred', 1, 0)


#selectedDF_dup <- selectedDF_dup %>%
#  dplyr::filter(GD != 'Clustered')



#t.test(test1, test2)
fit <- survfit(Surv(Survival, Status) ~ GD, data = selectedDF_dup)
#fit <- survfit(Surv(`Relapse Free Status (Months)`, relapse_status) ~ GD, data = selectedDF_dup)

p <- ggsurvplot(fit, data = selectedDF_dup, 
                palette = c('#e74d4d', 'grey', '#7688ab'),
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
                #pval = TRUE,
                #break.time.by = 2000
) +
  xlab('Time, (Months)') +
  ylab('Overall survival')
p
ggsave(print(p), file=paste0("Figures/Clustered_OS.png"), width = 7, height = 6, units = "in", dpi = 300, bg="black")



#-------------- Proportion of long versus short-term survivor ------------------#
selectedDF_dup <- selectedDF
mids_signif <- pvals[pvals$`mantel_test$signif` < 0.05, 'mid']

vasculature_DF$Survival <- ifelse(vasculature_DF$Survival <= 66.3, 'Survival low', 'Survival high')





signif_patients_data <- vasculature_DF[vasculature_DF$metabric_id %in% mids_signif,]
unsignif_patients_data <- vasculature_DF[!(vasculature_DF$metabric_id %in% mids_signif),]


perc_signif_high <- nrow(signif_patients_data[signif_patients_data$vd_level == 'high' , ]) / nrow(signif_patients_data)
perc_signif_low <- nrow(signif_patients_data[signif_patients_data$vd_level == 'low', ]) / nrow(signif_patients_data)

perc_unsignif_high <- nrow(unsignif_patients_data[unsignif_patients_data$vd_level == 'high', ]) / nrow(unsignif_patients_data)
perc_unsignif_low <- nrow(unsignif_patients_data[unsignif_patients_data$vd_level == 'low', ]) / nrow(unsignif_patients_data)


M <- as.table(rbind(c(11, 7), c(6, 34)))
colnames(M) <- c('Density high', 'Density low')
rownames(M) <- c('Clustered', 'Random')

chisq.test(M)

#-----------------------------------------------------------------#
#----- This section is to plot three vessel related figures ------#
#-----------------------------------------------------------------#




for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  #mid <- 'MB-3277'
  # get the single cell data for the current metabric object
  sg_metabricobj <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid)
  
  
  
  # get the serial number to match the boundary data
  serial <- paste0("MB", strsplit(mid, '-')[[1]][2], '_', unique(sg_metabricobj$ImageNumber))
  
  
  # read vessel data
  sg_vessel <- read.csv(paste0('./Vessel_Boundary/', 
                               vesselFiles[vesselFiles$filename %like% serial,]))
  
  
  # read single cell boundary data
  sg_boundary <- read.csv(paste0('./SingleCell_Boundary/', 
                                 boundaryFiles[boundaryFiles$filename %like% serial,]))
  
  pos_vessel <- holePolygon(sg_vessel)
  
  #---------- REAL ID ----------#
  
  
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
  
  
  # assign test result as single-cell attributes
  #sg_metabricobj <- sg_metabricobj %>%
  #  mutate(in_vasculature = result)
  
  #---------- get the filtered vessel mask --------#
  vesselmaskTRUE <- filter_vessel(pos_vessel, area_thresh = 300)
  
  
  
  
  #check_orientation <- function(x, y) {
  #  n <- length(x)
   # area <- sum(x[-n] * y[-1] - x[-1] * y[-n])
  #  area <- area + x[n] * y[1] - x[1] * y[n]
  #  area <- 0.5 * area
    
  #  if(area > 0) {
  #    return("Counter-Clockwise")
  #  } else {
  #    return("Clockwise")
  #  }
  #}
  
  # Test
  test <- sg_vessel[sg_vessel$ring == 2, ]
  x <- test$x
  y <- test$y
  
  
  ggplot()+
    geom_polygon(data =test, aes(x, y))
  print(check_orientation(x, y))
  
  
  
  pos_vesselTRUE <- pos_vessel %>%
    dplyr::filter(id %in% vesselmaskTRUE)
  
  
  pos_vesselFALSE <- pos_vessel %>%
    dplyr::filter(!(id %in% vesselmaskTRUE))
  
  p <- ggplot() +
    geom_polygon(data = pos_vesselTRUE, aes(x,y, group = id), fill = '#709fce') +
    geom_polygon(data = pos_vesselFALSE, aes(x,y, group = id), fill = 'grey') +
    theme_void()
  p  
  
  ggsave(p, file=paste0("Figures/Vessel_Mask_filtered/",mid, "_filteredMask.jpg"), width = 10, height = 10, units = "in", dpi = 300, bg="black")
  
  
  
  #---------- Point in polygon --------#
  # point in polygon test
  result <- apply(sg_metabricobj[, c('Location_Center_X', 'Location_Center_Y')], 1, function(x) 
    is_point_in_polygons(x, pos_vesselTRUE))
  
  sg_metabricobj <- sg_metabricobj %>%
    mutate(in_vasculature = result)
  
  
  in_out_vessel <- merge(sg_metabricobj, sg_boundary_w_ctype, by = 'CellID')
  p <- ggplot() +
    #geom_line() +
    #geom_point() +
    theme_void() +
    theme(legend.position = 'none') +
    geom_polygon(data = in_out_vessel, aes(X, Y, fill = in_vasculature, group = CellID)) +
    #geom_polygon(data = pos_vesselTRUE, aes(x,y, group = id), fill = '#709fce', alpha = 0.4) +
    #geom_polygon(data = pos_vesselFALSE, aes(x,y, group = id), fill = '#d3d3d3', alpha = 0.8) +
    
    scale_fill_manual(values = c('FALSE' = '#d3d3d3','TRUE' = '#709fce')) 
  p
  
  ggsave(p, file=paste0("Figures/Vessel_Point_in_Vessel/",mid, ".jpg"), width = 10, height = 10, units = "in", dpi = 300, bg="black")
  
  
  
  
  #----------- Vasculature density ----------------#
  
  #----------- Cell infiltration profile w.r.t vasculature ----------------#
  sg_metabricobj <- sg_metabricobj %>%
    mutate(in_vasculature = result)# %>%
  # dplyr::filter(in_vasculature == FALSE)
  
  try({
    nn_dist <- nn2(pos_vesselTRUE[, c('x','y')], query = sg_metabricobj[, c('Location_Center_X', 'Location_Center_Y')], k = nrow(pos_vesselTRUE),
                   treetype = 'kd', searchtype = 'priority') %>%
      `[[`('nn.dists') %>%
      as.data.frame()
    
    sg_metabricobj$nn_dist <- nn_dist[,1]
    
    # merge the distance to boundary
    
    cell_poly_with_distance <- merge(sg_metabricobj, sg_boundary_w_ctype, by = 'CellID')
    
    
    cell_poly_with_distance_filter <- cell_poly_with_distance %>%
      dplyr::filter(in_vasculature == FALSE)
    
    custom_gradient <- colorRampPalette(c("#4991b9", "#fdfebb", "#d7514d"))
    
    gg3 <- ggplot(cell_poly_with_distance_filter, aes(x=X, y=Y)) + 
      theme_void() +
      geom_polygon(aes(group=CellID, fill = -nn_dist)) +
      theme(legend.position = 'none') +
      #geom_point() +
      #geom_path(data = pos_vessel, aes(x=x, y=y, group=draw), color = 'black', size=1) +  #draw is the original ring
      #scale_fill_distiller(palette = 'Spectral', limits = c(-50, 0))
      scale_fill_gradientn(colours = custom_gradient(3),
                           values = scales::rescale(c(-160, -80, 0)),
                           breaks = c(-160, -80, 0),
                           labels = c("-160", "-80", "0"))# +
    #geom_polygon(data =pos_vesselTRUE, aes(x, y, group=id), fill = 'grey') 
    #geom_path(data = pos_vessel, aes(x=x, y=y, group=draw), color = 'black', size=1) +  #draw is the original ring
    gg3
    ggsave(gg3, file=paste0("Figures/Vessel_Mask_Distance/", mid, "_Vessel.png"), width = 10, height = 10, units = "in", dpi = 300, bg="black")
    
  })
  
}


#-----------------------------------------------------------------#
# Using Hierarchical Clustering to compare clustered versus random vasculature profiles
pe_Clustered <- pe_clustered %>%
  mutate(type = 'Clustered')

pe_Random <- pe_random %>%
  mutate(type = 'Random')

vasDF <- rbind.data.frame(pe_Clustered, pe_Random)

# Separate data and type columns
data_to_cluster <- vasDF[, !names(vasDF) %in% "type"]
type_column <- vasDF$type

# Define the column annotation
col_anno <- anno_block(gp = gpar(fill = c("red", "blue", "green")), 
                       labels = unique(type_column), 
                       which = "column", 
                       labels_rot = 90)

color_mapping <- colorRamp2(c(0, 0.5, 1), c("white", "#fb3640", "#0047ab"))


test <- numeric_matrix <- as.matrix(apply(vasDF[,1:39], 2, as.numeric))


p <- Heatmap(test, 
             col = color_mapping,
             #cluster_rows = FALSE,
             #cluster_columns = FALSE,
             column_dend_height = unit(1, "cm"),
             show_heatmap_legend = FALSE,
            # top_annotation = col_anno,
             row_split = type_column,
             heatmap_legend_param = list(legend_direction = "horizontal", 
                                         legend_width = unit(4, "cm"), 
                                         at = seq(-1, 1, by = 1),
                                         labels_gp = gpar(fontsize = 15)),
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 16)
)
p




