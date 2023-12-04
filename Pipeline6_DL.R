
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

cell_data <- read.csv('../TNBC_ClinicalTrial_Validation_Cohort/data/cells.csv') 


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
  
  return(list(area_filtered, areas))
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
    BiopsyPhase == 'Baseline',
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


#---------------------------------------#
# Create a ImageID - PatientID dictionary
ImageID_PatientID_Dict <- unique(val_cell_USETHIS[, c('ImageID', 'PatientID')])
saveRDS(ImageID_PatientID_Dict, 'ImageID_PatientID_dict.rds')

###!!!!
#----------------------Data Preparation Finish ---------------------------#
#---------------------------Feature Set 1 --------------------------------#


#---------------- Section 1: read vessel masks and do the processing -----#
# Plot Cell Boundaries for illustrative purposes
vesselFiles <- list.files('../TNBC_ClinicalTrial_Validation_Cohort/data/Vessel_Boundary') %>%
  data.frame() %>%
  `colnames<-` (c('filename'))






vasculature_features <- data.frame(matrix(nrow = 0, ncol = 0))
# parameters
area_thresh <- 500 #500, 700, 800
core_area <- pi * (0.5)^2
# factorize cell type
val_cell_USETHIS$Phenotype <- as.factor(val_cell_USETHIS$Phenotype)

for(pid in unique(val_cell_USETHIS$ImageID)){
  
  print(pid)
  try({
    #pid <- 'NTImg0011'
    # get the single cell data for the current patient object
    sg_obj <- val_cell_USETHIS %>%
      dplyr::filter(ImageID == pid)
    
    
    
    # get the serial number to match the boundary data
    imagnum <- unique(sg_obj$ImageID)
    
    # read vessel data
    sg_vessel <- read.csv(paste0('../TNBC_ClinicalTrial_Validation_Cohort/data/Vessel_Boundary/', 
                                 vesselFiles[vesselFiles$filename %like% imagnum,]))
    
    
    # read single cell boundary data
    #sg_boundary <- read.csv(paste0('./SingleCell_Boundary/', 
    #                              boundaryFiles[boundaryFiles$filename %like% serial,]))
    
    
    pos_vessel <- holePolygon(sg_vessel)
    
    
    #---------------------------------------------#
    #------- Plot Vasculature for Validation -----#
    #---------------------------------------------#
    
    
    #---- remove vessels with very small sizes
    vesselmaskTRUE <- filter_vessel(pos_vessel, area_thresh)
    
    
    
    pos_vesselTRUE <- pos_vessel %>%
      dplyr::filter(id %in% vesselmaskTRUE[[1]])
    
    
    pos_vesselFALSE <- pos_vessel %>%
      dplyr::filter(!(id %in% vesselmaskTRUE[[1]]))
    

    
    #----------- STAMP ----------------#
    
    
    #---------- Point in polygon --------#
    # point in polygon test
    result <- apply(sg_obj[, c('Location_Center_X', 'Location_Center_Y')], 1, function(x) 
      is_point_in_polygons(x, pos_vesselTRUE))
    
    sg_obj <- sg_obj %>%
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
      
      
      
      
      # molecular distance
      pe_mean <- sapply(ids, function(m) {
        # get df for the vessel
        #m <- 65
        df <- pos_vesselTRUE %>%
          dplyr::filter(id == m) %>%
          select(x, y)
        
        # find nearest neighbors
        
        nn_idx <- nn2(df, query = sg_obj[, c('Location_Center_X', 'Location_Center_Y')], k = nrow(df),
                      treetype = 'kd', searchtype = 'radius', radius = 20) %>%
          .$nn.idx
        
        non_zero_rows <- which(apply(nn_idx, 1, function(x) all(x == 0), simplify = TRUE) == FALSE)
        
        # use the index to get the protein expression profile
        
        #----------------- For molecular phenotyping ---------------#
        pe <- sg_obj[non_zero_rows,] %>%
          select(15:60) %>%
          dplyr::summarise(across(everything(), mean, na.rm = TRUE)) %>%
          mutate(across(everything(), ~ifelse(is.na(.), 0, .)))
        
        
        
      }) %>%
        t() %>%
        `rownames<-` (c(ids))
      
      # identify which row contains all NaN
      na_rows <- which(apply(pe_mean, 1, function(row) all(is.na(row))))
      # remove that ID as well
      
      # physical distance
      geographic_distance <- dist(centroids_df[, 2:3])
      
      if(length(na_rows) != 0){
        ids <- ids[-na_rows]     
        geographic_distance <- dist(centroids_df[-na_rows, 2:3])
      }
      
      
      pe_mean <- pe_mean[apply(pe_mean, 1, function(row) any(!is.na(row))), ]
      
      # first get the molecular profile 
      molecular_distance <- dist(pe_mean)
      molecular_distance <- dist_setNames(molecular_distance, ids)
      
      
      mantel_test <- mantel(geographic_distance, molecular_distance, permutations = 999)
      
      STAMP <- ifelse(mantel_test$signif < 0.05, 1, 2)
    }
    if(length(ids) <= 2){
      STAMP <- 0
    }
    
    
    
    
    
    #-----------------------------------------------#
    #--------------- Area metrics ------------------#
    #-----------------------------------------------#
    #STAMP:
    #      0 - Depleted
    #      1 - Clustered
    #      2 - Random
    
    vasculature_features <- rbind.data.frame(vasculature_features, cbind.data.frame(imagnum, length(vesselmaskTRUE[[1]]) / core_area, 
                                                                                  max(vesselmaskTRUE[[2]]),
                                                                                  min(vesselmaskTRUE[[2]]),
                                                                                  mean(vesselmaskTRUE[[2]]),
                                                                                  sd(vesselmaskTRUE[[2]]),
                                                                                  STAMP)
    )
    
  })
}


#colnames(vasculature_features)[1] <- 'ImageID'
# aggregate to patient level
pt_vasculature_features_mean <- merge(vasculature_features, ImageID_PatientID_Dict, by = 'ImageID') %>%
  group_by(PatientID) %>%
  select(-ImageID) %>%
  dplyr::summarise(across(everything(),mean),
                   .groups = 'drop') %>%
  as.data.frame()

colnames(pt_vasculature_features_mean) <- c('PatientID', 'Vasculature_density_mean', 'max_Vasculature area_mean', 'min_Vasculature area_mean',
                                           'mean_Vasculature area_mean', 'sd_Vasculature area_mean', 'STMAP_mean')

pt_vasculature_features_max <- merge(vasculature_features, ImageID_PatientID_Dict, by = 'ImageID') %>%
  group_by(PatientID) %>%
  select(-ImageID) %>%
  dplyr::summarise(across(everything(),max),
                   .groups = 'drop') %>%
  as.data.frame()


colnames(pt_vasculature_features_max) <- c('PatientID', 'Vasculature_density_min', 'max_Vasculature area_max', 'min_Vasculature area_max',
                                           'mean_Vasculature area_max', 'sd_Vasculature area_max', 'STAMP_max')
pt_vasculature_features_min <- merge(vasculature_features, ImageID_PatientID_Dict, by = 'ImageID') %>%
  group_by(PatientID) %>%
  select(-ImageID) %>%
  dplyr::summarise(across(everything(),min),
                   .groups = 'drop') %>%
  as.data.frame()

colnames(pt_vasculature_features_min) <- c('PatientID', 'Vasculature_density_min', 'max_Vasculature area_min', 'min_Vasculature area_min',
                                           'mean_Vasculature area_min', 'sd_Vasculature area_min', 'STAMP_min')


pt_vasculature_features <- cbind.data.frame(
  pt_vasculature_features_mean, pt_vasculature_features_max[, -1], pt_vasculature_features_min[,-1]
)


saveRDS(pt_vasculature_features, 'vasculature_features.rds')






###!!!!
#----------------------Data Preparation Finish ---------------------------#
#---------------------------Feature Set 2 --------------------------------#

library(foreach)
library(doParallel)
library(iterators)
num_cores <- detectCores() # Identify number of cores on your machine

registerDoParallel(cores = 20)



#nearest_neighbors_all <- data.frame(matrix(nrow = 0, ncol = 0))
#sg_data <- sg_metabricobj

neighbor_freq <- function(i, sg_data){
  #i <- 1
  #sg_data <- sg_obj
  nn_idx <- RANN::nn2(sg_data[, c('Location_Center_X', 'Location_Center_Y')], 
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



imagid_vector <- unique(val_cell_USETHIS$ImageID)
#iid <- 2

result <- foreach(iid = 1:length(imagid_vector), .combine = rbind, .packages = c('tidyr', 'dplyr')) %dopar% {
  
  #print(iid)
  pid <- imagid_vector[iid]
  ##mid <- 'MB-3277'
  # get the single cell data for the current metabric object
  sg_obj <- val_cell_USETHIS %>%
    dplyr::filter(ImageID == pid)
  
  
  
  # get the serial number to match the boundary data
  imagnum <- unique(sg_obj$ImageID)
  
  
  
  
  nearest_neighbors <- sapply(seq_len(nrow(sg_obj)),  neighbor_freq, sg_data = sg_obj) %>%
    t() %>%
    data.frame()
  
  nearest_neighbors <- sapply(nearest_neighbors, as.numeric)
  #nearest_neighbors <- nearest_neighbors/10 
  
  nearest_neighbors <- as.matrix(nearest_neighbors) %>%
    data.frame()
  
  # add cell id column
  nearest_neighbors$cellid <- sg_obj$ObjectNumber
  nearest_neighbors$imageid <- pid
  
  # the final object for rbind
  nearest_neighbors
  #nearest_neighbors_all <- rbind.data.frame(nearest_neighbors_all, cbind.data.frame(iid, sg_obj$ObjectNumber, nearest_neighbors))
  
}




#saveRDS(result, 'nearest_neighbors_all_validation.rds')




#-------------------------------------------------------------------------#
#----------------- Run this section to do unsupervised clustering --------#
#-------------------------------------------------------------------------#
library(ClusterR)
set.seed(2)
km <- MiniBatchKmeans(as.matrix(result[,1:21]), clusters = 10, batch_size = 100)
clusters <- predict_MBatchKMeans(as.matrix(result[,1:21]), km$centroids)


result_clus <- result

# assign clusters
result_clus$cluster <- clusters


result_clus <- result_clus %>%
  group_by(cluster) %>%
  select(1:21) %>%
  summarise_all(mean, na.rm = TRUE)

result_clus <- sapply(result_clus, as.numeric)


# scale clusters
result_clus[,2:22] <- scale(result_clus[,2:22])
result_clus[result_clus < -1] <- -1; result_clus[result_clus >  1] <- 1

# rename columns
colnames(result_clus)
colnames(result_clus) <- c("cluster", "CA9+", "CD20+ B cells", "CD4+PD1+ T cell", "CD4+TCF1+ T cell",
                                      "CD56+ NK", "CD79a+ Plasma", "CD8+ GZMB+ T cell", "CD8+PD1+ Tex", "CD8+ T", 'CD8+ TCF1+ T  cell', 'DCs',
                                      "Endothelial", "Firbroblasts", "M2 Macrophage", "Myofibroblasts", "Neutrophils",
                           "PDL1+ APCs", "PDL1+ IDO+ APCs", "PDPN+ Stromal", "Treg", "Tumor")

rownames(result_clus) <- c('1', '2', '3', '4', '5', '6', '7', '8', '9', '10')

library(ComplexHeatmap)
library(colorRamp2)

color_mapping <- colorRamp2(c(-1, 0, 1), c("#0047ab", "white", "#fb3640"))

png(file="Figures/Neighborhood_cluster_validation.png", width = 8.5, height = 9, units = "in", res = 300)
p <- Heatmap(result_clus[,2:22], col = color_mapping,
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
#-------------- Prepare CN related featuers -----------------------#

set.seed(2)
km <- MiniBatchKmeans(as.matrix(result[,1:21]), clusters = 10, batch_size = 100)
clusters <- predict_MBatchKMeans(as.matrix(result[,1:21]), km$centroids)


result_clus <- result

# assign clusters
result_clus$cluster <- clusters
result_clus$cluster <- as.factor(result_clus$cluster)


per_imageid_cluster_prop <- result_clus %>%
  group_by(imageid, cluster, .drop = FALSE) %>%
  tally() %>%
  pivot_wider(names_from = cluster, values_from = n) %>%
  rowwise() %>%
  mutate(across(1:10, ~ .x / sum(c_across(1:10)))) %>%
  as.data.frame()

colnames(per_imageid_cluster_prop) <- c('ImageID', 'percent_CN1', 'percent_CN2', 
                                        'percent_CN3', 'percent_CN4', 'percent_CN5',
                                        'percent_CN6', 'percent_CN7', 'percent_CN8',
                                        'percent_CN9', 'percent_CN10')

  
# aggregate to patient level
pt_cn_features_mean <- merge(per_imageid_cluster_prop, ImageID_PatientID_Dict, by = 'ImageID') %>%
  select(-ImageID) %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),mean, .names = "{col}_mean"),
            .groups = 'drop') %>%
  as.data.frame()


pt_cn_features_max <- merge(per_imageid_cluster_prop, ImageID_PatientID_Dict, by = 'ImageID') %>%
  select(-ImageID) %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),max, .names = "{col}_max"),
                   .groups = 'drop') %>%
  as.data.frame()

pt_cn_features_min <- merge(per_imageid_cluster_prop, ImageID_PatientID_Dict, by = 'ImageID') %>%
  select(-ImageID) %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),min, .names = "{col}_min"),
                   .groups = 'drop') %>%
  as.data.frame()
#colnames(pt_cn_features) <- NA
pt_cn_features <- cbind.data.frame(
  pt_cn_features_mean, pt_cn_features_max[, -1], pt_cn_features_min[,-1]
)



saveRDS(pt_cn_features, 'CN_features.rds')





#---------------------------------------------------------------------#
#-------------- Prepare STAIN related featuers -----------------------#
#---------------------------------------------------------------------#

ratios_all <- data.frame(matrix(nrow = 0, ncol = 0))


imageid_vector <- unique(val_cell_USETHIS$ImageID)
#iid <- 2

result <- foreach(iid = 1:length(imageid_vector), .combine = rbind, .packages = c('tidyr', 'dplyr', 'RANN')) %dopar% {
  
  #print(iid)
  #iid <- 849
  #mid <- 'MB-0316'
  # get the single cell data for the current metabric object
  pid <- imagid_vector[iid]
  
  
  sg_obj <- val_cell_USETHIS %>%
    dplyr::filter(ImageID == pid)
  
  
  Biopsy <- unique(sg_obj$BiopsyPhase)
  PatientID <- unique(sg_obj$PatientID)
  

  # get CD8 T cells
  Tumor <- val_cell_USETHIS %>%
    dplyr::filter(ImageID == pid) %>%
    dplyr::filter(Phenotype == 'Tumor') %>%
    select(Location_Center_X, Location_Center_Y)
  
  
  # Endothelial cells
  Treg <- val_cell_USETHIS %>%
    dplyr::filter(ImageID == pid) %>%
    dplyr::filter(Phenotype == "Treg" | Phenotype == "CD8^+PD1^+T_{Ex}") %>%
    select(Location_Center_X, Location_Center_Y)
  
  # Endothelial cells
  Myo <- val_cell_USETHIS %>%
    dplyr::filter(ImageID == pid) %>%
    dplyr::filter(Phenotype == "Myofibroblasts") %>%
    select(Location_Center_X, Location_Center_Y)
  
  
  if(nrow(Tumor) > 5 & nrow(Treg) > 5 & nrow(Myo) > 5){
    Tumor_Treg <- nn2(Treg, Tumor, k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.dists %>%
      as.vector() 
    
    Tumor_Myo <- nn2(Myo, Tumor, k = 1, searchtype = 'priority', treetype = 'kd') %>%
      .$nn.dists %>%
      as.vector()
    
    ratio <- Tumor_Treg / (Tumor_Treg + Tumor_Myo)  
    #ratio <- sqrt(CD8T_Tumor * CD8T_FoxP3)
    ratios_all <- cbind.data.frame(PatientID, max(ratio), mean(ratio), min(ratio), sd(ratio))
  }
  
  ratios_all
}

pt_stain_features_mean <- result %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(), mean, .names = "{col}_mean"), .groups = 'drop') %>%
  as.data.frame()
  
pt_stain_features_max <- result %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(), max, .names = "{col}_max"), .groups = 'drop') %>%
  as.data.frame()

pt_stain_features_mean <- result %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),min, .names = "{col}_min"), .groups = 'drop') %>%
  as.data.frame()


#colnames(pt_cn_features) <- NA
pt_stain_features <- cbind.data.frame(
  pt_stain_features_mean, pt_stain_features_max[, -1], pt_stain_features_mean[,-1]
)

saveRDS(pt_stain_features, 'STAIN_features.rds')




#------------------------------------------------------------------------#
#---------------------Prepare Graph data for deep learning --------------#
#-------------------          Cell Graph   ------------------------------#
library(jsonlite)
library(tripack)
euclidean_distance <- function(node1, node2, Dual_NodeList) {
  
  #node1 <- 1
  
  x1 <- Dual_NodeList$x[Dual_NodeList$nodes == node1]
  y1 <- Dual_NodeList$y[Dual_NodeList$nodes == node1]
  x2 <- Dual_NodeList$x[Dual_NodeList$nodes == node2]
  y2 <- Dual_NodeList$y[Dual_NodeList$nodes == node2]
  return(sqrt((x1 - x2)^2 + (y1 - y2)^2))
}

colnames(result_clus)[22] <- 'ObjectNumber'
colnames(result_clus)[23] <- 'ImageID'



for(pid in unique(val_cell_USETHIS$ImageID)){
  
  print(pid)
  #dir.create(paste0('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/networks/', pid))
  
  #pid <- 'NTImg0039'
  # get the single cell data for the current metabric object

  # corePts cluster
  corePts_cn <-result_clus %>%
    dplyr::filter(ImageID == pid)
  
  corePts <- val_cell_USETHIS %>%
    dplyr::filter(ImageID == pid) %>%
    merge(corePts_cn, by = c('ImageID', 'ObjectNumber'))
  
  
  
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
  
  Dual_NodeList <- cbind(seq(1, nrow(corePts)),  corePts) %>%
    select(c(1:3, 8, 14:15, 63:84))
  colnames(Dual_NodeList)[1] <- 'nodes'
  colnames(Dual_NodeList)[5] <- 'x'
  colnames(Dual_NodeList)[6] <- 'y'
  
  
  
  Dual_EdgeList <- data.frame(cbind(i, j))
  colnames(Dual_EdgeList) <- c('col1', 'col2')
  
  Dual_EdgeList <- Dual_EdgeList %>%
    mutate(from = pmin(col1, col2),
           to = pmax(col1, col2)) %>%
    distinct(from, to)
  
  
  #Dual_EdgeList$distance <- mapply(euclidean_distance, Dual_EdgeList$from, Dual_EdgeList$to, MoreArgs = list(Dual_NodeList = Dual_NodeList))
  
  # remove isolated nodes
  Dual_NodeList <- Dual_NodeList %>%
    dplyr::filter(nodes %in% Dual_EdgeList$from | nodes %in% Dual_EdgeList$to)
  
  
  
  #-------------------------------------------------------------------------#
  #------ This section make edge and node list id consecutive --------------#
  
  Dual_NodeList$nodes_new <- seq(nrow(Dual_NodeList))
  
  
  Dual_EdgeList$from_new <- Dual_NodeList[match(Dual_EdgeList$from, Dual_NodeList$nodes), 'nodes_new']
  Dual_EdgeList$to_new <- Dual_NodeList[match(Dual_EdgeList$to, Dual_NodeList$nodes), 'nodes_new']
  
  

  #-------------------------------------------------------------------------#
  
  p <- ggplot() +
    theme_bw() +
    theme_void() +
    geom_segment(aes(x = Dual_NodeList[match(Dual_EdgeList$node_1, Dual_NodeList$nodes), 5],
                     xend = Dual_NodeList[match(Dual_EdgeList$node_2, Dual_NodeList$nodes), 5],
                     y = Dual_NodeList[match(Dual_EdgeList$node_1, Dual_NodeList$nodes), 6],
                     yend = Dual_NodeList[match(Dual_EdgeList$node_2, Dual_NodeList$nodes), 6]), size = 0.3) +
    geom_point(data = Dual_NodeList, aes(x, y), size = 3, shape = 21, fill = 'grey', stroke = 0.5) +
    theme(legend.position = 'none') +
    coord_fixed(ratio = 1) 
  p
  
  
  
  # save edge list to file
  Dual_EdgeList <- Dual_EdgeList[, c('from_new', 'to_new')] %>%
    `colnames<-` (c('node_1', 'node_2'))
  
  ### Because handling in python, idx should start from 0
  Dual_EdgeList$node_1 <- Dual_EdgeList$node_1 - 1
  Dual_EdgeList$node_2 <- Dual_EdgeList$node_2 - 1
  
  
  write.csv(Dual_EdgeList, paste0('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/networks/', pid, '/cell_edges.csv'), row.names = FALSE)
  
  
  # save node features
  node_features <- Dual_NodeList %>%
    select(-c('ImageID', 'ObjectNumber', 'nodes'))
  
  colnames(node_features) <- paste0("x_", 0:(ncol(node_features) - 1))
  node_features$x_3 <- as.numeric(as.factor(node_features$x_3))
  
  write.csv(node_features, paste0('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/networks/', pid, '/cell_features.csv'), row.names = FALSE)
  
  
  # save graph
  #json_data <- list(
  #  "edges" = lapply(1:nrow(Dual_EdgeList), function(i) c(Dual_EdgeList$node_1[i], Dual_EdgeList$node_2[i])),
  #  'features' = setNames(as.list(node_features$Phenotype), as.character(node_features$nodes_new - 1))
  #)
  json_data <- list(
    "0" = lapply(1:nrow(Dual_EdgeList), function(i) c(Dual_EdgeList$node_1[i], Dual_EdgeList$node_2[i]))
  )
  json_text <- toJSON(json_data, auto_unbox = TRUE, pretty = TRUE)
  writeLines(json_text, paste0('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/networks/', pid, '/cell_grpah.json'))
  
}




#--------------------------------------------------------#
#-------------------Graph Embedding ---------------------#


allimgs <- list.dirs('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/networks/', recursive = FALSE)

embedding_all_cell <- data.frame(matrix(nrow = 0, ncol = 0))
for(img in allimgs){
  
  
  # get the image id
  image_id <- strsplit(img, '//')[[1]][2]
  print(image_id)
  # patient id
  PatientID <- ImageID_PatientID_Dict[ImageID_PatientID_Dict$ImageID == image_id, ]$PatientID %>%
    unique()
  
  # read embedding file
  embedding <- read.csv(paste0(img, '/cell_embedding.csv')) %>%
    dplyr::filter(id == '0')
  embedding$id <- PatientID
  
  
  embedding_all_cell <- rbind.data.frame(embedding_all_cell, embedding)
  
  
  
  }

colnames(embedding_all_cell)[1] <- 'PatientID'

pt_embeddings_mean <- embedding_all_cell %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),mean, .names = "{col}_mean"), .groups = 'drop') %>%
  as.data.frame()

pt_embeddings_min <- embedding_all_cell %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),min, .names = "{col}_min"), .groups = 'drop') %>%
  as.data.frame() %>%
  select(-PatientID)

pt_embeddings_max <- embedding_all_cell %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),max, .names = "{col}_max"), .groups = 'drop') %>%
  as.data.frame() %>%
  select(-PatientID)



pt_embeddings_cell <- cbind.data.frame(pt_embeddings_mean, pt_embeddings_min, pt_embeddings_max)

pca_res <- prcomp(pt_embeddings_cell[,2:751], scale. = TRUE)


label <- merge(pt_embeddings_cell, clinical_data, by ='PatientID')



p <- autoplot(pca_res, data = label, fill = 'pCR', size = 4, shape = 21, color = 'black') +
  scale_fill_manual(values = c('pCR' = '#3a569e','RD' = '#ae2024')) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.position = 'none') +
  xlab('PC1') +
  ylab('PC2')
ggsave(p, file="Figures/PCA_cell_features.png", width = 5, height = 5, units = "in", dpi = 300, bg="transparent")

#saveRDS(pt_embeddings_cell, 'cell_embedding.rds')



#------------------------------------------------------------------------#
#---------------------Prepare Graph data for deep learning --------------#
#-------------------          Vessel Graph   ----------------------------#
library(jsonlite)
vesselFiles <- list.files('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/Vessel_Boundary') %>%
  data.frame() %>%
  `colnames<-` (c('filename'))



for(pid in unique(val_cell_USETHIS$ImageID)){
  
  pid <- 'NTImg0222'
  print(pid)
  try({
    #pid <- 'NT181'
    # get the single cell data for the current patient object
    sg_obj <- val_cell_USETHIS %>%
      dplyr::filter(ImageID == pid)
    
    
    
    # get the serial number to match the boundary data
    imagnum <- unique(sg_obj$ImageID)
    
    # read vessel data
    sg_vessel <- read.csv(paste0('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/Vessel_Boundary/', 
                                 vesselFiles[vesselFiles$filename %like% imagnum,]))
    
    
    # read single cell boundary data
    #sg_boundary <- read.csv(paste0('./SingleCell_Boundary/', 
    #                              boundaryFiles[boundaryFiles$filename %like% serial,]))
    
    
    pos_vessel <- holePolygon(sg_vessel)
    
    
    #---------------------------------------------#
    #------- Plot Vasculature for Validation -----#
    #---------------------------------------------#
    
    #---------- Point in polygon --------#
    # point in polygon test
    result <- apply(sg_obj[, c('Location_Center_X', 'Location_Center_Y')], 1, function(x) 
      is_point_in_polygons(x, pos_vessel))
    
    sg_obj <- sg_obj %>%
      mutate(in_vasculature = result) %>%
      dplyr::filter(in_vasculature == FALSE)    

    
    
    # get the centroid of each vessel
    ids <- unique(pos_vessel$id)
    my_polygons_sp <- lapply(ids, function(m) {
      df <- pos_vessel %>%
        dplyr::filter(id == m)
      p <- Polygon(cbind(df$x, df$y))
      Polygons(list(p), ID = as.character(m))
    }) %>% SpatialPolygons()
    
    centroids <- gCentroid(my_polygons_sp, byid = TRUE)
    centroids_df <- data.frame(id = ids,
                               x = coordinates(centroids)[, 1],
                               y = coordinates(centroids)[, 2]) 
    
    
    
    # find the nearest neighbor of each vessel
    nn_vessel <- nn2(data = centroids_df[, c('x', 'y')], query = centroids_df[, c('x', 'y')],
                     k = 2, treetype = 'kd', searchtype = 'priority') %>%
      .$nn.idx %>%
      as.data.frame() %>%
      select(V2)
    
    
    
    Dual_EdgeList <- cbind.data.frame(centroids_df$id, nn_vessel) %>%
      `colnames<-` (c('col1', 'col2'))
    
    Dual_EdgeList <- Dual_EdgeList %>%
      mutate(from = pmin(col1, col2),
             to = pmax(col1, col2)) %>%
      distinct(from, to)
    
    
    p <- ggplot() + 
      theme_void() +
      geom_polygon(data = pos_vessel, aes(x=x, y=y, group=id), fill = 'grey') +
      #geom_path(data=pos_vessel, aes(x=x, y=y, group=draw), color = 'black', size=1) +  #draw is the original ring
      #scale_fill_discrete("ID") +
      #geom_segment(aes(x = centroids_df[match(Dual_EdgeList$from, centroids_df$id), 2],
      #                 xend = centroids_df[match(Dual_EdgeList$to, centroids_df$id), 2],
      #                 y = centroids_df[match(Dual_EdgeList$from, centroids_df$id), 3],
      #                 yend = centroids_df[match(Dual_EdgeList$to, centroids_df$id), 3]), linewidth = 0.3) +
      #geom_point(data = centroids_df, aes(x, y), size = 3, shape = 21, fill = 'grey', stroke = 0.5) +
      #theme(legend.position = 'none') +
      coord_fixed(ratio = 1) 
    p
    ggsave(p, file=paste0("../TNBC_ClinicalTrial_Validation_Cohort/Vessel_Graph/", pid, ".png"), width = 7, height = 7, units = "in", dpi = 300, bg="black")
    
    
    # node features
    pe_mean <- sapply(ids, function(m) {
      # get df for the vessel
      #m <- 1
      df <- pos_vessel %>%
        dplyr::filter(id == m) %>%
        select(x, y)
      
      # find nearest neighbors
      
      nn_idx <- nn2(df, query = sg_obj[, c('Location_Center_X', 'Location_Center_Y')], k = nrow(df),
                    treetype = 'kd', searchtype = 'radius', radius = 20) %>%
        .$nn.idx
      
      non_zero_rows <- which(apply(nn_idx, 1, function(x) all(x == 0), simplify = TRUE) == FALSE)
      
      # use the index to get the protein expression profile
      
      #----------------- For molecular phenotyping ---------------#
      pe <- sg_obj[non_zero_rows,] %>%
        select(15:60) %>%
        dplyr::summarise(across(everything(), mean, na.rm = TRUE)) %>%
        mutate(across(everything(), ~ifelse(is.na(.), 0, .)))
      
      
    }) %>%
      t() %>%
      `rownames<-` (c(ids))
    
    
    
 
    #--------------------------------------------------------------------------#
    #----------------------- Prepare data for save ----------------------------#
    colnames(Dual_EdgeList)[1] <- 'node_1'
    colnames(Dual_EdgeList)[2] <- 'node_2'
    
    Dual_EdgeList$node_1 <- Dual_EdgeList$node_1 - 1
    Dual_EdgeList$node_2 <- Dual_EdgeList$node_2 - 1
    #write.csv(Dual_EdgeList, paste0('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/networks/', pid, '/vessel_edges.csv'), row.names = FALSE)
    
    
    # save graph
    json_data <- list(
      "0" = lapply(1:nrow(Dual_EdgeList), function(i) c(Dual_EdgeList$node_1[i], Dual_EdgeList$node_2[i]))
    )
    json_text <- toJSON(json_data, auto_unbox = TRUE, pretty = TRUE)
    #writeLines(json_text, paste0('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/networks/', pid, '/vessel_graph.json'))

    
    
    
    # save node features
    pe_mean <- pe_mean %>%
      as.data.frame() %>%
      `colnames<-` (paste0("x_", 0:(ncol(pe_mean) - 1)))
    
    pe_mean <- data.frame(lapply(pe_mean, as.numeric))
    
    #write.csv(pe_mean, paste0('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/networks/', pid, '/vessel_features.csv'), row.names = FALSE)
    
  })
  
}



#--------------------------------------------------------#
#-------------------Graph Embedding ---------------------#


allimgs <- list.dirs('../TNBC_ClinicalTrial_Validation_Cohort/Codes/data/networks/', recursive = FALSE)

embedding_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(img in allimgs){
  
  print(img)
  # get the image id
  try({
    image_id <- strsplit(img, '//')[[1]][2]
    
    # patient id
    PatientID <- ImageID_PatientID_Dict[ImageID_PatientID_Dict$ImageID == image_id, ]$PatientID %>%
      unique()
    
    # read embedding file
    embedding <- read.csv(paste0(img, '/vessel_embedding.csv'))
    embedding$id <- PatientID
    
    
    embedding_all <- rbind.data.frame(embedding_all, embedding)
    
    
  })
  
  
}

library(ggfortify)



vessel_graph_features <- embedding_all
colnames(vessel_graph_features) <- 'PatientID'
colnames(vessel_graph_features)[2:251] <- (paste0("x_", 250:(ncol(vessel_graph_features) + 248)))
 

vessel_graph_features_mean <- vessel_graph_features %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),mean, .names = "{col}_mean"), .groups = 'drop') %>%
  as.data.frame()


vessel_graph_features_max <- vessel_graph_features %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),max, .names = "{col}_max"), .groups = 'drop') %>%
  as.data.frame() %>%
  select(-PatientID)

vessel_graph_features_min <- vessel_graph_features %>%
  group_by(PatientID) %>%
  dplyr::summarise(across(everything(),min, .names = "{col}_min"), .groups = 'drop') %>%
  as.data.frame() %>%
  select(-PatientID)

pt_embeddings_vessel <- cbind.data.frame(vessel_graph_features_mean, vessel_graph_features_min, vessel_graph_features_max)


pca_res <- prcomp(pt_embeddings_vessel[,2:751], scale. = TRUE)


label <- merge(pt_embeddings_vessel, clinical_data, by ='PatientID')



p <- autoplot(pca_res, data = label, fill = 'pCR', size = 4, shape = 21, color = 'black') +
  scale_fill_manual(values = c('pCR' = '#3a569e','RD' = '#ae2024')) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.position = 'none') +
  xlab('PC1') +
  ylab('PC2')
p
ggsave(p, file="Figures/PCA_vessel_features.png", width = 5, height = 5, units = "in", dpi = 300, bg="transparent")






#-----------------------------------------------------------#
#-------------------Combine PCA ----------------------------#


pt_embeddings <- cbind.data.frame(pt_embeddings_cell, pt_embeddings_vessel[,-1])


pca_res <- prcomp(pt_embeddings[,-1], scale. = TRUE)


label <- merge(pt_embeddings, clinical_data, by ='PatientID')



p <- autoplot(pca_res, data = label, fill = 'pCR', size = 4, shape = 21, color = 'black') +
  scale_fill_manual(values = c('pCR' = '#3a569e','RD' = '#ae2024')) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.position = 'none') +
  xlab('PC1') +
  ylab('PC2')
p
ggsave(p, file="Figures/PCA_all_features.png", width = 5, height = 5, units = "in", dpi = 300, bg="transparent")

