#-------------------------------------------#
# METABRIC TNBC validation process using Ali's new dataset-------#
# @Author: Haoyang Mi ----------------------#
# Date: Sep 26th 2023------------------------#
# REF: 


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



#----------------------Data Preparation Finish -------------------------------#


#------------------------ Analysis Starts --------------------------------#


#---------------- Section 1: read vessel masks and do the processing -----#
# Plot Cell Boundaries for illustrative purposes
vesselFiles <- list.files('../TNBC_ClinicalTrial_Validation_Cohort/data/Vessel_Boundary') %>%
  data.frame() %>%
  `colnames<-` (c('filename'))






vasculature_density <- data.frame(matrix(nrow = 0, ncol = 0))
# parameters
area_thresh <- 500 #500, 700, 800
core_area <- pi * (0.5)^2
# factorize cell type
val_cell_USETHIS$Phenotype <- as.factor(val_cell_USETHIS$Phenotype)

for(pid in unique(val_cell_USETHIS$ImageID)){
  
  print(pid)
  try({
    #pid <- 'NT181'
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
      dplyr::filter(id %in% vesselmaskTRUE)
    
    
    pos_vesselFALSE <- pos_vessel %>%
      dplyr::filter(!(id %in% vesselmaskTRUE))
    
    gg3 <- ggplot(pos_vesselTRUE, aes(x=x, y=y)) + 
      theme_void() +
      geom_polygon(aes(group=id), fill = 'grey') +
      geom_path(data=pos_vesselTRUE, aes(x=x, y=y, group=draw), color = 'black', size=1) +  #draw is the original ring
      scale_fill_discrete("ID") #+
    gg3
    
    
    #----------- Vasculature density ----------------#
    #mid <- 'MB-0399'
    #vasculature_area <- compute_polygons_area(pos_vessel, area_thresh)
    
    # vasculature density
    #vasculature_density <- rbind.data.frame(vasculature_density, cbind.data.frame(mid, vasculature_area / core_area))
    vasculature_density <- rbind.data.frame(vasculature_density, cbind.data.frame(imagnum, length(vesselmaskTRUE) / core_area))
    #----------- Cell infiltration profile w.r.t vasculature ----------------#
    
    
  })
  
}

colnames(vasculature_density) <- c('ImageID', 'density')

image_dict <- unique(val_cell_USETHIS[, c('PatientID', 'ImageID')])

test <- merge(vasculature_density, image_dict, by = 'ImageID')
test <- merge(test, clinical_data, by = 'PatientID')


test1 <- test[test$pCR == 'pCR' & test$BiopsyPhase == 'Post-treatment', 'density']
test2 <- test[test$pCR == 'RD' & test$BiopsyPhase == 'Post-treatment', 'density']

mean(test1)
mean(test2, na.rm = TRUE)
t.test(test1, test2)



boxplot(test1, test2)



# STMAP^2 to evaluate 
pvals <- data.frame(matrix(nrow = 0, ncol = 0))
pe_clustered <- data.frame(matrix(nrow = 0, ncol = 0))
pe_random <- data.frame(matrix(nrow = 0, ncol = 0))
area_thresh <- 442.3


for(pid in unique(val_cell_USETHIS$ImageID)){
  
  print(pid)
 # mid <- 'MB-0350'
  # get the single cell data for the current metabric object
  sg_obj <- val_cell_USETHIS %>%
    dplyr::filter(ImageID == pid)
  
  
  
  # get the image number to match the boundary data
  imagnum <- unique(sg_obj$ImageID)
  
  # read vessel data
  sg_vessel <- read.csv(paste0('../TNBC_ClinicalTrial_Validation_Cohort/data/Vessel_Boundary/', 
                               vesselFiles[vesselFiles$filename %like% imagnum,]))
  
  
  pos_vessel <- holePolygon(sg_vessel)
  
  #---- remove vessels with very small sizes
  vesselmaskTRUE <- filter_vessel(pos_vessel, area_thresh)
  
  pos_vesselTRUE <- pos_vessel %>%
    dplyr::filter(id %in% vesselmaskTRUE)
  pos_vesselFALSE <- pos_vessel %>%
    dplyr::filter(!(id %in% vesselmaskTRUE))
  
  
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
      
      
      
      #pe[is.na(pe)] <- 0
      #----------------- For cellular phenotyping ---------------#
      #pe <- sg_metabricobj[non_zero_rows,] %>%
      #  group_by(Phenotype, .drop = FALSE) %>%
      #  tally() %>%
      #  spread(key = Phenotype, value = n, fill = 0)
      
      
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
    
    
    #pe_mean[is.na(pe_mean)] <- 0
    pe_mean <- pe_mean[apply(pe_mean, 1, function(row) any(!is.na(row))), ]
    
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
      pe_clustered <- rbind.data.frame(pe_clustered, cbind.data.frame(pe_mean, imagnum))
    }
    if(mantel_test$signif >= 0.05){
      # these are the proteomic profiles stacked for all valid vessels
      pe_random <- rbind.data.frame(pe_random, cbind.data.frame(pe_mean, imagnum))
    }
    
    #colorcode <- 'red'
    
    
    #------------------- Plot --------------------#
    
    #       PLACE HOLDER FOR TISSUE PLOT          #
    
    #---------------------------------------------#
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
    #p
    #ggsave(p, file=paste0("Figures/Random_correlation.png"), width = 6.5, height = 6, units = "in", dpi = 300)
    
    
    
    pvals <- rbind.data.frame(pvals, cbind.data.frame(imagnum, mantel_test$signif, length(ids)))
  }
  
  
}

colnames(pvals) <- c('ImageID', 'Signif', 'Num')

image_dict <- unique(val_cell_USETHIS[, c('PatientID', 'ImageID')])

test <- merge(pvals, image_dict, by = 'ImageID')
test <- merge(test, clinical_data, by = 'PatientID')
test$Label <- ifelse(test$Signif < 0.05 , 'Clustered', 'Random')



test5 <- test %>%
  dplyr::filter(BiopsyPhase == 'On-treatment')

pt_DF <- data.frame(matrix(nrow = 0, ncol = 0))
for(pid in unique(test5$PatientID)){
  
  
  n1 <- test5 %>%
    dplyr::filter(PatientID == pid) %>%
    dplyr::filter(Label == 'Random') %>%
    nrow()
  
  n2 <- test5 %>%
    dplyr::filter(PatientID == pid) %>%
    dplyr::filter(Label == 'Clustered') %>%
    nrow()
  
  
  if(n1 > n2){
    pt_label <- 'Random'
  }
  if(n2 > n1){
    pt_label <- 'Clustered'
  }
  
  if(n2 == n1){
    pt_label <- 'Undefined'
  }
  
  pt_DF <- rbind.data.frame(pt_DF, cbind.data.frame(pid, pt_label))
  
  
  
  
}

colnames(pt_DF) <- c('PatientID', 'Label')

dict <- unique(clinical_data[, c('PatientID', 'pCR')])


merge(pt_DF, dict, by = 'PatientID') %>%
  group_by(pCR, Label) %>%
  tally()


test_pos <- val_cell_USETHIS %>%
  dplyr::filter(BiopsyPhase == 'Baseline') 

test_pos <- test_pos[1:100000,]

plot(test_pos$Location_Center_X, test_pos$Location_Center_Y)


test1 <- image_dict[!(image_dict$ImageID %in% test$ImageID), ]
test1 <- merge(test1, clinical_data, by = 'PatientID')

test1 %>%
  group_by(BiopsyPhase, pCR) %>%
  tally()

length(unique(test1[test1$BiopsyPhase.x == 'Post-treatment','ImageID']))
length(unique(test1[test1$BiopsyPhase.y == 'Baseline' & test1$pCR == 'pCR','ImageID']))
length(unique(test1[test1$BiopsyPhase.y == 'Post-treatment' & test1$pCR == 'pCR','ImageID']))


library(knitr)


observed <- matrix(c(461, 432, 46, 46), nrow = 2)
rownames(observed) <- c("Baseline", "Post-treatment")
colnames(observed) <- c("Depleted", "Random")
observed

chi_result <- chisq.test(observed)

print(chi_result)
















# Deprecated
#-----------------------------------------------------------------------#
#-------- CD8+ T, Granulocytes, Endothelial cells interactions ---------#
#-----------------------------------------------------------------------#

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


#boxplot(test1$ratio, test2$ratio)
#-----------------------------------------------#
# Are there patient-wise increase or decrease ? #
#-----------------------------------------------#
colnames(ratios_all)[2] <- 'ImageNumber'
iid_patient_dict <- unique(val_cell_USETHIS[, c('ImageNumber', 'PatientID')])

testttt <- ratios_all %>%
  merge(iid_patient_dict, by = 'ImageNumber') %>%
  group_by(PatientID, pCR, Biopsy) %>%
  dplyr::summarise(mean_ratio = mean(ratio)) %>%
  dplyr::filter(pCR == 'RD')


t1 <- testttt[testttt$Biopsy == 'Post-treatment', ]
t2 <- testttt[testttt$Biopsy == 'Baseline', ]

t3 <- rbind.data.frame(t1, t2)
# Identify IDs that are present in all unique groups
valid_ids <- t3 %>%
  group_by(PatientID) %>%
  dplyr::summarise(n_unique_groups = n_distinct(Biopsy)) %>%
  filter(n_unique_groups == length(unique(t3$Biopsy))) %>%
  select(PatientID)

# Filter the original data frame to only include valid IDs
filtered_df <- t3 %>%
  data.frame() %>%
  dplyr::filter(PatientID %in% valid_ids$PatientID)


wilcox.test(filtered_df[filtered_df$Biopsy == 'Baseline', 'mean_ratio'], filtered_df[filtered_df$Biopsy == 'Post-treatment', 'mean_ratio'], paired = TRUE)
colnames(t3)[2] <- 'group'
# create a dictionary
pd <- position_dodge(0.1) # move them .05 to the left and right

p <- ggplot(filtered_df, aes(Biopsy,mean_ratio)) +
  theme_bw() +
  geom_boxplot(aes(fill=Biopsy),outlier.shape = NA,alpha=0.6, size = 1) +
  geom_line(aes(group=PatientID), position = pd, size = 1) +
  geom_point(aes(fill=Biopsy,group=PatientID),size=5,shape=21, position = pd) +
  #facet_grid(~Side,scales = "free") +
  #scale_color_manual(values = c('Baseline' = '#5ba6ea', 'Post-treatment' = "#312cb4")) +
  #scale_fill_manual(values = c('Baseline' = '#5ba6ea', 'Post-treatment' = "#312cb4")) +
  scale_color_manual(values = c('Baseline' = '#848484', 'Post-treatment' = "black")) +
  scale_fill_manual(values = c('Baseline' = '#848484', 'Post-treatment' = "black")) +
  
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26),
        legend.position = 'none')  +
  ylab('STAIN') +
  geom_signif(comparisons = list(c("Baseline", "Post-treatment")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 0.9, annotations = '*'
  ) +
  ylim(0, 1) 
p
ggsave(p, file="Figures/RD_STAIN_Base_Post.png", width = 6, height = 8, units = "in", dpi = 300)














#

library(Rmisc)
library(ggpubr)
library(ggsignif)
tgc <- summarySE(ratios_all, measurevar= "ratio", groupvars=c("pCR"))
mycomparison <- compare_means(ratio ~ pCR,  data = ratios_all)
mycomparison

pd <- position_dodge(0.1) # move them .05 to the left and right

p <- ggplot(tgc, aes(x= pCR, y = ratio, colour=pCR, group=pCR)) + 
  theme_bw() +
  geom_errorbar(aes(ymin=ratio - 10*se, ymax=ratio + 10*se, color = pCR), width=.05, size = 1) +
  geom_point(aes(fill = pCR, color = 'black'), position=pd, size=5, shape=21) + # 21 is filled circle
  
  #scale_color_manual(values = c('pCR' = '#9dc9e6', 'RD' = "#787878")) +
  #scale_fill_manual(values = c('pCR' = '#9dc9e6', 'RD' = "#787878")) +
  scale_color_manual(values = c('pCR' = '#1e19b3', 'RD' = "black")) +
  scale_fill_manual(values = c('pCR' = '#1e19b3', 'RD' = "black")) +
  expand_limits(y=0) +                        # Expand y range
  ylim(0, 1) +
  geom_signif(comparisons = list(c("pCR", "RD")), 
              map_signif_level=TRUE, size = 0.5, textsize = 10, color = 'black', 
              y_position = 0.9, annotations = '****'
  ) +
  theme_prism() +
  theme(axis.title.x = element_blank(),
        #axis.text.x = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.text = element_text(size = 24),
        axis.title = element_text(size = 26),
        legend.position = 'none')  +

  ylab('STAIN')
p

ggsave(p, file="Figures/RISC_Post.png", width = 5, height = 7, units = "in", dpi = 300)












