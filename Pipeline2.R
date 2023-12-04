
#-------------------------------------------#
# METABRIC TNBC processing pipeline 2-------#
# @Author: Haoyang Mi ----------------------#
# Date: April 27th 2023------- --------------#
# REF: 


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
  # 3. 归一化至 0-1 范围
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

testt <- sg_cell_expr[sg_cell_expr$cellPhenotype == 'Fibroblasts',]
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










#---------------Comparison of cell types between different clinical groups ---------------#

selectedDF_dup <- selectedDF
# CD8 / Tumor cell ratio
colnames(selectedDF_dup)[2] <- 'metabric_id'



selectedDF_dup$Stage <- ifelse(selectedDF_dup$Stage <= 2, 'Stage I-II', 'Stage III-IV')
selectedDF_dup$`Age at Diagnosis` <- ifelse(selectedDF_dup$`Age at Diagnosis` <= 55, '<= 55 years of age', '>55 years of age')
selectedDF_dup$Grade <- ifelse(selectedDF_dup$Grade <= 2, 'Grade I-II', 'Grade III')
selectedDF_dup$`Relapse Free Status` <- ifelse(selectedDF_dup$`Relapse Free Status` == '1:Recurred', 'Recurred', 'Not recurred')
selectedDF_dup$Survival <- ifelse(selectedDF_dup$Survival < 66, 'Survival low', 'Survival high')
selectedDF_dup$`TMB (nonsynonymous)` <- ifelse(selectedDF_dup$`TMB (nonsynonymous)` < 5, 'TMB low', 'TMB high')
selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'High'] <- 'Cellularity high'
selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'Moderate'] <- 'Cellularity moderate'
selectedDF_dup['Cellularity'][selectedDF_dup['Cellularity'] == 'Low'] <- 'Cellularity low'

table(selectedDF_dup$`Age at Diagnosis`)

cd8_tumor_ratio <- sg_cell_USETHIS %>%
  group_by(metabric_id, Phenotype) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  dplyr::select(-metabric_id) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate(across(everything(), ~ . * 100/ row_sum)) %>%
  replace(is.na(.), 0) %>%
  mutate(metabric_id = selectedDF_dup$metabric_id) %>%
  merge(selectedDF_dup) %>%
  as.data.frame() %>%
  #mutate(ratio = (`CD8+ T cells` + `CD4+ T cells` + `B cells` + `CD38+ lymphocytes` + `CD4+ T cells & APCs`
  #                + `Tregs and Tex`) / (`Macrophages` + `Granulocytes` + `Macrophages & granulocytes`)) %>%
  #mutate(ratio = (`CD8+ T cells` + `CD4+ T cells` +  `CD4+ T cells & APCs`) / `Fibroblasts`) %>%
  mutate(ratio =  `CD4+ T cells` / `Macrophages`) %>%
  dplyr::select(metabric_id, ratio , `Age at Diagnosis`, Survival, Grade, Cellularity, `Primary Tumor Laterality`, `TMB (nonsynonymous)`, `Relapse Free Status`, Stage) %>%
  #dplyr::select(metabric_id, `B cells` , `Age at Diagnosis`, Survival, Grade, Cellularity, `Primary Tumor Laterality`, `TMB (nonsynonymous)`, `Relapse Free Status`, Stage) %>%
  pivot_longer(cols = `Age at Diagnosis`:Stage, names_to = 'Category') %>%  
  filter(complete.cases(.)) %>%
  filter(value != 'Normal') %>%
  filter(ratio != 'Inf')





#ggboxplot(longTable, x = 'response', y = 'value', palette = 'jco', color = 'response', facet.by = c('variable')) +
cd8_tumor_ratio$value <- factor(cd8_tumor_ratio$value, levels = c('<= 55 years of age',  '>55 years of age', 'Cellularity low', 
                                                                  'Cellularity moderate', 'Cellularity high',
                                                                  'Basal', 'claudin-low', 'Her2', 'Grade I-II', 'Grade III',
                                                                  'Stage I-II', 'Stage III-IV', 'Left', 'Right', 'Not recurred', 'Recurred',
                                                                  'TMB low', 'TMB high', 'Survival low', 'Survival high'))

test1 <- cd8_tumor_ratio[cd8_tumor_ratio$value == 'Stage I-II',]
test2 <- cd8_tumor_ratio[cd8_tumor_ratio$value == 'Stage III-IV',]
wilcox.test(test1$ratio, test2$ratio)

test1 <- cd8_tumor_ratio[cd8_tumor_ratio$value == 'TMB low',]
test2 <- cd8_tumor_ratio[cd8_tumor_ratio$value == 'TMB high',]
wilcox.test(test1$ratio, test2$ratio)

test1 <- cd8_tumor_ratio[cd8_tumor_ratio$value == 'Survival high',]
test2 <- cd8_tumor_ratio[cd8_tumor_ratio$value == 'Survival low',]
wilcox.test(test1$ratio, test2$ratio)


test1 <- cd8_tumor_ratio[cd8_tumor_ratio$value == 'Recurred',]
test2 <- cd8_tumor_ratio[cd8_tumor_ratio$value == 'Not recurred',]
wilcox.test(test1$ratio, test2$ratio)



p <- ggplot(cd8_tumor_ratio, aes(x = value, y = ratio), fill = 'black', fill = 'transparent') +
  stat_boxplot( aes(value, ratio), 
                geom='errorbar', linetype=1, width=0.3)+  #whiskers  #stat_pvalue_manual(
  #  size = 6, textsize = 12, color = 'black', 
  #  stat.test, bracket.nudge.y = -2, hide.ns = FALSE,
  #  label = "{p.adj.signif}"
  #) +
  geom_boxplot(aes(value, ratio),  outlier.size = -1) +
  geom_jitter(shape = 21, size = 2, aes(fill = value), alpha = 0.6)+
  
  ylab("Cell abundance" ) +
  #scale_fill_manual(values = c('LumA' = '#2199bd', 'Her2' = '#ec6d6d')) +
  facet_wrap(~value, ncol = 17, scales = 'free_x') +
  
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        strip.text = element_blank(),
        axis.text.x = element_text(angle = 90),
        panel.spacing.x = unit(0, "pt"),
        legend.position = 'none',
        axis.title = element_text(size = 20)) +
  xlab('') +
  ylab(expression(CD4^'+'~'T cell-macrophage ratio')) +
  #ylab(expression('T cell'-'Fibroblasts ratio')) +
  ylim(0, 3)

p
ggsave(p, file=paste0("Figures/CD4T cell-macrophage_ratio.png"), width = 10, height = 8, units = "in", dpi = 300)


#---------------------------------------------------#
#--------------- Cell-cell interaction -------------#
# --------------------------------------------------#
aType <- unique(sg_cell_USETHIS$Phenotype)

compHeatmap <- data.frame(matrix(nrow = 9, ncol = 15))


colorHeatmap <- data.frame(matrix(nrow = 9, ncol = 15))
colorDictionary <- data.frame(`A term` = c('#0868ac', '#ece7f2', '#ece7f2', '#a6bedb', '#9ebedb', '#8c96c6', '#4d0149', '#8c6bb2', '#31a353'),
                              `B term` = c('#a8deb5', '#3690c0', '#4db4d4', '#ffebe3', '#ffffcc', '#fdc5c0', '#7bccc4', '#fa9fb5', '#c2e699'))

sizeHeatmap <- data.frame(matrix(nrow = 9, ncol = 15))

colnames(compHeatmap) <- aType
colnames(colorHeatmap) <- aType
colnames(sizeHeatmap) <- aType


rowVariables <- c('<= 55 years of age',  '>55 years of age', 'Cellularity low', 
                  'Cellularity moderate', 'Cellularity high', 'skip', 'Grade I-II', 'Grade III',
                  'Stage I-II', 'Stage III-IV', 'Left', 'Right', 'Not recurred', 'Recurred',
                  'TMB low', 'TMB high', 'Survival low', 'Survival high')

DF <- sg_cell_USETHIS %>%
  group_by(metabric_id, Phenotype) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  dplyr::select(-metabric_id) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate(across(everything(), ~ . * 100/ row_sum)) %>%
  replace(is.na(.), 0) %>%
  mutate(metabric_id = selectedDF_dup$metabric_id) %>%
  merge(selectedDF_dup) %>%
  as.data.frame()


test1 <- DF[DF$Survival == 'Survival low',]
test2 <- DF[DF$Survival == 'Survival high',]

wilcox.test(test1$`CD8+ T cells`, test2$`CD8+ T cells`)
#DF <- sg_cell_USETHIS %>%
#  group_by(metabric_id, Phenotype) %>%
#  tally() %>%
#  pivot_wider(names_from = Phenotype, values_from = n) %>%
#  replace(is.na(.), 0) %>%
#  mutate(ratio = `CD8+ T cells` / Tumor) %>%
#  merge(selectedDF_dup) %>%
#  as.data.frame() 

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
      
      
      rv_2_data <- DF %>%
        filter(apply(., 1, function(row) any(sapply(rv_2, grepl, row)))) %>%
        dplyr::select(all_of(ct)) %>%
        as.matrix()
      
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
      
      
      rv_2_data <- DF %>%
        filter(apply(., 1, function(row) any(grepl(all_of(rv_2), row))))%>%
        dplyr::select(all_of(ct)) %>%
        as.matrix()
      
    }
    
    # statistical test
    
    
    pval <- wilcox.test(rv_1_data, rv_2_data)$p.value
    
    
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

p <- ggplot(melted_matrix, aes(x = Var2, y = Var1, fill = custom_color, size = custom_size)) +
  geom_point(shape = 21) +
  scale_fill_identity() +
  scale_size_identity() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 24),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) 
p

ggsave(p, file=paste0("Figures/Cell_Freq_Comparison.png"), width = 17, height = 12, units = "in", dpi = 300)









#------------- Survival Analysis on Clustered Tissue Architectures ------------#

selectedDF_dup2 <- selectedDF %>%
  mutate(group = ifelse(Survival > survival_thresh, 'long', 'short')) %>%
  data.frame()

fit <- survfit(Surv(Survival, Status) ~ group, data = selectedDF_dup2)
fit


p <- ggsurvplot(fit, data = selectedDF_dup2, 
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
                #pval = TRUE
) +
  xlab('Time, (months)')
p
ggsave(file= paste("./Figures/Survival" , '.pdf'), plot = p$plot, width = 9, height = 6, dpi = 300)


colnames(selectedDF)[2] <- 'metabric_id'

DF <- sg_cell_USETHIS %>%
  group_by(metabric_id, Phenotype) %>%
  tally() %>%
  pivot_wider(names_from = Phenotype, values_from = n) %>%
  replace(is.na(.), 0) %>%
  as.data.frame() %>%
  dplyr::select(-metabric_id) %>%
  mutate(row_sum = rowSums(.)) %>%
  mutate(across(everything(), ~ . * 100/ row_sum)) %>%
  replace(is.na(.), 0) %>%
  mutate(metabric_id = selectedDF_dup$metabric_id) %>%
  merge(selectedDF_dup) %>%
  as.data.frame()


DF$ratio <- DF$`CD8+ T cells` / DF$Tumor


test1 <- DF[DF$`Age at Diagnosis` == '<= 55 years of age',]
test2 <- DF[DF$`Age at Diagnosis` == '>55 years of age',]
test1 <- DF[DF$`Relapse Free Status` == 'Not recurred',]
test2 <- DF[DF$`Relapse Free Status` == 'Recurred',]

wilcox.test(test1$ratio, test2$ratio)
wilcox.test(test1$`CD8+ T cells`, test2$`CD8+ T cells`)

mycomparison <- compare_means(`CD4+ T cells` ~ Survival,  data = DF)
mycomparison

my_comparisons <- list( c("Survival low", "Survival high"))

p <- ggplot(DF, aes(Survival, Fibroblasts, color = Survival)) +
  stat_boxplot( aes(Survival, Fibroblasts), 
                geom='errorbar', linetype=1, width= 0.2, size = 2)+  #whiskers
  geom_boxplot(size = 2) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width = 0.3), 
             aes(fill=factor(Survival)), show.legend = F, size = 4, shape = 21, color = 'black') +
  theme_foundation() +
  scale_color_manual(values = c('Survival low' = '#31a353', 'Survival high' = '#c2e699')) +
  scale_fill_manual(values = c('Survival low' = '#31a353', 'Survival high' = '#c2e699')) +
  theme_prism() +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 1)) +
  geom_signif(comparisons = list(c("Survival low", "Survival high")), 
              map_signif_level=TRUE, size = 1, textsize = 20, color = 'black', 
              y_position = max(DF$Fibroblasts) + 1 , annotations = mycomparison$p.signif
  ) +
  scale_y_continuous(limits = c(min(DF$Fibroblasts)-1, 60), expand = expansion(mult = c(0, 0))) +
  xlab('') +
  ylab(expression(Percentage~of~'Fibroblasts'~'(%)')) 
p
ggsave(p, file=paste0("Figures/Fibroblasts_Survival.png"), width = 4, height = 7, units = "in", dpi = 300)



#------- Laterality  ---------##
colnames(DF)[37] <- 'Laterality'
DF_laterality <- DF %>%
  dplyr::filter(Laterality != 'NA')

mycomparison <- compare_means(Granulocytes ~ Laterality,  data = DF_laterality)
mycomparison

colorDictionary <- data.frame(`A term` = c('#0868ac', '#ece7f2', '#ece7f2', '#a6bedb', '#9ebedb', '#8c96c6', '#4d0149', '#8c6bb2', '#31a353'),
                              `B term` = c('#a8deb5', '#3690c0', '#4db4d4', '#ffebe3', '#ffffcc', '#fdc5c0', '#7bccc4', '#fa9fb5', '#c2e699'))
p <- ggplot(DF_laterality, aes(Laterality, Granulocytes, color = Laterality)) +
  stat_boxplot( aes(Laterality, Granulocytes), 
                geom='errorbar', linetype=1, width= 0.2, size = 2)+  #whiskers
  geom_boxplot(size = 2) +
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width = 0.3), 
             aes(fill=factor(Laterality)), show.legend = F, size = 4, shape = 21, color = 'black') +
  theme_foundation() +
  scale_color_manual(values = c('Left' = '#8c96c6', 'Right' = '#fdc5c0')) +
  scale_fill_manual(values = c('Left' = '#8c96c6', 'Right' = '#fdc5c0')) +
  theme_prism() +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 1)) +
  geom_signif(comparisons = list(c("Left", "Right")), 
              map_signif_level=TRUE, size = 1, textsize = 20, color = 'black', 
              y_position = max(DF_laterality$Granulocytes) + 2, annotations = mycomparison$p.signif
  ) +
  scale_y_continuous(limits = c(min(DF_laterality$Granulocytes)-1, 15), expand = expansion(mult = c(0, 0))) +
  xlab('') +
  ylab(expression(Percentage~of~'Granulocytes'~'(%)')) 
p
ggsave(p, file=paste0("Figures/Granulocytes_Laterality.png"), width = 4, height = 7, units = "in", dpi = 300)







#----------- Age --------------#
colnames(DF)[20] <- 'Age'

DF_age <- DF
colnames(DF_age)[5] <- 'CD4_T'
colorDictionary <- data.frame(`A term` = c('#0868ac', '#ece7f2', '#ece7f2', '#a6bedb', '#9ebedb', '#8c96c6', '#4d0149', '#8c6bb2', '#31a353'),
                              `B term` = c('#a8deb5', '#3690c0', '#4db4d4', '#ffebe3', '#ffffcc', '#fdc5c0', '#7bccc4', '#fa9fb5', '#c2e699'))

mycomparison <- compare_means(`CD4_T` ~ Age,  data = DF_age)
mycomparison

dp <- ggplot(DF_age, aes(Age, `B_cells`, color = Age)) +
  stat_boxplot( aes(Age, `CD8+ T cells`), 
                geom='errorbar', linetype=1, width= 0.2, size = 2)+  #whiskers
  geom_boxplot(size = 2, outlier.color = 'white') +
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width = 0.3), 
             aes(fill=factor(Age)), show.legend = F, size = 4, shape = 21, color = 'black') +
  theme_foundation() +
  scale_color_manual(values = c('<= 55 years of age' = '#0868ac', '>55 years of age' = '#a8deb5')) +
  scale_fill_manual(values = c('<= 55 years of age' = '#0868ac', '>55 years of age' = '#a8deb5')) +
  theme_prism() +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 1)) +
  geom_signif(comparisons = list(c("<= 55 years of age", ">55 years of age")), 
              map_signif_level=TRUE, size = 1, textsize = 20, color = 'black', 
              y_position = max(DF_age$`B_cells`) + 2, annotations = mycomparison$p.signif
  ) +
  scale_y_continuous(limits = c(min(DF_age$`B_cells`)-1, 50), expand = expansion(mult = c(0, 0))) +
  xlab('') +
  #ylab(expression(Percentage~of~'CD8'^'+'~'T cells (%)')) 
  ylab(expression(Percentage~of~'B '~'cells (%)')) 
p
ggsave(p, file=paste0("Figures/B cells_Age.png"), width = 4, height = 7, units = "in", dpi = 300)



#------- Grade  ---------##
DF_grade <- DF %>%
  dplyr::filter(Grade != 'NA')
colnames(DF_grade)[3] <- 'CD4_T'

mycomparison <- compare_means(CD4_T ~ Grade,  data = DF_grade)
mycomparison

colorDictionary <- data.frame(`A term` = c('#0868ac', '#ece7f2', '#ece7f2', '#a6bedb', '#9ebedb', '#8c96c6', '#4d0149', '#8c6bb2', '#31a353'),
                              `B term` = c('#a8deb5', '#3690c0', '#4db4d4', '#ffebe3', '#ffffcc', '#fdc5c0', '#7bccc4', '#fa9fb5', '#c2e699'))

p <- ggplot(DF_grade, aes(Grade, CD4_T, color = Grade)) +
  stat_boxplot( aes(Grade, CD4_T), 
                geom='errorbar', linetype=1, width= 0.2, size = 2)+  #whiskers
  geom_boxplot(size = 2, outlier.color = 'white') +
  geom_point(position=position_jitterdodge(jitter.width=0.2, dodge.width = 0.3), 
             aes(fill=factor(Grade)), show.legend = F, size = 4, shape = 21, color = 'black') +
  theme_foundation() +
  scale_color_manual(values = c('Grade I-II' = '#a6bedb', 'Grade III' = '#ffebe3')) +
  scale_fill_manual(values = c('Grade I-II' = '#a6bedb', 'Grade III' = '#ffebe3')) +
  theme_prism() +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 25),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 1)) +
  geom_signif(comparisons = list(c("Grade I-II", "Grade III")), 
              map_signif_level=TRUE, size = 1, textsize = 20, color = 'black', 
              y_position = max(DF_grade$CD4_T) + 2, annotations = mycomparison$p.signif
  ) +
  scale_y_continuous(limits = c(min(DF_grade$CD4_T)-1, 12), expand = expansion(mult = c(0, 0))) +
  xlab('') +
  ylab(expression(Percentage~of~'CD4'^'+'~'T cells (%)')) 
p
ggsave(p, file=paste0("Figures/CD4_Grade.png"), width = 4, height = 7, units = "in", dpi = 300)



#----------------------------------------------------#
#----------- Fibroblasts heterogeneity --------------#
#----------------------------------------------------#
library(ClusterR)

Fibroblasts <- sg_cell_USETHIS %>%
  dplyr::filter(Phenotype == 'Fibroblasts') %>%
  select(metabric_id, SMA, FSP1, PDGFRB, `Caveolin-1`, Location_Center_X, Location_Center_Y)

nFeature <- 5

set.seed(2)
km <- MiniBatchKmeans(as.matrix(Fibroblasts[,2:nFeature]), clusters = 4, batch_size = 100)
clusters <- predict_MBatchKMeans(as.matrix(Fibroblasts[,2:nFeature]), km$centroids)



Fibroblasts_clus <- Fibroblasts

# assign clusters
Fibroblasts_clus$cluster <- as.factor(clusters)

# read tissue area
tissue_area <- read.csv('TA_area.csv')


# once assigned, compared between long and short-term survivors.
Fibroblasts_clus_prop <- Fibroblasts_clus %>%
  group_by(metabric_id, cluster) %>%
  tally() %>%
  pivot_wider(names_from = cluster, values_from = n) %>%
  mutate_all(~replace_na(., 0)) %>%
  #rowwise() %>%
  #mutate(across(1:7, ~ . / sum(c_across(1:7)))) %>%
  merge(selectedDF_dup, by = 'metabric_id') %>%
  merge(tissue_area, by = 'metabric_id') %>%
  mutate(density_1 = `1` / TA,
         density_2 = `2` / TA,
         density_3 = `3` / TA,
         density_4 = `4` / TA)


test1 <- Fibroblasts_clus_prop[Fibroblasts_clus_prop$Survival == 'Survival high', 'density_4']
test2 <- Fibroblasts_clus_prop[Fibroblasts_clus_prop$Survival == 'Survival low', 'density_4'] 


wilcox.test(test1, test2)


# Cluster 2: CSF-S4

Fibroblasts_clus <- Fibroblasts_clus %>%
  group_by(cluster) %>%
  select(2:nFeature) %>%
  summarise_all(mean, na.rm = TRUE)

Fibroblasts_clus <- sapply(Fibroblasts_clus, as.numeric)


# scale clusters
Fibroblasts_clus[,2:nFeature] <- scale(Fibroblasts_clus[,2:nFeature])
Fibroblasts_clus[Fibroblasts_clus < -1] <- -1; Fibroblasts_clus[Fibroblasts_clus >  1] <- 1

# rename columns
rownames(Fibroblasts_clus) <- c('1', '2', '3', '4')

library(ComplexHeatmap)
library(colorRamp2)

color_mapping <- colorRamp2(c(-1, 0, 1), c("#0047ab", "white", "#fb3640"))

png(file="Figures/Fibroblast_cluster.png", width = 3, height = 4, units = "in", res = 300)
p <- Heatmap(Fibroblasts_clus[,2:nFeature], col = color_mapping,
             #cluster_rows = FALSE,
             #cluster_columns = FALSE,
             column_dend_height = unit(1, "cm"),
             show_heatmap_legend = FALSE,
             heatmap_legend_param = list(legend_direction = "horizontal", 
                                         legend_width = unit(4, "cm"), 
                                         at = seq(-1, 1, by = 1),
                                         labels_gp = gpar(fontsize = 15)),
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 16)
)
p
dev.off()



#-----------------------------------------------#
# Composition as bar plots ---------------------#
#-----------------------------------------------#
# I want to compare frequencies
allCell_counts <- sg_cell_USETHIS %>%
  group_by(metabric_id) %>%
  tally()

df_long <- Fibroblasts_clus %>%
  #merge(selectedDF_dup, by = 'metabric_id') %>%
  group_by(metabric_id, cluster) %>%
  tally() %>%
  pivot_wider(names_from = cluster, values_from = n) %>%
  mutate_all(~replace_na(., 0)) %>%
  merge(allCell_counts, by = 'metabric_id') %>%
  mutate_at(vars(2:5), ~ . / n) %>%
  merge(selectedDF_dup, by = 'metabric_id')
  #rowwise() %>%
  #mutate(across(1:4, ~ . / n))# %>%
  #pivot_longer(cols = c('1', '2', '3', '4'), names_to = "variable", values_to = "value")

test1 <- df_long[df_long$Survival == 'Survival low', 2]
test2 <- df_long[df_long$Survival == 'Survival high', 2]

wilcox.test(test1, test2)

boxplot(test1, test2)

#------------------#


df_long <- Fibroblasts_clus %>%
  merge(selectedDF_dup, by = 'metabric_id') %>%
  group_by(Survival, cluster) %>%
  tally() %>%
  pivot_wider(names_from = cluster, values_from = n) %>%
  mutate_all(~replace_na(., 0)) %>%
  rowwise() %>%
  mutate(across(1:4, ~ . / sum(c_across(1:4)))) %>%
  pivot_longer(cols = c('1', '2', '3', '4'), names_to = "variable", values_to = "value")



p <- ggplot(df_long, aes(x = Survival, y = value * 100, fill = variable)) +
  theme_classic() +
  geom_bar(stat = "identity", position = "stack") +
  labs(y = "Percentage (%)", x = "") +
  scale_fill_manual(values = c('1' = '#67a61f', '3' = '#e4191c',
                               '2' = '#1c79b5', '4' = '#ec8a48')) +
  scale_y_continuous(limits = c(0, 100), expand = expansion(mult = c(0, 0))) +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 20),
        axis.title = element_text(size = 22),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 1)) 
p
ggsave(p, file=paste0("Figures/Fibroblasts_Percentage_StackedBar.png"), width = 3.5, height = 5.5, units = "in", dpi = 300)




# Create contingency table
contingency_table <- table(df_long$Survival, df_long$cluster)

# Perform chi-square test
chi_sq_test <- chisq.test(contingency_table)
print(chi_sq_test)




mycomparison <- compare_means(`density_4` ~ Survival,  data = Fibroblasts_clus_prop)
mycomparison

p <- ggplot(Fibroblasts_clus_prop, aes(Survival, `density_4`, color = Survival)) +
  stat_boxplot( aes(Survival, `density_4`), 
                geom='errorbar', linetype=1, width= 0.2, size = 2)+  #whiskers
  geom_boxplot(size = 2, outlier.color = 'white') +
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0.3), 
             aes(fill=factor(Survival)), show.legend = F, size = 4, shape = 21, color = 'black') +
  #theme_foundation() +
  scale_color_manual(values = c('Survival high' = '#c2e699', 'Survival low' = '#31a353')) +
  scale_fill_manual(values = c('Survival high' = '#c2e699', 'Survival low' = '#31a353')) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 22),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 1)) +
  geom_signif(comparisons = list(c("Survival high", "Survival low")), 
              map_signif_level=TRUE, size = 1, textsize = 20, color = 'black', 
              y_position = max(Fibroblasts_clus_prop$`2`) + 10, annotations = mycomparison$p.signif
  ) +
  scale_y_continuous(limits = c(-20, 1200), expand = expansion(mult = c(0, 0))) +
  xlab('') +
  ylab(expression('Densities of Fibro-C3, mm'^-2))
p

#----------------#
# Note that cluster 2 here is cluster 3 in the figure!!!!!!
#########################################################

ggsave(p, file=paste0("Figures/Fibroblast_Cluster2_Compare.png"), width = 4, height = 7, units = "in", dpi = 300)





#-----------------------------------------------#
# Relative distance to CD8 T cells -------------#
#-----------------------------------------------#
nnDist_all <- data.frame(matrix(nrow = 0, ncol = 0))

for(mid in unique(sg_cell_USETHIS$metabric_id)){
  
  
  # Fibroblasts:
  pt_FibroC2 <- Fibroblasts_clus[Fibroblasts_clus$metabric_id == mid & Fibroblasts_clus$cluster == 2,]
  
  
  # CD8+ T cells
  pt_CD8T <- sg_cell_USETHIS[sg_cell_USETHIS$Phenotype == 'Tregs and Tex' & sg_cell_USETHIS$metabric_id == mid,]

  
  tryCatch({
    Fibro_CD8T <- RANN::nn2(data = pt_CD8T[, c('Location_Center_X', 'Location_Center_Y')], query = pt_FibroC2[, c('Location_Center_X', 'Location_Center_Y')],
                                 k = nrow(pt_CD8T), searchtype = 'radius', radius = 10, treetype = 'kd')  
    
    Fibro_Firbo <- RANN::nn2(data = pt_FibroC2[, c('Location_Center_X', 'Location_Center_Y')], query = pt_FibroC2[, c('Location_Center_X', 'Location_Center_Y')],
                                  k = nrow(pt_FibroC2), searchtype = 'radius', radius = 10, treetype = 'kd')  
    
    interactions1 <- sum(Fibro_CD8T$nn.idx != 0)
    interactions2 <- Fibro_Firbo$nn.idx %>%
      as.data.frame() %>% 
      pivot_longer(cols = -1, names_to = "Variable", values_to = "Value") %>%
      dplyr::filter(Value != 0) %>%
      dplyr::filter() %>%
      rowwise() %>%
      mutate(V1 = min(V1, Value), V2 = max(V1, Value)) %>%
      ungroup() %>%
      distinct(V1, Value) %>%
      dplyr::filter(Value != V1) %>%
      nrow()

    
    score <- interactions1 / interactions2
    
    
    
    # mixing score
    
    
    nnDist_all <- rbind.data.frame(nnDist_all, cbind.data.frame(mid, score))
  }, error = function(e) {
    # This function gets executed when there's an error
  })
  
  
}

colnames(nnDist_all)[1] <- 'metabric_id'
colnames(nnDist_all)[2] <- 'dist'

nrow(selectedDF_dup[selectedDF_dup$Survival == 'Survival low',])

test <- merge(nnDist_all, selectedDF_dup, by = 'metabric_id')
test1 <- test[test$Survival == 'Survival low', 'dist']
test2 <- test[test$Survival == 'Survival high', 'dist']
  
test_val <- wilcox.test(test1, test2)
test_val$p.value


Gcross_AUC <- function(pt_Fibro, pt_FoxP3, Fibro_type){
  #pt_Fibro <- pt_FibroC2
  #Fibro_type <- 'Fibro-C2'
  
  n1 <- nrow(pt_Fibro)
  n2 <- nrow(pt_FoxP3)
  
  if((n1 + n2) != 0){
    #pts_OI <- rbind(pt_Fibro[,c('Location_Center_X', 'Location_Center_Y', 'Phenotype')], pt_FoxP3[,c('Location_Center_X', 'Location_Center_Y', 'Phenotype')])
    
    # define the type
    #species <- factor(pts_OI$Phenotype)
   # multitype_ppp <- ppp(pts_OI$Location_Center_X, pts_OI$Location_Center_Y, marks = species, window)
    
    #ctype2 <- 'CD8+ T cells'
  #  Gihc <- data.frame(Gcross(multitype_ppp, i = Fibro_type, j = ctype2, r = seq(0, 100, 1)))
  #  Gihc <- Gihc[, c(1,2,4,5)]
    
    #Gihc <- Gihc[complete.cases(Gihc),]
    
    #auc <- trapz(Gihc$r, Gihc$rs)

    MH = 2 * n1 * n2 / (n1^2 +n2^2)
  }
  
  if(n1 == 0 & n2 == 0){
    #auc <- 0
    MH <- 0
  }
  
  return(MH)
}
  
  
  
  
  
for(ctype1 in unique(sg_cell_USETHIS$Phenotype)){
  for(ctype2 in unique(sg_cell_USETHIS$Phenotype)){
    
    
    nnDist_all <- data.frame(matrix(nrow = 0, ncol = 0))
    
    for(mid in unique(sg_cell_USETHIS$metabric_id)){
      
      allCells <- sg_cell_USETHIS[sg_cell_USETHIS$metabric_id == mid,]
      #ctype1 <- 'Fibro-C2'
      #ctype2 <- 'CD4+ T cells'
      # Fibroblasts:
      pt_FibroC1 <- Fibroblasts_clus[Fibroblasts_clus$metabric_id == mid & Fibroblasts_clus$cluster == 1,]
      pt_FibroC1$Phenotype <- 'Fibro-C1'
      
      pt_FibroC2 <- Fibroblasts_clus[Fibroblasts_clus$metabric_id == mid & Fibroblasts_clus$cluster == 2,]
      pt_FibroC2$Phenotype <- 'Fibro-C2'
      
      pt_FibroC3 <- Fibroblasts_clus[Fibroblasts_clus$metabric_id == mid & Fibroblasts_clus$cluster == 3,]
      pt_FibroC3$Phenotype <- 'Fibro-C3'
      
      pt_FibroC4 <- Fibroblasts_clus[Fibroblasts_clus$metabric_id == mid & Fibroblasts_clus$cluster == 4,]
      pt_FibroC4$Phenotype <- 'Fibro-C4'
      # CD8+ T cells
      #pt_CD8T <- sg_cell_USETHIS[sg_cell_USETHIS$Phenotype == ctype1 & sg_cell_USETHIS$metabric_id == mid,]
      pt_FoxP3 <- sg_cell_USETHIS[sg_cell_USETHIS$Phenotype == ctype2 & sg_cell_USETHIS$metabric_id == mid,]
      
      
      tryCatch({
        #Fibro_CD8T_dist <- RANN::nn2(data = pt_CD8T[, c('Location_Center_X', 'Location_Center_Y')], query = pt_FibroC2[, c('Location_Center_X', 'Location_Center_Y')],
        #                             k = nrow(pt_CD8T), searchtype = 'priority', treetype = 'kd')  
        
        #Fibro_Tumor_dist <- RANN::nn2(data = pt_FoxP3[, c('Location_Center_X', 'Location_Center_Y')], query = pt_FibroC2[, c('Location_Center_X', 'Location_Center_Y')],
        #                              k = nrow(pt_FoxP3), searchtype = 'priority', treetype = 'kd')  
        
        #nnDist <- Fibro_CD8T_dist$nn.dists %>%
        #  as.data.frame() %>%
        #  select(1) %>%
        #  as.matrix() %>%
        #  mean()
        
        #nnDist2 <- Fibro_Tumor_dist$nn.dists %>%
        #  as.data.frame() %>%
        #  select(1) %>%
        #  as.matrix() %>%
        #  mean()
        
        
        #score <- nnDist / (nnDist + nnDist2)
        
        
        
        # spat kcross
        library(spatstat)
        library(pracma)
        
        window <- owin(c(0, max(allCells$Location_Center_X)), c(0, max(allCells$Location_Center_Y)))
        
        # create multitype df
        auc <- Gcross_AUC(pt_FibroC1, pt_FoxP3, 'Fibro-C1')  
        auc2 <- Gcross_AUC(pt_FibroC2, pt_FoxP3, 'Fibro-C2')  
        auc3 <- Gcross_AUC(pt_FibroC3, pt_FoxP3, 'Fibro-C3')  
        auc4 <- Gcross_AUC(pt_FibroC4, pt_FoxP3, 'Fibro-C4')  
        
        
        
        # get the 'dense distance'
       # diff <- diff(Gihc$theo)/diff(Gihc$r)
        #common.distance1 <- Gihc$r[which.max(diff)]
        
        #plot(multitype_ppp)
        # calculate the area (positive - negative )  
        #Gihc$y <- Gihc$km - Gihc$theo
        
        
        
        #zero_crossings <- which(diff(sign(Gihc$y)) != 0)
        
        # Split data at zero crossings and compute AUC for each segment
        #auc_total <- 0
        #start_idx <- 1
        #for (end_idx in zero_crossings) {
        #  auc_total <- auc_total + trapz(Gihc$r[start_idx:end_idx], Gihc$y[start_idx:end_idx])
        #  start_idx <- end_idx + 1
        #}
        #auc_total <- auc_total + trapz(Gihc$r[start_idx:length(Gihc$r)], Gihc$y[start_idx:length(Gihc$r)])
        
        
        
        n1 <- nrow(pt_FibroC2)
        n2 <- nrow(pt_FoxP3)
        MH = 2 * n1 * n2 / (n1^2 +n2^2)
        
        
        #apprx_index1 <- trapz(Gihc$r, Gihc$rs) 
        #apprx_index2 <- trapz(Gihc2$r, Gihc2$rs) 
       # apprx_index3 <- trapz(Gihc3$r, Gihc3$rs) 
        #apprx_index4 <- trapz(Gihc4$r, Gihc4$rs) 
        #plot(Gihc$r, Gihc$rs)
        #points(Gihc$r, Gihc$theo, col ='red')
        
        nnDist_all <- rbind.data.frame(nnDist_all, cbind.data.frame(mid, auc, auc2, auc3, auc4))
        
      }, error = function(e) {
        # This function gets executed when there's an error
      })
      
      
    }
    
    colnames(nnDist_all)[1] <- 'metabric_id'
    #colnames(nnDist_all)[2] <- 'dist'
    
    
    
    tryCatch({
      test <- merge(nnDist_all, selectedDF_dup, by = 'metabric_id')
      
      test1 <- test[, 'auc']
      test2 <- test[, 'auc2']
      test3 <- test[, 'auc3']
      test4 <- test[, 'auc4']
      
      #mean(test1)
      #mean(test2)
      #boxplot(test1, test2)
      test_val1 <- wilcox.test(test2, test1)
      test_val3 <- wilcox.test(test2, test3)
      test_val4 <- wilcox.test(test2, test4)
      
      if(test_val1$p.value < 0.05 & test_val3$p.value < 0.05 & test_val4$p.value < 0.05){
        #print(test_val$p.value)
        print(ctype1)
        print(ctype2)
      }
    }, error = function(e){
      
    })
    
    
  }
  
  
}

mycomparison <- compare_means(dist ~ Survival,  data = test)
mycomparison

p <- ggplot(test, aes(x=Survival, y=dist)) +
  stat_boxplot(aes(Survival, dist, color = Survival),
               geom='errorbar', linetype=1, width= 0.2, size = 2)+  #whiskers
  geom_boxplot(size = 2, outlier.color = 'white', aes(color = Survival)) +
  geom_point(position=position_jitter(width=0), aes(fill=Survival), size = 4, shape = 21) + # Individual data points with some jitter for better visualization

  #stat_summary(fun=mean, geom="crossbar", width=0.05, color="black", fatten=1) +      # Horizontal line at mean
  #stat_summary(fun.data=mean_se, geom="errorbar", width=0.05, color="black") +     # Error bars
  theme_classic() +
  scale_color_manual(values = c('Survival high' = '#c2e699', 'Survival low' = '#31a353')) +
  scale_fill_manual(values = c('Survival high' = '#c2e699', 'Survival low' = '#31a353')) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.text = element_text(size = 22),
        axis.title = element_text(size = 22),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        axis.line = element_line(size = 1)) +
  geom_signif(comparisons = list(c("Survival high", "Survival low")), 
              map_signif_level=TRUE, size = 1, textsize = 15, color = 'black', 
              y_position = max(test$dist) + 0.2, annotations = mycomparison$p.signif
  ) +
  scale_y_continuous(limits = c(-1, 6), expand = expansion(mult = c(0, 0))) +
  xlab('') +
  ylab('Area between G(r) curves') 

p
ggsave(p, file=paste0("Figures/ABC_Gcross.png"), width = 3, height = 5, units = "in", dpi = 300)


#---------------------------------------------------#
#--------------- Cell-cell interaction -------------#
# --------------------------------------------------#

library(ClusterR)

Fibroblasts <- sg_cell_USETHIS %>%
  dplyr::filter(Phenotype == 'Fibroblasts') %>%
  select(metabric_id, SMA, FSP1, PDGFRB, `Caveolin-1`, Location_Center_X, Location_Center_Y)

nFeature <- 5

set.seed(2)
km <- MiniBatchKmeans(as.matrix(Fibroblasts[,2:nFeature]), clusters = 4, batch_size = 100)
clusters <- predict_MBatchKMeans(as.matrix(Fibroblasts[,2:nFeature]), km$centroids)



Fibroblasts_clus <- Fibroblasts

# assign clusters
Fibroblasts_clus$cluster <- as.factor(clusters)
sg_cell_USETHIS$Phenotype <- as.factor(sg_cell_USETHIS$Phenotype)


nnFibro_all <- data.frame(matrix(nrow = 0, ncol = 0))
for(mid in unique(Fibroblasts_clus$metabric_id)){
  
  
  # get the fibroblasts of that mid
  fibro_mid <- Fibroblasts_clus %>%
    dplyr::filter(metabric_id == mid)
  
  # get other cells
  other_mid <- sg_cell_USETHIS %>%
    dplyr::filter(metabric_id == mid) %>%
    dplyr::filter(Phenotype != 'Fibroblasts')
  
  
  # get the nearest neighbors of fibroblasts
  library(purrr)
  nnFibro <- nn2(data = other_mid[, c('Location_Center_X', 'Location_Center_Y')], 
                 query = fibro_mid[, c('Location_Center_X', 'Location_Center_Y')],
                 k = nrow(other_mid),
                 treetype = 'kd',
                 searchtype = 'radius',
                 radius = 30) %>%
    pluck('nn.idx') %>%
    data.frame() %>%
    mutate(Cluster = fibro_mid$cluster) %>%
    gather(key = "variable", value = "value", -Cluster) %>%
    dplyr::filter(value != 0) %>%
    mutate(Phenotype = other_mid$Phenotype[value]) %>%
    group_by(Cluster, Phenotype, .drop = FALSE) %>%
    tally() %>%
    dplyr::filter(Phenotype != 'Fibroblasts') %>%
    pivot_wider(names_from = Phenotype, values_from = n) %>%
    data.frame() %>%
    mutate(sum = rowSums(.[,-1])) %>%
    dplyr::filter(sum != 0) %>%
    mutate_at(vars(-Cluster, -sum), ~ ./sum) %>%
    select(-sum)
  
    

  nnFibro_all <- rbind.data.frame(nnFibro_all, cbind.data.frame(mid, nnFibro))

  
}



colnames(nnFibro_all)[1] <- 'metabric_id'

test1 <- nnFibro_all[nnFibro_all$Cluster == 1, 'Tregs.and.Tex']
test2 <- nnFibro_all[nnFibro_all$Cluster == 2, 'Tregs.and.Tex']
test3 <- nnFibro_all[nnFibro_all$Cluster == 3, 'Tregs.and.Tex']
test4 <- nnFibro_all[nnFibro_all$Cluster == 4, 'Tregs.and.Tex']

mean(test4)

wilcox.test(test4, test2)


nnFibro_all_forHierClus <- nnFibro_all

nnFibro_all_forHierClus <- nnFibro_all_forHierClus %>%
  group_by(Cluster) %>%
  select(3:16) %>%
  summarise_all(mean, na.rm = TRUE)



nnFibro_all_forHierClus[,2:15] <- scale(nnFibro_all_forHierClus[,2:15])
nnFibro_all_forHierClus[nnFibro_all_forHierClus < -1] <- -1; nnFibro_all_forHierClus[nnFibro_all_forHierClus >  1] <- 1

# rename columns
nnFibro_all_forHierClus <- data.frame(nnFibro_all_forHierClus)
rownames(nnFibro_all_forHierClus) <- c('1', '2', '3', '4')

library(ComplexHeatmap)
library(colorRamp2)

color_mapping <- colorRamp2(c(-1, 0, 1), c("#0047ab", "white", "#fb3640"))

#png(file="Figures/Fibroblast_cluster.png", width = 3, height = 4, units = "in", res = 300)
p <- Heatmap(nnFibro_all_forHierClus[,2:15], col = color_mapping,
             #cluster_rows = FALSE,
             #cluster_columns = FALSE,
             show_row_names = TRUE,
             column_dend_height = unit(1, "cm"),
             show_heatmap_legend = FALSE,
             heatmap_legend_param = list(legend_direction = "horizontal", 
                                         legend_width = unit(4, "cm"), 
                                         at = seq(-1, 1, by = 1),
                                         labels_gp = gpar(fontsize = 15)),
             row_names_gp = gpar(fontsize = 12),
             column_names_gp = gpar(fontsize = 16)
)
p
#dev.off()

