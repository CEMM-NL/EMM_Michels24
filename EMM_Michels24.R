## R-code
# Script for manuscript michels et al., 2024 neutrophil proteomics
### load packages
pacman::p_load(pacman, tidyverse, readr, dplyr, tidyr, rio, ggplot2, naniar,broom,
               nlme, lcmm, tableone, lattice, reshape2, data.table, scales, plyr, Hmisc,
               tableone, kmed, survminer, skimr, lattice, reshape2, grDevices, ggrepel,
               org.Hs.eg.db, circlize, limma, ReactomePA, clusterProfiler, enrichplot,
               factoextra, ggthemes, devtools, ggbiplot, WGCNA,ComplexHeatmap,corrplot,
               ggvenn, ReactomeGSA)

# Human
org <- "Homo sapiens"  
p_load("org.Hs.eg.db")  

######################################################
## PREPROCESS,  IMPUTING, logFC between conditions  ##
######################################################

### Load data, filtering for proteotypic and with minimal 2 precursors
{
  # prefiltering of precursors
  dia_proteotypic_2pr <- read.delim('PMN_CAP_AMC_output.pr_matrix.tsv') %>%
    dplyr::select(Protein.Group, Genes, Proteotypic, Precursor.Id, contains('raw')) %>%
    gather('Sample', 'Intensity', contains('raw')) %>%
    filter(!is.na(Intensity)) %>%
    filter(Proteotypic == T) %>%
    group_by(Sample, Protein.Group) %>%
    distinct(Precursor.Id, .keep_all = T) %>%
    summarise(count = n()) %>%
    filter(count > 1)  %>%
    dplyr::select(-count)
  
  # filtering of proteins
  prot_filtered <- read.delim('PMN_CAP_AMC_output.pg_matrix.tsv') %>%
    gather('Sample', 'Intensity', -Protein.Group, -Protein.Ids, -Protein.Names, -Genes, -First.Protein.Description) %>%
    filter(!is.na(Intensity)) %>%
    inner_join(dia_proteotypic_2pr, c('Protein.Group', 'Sample')) %>%
    spread('Sample', 'Intensity') %>% 
    rownames_to_column('UID') %>% 
    mutate(UID = paste0('UID', UID)) %>% 
    mutate(across(contains('raw'),log2)) %>% 
    set_names(~gsub('(.*)N_|(.*)l_', '', .x))
}  

### load clinical data and condition
## load condition
pheno <- read.csv2('sample_info.csv', stringsAsFactors = F) %>% 
  filter(Sample %in% gsub('.raw','',colnames(prot_filtered))) %>% 
  mutate(ct = factor(Type,levels = unique(Type)[2:1])) %>% 
  mutate(Type = factor(Type, levels = unique(Type)[2:1])) 

### load clinical information
clin_dat <- read.csv2('../AMC_clin_dat/sanquin_patient_adj_20232414.csv') %>% 
  mutate(Record.Id = as.character(Record.Id))

### extend pheno with clin data
pheno <- pheno %>% 
  left_join(clin_dat, by = c('Sample'='Record.Id'))

###valid values, >75% in at least one (IALO)
vvs <- prot_filtered %>% 
  dplyr::select(contains('raw'),UID) %>% 
  column_to_rownames('UID') %>% 
  mutate(across(.cols = everything(), ~ !is.na(.))) %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column('Sample') %>% 
  mutate(Sample = str_remove(Sample,'\\.raw')) %>% 
  left_join(pheno %>% dplyr::select(Sample,ct), c('Sample')) %>% 
  column_to_rownames('Sample') %>% 
  filter(!is.na(ct)) %>% 
  group_by(ct) %>% 
  summarise_all(sum) %>% 
  column_to_rownames('ct') %>% 
  mutate(across(.cols = everything(), ~ . / as.numeric(table(list(pheno$ct))))*100) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(No_pass = rowSums(across(.cols = everything(), ~ . >= 75 ))) %>% # valid value threshold == 75%
  rownames_to_column('UID')

  # filter data
  prot_valid <- prot_filtered %>% 
                  subset(UID %in% vvs$UID[vvs$No_pass >=1])
###


###

### impute missing values
# define parameters
mann_downshift <- 1.8
mann_width <- 0.3

# create imputed values
set.seed(2112)
imputed_values <- prot_valid %>% 
  dplyr::select(paste0(pheno$Sample,'.raw')) %>% 
  mutate(across(.cols = everything(), ~ rnorm(length(.), (mean(na.omit(.) -  mann_downshift*sd(na.omit(.)))), (sd(na.omit(.))*mann_width))))

# replace NA's in DF with imputed
prot_imp <- prot_valid %>% 
  dplyr::select(paste0(pheno$Sample,'.raw')) %>% 
  coalesce(., imputed_values) %>% 
  bind_cols(Genes = prot_valid$Genes, UID = prot_valid$UID)

# check
prot_imp %>% 
  dplyr::select(paste0(pheno$Sample,'.raw')) %>% 
  boxplot()
###

### Calculate logFC
# make a design
design <- model.matrix(~0+as.factor(as.numeric(pheno$ct)))
colnames(design) <- levels(pheno$ct)

# fit in limma
fit <- prot_imp %>% 
  dplyr::select(paste0(pheno$Sample,'.raw')) %>% 
  lmFit(.,design)

contrast_matrix <- makeContrasts(CAP-HC, levels =design)

# fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
# ebayes
ebayes <- eBayes(fit2)

#all results
results = decideTests(ebayes, p.value = 0.05, adjust.method = 'BH', lfc = 1)
summary(results)

# get significance table
results_df <- ebayes %>% 
  topTable(number = Inf, coef = 1, sort.by = 'none') %>% 
  bind_cols(prot_valid)

sig_df <- ebayes %>% 
  topTable(number = Inf, coef = 1, sort.by = 'none') %>% 
  bind_cols(prot_valid) %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC)  >= 1)


# Make table for export
export_df <- prot_filtered %>%
  set_names(~gsub('\\.raw', '_intensity', .x)) %>% 
  left_join(prot_imp %>% 
              dplyr::select(UID,contains('.raw')) %>% 
              set_names(~gsub('\\.raw', '_imputed_intensity', .x)) %>% 
              add_column(Quantified = 'True', .before = 1) %>% 
              bind_cols(ebayes %>% 
                          topTable(number = Inf, coef = 1, sort.by = 'none'))
            ,
            c('UID'))

write.csv2(export_df, file ='proteomics_AMC_CAP_PMN_export2.csv"')


############################################
### Clean and format data for analyses #####
############################################

### Load clinical data 
Clinical_data <- import("all_clinical_data.xlsx")

## load PMN protein data
PMN <- read.csv2("proteomics_AMC_CAP_PMN_export2.csv")

#############################################
##       Figure 1a quantified proteins     ##
##############################################

## patients columns, measured proteins rows
## get columns of raw data => remove surfix => transpose (1 patient is 1 row)
raw_data <- PMN[ , 8:90]
colnames(raw_data) <- gsub("^X|_intensity$", "", colnames(raw_data))

## filter on CAP/control
CAP <- filter(Clinical_data, Clinical_data$patient_volunteer == "Pneumonia patient")
healthy <- filter(Clinical_data, Clinical_data$patient_volunteer == "Healthy volunteer")

##
PMN_CAP <- raw_data[, colnames(raw_data) %in% CAP$EB_id]
PMN_healthy <- raw_data[, colnames(raw_data) %in% healthy$EB_id, ]

## count missing 
PMN_CAP$percentage_missing <- rowSums(is.na(PMN_CAP))/length(PMN_CAP) * 100
sum(PMN_CAP$percentage_missing <25) ## 3433

PMN_healthy$percentage_missing <- rowSums(is.na(PMN_healthy))/length(PMN_healthy) * 100
sum(PMN_healthy$percentage_missing <25) ## 3168

##
categories <- c("CAP", "Healthy")
quantified_proteins <- c(3433, 3168)

# Create a data frame
data <- data.frame(Category = categories, Quantified_Proteins = quantified_proteins)
data$Category <- factor(data$Category, levels = c("Healthy", "CAP"))

# Plot
ggplot(data, aes(x = Quantified_Proteins, y = Category)) +
  geom_bar(stat = "identity", fill =  c("#7e53a1", "#1b9e77")) +
  geom_text(aes(label = Quantified_Proteins), hjust = -0.2, color = "white") +
  labs(x = "Quantified Proteins", y = "") +
  coord_cartesian(xlim = c(0, 4000)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line())

##################################################
##        Figure 1b PCA analysis             #####
##################################################

## remove unquantified proteins
PMN <- filter(PMN, PMN$Quantified == "True")

## manually change duplicated gene names and empty gene names
x <- PMN[duplicated(PMN$Genes), ] ## TMPO duplicated
PMN$Genes <- ifelse(PMN$Protein.Ids == "P42166", "TMPOa", PMN$Genes)
PMN$Genes[PMN$Genes == ""] <- "FLJ45252"

## remove raw values and other for PCA abundant columns
PMN_clean <- PMN[ ,c(6, 92:174)]

## modify rownames and columns names and transpose
## use gene names
row.names(PMN_clean) <- PMN_clean$Genes
PMN_clean$Genes <- NULL
t_PMN_clean <- t(PMN_clean) %>% as.data.frame()
new_row_names <- gsub("^X|_imputed_intensity$", "", rownames(t_PMN_clean))
rownames(t_PMN_clean) <- new_row_names

## save intensities and gene ID
PMN_clean_save <- PMN_clean
write.csv(t(PMN_clean), "PMN_clean.csv")

## dataset with rows as patients and columns as proteins (t_PMN_clean)
## get condition
t_PMN_clean$CAP <- as.factor(ifelse(row.names(t_PMN_clean) %in% CAP$EB_id, 1, 0))

## make centered and scaled PCA
Groups <- t_PMN_clean[,"CAP"]
e.pca <- prcomp(t_PMN_clean[ ,1:3482], center = T, scale = T) 
summary(e.pca)
str(e.pca)

##
colors <- c("#1b9e77","#7e53a1")
e.plot1 <- ggbiplot(e.pca, ellipse = FALSE, obs.scale = 1, var.scale = 1, var.axes = FALSE,
                    group = Groups, circle = FALSE, varname.size = 0, alpha = 0,
                    varname.adjust = c(1), ellipse.prob = 0.50) +
  scale_color_manual(name = "Groups", values = colors) +
  geom_point(aes(colour = Groups), size = 0.8, alpha = 1) +
  stat_ellipse(aes(colour = Groups), size = 1, type = "norm", level = 0.50, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "right") +
  ggtitle("Principal component analysis")
e.plot1

## statistically test components
pcdat <- as.data.frame(e.pca[["x"]])
pcdat$group <- Groups
summary(aov(pcdat$PC1 ~ pcdat$group))
summary(aov(pcdat$PC2 ~ pcdat$group))

##
PMN_clean <- t_PMN_clean

##################################################
##        Figure 1c Volcano                  #####
##################################################

## PMN from previous section contains limma results (logFC and t-stastic)
#color for volcanoplot
PMN$col <- "A"
PMN$col <- ifelse(PMN$logFC < 0 & PMN$adj.P.Val <=0.05, "B", PMN$col)
PMN$col <- ifelse(PMN$logFC > 0 & PMN$adj.P.Val <= 0.05, "C", PMN$col)
PMN$col <- ifelse(PMN$logFC > 0 &  PMN$logFC <1 & PMN$adj.P.Val<=0.05, "D", PMN$col)
PMN$col <- ifelse(PMN$logFC < 0 &  PMN$logFC > -1 & PMN$adj.P.Val<=0.05, "E", PMN$col)
table(PMN$col)

#label top 10 each side
top10_high <- filter(PMN, PMN$adj.P.Val <= 0.05 & PMN$logFC > 1)
top10_high <- top10_high[order(top10_high$adj.P.Val), ]
top10_high <- top10_high[1:10, ]

top10_low <- filter(PMN, PMN$adj.P.Val <= 0.05 & PMN$logFC < -1)
top10_low <- top10_low[order(top10_low$adj.P.Val), ]
top10_low <- top10_low[1:10, ]

PMN <- PMN[order(PMN$adj.P.Val), ]
PMN$label <- NA
PMN$label <- ifelse(PMN$Genes %in% top10_high$Genes, PMN$Genes, 
                    ifelse(PMN$Genes %in% top10_low$Genes, PMN$Genes, NA))
primary_results <- PMN

# Modify plot aesthetics
p <- ggplot(PMN, aes(x = logFC, y = -log10(adj.P.Val), label = label)) +
  geom_point(aes(col = col), size = 0.7, alpha = 1) +
  scale_color_manual(values = c(A = "gray", B= "Blue", C = "red", D = "#FFB2B2", E = "#B2B2FF")) +
  scale_shape_manual(values = 1) +
  geom_text_repel(size = 4, nudge_x = 0.1, nudge_y = 0.1, max.overlaps = Inf, segment.color = "black") +
  xlab(expression(paste("logFC", " Fold Change"))) +
  ylab(expression(paste("-log"[10], " (BH adjusted P)"))) +
  geom_hline(yintercept = 1.30103, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
  geom_vline(xintercept = -1, linetype = "dotted", color = "black") +
  ggtitle("Community-Acquired Pneumonia vs. Healthy Control") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

##
# Create pie chart
table(PMN$col)
group_counts <- c( 943,  77,25, 502, 1935)
group_colors <- c("#FFB2B2","red", "blue","#B2B2FF", "gray")
pie(group_counts, labels = NA, col = group_colors)

##################################
##### fig 1d pathway analysis  ##
###################################

## Make genelist using t-statistis
set.seed(5)
genelist <- PMN$t
names(genelist) <- PMN$Genes
genelist <- genelist[order(-genelist)]
head(genelist)

# Perform GSEA with the GO databases, first BP
result_BP_CAP_healthy <- gseGO(geneList = genelist,
                               OrgDb = org.Hs.eg.db,
                               keyType = "SYMBOL",
                               ont = "BP",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               minGSSize = 10,
                               maxGSSize = 500,
                               eps = 0,
                               nPermSimple = 10000)

## make dataframe
result_BP_CAP_healthy_df <- result_BP_CAP_healthy %>% as.data.frame()
openxlsx::write.xlsx(result_BP_CAP_healthy_df, "result_BP_CAP_healthy_df.xlsx")

## make enrichplot tree of pathways => cluster similar pathways using ward.D2
BP <- pairwise_termsim(result_BP_CAP_healthy)
BP_tree <- treeplot(BP, nCluster = 3, hclust_method = "ward.D2",offset = rel(4), color = "NES",cex_category = 0.9,   showCategory = 20, group_color = c( "#00C957","#FFD700", "#008080")) +
  scale_color_gradient(low = "blue", high = "red", limits = c(-3, 3))
BP_tree

# Perform GSEA with the GO databases, CC
result_CC_CAP_healthy <- gseGO(geneList = genelist,
                               OrgDb = org.Hs.eg.db,
                               keyType = "SYMBOL",
                               ont = "CC",
                               pvalueCutoff = 0.05,
                               pAdjustMethod = "BH",
                               minGSSize = 10,
                               maxGSSize = 500,
                               eps = 0,
                               nPermSimple = 10000)

result_CC_CAP_healthy_df <- result_CC_CAP_healthy %>% as.data.frame()
openxlsx::write.xlsx(result_CC_CAP_healthy_df, "result_CC_CAP_healthy_df.xlsx")

## CC
CC <- pairwise_termsim(result_CC_CAP_healthy)
CC_tree <- treeplot(CC, nCluster = 3, hclust_method = "ward.D2", color = "NES",offset = rel(4),cex_category = 0.9,
                    showCategory = 20, group_color = c( "#00C957","#FFD700", "#008080")) +
  scale_color_gradient(low = "blue", high = "red", limits = c(-3, 3))
CC_tree

#################################################################
## Figure 1e heatmap per most substantially alter pathways #####
################################################################

## identify pathway with highest absolute NES per pathway cluster using the BP database
## cellular macromolecule biosynthetic process, protein folding, actin filament-based process
result_BP_CAP_healthy_df <- result_BP_CAP_healthy_df[order(-abs(result_BP_CAP_healthy_df$NES)), ]
top10_BP <- result_BP_CAP_healthy[1:10, ]
top10_BP$Description

## cellular macromolecule biosynthetic process - acquired involved proteins
## select top 20
top_proteins_biosynthetic <- filter(result_BP_CAP_healthy_df, result_BP_CAP_healthy_df$Description == "cellular macromolecule biosynthetic process")
top_proteins_biosynthetic <- top_proteins_biosynthetic$core_enrichment
top_proteins_biosynthetic <-  unlist(strsplit(top_proteins_biosynthetic , "/"))
heatmap_biosynthetic <- PMN_clean[ , colnames(PMN_clean) %in% top_proteins_biosynthetic[ 1:20]]

## format for heatmap check order
data_matrix1 <- t(scale(as.matrix(heatmap_biosynthetic)))
colnames(data_matrix1) == row.names(PMN_clean)

# Create a color scale
row_ha <- rowAnnotation(cap = PMN_clean$CAP, 
                        col = list(cap = c("1" = "#7E53A1", "0" = "#1B9E77")))

hist_major <- Heatmap(as.matrix(t(data_matrix1)), 
                      split = PMN_clean$CAP, 
                      column_names_gp = gpar(fontsize = 10),
                      name = "Scaled\nvalue",
                      show_column_dend = FALSE,
                      show_row_names = T,
                      show_column_names = T,
                      clustering_distance_rows = "euclidean",
                      clustering_distance_columns = "euclidean",
                      clustering_method_rows = "ward.D2",
                      clustering_method_columns = "ward.D2",
                      use_raster= TRUE,
                      raster_resize_mat = max,
                      border_gp = gpar(col = "black", lty = 2)) + row_ha
hist_major


## repeat for second pathway => protein folding
top_proteins_folding <- filter(result_BP_CAP_healthy_df, result_BP_CAP_healthy_df$Description == "protein folding")
top_proteins_folding <- top_proteins_folding$core_enrichment
top_proteins_folding <-  unlist(strsplit(top_proteins_folding  , "/"))
heatmap_folding <- PMN_clean[ , colnames(PMN_clean) %in% top_proteins_folding[1:20]]

## make heatmap
data_matrix2 <- t(scale(as.matrix(heatmap_folding)))
colnames(data_matrix2) == row.names(PMN_clean)

row_ha <- rowAnnotation(cap = PMN_clean$CAP, 
                        col = list(cap = c("1" = "#7E53A1", "0" = "#1B9E77")))

hist_major <- Heatmap(as.matrix(t(data_matrix2)), 
                      split = PMN_clean$CAP, 
                      column_names_gp = gpar(fontsize = 10),
                      name = "Scaled\nvalue",
                      show_column_dend = FALSE,
                      show_row_names = T,
                      show_column_names = T,
                      clustering_distance_rows = "euclidean",
                      clustering_distance_columns = "euclidean",
                      clustering_method_rows = "ward.D2",
                      clustering_method_columns = "ward.D2",
                      use_raster= TRUE,
                      raster_resize_mat = max,
                      border_gp = gpar(col = "black", lty = 2)) + row_ha
hist_major


## repeat for last pathway => actin filament-based process
top_proteins_actin <- filter(result_BP_CAP_healthy_df, result_BP_CAP_healthy_df$Description == "actin filament-based process")
top_proteins_actin  <- top_proteins_actin$core_enrichment
top_proteins_actin <-   unlist(strsplit(top_proteins_actin   , "/"))
heatmap_actin <- PMN_clean[ , colnames(PMN_clean) %in% top_proteins_actin[1:20]]

## make heatmap
data_matrix3 <- t(scale(as.matrix(heatmap_actin)))
colnames(data_matrix3) == row.names(PMN_clean)

row_ha <- rowAnnotation(cap = PMN_clean$CAP, 
                        col = list(cap = c("1" = "#7E53A1", "0" = "#1B9E77")))

hist_major <- Heatmap(as.matrix(t(data_matrix3)), 
                      split = PMN_clean$CAP, 
                      column_names_gp = gpar(fontsize = 10),
                      name = "Scaled\nvalue",
                      show_column_dend = FALSE,
                      show_row_names = T,
                      show_column_names = T,
                      clustering_distance_rows = "euclidean",
                      clustering_distance_columns = "euclidean",
                      clustering_method_rows = "ward.D2",
                      clustering_method_columns = "ward.D2",
                      use_raster= TRUE,
                      raster_resize_mat = max,
                      border_gp = gpar(col = "black", lty = 2)) + row_ha
hist_major

##################################
#### FIGURE 2a TCS association ###
#################################

## merge TCS and CAP with proteins
CAP_only <- filter(PMN_clean, PMN_clean$CAP == 1)
TTCS <- select(Clinical_data, "EB_id", "clin_stab_hosp_old")
CAP_only <- merge(CAP_only, TTCS, by.y = "EB_id", by.x = 0)
CAP_only$CAP <- NULL

## Spearman loop
correlation_spear <- function(x){
  r <- cor.test(CAP_only[, x], CAP_only[ , "clin_stab_hosp_old"], method = "spearman")$estimate
  p <- cor.test(CAP_only[, x], CAP_only[ , "clin_stab_hosp_old"], method = "spearman")$p.value
  vec <- c(r, p)
  names(vec) <- c("r", "p")
  return(vec)
}


## apply and format
proteins <- colnames(CAP_only[,c(2:3482)])
cor.results_spear <- sapply(proteins, correlation_spear)
cor.results_spear <- as.data.frame(t(cor.results_spear))
cor.results_spear$marker <- row.names(cor.results_spear)
cor.results_spear$abs_r<- abs(cor.results_spear$r)
cor.results_spear$adj.p <- p.adjust(cor.results_spear$p, method = "BH")
cor.results_spear <- cor.results_spear[order(-cor.results_spear$abs_r), ]

## retrieve logFC and adj.P.Val from CAP vs control comparison for visualisation and interpretation
logfc <- select(PMN, "Genes", "logFC", "adj.P.Val")
cor.results_spear <- merge(cor.results_spear, logfc, by.x = "marker", by.y = "Genes")

#color for volcanoplot
cor.results_spear$col <- "A"
cor.results_spear$col <- ifelse(cor.results_spear$r >0.2, "D", cor.results_spear$col)
cor.results_spear$col <- ifelse(cor.results_spear$r < -0.2, "E", cor.results_spear$col)
cor.results_spear$col <- ifelse(cor.results_spear$r >0.4, "B", cor.results_spear$col)
cor.results_spear$col <- ifelse(cor.results_spear$r < -0.4, "C", cor.results_spear$col)
table(cor.results_spear$col)

## label
cor.results_spear$label <- NA
cor.results_spear$label <- ifelse(cor.results_spear$col == "B" | cor.results_spear$col == "C", cor.results_spear$marker, NA)

##
p2 <- ggplot(cor.results_spear, aes(x=r, y= -log10(adj.P.Val), label=label )) +
  geom_point(aes(col=col), size = 1) +
  scale_color_manual(values = c("A" = "grey", "B" = "darkred", "C" = "blue", "D" = "#FFB3B3", "E" = "#ADD8E6", "F" = "darkred")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_text_repel(max.overlaps = 9999) +
  xlab("Correlation with TTCS in CAP") +
  ylab("-log10(Adjusted P-value) comparing CAP vs controls") +
  scale_x_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +
  coord_cartesian(xlim = c(-0.6, 0.6)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted")
p2

## established proteins from literature
primary_granules <- c("MPO", "ELANE", "PRTN3", "AZU1", "CTSG")
secondary_granules <- c("PTX3", "ORM1", "CRISP3", "LCN2", "LTF")
tertiary_granules <- c("MMP8", "MMP9")
secretory_granules <- c("CR1", "FCGR3B")
lysosomal <- c("LAMP1", "LAMP2", "GUSB", "CTSC", "CTSD")
antimicrobial <- c("DEFA1", "DEFA3", "DEFA4", "DEF6", "CAMP", 
                   "BPI") 
ROS <- c("CYBC1", "CYBA", "CYBB", "NCF2", "NCF1", "NCF4", "RAC1", 
         "RAC2", "SOD1", "S100A8", "S100A9")
apoptosis <- c("BAX", "BID", "PCNA", "BCL10")
other <- c("ITGAM", "ITGB2", "ARG1", "CXCR2", "LAIR1", "CD200R1")
all <- c(primary_granules, secondary_granules, tertiary_granules, secretory_granules, antimicrobial, lysosomal,  ROS, apoptosis, other)


##
target <- cor.results_spear[cor.results_spear$marker %in% all, ]
target$label <- ifelse(!target$col == "A", target$marker, NA)

##
p2 <- ggplot(target, aes(x=r, y= -log10(adj.P.Val), label=label )) +
  geom_point(aes(col=col), size = 1) +
  scale_color_manual(values = c("A" = "grey", "B" = "darkred", "C" = "blue", "D" = "#FFB3B3", "E" = "#ADD8E6", "F" = "darkred")) +
  theme_bw() +
  theme(legend.position = "none",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_text_repel(max.overlaps = 9999) +
  xlab("Correlation with TTCS in CAP") +
  ylab("-log10(Adjusted P-value) comparing CAP vs controls") +
  scale_x_continuous(breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6)) +
  coord_cartesian(xlim = c(-0.6, 0.6)) +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted")
p2

##################################
#### FIGURE 2b TCS pathways ###
#################################

## Make genelist based on correlation to TCS
set.seed(10)
genelist <- cor.results_spear$r
names(genelist) <- cor.results_spear$marker
genelist <- genelist[order(-genelist)]
head(genelist)

# Perform GSEA with the GO database
result_BP_TTCS <- gseGO(geneList = genelist,
                        OrgDb = org.Hs.eg.db,
                        keyType = "SYMBOL",
                        ont = "BP",
                        pvalueCutoff = 0.05,
                        pAdjustMethod = "BH",
                        minGSSize = 10,
                        maxGSSize = 500,
                        eps = 0,
                        nPermSimple = 10000)

result_BP_TTCS_df  <- result_BP_TTCS %>% as.data.frame() 

## make enrichplot tree of pathways => cluster similar pathways using ward.D2
## simplify results for better overview of pathways related to TCS
result_BP_TTCS_simplified <- simplify(result_BP_TTCS, cutoff = 0.7, by = "p.adjust", select_fun = min, measure = "Wang")
result_BP_TTCS_simplified_df <- result_BP_TTCS_simplified %>% as.data.frame()

BP <- pairwise_termsim(result_BP_TTCS_simplified)
BP_tree <- treeplot(BP, nCluster = 5, hclust_method = "ward.D2",offset = rel(4), color = "NES",legend_n = 1000, showCategory = 20) +
  scale_color_gradient(low = "blue", high = "red", limits = c(-3, 3))
BP_tree

######################
## Figure 2c ssPSEA ##
#####################

## REACTOME plugin single sample protein set enrichemnt available at https://reactome.org/gsa/home
ssGSEA <- import("ssPSEA_CAP.rds")

# add pathway names 
all_pathways <- pathways(ssGSEA) %>% as.data.frame()
names <- row.names(all_pathways)
names2 <- select(all_pathways, "Name")
names <- cbind(names, names2)

all_pathways$Name <- NULL
all_pathways <- t(all_pathways) %>% as.data.frame()
row.names(all_pathways) <- gsub("X", "", row.names(all_pathways))
row.names(all_pathways) <- gsub(".PMN_clean", "", row.names(all_pathways))
row.names(all_pathways) <- gsub(".PROTEOMICS_INT_1", "", row.names(all_pathways))

## merge with TCS per CAP patient
TTCS_long <- select(Clinical_data, "EB_id", "clin_stab_hosp_old")
TTCS_long <- TTCS_long[!is.na(TTCS_long$clin_stab_hosp_old), ]
all_pathways <- merge(all_pathways, TTCS_long, by.x = 0, by.y = "EB_id")

## correlate pathways scores with TCS
correlation_spear <- function(x){
  r <- cor.test(all_pathways[, x], all_pathways[ , "clin_stab_hosp_old"], method = "spearman", exact = T)$estimate
  p <- cor.test(all_pathways[, x], all_pathways[ , "clin_stab_hosp_old"], method = "spearman", exact = T)$p.value
  vec <- c(r, p)
  names(vec) <- c("r", "p")
  return(vec)
}

## apply and format
proteins <- colnames(all_pathways[,c(2:1637)])
cor.results_spear <- sapply(proteins, correlation_spear)
cor.results_spear <- as.data.frame(t(cor.results_spear))
cor.results_spear$marker <- row.names(cor.results_spear)

## merge pathway names
cor.results_spear <- cor.results_spear[order(-cor.results_spear$r), ]
cor.results_spear <- merge(cor.results_spear, names, by.x = "marker", by.y = "names")

### open reactome structure
RPR <- import("reactomepathways_19sep2024.xlsx")
names(RPR) <- c("parent", "child")

## Filter on immune pathways and it's children
## Innate Immune System (R-HSA-168249,  Metabolism (R-HSA-1430728), Programmed Cell Death (R-HSA-5357801), 
## Death Receptor Signalling (R-HSA-73887), and Cellular Response to Stress (R-HSA-2262752).
parents <-  c("R-HSA-168249", "R-HSA-1430728", "R-HSA-5357801", "R-HSA-73887", "R-HSA-2262752")
level2 <- subset(RPR, parent %in% parents)
level2$parent <- NULL
immune_pathways <- c(level2$child, parents)

#color and label for volcanoplot
cor.results_spear2 <- cor.results_spear[cor.results_spear$marker %in% immune_pathways, ]
cor.results_spear2$col <- "A"
cor.results_spear2$col <- ifelse(cor.results_spear2$r < -0.3 , "B", cor.results_spear2$col)
cor.results_spear2$col <- ifelse(cor.results_spear2$r > 0.3 , "C", cor.results_spear2$col)
cor.results_spear2$label <- NA
cor.results_spear2$label <- ifelse(!cor.results_spear2$col == "A", cor.results_spear2$Name, NA)

# ggplot volcano (note, unadjusted p)
p2 <- ggplot(cor.results_spear2, aes(x=r, y= -log10(p), label=label )) +
  geom_point(aes(col=col), size = 1) +
  scale_color_manual(values = c("A" = "grey", "B" = "blue", "C" = "red")) +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(face="bold", size = 15),
        axis.title.y = element_text(face="bold",size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(#face="bold",
          size=13, color = "black"),
        axis.text.y = element_text(#face="bold", 
          size=13, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_text_repel(max.overlaps = 9999) +
  xlab("Correlation (Rho)") +
  ylab(expression('-log'[10]*'(p-value)')) +
  geom_hline(yintercept = 1.30103, linetype = "dotted", color = "black") +
  ggtitle("Correlation with TTCS") + 
  scale_x_continuous(breaks = c(-1, -0.8, -0.6, -0.4, -0.2, -0.4, 0, 0.2, 0.4, 0.6, 0.8, 1)) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(0, 3))
p2

################################################
################ Figure 3 WGCNA ################
###############################################

## take neutrophil data
## samples in rows, proteins in columns
tree <- t_PMN_clean
tree$CAP <- NULL

#===============================================================================
#  Prepare/thresholds
#===============================================================================

## make sample tree to look for outliers => no big outliers
dev.off()
sampleTree<- hclust(dist(tree), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1,cex.axis = 1, cex.main = 2)

## decide which is the best beta to transform matrix to scale-free network
tree <- tree %>% as.data.frame()
powers <- c(c(1:20), seq(from=22, to=30, by=2))
sft <- pickSoftThreshold(tree, powerVector = powers, verbose = 5, networkType = "signed")
sft

## plot scale independence and mean connectivity 
par(mfrow = c(1,2));cex1 = 0.8;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.9,col="red") 

# Mean connectivity as a function of the soft-threshold power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

# select best power estimate (plot and function align)
power <- 15
sft$powerEstimate

#===============================================================================
#  Creating the similarity matrix => SIGNED (biological meaning)
#===============================================================================

## calculate similarity and dissimilarity matrix 
TOM <- TOMsimilarityFromExpr(tree, power = power, networkType = "signed",corType = "pearson")
dissTOM <- 1-TOM 

# Plot gene tree
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# Module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 3, 
                             pamRespectsDendro = FALSE, minClusterSize = 40)
table(dynamicMods)
length(table(dynamicMods)) 

# Convert numeric labels into colors
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

#===============================================================================
#  Merge modules if correlation coefficient is greater than 0.75
#===============================================================================
# Calculate eigengenes
MEList <- moduleEigengenes(tree, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs);

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average");
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres=0.25
abline(h=MEDissThres, col = "red")
merge <- mergeCloseModules(tree, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors <- merge$colors
table(mergedColors)
mergedMEs <- merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
names(mergedColors) <- colnames(tree)

## plot dimensions => nice visual split
cmd1 <- cmdscale(as.dist(dissTOM),2) 
par(mfrow=c(1,1)) 
plot(cmd1, col=as.character(mergedColors), main="MDS plot",xlab="Scaling Dimension 1",ylab="Scaling Dimension 2", 
     cex.axis=1.5,cex.lab=1.5, cex.main=1.5)

#===============================================================================
#  Evaluate dimensions and calculate module membership
#===============================================================================

## inspect intramodular connectivity of genes in module based on module membership (KME)
## MMPvalue = how significant a gene is assigned to a module
datKME<- signedKME(tree, mergedMEs, outputColumnName = "kME",corFnc = "cor", corOptions = "use = 'p', method='spearman'")
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(datKME), 83))

# check labeling => should be blue => correct
color_df <- mergedColors %>% as.data.frame()
x <- color_df[row.names(color_df) == "ARHGEF6", ]

#=================================================================================================================================
#  We aim to focus on the strongly assinged proteins per module
#  Remove noisy proteins that do not have a high KME to any of the modules  = proteins who are not strongly assigned to any module
#  Given that we used signed => genes are only stongly assigned when having a strong positive correlation => no need to take absolute values
#=================================================================================================================================

# Filter out genes with low KME for all modules
## check if correct
high_MM_genes <- rownames(datKME) %>% 
  setdiff(rownames(filter(datKME, if_all(starts_with("kME"), ~ . < 0.7))))
datKME_highMM <- datKME[row.names(datKME) %in% high_MM_genes, ]  
datKME_highMM$max <- apply(datKME_highMM, 1, max)
min(datKME_highMM$max ) ## correct

## 
tree_highMM <- tree[ , colnames(tree) %in% high_MM_genes]

## filter TOM 
TOM_highMM <- TOM
colnames(TOM_highMM) <- colnames(tree)
row.names(TOM_highMM) <- colnames(tree)
TOM_highMM <- TOM_highMM[rownames(TOM_highMM) %in% colnames(tree_highMM), colnames(TOM_highMM) %in% colnames(tree_highMM)]
dim(TOM_highMM)

## filter color labels => check if correct
merged_colors_no_low <- mergedColors[names(mergedColors) %in% high_MM_genes]
table(merged_colors_no_low)
table(colnames(tree_highMM) == names(merged_colors_no_low))

## still blue? => yes
color_df <- merged_colors_no_low %>% as.data.frame()
x <- color_df[row.names(color_df) == "ARHGEF6", ]

# Calculate slightly new eigengene using "core" proteins
MEList <- moduleEigengenes(tree_highMM, colors = merged_colors_no_low)
MEs <- MEList$eigengenes
MEDiss = 1-cor(MEs);
METree = hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

## see if modules are now too correlated => no
merge <- mergeCloseModules(tree_highMM, merged_colors_no_low, cutHeight = 0.25, verbose = 3) 
merged_colors_no_low <- merge$colors
table(merged_colors_no_low)
mergedMEs2 <- merge$newMEs

#=================================================================================================================================
#  Merged with clinical data
#=================================================================================================================================
## merge phenotype data
ME <- merge(mergedMEs2, Clinical_data, by.x = 0, by.y = "EB_id")
row.names(ME) <- ME$Row.names
ME$Row.names <- NULL
ME$EB_id <- row.names(ME)

## calculate gene trait significance 
ME <- ME[order(match(row.names(ME), row.names(tree_highMM))), ] ## same order
row.names(ME) == row.names(tree_highMM) ## correct

## relate modules to CAP
ME$CAP <- ifelse(ME$patient_volunteer == "Pneumonia patient", 1, 0)
CAP <- select(ME, "CAP")
CAP$CAP<- as.numeric(CAP$CAP)

geneTraitSignificance <- as.data.frame(cor(tree_highMM, CAP, use = "p", method = "spearman")) 
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 83)) 
GSPvalue$CAP <- p.adjust(GSPvalue$CAP, method = "BH")
GSPvalue <- -log10(GSPvalue$CAP)
plotModuleSignificance(GSPvalue, merged_colors_no_low, boxplot = F, ylim = c(0.0,10))
abline(h=-log10(0.05))

## confirm borderline with wilcox 
wilcox_loop <- function(x){
  p_value <- wilcox.test(ME[, x] ~ ME[ , "CAP"])[["p.value"]]
  vec <- p_value
  names(vec) <- "p_value"
  return(vec)
}

wilcox <- sapply(colnames(ME[1:9]), wilcox_loop) %>% as.data.frame() 
wilcox$p.adjust <- p.adjust(wilcox$., method = "BH")


## Evaluate intermodule correlation
ME_matrix <- select(ME,  "MEbrown",  "MEgreen", "MEpurple", "MEmagenta", "MEturquoise")
signif(cor(mergedMEs2, use = "p", method = "spearman"), 2)
correlation_matrix <- cor(ME_matrix)

### plot eigengene network
heat <- select(ME, "MEbrown", "MEgreen" ,  "MEpurple"    ,  "MEturquoise"   ,   "MEmagenta"   , "CAP")
plotEigengeneNetworks(heat, setLabels = NA, signed = T, plotAdjacency = T, colorLabels = T, coloredBarplot = T)

##
getwd()
ME_save <- select(ME, "EB_id", "MEbrown",  "MEgreen", "MEpurple", "MEmagenta", "MEturquoise")
openxlsx::write.xlsx(ME_save, "ME_matrix.PMN_final.xlsx")

## make corplot
rd_bu_palette <- COL2('RdBu', 200)
rd_bu_palette <- rev(rd_bu_palette)
corrplot(correlation_matrix, type = "upper", col = rd_bu_palette)


#####################################################################################
## Figure 3 - Extract proteins and plot enrichment results of CAP related modules ###
#####################################################################################

## green module
green <- names(merged_colors_no_low)[merged_colors_no_low == "green"]

## plot green module
ME$CAP <- factor(ME$CAP, levels = c("1", "0"))
colors <- c("#7e53a1","#1b9e77")
ggplot(ME, aes(x = CAP, y = MEgreen, fill = CAP)) +
  geom_boxplot() +
  scale_fill_manual(name = "Groups", values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")) +
  coord_cartesian(ylim = c(-0.4, 0.2))

##
result_green_bp <- enrichGO(gene          = green,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = "SYMBOL",
                            ont           = "BP",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            minGSSize     = 10,
                            maxGSSize     = 500)
result_green_bp_df <- result_green_bp %>% as.data.frame()

## change
p1_green <- cnetplot(result_green_bp, showCategory = 3, color_category='darkgreen', color_gene = "#AFE1AF", cex_label_gene = 0.6,  
                     cex_label_category = 0.7)

##
result_green_cc <- enrichGO(gene          = green,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = "SYMBOL",
                            ont           = "CC",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            minGSSize     = 10,
                            maxGSSize     = 500)
result_green_cc_df <- result_green_cc %>% as.data.frame()
top_terms <- head(result_green_cc_df, 3)

# Create a horizontal barplot
ggplot(top_terms, aes(x = -log10(p.adjust), y = reorder(Description, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(title = "Cellular component", x = "-log10(p.adjust)", y = "GO Term") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) + coord_flip()

###############
## brown module
################
brown <- names(merged_colors_no_low)[merged_colors_no_low == "brown"]

## plot brown module
boxplotbrown <- ggplot(ME, aes(x = CAP, y = MEbrown, fill = CAP)) +
  geom_boxplot() +
  scale_fill_manual(name = "Groups", values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )

##
result_brown_bp <- enrichGO(gene          = brown,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = "SYMBOL",
                            ont           = "BP",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            minGSSize     = 10,
                            maxGSSize     = 500)
result_brown_bp_df <- result_brown_bp %>% as.data.frame()

## cnetplot
p1_brown <- cnetplot(result_brown_bp, showCategory = 3, color_category='brown', 
                     color_gene = "#AFE1AF", cex_label_gene = 0.6,   cex_label_category = 0.7)

##
result_brown_cc <- enrichGO(gene          = brown,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = "SYMBOL",
                            ont           = "CC",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            minGSSize     = 10,
                            maxGSSize     = 500)
result_brown_cc_df <- result_brown_cc %>% as.data.frame()
top_terms <- head(result_brown_cc_df, 3)

# Create CC barplot
ggplot(top_terms, aes(x = -log10(p.adjust), y = reorder(Description, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(title = "Cellular component", x = "-log10(p.adjust)", y = "GO Term") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )  +
  coord_flip()


## magenta
magenta <- names(merged_colors_no_low)[merged_colors_no_low == "magenta"]

## plot magenta module
boxplotmagenta <- ggplot(ME, aes(x = CAP, y = MEmagenta, fill = CAP)) +
  geom_boxplot() +
  scale_fill_manual(name = "Groups", values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )

##
result_magenta_bp <- enrichGO(gene          = magenta,
                              OrgDb         = org.Hs.eg.db,
                              keyType       = "SYMBOL",
                              ont           = "BP",
                              pvalueCutoff  = 0.05,
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500)
result_magenta_bp_df <- result_magenta_bp %>% as.data.frame()
p1_magenta <- cnetplot(result_magenta_bp, showCategory = 3, color_category='magenta', color_gene = "#AFE1AF", cex_label_gene = 0.6,   cex_label_category = 0.7)

##
result_magenta_cc <- enrichGO(gene          = magenta,
                              OrgDb         = org.Hs.eg.db,
                              keyType       = "SYMBOL",
                              ont           = "CC",
                              pvalueCutoff  = 0.05,
                              pAdjustMethod = "BH",
                              minGSSize     = 10,
                              maxGSSize     = 500)
result_magenta_cc_df <- result_magenta_cc %>% as.data.frame()
top_terms <- head(result_magenta_cc_df, 3)

# Create CC barplot
ggplot(top_terms, aes(x = -log10(p.adjust), y = reorder(Description, p.adjust))) +
  geom_bar(stat = "identity") +
  labs(title = "Cellular component", x = "-log10(p.adjust)", y = "GO Term") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )  +
  coord_flip()

## purple
purple <- names(merged_colors_no_low)[merged_colors_no_low == "purple"]

## plot purple module
boxplotpurple <- ggplot(ME, aes(x = CAP, y = MEpurple, fill = CAP)) +
  geom_boxplot() +
  scale_fill_manual(name = "Groups", values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )

##
result_purple_bp <- enrichGO(gene          = purple,
                             OrgDb         = org.Hs.eg.db,
                             keyType       = "SYMBOL",
                             ont           = "BP",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH",
                             minGSSize     = 10,
                             maxGSSize     = 500)
result_purple_bp_df <- result_purple_bp %>% as.data.frame()
p1_purple <- cnetplot(result_purple_bp, showCategory = 3, color_category='purple', 
                      color_gene = "#AFE1AF", cex_label_gene = 0.6,   cex_label_category = 0.7)


##
result_purple_cc <- enrichGO(gene          = purple,
                             OrgDb         = org.Hs.eg.db,
                             keyType       = "SYMBOL",
                             ont           = "CC",
                             pvalueCutoff  = 0.05,
                             pAdjustMethod = "BH",
                             minGSSize     = 10,
                             maxGSSize     = 500)
result_purple_cc_df <- result_purple_cc %>% as.data.frame()
top_terms <- head(result_purple_cc_df, 3)

# Create CC barplot
ggplot(top_terms, aes(x = -log10(p.adjust), y = reorder(Description, -p.adjust))) +
  geom_bar(stat = "identity") +
  labs(title = "Cellular component", x = "-log10(p.adjust)", y = "GO Term") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  ) 
## turquoise
turquoise <- names(merged_colors_no_low)[merged_colors_no_low == "turquoise"]

## plot turquoise module
boxplotturquoise <- ggplot(ME, aes(x = CAP, y = MEturquoise, fill = CAP)) +
  geom_boxplot() +
  scale_fill_manual(name = "Groups", values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )

##
result_turquoise_bp <- enrichGO(gene          = turquoise,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = "SYMBOL",
                                ont           = "BP",
                                pvalueCutoff  = 0.05,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)
result_turquoise_bp_df <- result_turquoise_bp %>% as.data.frame()

p1_turquoise <- cnetplot(result_turquoise_bp, showCategory = 3, color_category='turquoise', 
                         color_gene = "#AFE1AF", cex_label_gene = 0.6,   cex_label_category = 0.7)
# Create CC barplot
result_turquoise_cc <- enrichGO(gene          = turquoise,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = "SYMBOL",
                                ont           = "CC",
                                pvalueCutoff  = 0.05,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)
result_turquoise_cc_df <- result_turquoise_cc %>% as.data.frame()
top_terms <- head(result_turquoise_cc_df, 3)

# Create a horizontal barplot
ggplot(top_terms, aes(x = -log10(p.adjust), y = reorder(Description, -p.adjust))) +
  geom_bar(stat = "identity") +
  labs(title = "Cellular component", x = "-log10(p.adjust)", y = "GO Term") +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "right",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )  

##
ggarrange(p1_green,p1_brown, p1_magenta, p1_purple, p1_turquoise, common.legend = T)

#########################################################
### Figure 4 - correlate with characteristics in CAP ###
########################################################

## filter module scores of CAP 
ME_pneumonia <- filter(ME, ME$CAP == 1)
mergedMEs2_pneumonia <- mergedMEs2[row.names(mergedMEs2) %in% ME_pneumonia$EB_id, ]
ME_pneumonia$EB_id == row.names(mergedMEs2_pneumonia) ## double check order => correct

## filter out non-significant modules
mergedMEs2_pneumonia <- mergedMEs2_pneumonia[ , !colnames(mergedMEs2_pneumonia) == "MEblue"]
mergedMEs2_pneumonia <- mergedMEs2_pneumonia[ , !colnames(mergedMEs2_pneumonia) == "MEblack"]
mergedMEs2_pneumonia <- mergedMEs2_pneumonia[ , !colnames(mergedMEs2_pneumonia) == "MEgreenyellow"]
mergedMEs2_pneumonia <- mergedMEs2_pneumonia[ , !colnames(mergedMEs2_pneumonia) == "MEtan"]

datTraits <- select(ME_pneumonia,	"age_yrs",	"gender",	
                    "BMI",	"symptom_days",
                    "CCI","ccd",	"hypertension",	"cpd",	"COPD", "Asthma",	"diabetes",	
                    "ckd",	"mneoplasm",	
                    "qSOFA_score",	"MEWS_score",	"CURB_score",	
                    "crp_1_1",	 "WBC_2_1",	"Neutrophil_unit_1",	"Lymphocyte_1_1",	
                    "Creatinine_value_1",
                    "Platelets_value_1",	
                    "Blood_Urea_Nitrogen_value_1",	
                    "Gram_positive", "Gram_negative", "Virus", "unknown_pathogen",
                    "clin_stab_hosp_old")

## rename variables
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "age_yrs"] <- "Age"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "gender"] <- "Sex"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "hypertension"] <- "Hypertension"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "ckd"] <- "CKD"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "mneoplasm"] <- "Malignancy"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "CURB_score"] <- "CURB-65 score"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "crp_1_1"] <- "CRP"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "WBC_2_1"] <- "Leukocyte count"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "Neutrophil_unit_1"] <- "Neutrophil count"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "Lymphocyte_1_1"] <- "Lymphocyte count"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "Platelets_value_1"] <- "Platelet count"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "Creatinine_value_1"] <- "Creatinin"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "Virus"] <- "Viral infection"
colnames(ME_pneumonia)[colnames(ME_pneumonia) == "clin_stab_hosp_old"] <- "TTCS"

## correct format
datTraits$gender <- ifelse(datTraits$gender == "Male", 1, 0)
datTraits <- data.frame(lapply(datTraits, function(x) ifelse(x == "Yes", 1, ifelse(x == "No", 0, x))))
datTraits <- data.frame(sapply(datTraits, function(x) as.numeric(as.character(x))))

##
moduleTraitCor <- round(cor(mergedMEs2_pneumonia, datTraits, use = "p", method = "spearman"), digits =2)
moduleTraitPvalue <- round(corPvalueStudent(moduleTraitCor, 57), digits= 3)

# Create a dataframe to store correlation results as excel
combined <- cbind(moduleTraitCor, moduleTraitPvalue)
colnames(combined) <- c(colnames(moduleTraitCor), paste0(colnames(moduleTraitCor), ".p"))
combined <- combined[ , order(colnames(combined))] %>% as.data.frame
combined <- t(round(combined, digits = 5))

## format dataframe names
for (i in seq(2, nrow(combined), by = 2)) {
  combined[(i-1), ] <- paste(combined[(i-1), ], combined[i, ], sep = " ")
}

# Iterate through each row and column and modify the values
for (i in 1:nrow(combined)) {
  for (j in 1:ncol(combined)) {
    cell_value <- combined[i, j]
    components <- strsplit(as.character(cell_value), " ")
    rho <- components[[1]][1]
    p <- components[[1]][2]
    formatted <- paste("Rho=", rho, ", p=", p)
    combined[i, j] <- formatted
  }
}

combined <- combined[-seq(2, nrow(combined), by = 2), ] %>% as.data.frame()
combined$variable <- row.names(combined)
openxlsx::write.xlsx(combined, "module_clinical_association_neutro.xlsx")

## make heatmap of significant onces
datTraits <- select(ME_pneumonia,	"Age",
                    "CCI", "Hypertension","Asthma",
                    "CKD",	
                    "CRP",	 "Leukocyte count",	"Lymphocyte count", "Neutrophil count",
                    "Viral infection", 
                    "TTCS")

datTraits <- data.frame(lapply(datTraits, function(x) ifelse(x == "Yes", 1, ifelse(x == "No", 0, x))))
datTraits <- data.frame(sapply(datTraits, function(x) as.numeric(as.character(x))))

moduleTraitCor <- cor(mergedMEs2_pneumonia, datTraits, use = "p", method = "spearman")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, 57)
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) == dim(moduleTraitCor)

dev.off()
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(mergedMEs2_pneumonia),
               colorLabels = FALSE,
               colors = blueWhiteRed(10),
               textMatrix = textMatrix,
               setStdMargins = T,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

#########################################################
### Figure 5 - plasma protein analysis               ###
########################################################

######################################################
## PREPROCESS,  IMPUTING, logFC between conditions  ##
######################################################

  ### Load data, filtering for proteotypic and with minimal 2 precursors
  {
    # prefiltering
    dia_proteotypic_2pr <- read.delim('Plasma_CAP_AMC_output.pr_matrix.tsv') %>%
      dplyr::select(Protein.Group, Genes, Proteotypic, Precursor.Id, contains('raw')) %>%
      gather('Sample', 'Intensity', contains('raw')) %>%
      filter(!is.na(Intensity)) %>%
      filter(Proteotypic == T) %>%
      group_by(Sample, Protein.Group) %>%
      distinct(Precursor.Id, .keep_all = T) %>%
      summarise(count = n()) %>%
      filter(count > 1)  %>%
      dplyr::select(-count)
    
    prot_filtered <- read.delim('Plasma_CAP_AMC_output.pg_matrix.tsv') %>%
      gather('Sample', 'Intensity', -Protein.Group, -Protein.Ids, -Protein.Names, -Genes, -First.Protein.Description) %>%
      filter(!is.na(Intensity)) %>%
      inner_join(dia_proteotypic_2pr, c('Protein.Group', 'Sample')) %>%
      # mutate(Sample = str_remove(Sample,'C..Users.massspecuser.Desktop.AjH.test.CAP.')) %>% 
      spread('Sample', 'Intensity') %>% 
      rownames_to_column('UID') %>% # rownames to UID
      mutate(UID = paste0('UID', UID)) %>% # paste UID
      mutate(across(contains('raw'),log2)) %>% 
      set_names(~gsub('(.*)a_|(.*)l_', '', .x))
  }  


### load condition
pheno <- read.csv2('sample_info.csv', stringsAsFactors = F) %>% 
  filter(Sample %in% gsub('.raw','',colnames(prot_filtered))) %>% 
  mutate(ct = factor(Type,levels = unique(Type)[2:1])) %>% 
  mutate(Type = factor(Type, levels = unique(Type)[2:1]))

###


###valid values, >75% in at least one (IALO)
vvs <- prot_filtered %>% 
  dplyr::select(contains('raw'),UID) %>% 
  column_to_rownames('UID') %>% 
  mutate(across(.cols = everything(), ~ !is.na(.))) %>%
  t() %>% 
  as.data.frame() %>%
  rownames_to_column('Sample') %>% 
  mutate(Sample = str_remove(Sample,'\\.raw')) %>% 
  left_join(pheno %>% dplyr::select(Sample,ct), c('Sample')) %>% 
  column_to_rownames('Sample') %>% 
  filter(!is.na(ct)) %>% 
  group_by(ct) %>% 
  summarise_all(sum) %>% 
  column_to_rownames('ct') %>% 
  mutate(across(.cols = everything(), ~ . / as.numeric(table(list(pheno$ct))))*100) %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(No_pass = rowSums(across(.cols = everything(), ~ . >= 75 ))) %>% # valid value threshold == 75%
  rownames_to_column('UID')
###

### IALO # percent
prot_valid <- prot_filtered %>% 
  subset(UID %in% vvs$UID[vvs$No_pass >=1])
###


### impute missing values
# define parameters
mann_downshift <- 1.8
mann_width <- 0.3

# create imputed values
set.seed(2112)
imputed_values <- prot_valid %>% 
  dplyr::select(paste0(pheno$Sample,'.raw')) %>% 
  mutate(across(.cols = everything(), ~ rnorm(length(.), (mean(na.omit(.) -  mann_downshift*sd(na.omit(.)))), (sd(na.omit(.))*mann_width))))

# replace NA's in DF with imputed
prot_imp <- prot_valid %>% 
  dplyr::select(paste0(pheno$Sample,'.raw')) %>% 
  map2_df(., imputed_values, coalesce) %>% # new dplyr
  bind_cols(Genes = prot_valid$Genes, UID = prot_valid$UID)

# check
prot_imp %>% 
  dplyr::select(paste0(pheno$Sample,'.raw')) %>% 
  boxplot()
###

### Statistical analysis
# make a design
design <- model.matrix(~0+as.factor(as.numeric(pheno$ct)))
colnames(design) = levels(pheno$ct)

# fit in limma
fit <- prot_imp %>% 
  dplyr::select(paste0(pheno$Sample,'.raw')) %>% 
  lmFit(.,design)

contrast_matrix <- makeContrasts(CAP-HC, levels =design)


# fit contrasts
fit2 <- contrasts.fit(fit, contrast_matrix)
# ebayes
ebayes <- eBayes(fit2)

#all results
results <- decideTests(ebayes, p.value = 0.05, adjust.method = 'BH', lfc = 1)
summary(results)


# get significance table
results_df <- ebayes %>% 
  topTable(number = Inf, coef = 1, sort.by = 'none') %>% 
  bind_cols(prot_valid)

sig_df <- ebayes %>% 
  topTable(number = Inf, coef = 1, sort.by = 'none') %>% 
  bind_cols(prot_valid) %>% 
  filter(adj.P.Val < 0.05) %>% 
  filter(abs(logFC)  > 1)


# Make table for export
export_df <- prot_filtered %>%
  set_names(~gsub('\\.raw', '_intensity', .x)) %>% 
  left_join(prot_imp %>% 
              dplyr::select(UID,contains('.raw')) %>% 
              set_names(~gsub('\\.raw', '_imputed_intensity', .x)) %>% 
              add_column(Quantified = 'True', .before = 1) %>% 
              bind_cols(ebayes %>% 
                          topTable(number = Inf, coef = 1, sort.by = 'none'))
            ,
            c('UID'))
# write out
write.csv2(export_df, file ='proteomics_AMC_CAP_plasma_export.csv')

## proteins per condition
plasma <- read.csv2("proteomics_AMC_CAP_plasma_export.csv")

#############################################
##       Figure 5a quantified proteins     ##
##############################################

## patients columns, measured proteins rows
## get columns of raw data => remove surfix => transpose (1 patient is 1 row)
raw_data <- plasma[ , 8:89]
colnames(raw_data) <- gsub("^X|_intensity$", "", colnames(raw_data))

## filter on CAP/healthy
CAP <- filter(Clinical_data, Clinical_data$patient_volunteer == "Pneumonia patient")
healthy <- filter(Clinical_data, Clinical_data$patient_volunteer == "Healthy volunteer")

##
plasma_CAP <- raw_data[, colnames(raw_data) %in% CAP$EB_id]
plasma_healthy <- raw_data[, colnames(raw_data) %in% healthy$EB_id, ]

## count missing 
plasma_CAP$percentage_missing <- rowSums(is.na(plasma_CAP))/length(plasma_CAP) * 100
sum(plasma_CAP$percentage_missing <25) ## 386

plasma_healthy$percentage_missing <- rowSums(is.na(plasma_healthy))/length(plasma_healthy) * 100
sum(plasma_healthy$percentage_missing <25) ## 367

##
categories <- c("CAP", "Healthy")
quantified_proteins <- c(386, 367)

# Create a data frame
data <- data.frame(Category = categories, Quantified_Proteins = quantified_proteins)
data$Category <- factor(data$Category, levels = c("Healthy", "CAP"))

# Plot
ggplot(data, aes(x = Quantified_Proteins, y = Category)) +
  geom_bar(stat = "identity", fill =  c("#800080", "#008000")) +
  geom_text(aes(label = Quantified_Proteins), hjust = -0.2, color = "white") +
  labs(x = "Quantified Proteins", y = "") +
  coord_cartesian(xlim = c(0, 400)) +
  theme_minimal() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        axis.ticks = element_line())

##################################################
##        Figure 1b PCA analysis             #####
##################################################

## take only quantified proteins
plasma <- filter(plasma, plasma$Quantified == "True")
x <- plasma[duplicated(plasma$Genes), ]
x <- plasma[plasma$Genes == "", ]

## modify empty gene names
plasma$Genes[plasma$Protein.Ids == "P0DOX2"] <- "IGHA2"
plasma$Genes[plasma$Protein.Ids == "P0DOX3"] <- "IGHD"
plasma$Genes[plasma$Protein.Ids == "P0DOX5"] <- "IGHG1"
plasma$Genes[plasma$Protein.Ids == "P0DOX7"] <- "IGKV"

## remove raw values and other abundant columns
plasma_clean <- plasma[ ,c(6, 91:172)]

## use gene names
row.names(plasma_clean) <- plasma_clean$Genes
plasma_clean$Genes <- NULL
t_plasma_clean <- t(plasma_clean) %>% as.data.frame()
new_row_names <- gsub("^X|_imputed_intensity$", "", rownames(t_plasma_clean))
rownames(t_plasma_clean) <- new_row_names

## save
plasma_clean_save <- plasma_clean
plasma_clean_save$CAP <- NULL
write.csv(t(plasma_clean), "plasma_clean.csv")

## dataset with rows as patients and columns as proteins (t_plasma_clean) => get condition
t_plasma_clean$CAP <- as.factor(ifelse(row.names(t_plasma_clean) %in% CAP$EB_id, 1, 0))

## make centered and scaled PCA
Groups <- t_plasma_clean[,"CAP"]
e.pca <- prcomp(t_plasma_clean[ ,1:395], center = T, scale = T) 
summary(e.pca)
str(e.pca)

##
colors <- c("#008000","#800080")
e.plot1 <- ggbiplot(e.pca, ellipse = FALSE, obs.scale = 1, var.scale = 1, var.axes = FALSE,
                    group = Groups, circle = FALSE, varname.size = 0, alpha = 0,
                    varname.adjust = c(1), ellipse.prob = 0.50) +
  scale_color_manual(name = "Groups", values = colors) +
  geom_point(aes(colour = Groups), size = 0.8, alpha = 1) +
  stat_ellipse(aes(colour = Groups), size = 1, type = "norm", level = 0.50, alpha = 0.7) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        text = element_text(size = 10),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        legend.position = "right") +
  ggtitle("Principal component analysis")
e.plot1

## statistically test components
pcdat <- as.data.frame(e.pca[["x"]])
pcdat$group <- Groups
summary(aov(pcdat$PC1 ~ pcdat$group))
summary(aov(pcdat$PC2 ~ pcdat$group))

#################################################
##        Figure 5c Volcano                  #####
##################################################

#color for volcanoplot
plasma$col <- "A"
plasma$col <- ifelse(plasma$logFC < 0 & plasma$adj.P.Val <=0.05, "B", plasma$col)
plasma$col <- ifelse(plasma$logFC > 0 & plasma$adj.P.Val <= 0.05, "C", plasma$col)
plasma$col <- ifelse(plasma$logFC > 0 &  plasma$logFC <1 & plasma$adj.P.Val<=0.05, "D", plasma$col)
plasma$col <- ifelse(plasma$logFC < 0 &  plasma$logFC > -1 & plasma$adj.P.Val<=0.05, "E", plasma$col)
table(plasma$col)

#label any sign markers
top10_high <- filter(plasma, plasma$adj.P.Val <= 0.05 & plasma$logFC > 1)
top10_high <- top10_high[order(top10_high$adj.P.Val), ]
top10_high <- top10_high[1:10, ]

top10_low <- filter(plasma, plasma$adj.P.Val <= 0.05 & plasma$logFC < -1)
top10_low <- top10_low[order(top10_low$adj.P.Val), ]
top10_low <- top10_low[1:10, ]

plasma$label <- NA
plasma <- plasma[order(plasma$adj.P.Val), ]
plasma$label <- ifelse(plasma$Genes %in% top10_high$Genes, plasma$Genes, 
                       ifelse(plasma$Genes %in% top10_low$Genes, plasma$Genes, NA))


# Modify plot aesthetics
p <- ggplot(plasma, aes(x = logFC, y = -log10(adj.P.Val), label = label)) +
  geom_point(aes(col = col), size = 0.7, alpha = 1) +
  scale_color_manual(values = c(A = "gray", B= "Blue", C = "red", D = "#FFB2B2", E = "#B2B2FF")) +
  scale_shape_manual(values = 1) +
  geom_text_repel(size = 4, nudge_x = 0.1, nudge_y = 0.1, max.overlaps = Inf, segment.color = "black") +
  xlab(expression(paste("logFC", " Fold Change"))) +
  ylab(expression(paste("-log"[10], " (BH adjusted P)"))) +
  geom_hline(yintercept = 1.30103, linetype = "dotted", color = "black") +
  geom_vline(xintercept = 1, linetype = "dotted", color = "black") +
  geom_vline(xintercept = -1, linetype = "dotted", color = "black") +
  ggtitle("Community-Acquired Pneumonia vs. Healthy Control") +
  theme_bw() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10),
    legend.title = element_blank(),
    legend.text = element_text(size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()) +
  coord_cartesian(xlim = c(-4, 10)) +
  scale_x_continuous(breaks = seq(-4, 10, by = 2))

##
# Create pie chart
table(plasma$col)
group_counts <- c( 92, 35, 13,  68,187)
group_colors <- c("#FFB2B2","red", "blue","#B2B2FF", "gray")
pie(group_counts, labels = NA, col = group_colors)

#################################################
##        Figure 5c WGCNA plasma            #####
##################################################

## samples in rows, proteins in columns
tree <- t_plasma_clean
dim(tree)
tree$CAP <- NULL

#===============================================================================
#  Prepare/thresholds
#===============================================================================

## make sample tree to look for outliers => no outliers
sampleTree<- hclust(dist(tree), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1,cex.axis = 1, cex.main = 2)
tree <- tree %>% as.data.frame()

## decide which is the best beta to transform matrix to scale-free network
powers <- c(c(1:20), seq(from=22, to=30, by=2))
sft <- pickSoftThreshold(tree, powerVector = powers, verbose = 5, networkType = "signed")

# plot scale independence and mean connectivity 
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") 

# Mean connectivity as a function of the soft-threshold power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity")) 
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

# select best power estimate (plot and function align) 
power <- 6
sft$powerEstimate

#===============================================================================
#  Creating the similarity matrix => SIGNED (biological meaning)
#===============================================================================

## calculate similarity and dissimilarity matrix 
TOM <- TOMsimilarityFromExpr(tree, power = power, networkType = "signed",corType = "pearson")
dissTOM <- 1-TOM 

# Plot gene tree
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# Module identification using dynamic tree cut
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 3, 
                             pamRespectsDendro = FALSE, minClusterSize = 40)
table(dynamicMods)
length(table(dynamicMods)) 

# Convert numeric labels into colors
## choose colors not used in significant WGCNA modules neutrophil
dynamicColors <- labels2colors(dynamicMods, colorSeq = c( "4" = "red","3" = "yellow", "2"= "blue","1"= "pink"))
table(dynamicColors)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")

#===============================================================================
#  Merge modules if correlation coefficient is greater than 0.75
#===============================================================================
# Calculate eigengenes
dev.off()
MEList <- moduleEigengenes(tree, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs);

# Cluster module eigengenes
METree <- hclust(as.dist(MEDiss), method = "average");
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres <- 0.25
abline(h=MEDissThres, col = "red")
merge <- mergeCloseModules(tree, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors <- merge$colors
table(mergedColors)
mergedMEs <- merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
names(mergedColors) <- colnames(tree)

## check => red
color_df <- mergedColors %>% as.data.frame()
x <- color_df[row.names(color_df) == "TAGLN2", ]

## plot dimensions => nice visual split
cmd1 <- cmdscale(as.dist(dissTOM),2) 
par(mfrow=c(1,1)) 
plot(cmd1, col=as.character(mergedColors), main="MDS plot",xlab="Scaling Dimension 1",ylab="Scaling Dimension 2", 
     cex.axis=1.5,cex.lab=1.5, cex.main=1.5)

#===============================================================================
#  Evaluate dimensions and calculate module membership
#===============================================================================

## inspect intramodular connectivity of genes in module based on module membership (KME)
## MMPvalue = how significant a gene is assigned to a module
datKME<- signedKME(tree, mergedMEs, outputColumnName = "kME",corFnc = "cor", corOptions = "use = 'p', method='spearman'")
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(datKME), 82)) ## Note; 1 patient less, no plasma

## red indeed highest MM
x <- datKME[row.names(datKME) == "TAGLN2", ]

#=================================================================================================================================
#  We aim to focus on the strongly assinged proteins per module
#  Remove noisy proteins that do not have a high KME to any of the modules  = genes who are not strongly assigned to any module
#  Given that we used signed => genes are only stongly assigned when having a strong positive correlation => no need to take absolute values
#=================================================================================================================================

# Filter out genes with low KME for all modules
## check if correct
high_MM_genes <- rownames(datKME) %>% 
  setdiff(rownames(filter(datKME, if_all(starts_with("kME"), ~ . < 0.7))))
datKME_highMM <- datKME[row.names(datKME) %in% high_MM_genes, ]  
datKME_highMM$max <- apply(datKME_highMM, 1, max) ## correct
min(datKME_highMM$max)

## 
tree_highMM <- tree[ , colnames(tree) %in% high_MM_genes]

## filter TOM 
TOM_highMM <- TOM
colnames(TOM_highMM) <- colnames(tree)
row.names(TOM_highMM) <- colnames(tree)
TOM_highMM <- TOM_highMM[rownames(TOM_highMM) %in% colnames(tree_highMM), colnames(TOM_highMM) %in% colnames(tree_highMM)]
dim(TOM_highMM)

## filter color labels
merged_colors_no_low <- mergedColors[names(mergedColors) %in% high_MM_genes]
table(merged_colors_no_low)
table(colnames(tree_highMM) == names(merged_colors_no_low)) ## correct

## check ≠> still red
color_df <- merged_colors_no_low %>% as.data.frame()
x <- color_df[row.names(color_df) == "TAGLN2", ]

# Calculate new eigengene using "core" proteins
MEList <- moduleEigengenes(tree_highMM, colors = merged_colors_no_low)
MEs <- MEList$eigengenes
MEDiss <- 1-cor(MEs);
METree <- hclust(as.dist(MEDiss), method = "average");
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

## see if modules are now too correlated => no
merge <- mergeCloseModules(tree_highMM, merged_colors_no_low, cutHeight = 0.25, verbose = 3) 
merged_colors_no_low <- merge$colors
table(merged_colors_no_low)
mergedMEs2 <- merge$newMEs

## check ≠> still red
color_df <- merged_colors_no_low %>% as.data.frame()
x <- color_df[row.names(color_df) == "TAGLN2", ]

#=================================================================================================================================
#  Merged with clinical data
#=================================================================================================================================
## merge phenotype data
ME <- merge(mergedMEs2, Clinical_data, by.x = 0, by.y = "EB_id")
row.names(ME) <- ME$Row.names
ME$Row.names <- NULL
ME$EB_id <- row.names(ME)

## calculate gene trait significance 
ME <- ME[order(match(row.names(ME), row.names(tree_highMM))), ] ## same order
row.names(ME) == row.names(tree_highMM) ## correct

## relate modules to CAP
ME$CAP <- ifelse(ME$patient_volunteer == "Pneumonia patient", 1, 0)
CAP <- select(ME, "CAP")
CAP$CAP<- as.numeric(CAP$CAP)

geneTraitSignificance <- as.data.frame(cor(tree_highMM, CAP, use = "p", method = "spearman")) 
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 82)) ## note 82 
GSPvalue$CAP <- p.adjust(GSPvalue$CAP, method = "BH")
GSPvalue <- -log10(GSPvalue$CAP)
plotModuleSignificance(GSPvalue, merged_colors_no_low, boxplot = F, ylim = c(0.0,20))
abline(h=-log10(0.05))

## check borderline with wilcox (pink borderline significant, seems borderline ns in plot)
wilcox_loop <- function(x){
  p_value <- wilcox.test(ME[, x] ~ ME[ , "CAP"])[["p.value"]]
  vec <- p_value
  names(vec) <- "p_value"
  return(vec)
}

wilcox <- sapply(colnames(ME[1:4]), wilcox_loop) %>% as.data.frame() 
wilcox$p.adjust <- p.adjust(wilcox$., method = "BH")

## check with WGCNA method => isolate p-value pink module => non significant
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), 82)) ## note 82 
GSPvalue$Module <- merged_colors_no_low  
pinkModulePValue <- GSPvalue[GSPvalue$Module == "pink", ]
median(pinkModulePValue$CAP)
mean(pinkModulePValue$CAP)

#####################################################################################
## Figure 3 - Extract proteins and plot enrichment results of CAP related modules ###
#####################################################################################

## blue module
blue <- names(merged_colors_no_low)[merged_colors_no_low == "blue"]

## plot blue module
ME$CAP <- factor(ME$CAP, levels = c("1", "0"))
colors <- c("#7e53a1","#1b9e77")
ggplot(ME, aes(x = CAP, y = MEblue, fill = CAP)) +
  geom_boxplot() +
  scale_fill_manual(name = "Groups", values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")) 

##
result_blue_bp <- enrichGO(gene          = blue,
                           OrgDb         = org.Hs.eg.db,
                           keyType       = "SYMBOL",
                           ont           = "BP",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH",
                           minGSSize     = 10,
                           maxGSSize     = 500)
result_blue_bp_df <- result_blue_bp %>% as.data.frame()
p1_blue <- cnetplot(result_blue_bp, showCategory = 3, color_category='darkblue', color_gene = "#AFE1AF", cex_label_gene = 0.6,  
                    cex_label_category = 0.7)

###############
## red module
################
red <- names(merged_colors_no_low)[merged_colors_no_low == "red"]

## plot red module
boxplotred <- ggplot(ME, aes(x = CAP, y = MEred, fill = CAP)) +
  geom_boxplot() +
  scale_fill_manual(name = "Groups", values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )

##
result_red_bp <- enrichGO(gene          = red,
                            OrgDb         = org.Hs.eg.db,
                            keyType       = "SYMBOL",
                            ont           = "BP",
                            pvalueCutoff  = 0.05,
                            pAdjustMethod = "BH",
                            minGSSize     = 10,
                            maxGSSize     = 500)
result_red_bp_df <- result_red_bp %>% as.data.frame()

## cnetplot
p1_red <- cnetplot(result_red_bp, showCategory = 3, color_category='red', 
                     color_gene = "#AFE1AF", cex_label_gene = 0.6,   cex_label_category = 0.7)


## yellow
yellow <- names(merged_colors_no_low)[merged_colors_no_low == "yellow"]

## plot yellow module
boxplotyellow <- ggplot(ME, aes(x = CAP, y = MEyellow, fill = CAP)) +
  geom_boxplot() +
  scale_fill_manual(name = "Groups", values = colors) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    legend.position = "none",
    axis.line.x = element_line(color = "black"),
    axis.line.y = element_line(color = "black")
  )

##
result_yellow_bp <- enrichGO(gene          = yellow,
                                OrgDb         = org.Hs.eg.db,
                                keyType       = "SYMBOL",
                                ont           = "BP",
                                pvalueCutoff  = 0.05,
                                pAdjustMethod = "BH",
                                minGSSize     = 10,
                                maxGSSize     = 500)
result_yellow_bp_df <- result_yellow_bp %>% as.data.frame()
p1_yellow <- cnetplot(result_yellow_bp, showCategory = 3, color_category='yellow', 
                         color_gene = "#AFE1AF", cex_label_gene = 0.6,   cex_label_category = 0.7)


### 
## plot eigengenenetwork
###
## Evaluate intermodule correlation
ME_matrix <- select(ME,  "MEred",  "MEblue", "MEyellow")
signif(cor(mergedMEs2, use = "p", method = "spearman"), 2)
correlation_matrix <- cor(ME_matrix)

# Get the "RdBu" palette with 200 colors
rd_bu_palette <- COL2('RdBu', 200)
rd_bu_palette <- rev(rd_bu_palette)
corrplot(correlation_matrix, type = "upper", col = rd_bu_palette)

### plot eigengene network
heat <- select(ME,  "MEred",  "MEblue", "MEyellow", "CAP")
plotEigengeneNetworks(heat, setLabels = NA, signed = T, plotAdjacency = T, colorLabels = T, coloredBarplot = T)

ME_save <- select(ME, "EB_id",  "MEred",  "MEblue", "MEyellow")
colnames(ME_save) <- c("EB_id", "PEred", "PEblue", "PEyellow")
openxlsx::write.xlsx(ME_save, "ME_matrix.plasma.final.xlsx")

#################################################
##        Figure 6 Interaction               #####
##################################################

## Figure 6a log FC plot       
## Create venn of overlap (right corner)
list_venn <- list(PMN = PMN$Genes, plasma = plasma$Genes)

dev.off()
venn_plot <- ggvenn(list_venn, c("PMN", "plasma")) +
  theme_minimal()

## overlapping proteins
PMN_similar <- PMN[PMN$Genes %in% plasma$Genes, ]
plasma_similar <- plasma[plasma$Genes %in% PMN$Genes, ]

logFC_PMN <- select(PMN_similar, "Genes", "logFC", "adj.P.Val")
colnames(logFC_PMN) <- c("Genes", "logFC_PMN", "adj.P.Val_PMN")
logFC_plasma <- select(plasma_similar, "Genes", "logFC", "adj.P.Val")
colnames(logFC_plasma) <- c("Genes", "logFC_plasma", "adj.P.Val_plasma")
logFC <- merge(logFC_PMN, logFC_plasma, by = "Genes")

logFC$difference <- abs(logFC$logFC_PMN - logFC$logFC_plasma)
logFC <- logFC[order(logFC$logFC_PMN), ]
logFC$col <- ifelse(logFC$adj.P.Val_PMN >= 0.05 & logFC$adj.P.Val_plasma >= 0.05, "A",
                    ifelse(logFC$adj.P.Val_PMN >= 0.05 & logFC$adj.P.Val_plasma <= 0.05, "B", 
                           ifelse(logFC$adj.P.Val_PMN <= 0.05 & logFC$adj.P.Val_plasma >= 0.05, "C", "D")))
table(logFC$col)

logFC <- logFC[order(logFC$adj.P.Val_PMN), ]
logFC$label <- NA
logFC$label[1:20] <- logFC$Genes[1:20]

logFC <- logFC[order(logFC$adj.P.Val_plasma), ]
logFC$label[1:20] <- logFC$Genes[1:20]

##
ggplot(logFC, aes(x=logFC_PMN, y= logFC_plasma, label=label )) +
  theme_bw() +
  geom_point(aes(col = col), size = 2, alpha = 1) +
  scale_color_manual(values = c("A" = "grey", "B" = "purple", "C" = "orange", "D" = "red")) +
  theme(legend.position = "none",
        plot.title = element_text(face="bold", size = 15),
        axis.title.y = element_text(face="bold",size=15),
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(#face="bold",
          size=13, color = "black"),
        axis.text.y = element_text(#face="bold", 
          size=13, color = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  geom_text_repel(max.overlaps = 9999) +
  xlab("LogFC Neutrophil") +
  ylab(expression("logFC plasma")) +
  ggtitle("Correlation PMN and plasma") + 
  geom_vline(xintercept = 0, linetype = "dotted", color = "black") +
  geom_hline(yintercept = 0, linetype = "dotted", color = "black") +
  coord_cartesian(ylim = c(-2, 3), xlim = c(-4,4))

# Create pie chart
table(logFC$col)
group_counts <- c(24, 40, 36, 58)
group_colors <- c("gray", "purple", "orange", "red")
pie(group_counts, labels = NA, col = group_colors)


#################################################
##        Figure 6b module correlation       #####
##################################################

plasma_matrix <- import("ME_matrix.plasma.final.xlsx")
PMN_matrix <- import("ME_matrix.PMN_final.xlsx")

#
both <- merge(plasma_matrix, PMN_matrix, by = "EB_id")
both_matrix <- select(both, "PEyellow","PEred", "PEblue", "MEbrown", "MEgreen", "MEpurple", "MEmagenta", "MEturquoise")

## corplot
dev.off()
correlation_matrix <- cor(both_matrix, method = "spearman")
rd_bu_palette <- COL2('RdBu', 200)
rd_bu_palette <- rev(rd_bu_palette)
corrplot(correlation_matrix, type = "upper", col = rd_bu_palette)


# Calculate correlation matrix and p-values
cor_results <- rcorr(as.matrix(both_matrix), type = "spearman")
correlation_matrix <- cor_results$r
p_values_matrix <- cor_results$P
print(p_values_matrix)
p_values_matrix <- ifelse(p_values_matrix  >0.05, "ns", p_values_matrix)


######################################################################
##     Figure 6c  correlation with inflamm response in plasma CAP    #####
#######################################################################

## load module data
plasma_module <- import("ME_matrix.plasma.final.xlsx")

##
CAP <- filter(Clinical_data, Clinical_data$patient_volunteer == "Pneumonia patient")
PMN_clean_CAP <- PMN_clean[row.names(PMN_clean) %in% CAP$EB_id, ]
PMN_clean_CAP <- merge(PMN_clean_CAP, plasma_module, by.x = 0 ,by.y = "EB_id")

#
correlation_spear <- function(x){
  test_result <- cor.test(PMN_clean_CAP [, x], PMN_clean_CAP [, "PEblue"], method = "spearman", exact= T)
  r <- test_result$estimate
  p <- test_result$p.value  
  vec <- c(r, p)
  names(vec) <- c("r", "p")
  return(vec)
}

## apply and format
proteins <- colnames(PMN_clean_CAP[,c(2:3483)])
cor.results_spear <- sapply(proteins, correlation_spear)
cor.results_spear <- as.data.frame(t(cor.results_spear))
cor.results_spear$marker <- row.names(cor.results_spear)
cor.results_spear$adj.p <- p.adjust(cor.results_spear$p, method = "BH")

#color for volcanoplot and label
cor.results_spear$col <- "A"
cor.results_spear$col <- ifelse(cor.results_spear$r<0 & cor.results_spear$adj.p<=0.05, "B", cor.results_spear$col)
cor.results_spear$col <- ifelse(cor.results_spear$r>0 & cor.results_spear$adj.p<=0.05, "C", cor.results_spear$col)

#label top markers
cor.results_spear$label <- NA
cor.results_spear <- cor.results_spear[order(-abs(cor.results_spear$r)), ]
cor.results_spear$label[1:30] <- cor.results_spear$marker[1:30]

# volcano
p2 <- ggplot(cor.results_spear, aes(x=r, y= -log10(adj.p), label=label)) +
  geom_point(aes(color=col), size=2) +
  scale_color_manual(values = c("A" = "grey", "B" = "blue", "C" = "red")) +
  theme_bw(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face="bold"),
    axis.title = element_text(face="bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill=NA, size=1)
  ) +
  geom_text_repel(max.overlaps = Inf, nudge_x = 0.05, nudge_y = 0.05, size = 3, color = "black") +
  xlab("Correlation (Rho)") +
  ylab(expression('-log'[10]*'(BH adjusted P)')) +
  geom_hline(yintercept = 1.30103, linetype = "dotted", color = "red") +
  coord_cartesian(xlim = c(-1, 1))
p2



