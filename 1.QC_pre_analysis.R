#
#1.QC_pre_analysis.R
#Script to build a PCA with RNA-seq data counts 
#paulinapglz.99@gmail.com

#I used the reference to build a PCA https://cran.r-project.org/web/packages/LearnPCA/vignettes/Vig_03_Step_By_Step_PCA.pdf

#The steps to carry out PCA are:
#1. Center the data
#2. Scale the data
#3. Carry out data reduction
#4.Plot! (Scree, PCA, etc.)

###Libraries --- ---

pacman::p_load("tidyverse", 
               "ggplot2", 
               "ggfortify", 
               "gridExtra")

#Get data --- --- 

counts <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/ROSMAP_RNAseq_rawcounts_DLPFC.txt")
dim(counts)
#[1] 60606   891 for DLFPC
#[1] 60606   750 for

#Obtain factors from metadata

metadata <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/RNA_seq_metadata_DLPFC.txt")
dim(metadata)
#[1] 1141   41

#Data handling --- ---

#If needed, trim first rows

counts <- counts[-c(1:3), ]

# Delete the "_PAR_Y" elements

counts <- counts[-grep("_PAR_Y", counts$feature),]
dim(counts)
# [1] 60558   891

#Prepare feature column

counts <- counts %>% mutate(feature = str_remove(feature, "\\..*$"), .before = 1)
dim(counts)
#[1] 60558   891

#Create sub metadata

sub_metadata <- metadata %>% filter(!is.na(libraryBatch))
dim(sub_metadata)
#[1] 637  41

sub_counts <- counts %>% select(c(1, any_of(sub_metadata$specimenID)))
dim(sub_counts)
#[1] 60558   637

#PCA --- ---

#Build matrix

pca_matrix <- sub_counts %>% 
  column_to_rownames("feature") %>% 
  as.matrix() %>% 
  t()  # transpose the matrix so that rows = samples and columns = variables, this because dots in the PCA scatterplot will be the ones in the rows

# Look at the first 10 rows and first 5 columns of the matrix
pca_matrix[1:10, 1:10]

# Perform the PCA

pca <- prcomp(pca_matrix, retx = TRUE, center = TRUE, scale. = FALSE)

#PCA to table

pca_df <- pca$x %>% as.data.frame() %>% rownames_to_column(var = 'specimenID')

# Create a data frame with PC number and percentage of variance

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

variance_table <- data.frame(
  PC = 1:length(pca$sdev),
  Variance_Percentage = pca$sdev^2 / sum(pca$sdev^2) * 100,
  cumulative_percentage = cumsum(pca$sdev^2 / sum(pca$sdev^2) * 100))

#Elbow (Scree plot)

ggplot(variance_table, aes(x = PC, y = Variance_Percentage)) +
  geom_bar(stat = 'identity', fill = variance_table$PC, color = 'black', position = 'dodge') +
  labs(title = 'Scree plot',
       subtitle = 'filtered data',
       x = 'Principal Components',
       y = 'Variance percentage') +
  scale_x_discrete() +
  theme_minimal()

#Table with the PCs explaining the 95% of data

variance_table_95 <- variance_table %>% 
  filter(cumulative_percentage >= 95)

#Plotting PCs explaining the 95% of data

ggplot(variance_table_95, aes(x = PC, y = Variance_Percentage)) +
  geom_bar(stat = 'identity', fill = '#95B8D1', color = 'black', position = 'dodge') +
  labs(title = 'Scree plot',
       subtitle = 'for the PCs explaining the 95% of data',
       x = 'Principal Components',
       y = 'Variance percentage') +
  theme_minimal()

#Plot PCA scatterplots --- --- 

#Option1: use autoplot

autoplot(pca, label = TRUE)

#Option 2: use ggplot 

#Colour by library batch --- ---

#By library batch

#Plot PC1 and PC2

PC1_PC2_librarybatch <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2, colour = as.factor(sub_metadata$libraryBatch)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by library batch",
       subtitle = paste("PC1 vs PC2"),
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) +
  theme_minimal()

#Plot PC1 and PC3

PC1_PC3_librarybatch <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC3, colour = as.factor(sub_metadata$libraryBatch)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by library batch",
       subtitle = "de PC1 vs PC3",
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC3 (",  sprintf("%.2f", variance_table$Variance_Percentage[3]), "%)"))
  theme_minimal()

#Plot PC2 and PC3

PC2_PC3_librarybatch <- pca_df %>% 
  ggplot() +
  aes(x = PC2, y = PC3, colour = as.factor(sub_metadata$libraryBatch)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID))  +
  labs(title = "PCA Scatterplot coloured by library batch",
       subtitle = "PC2 vs PC3", 
       x = paste("PC2 (", sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)"),
       y = paste("PC3 (",  sprintf("%.2f", variance_table$Variance_Percentage[3]), "%)")) +
  theme_minimal()

#By sequencing batch --- ---

#PC1 and PC2

PC1_PC2_seq_batch <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2, colour = as.factor(sub_metadata$sequencingBatch)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by sequencing batch",
       subtitle = "PC1 vs PC2", 
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) +
  theme_minimal()

#Plot PC1 and PC3

PC1_PC3_seq_batch <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC3, colour = as.factor(sub_metadata$sequencingBatch)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by sequencing batch",
       subtitle = "PC1 vs PC3",
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC3 (",  sprintf("%.2f", variance_table$Variance_Percentage[3]), "%)")) +
  theme_minimal()

#Plot PC2 and PC3

PC2_PC3_seq_batch <- pca_df %>% 
  # filter(PC2>-100000 & PC3 < 100000) %>%  #trampeando
  ggplot() +
  aes(x = PC2, y = PC3, colour = as.factor(sub_metadata$sequencingBatch)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID))  +
  labs(title = "PCA Scatterplot coloured by sequencing batch",
       subtitle = "PC2 vs PC3", 
       x = paste("PC2 (", sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)"),
       y = paste("PC3 (",  sprintf("%.2f", variance_table$Variance_Percentage[3]), "%)")) +
  theme_minimal()

# By study --- ---

#PC1 and PC2

PC1_PC2_study_batch <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2, colour = as.factor(sub_metadata$Study)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by study batch",
       subtitle = "PC1 vs PC2", 
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) +
  theme_minimal()

#Plot PC1 and PC3

PC1_PC3_study_batch <- pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC3, colour = as.factor(sub_metadata$Study)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by study batch",
       subtitle = "PC1 vs PC3", 
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC3 (",  sprintf("%.2f", variance_table$Variance_Percentage[3]), "%)")) +
  theme_minimal()

#Plot PC2 and PC3

PC2_PC3_study_batch<- pca_df %>% 
  ggplot() +
  aes(x = PC2, y = PC3, colour = as.factor(sub_metadata$Study)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID))  +
  labs(title = "PCA Scatterplot coloured by study batch",
       subtitle = "de PC2 vs PC3",
       x = paste("PC2 (", sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)"),
       y = paste("PC3 (",  sprintf("%.2f", variance_table$Variance_Percentage[3]), "%)")) +
  theme_minimal()

# By library prep --- ---

PC1_PC2_study_batch <- pca_df %>% 
  ggplot() +
    aes(x = PC1, y = PC2, colour = as.factor(sub_metadata$cogdx)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by study batch",
       subtitle = "PC1 vs PC2", 
       x = paste("PC1 (", sprintf("%.2f", variance_table$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table$Variance_Percentage[2]), "%)")) +
  theme_minimal()

#Grid plots

PC1_PC2 <- grid.arrange(PC1_PC2_librarybatch, PC1_PC2_seq_batch, PC1_PC2_study_batch)

PC1_PC3 <- grid.arrange(PC1_PC3_librarybatch, PC1_PC3_seq_batch, PC1_PC3_study_batch)

PC2_PC3 <- grid.arrange(PC2_PC3_librarybatch, PC2_PC3_seq_batch, PC2_PC3_study_batch)

#Save plots

ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/bias_plots/PC1_PC2_batches.png", PC1_PC2, height = 25,  dpi = 300)
ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/bias_plots/PC1_PC3_batches.png", PC1_PC3, height = 25,  dpi = 300)
ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/bias_plots/PC2_PC3_batches.png", PC2_PC3, height = 25,  dpi = 300)

#Loading plot --- ---

# Create a dataframe with the variable names and load values for PC1

loading_data_PC1 <- data.frame(feature = names(pca$rotation[, 1]),
                               PC1_Loadings = pca$rotation[, 1]) %>% 
                               mutate(loadings_Percentage = PC1_Loadings^2 / sum(PC1_Loadings^2) * 100, 
                               cumulative_percentages = cumsum(PC1_Loadings^2 / sum(PC1_Loadings^2) * 100)) %>%
  arrange(desc(loadings_Percentage))

loading_data_PC1 <- arrange(loading_data_PC1, PC1_Loadings)

# Extract the smallest value
min(loading_data_PC1$PC1_Loadings)
which.min(loading_data_PC1$PC1_Loadings)

#Plot contribution for PC1

#Option 1: by barplot
barplot(pca$rotation[,1])

#Option 2: ggplot2, slooow

ggplot(loading_data_PC1, aes(x = feature , y = loadings_Percentage)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(loadings_Percentage, 2)), vjust = -0.5, size = 3.5) +
  labs(title = "PC1 Loading plot",
    x = "Features",
    y = "Loadings percentage") +
  theme_minimal()

################################################################################

#There are genes that are presenting a very high variability

#loading_data_PC1_10th <- loading_data_PC1[1:20,]
#features_to_deplete <- loading_data_PC1_10th$feature

#sub_counts_depleted <- sub_counts %>% filter(!feature %in% features_to_deplete)
#dim(sub_counts_depleted)
#[1] 60548   638

#Rebuild matrix with depleted genes

#pca_matrix_depleted <- sub_counts %>% 
#  column_to_rownames("feature") %>% 
#  as.matrix() %>% 
#  t()  # transpose the matrix so that rows = samples and columns = variables, this because dots in the PCA scatterplot will be the ones in the rows

# Look at the first 10 rows and first 5 columns of the matrix
#pca_matrix_depleted[1:10, 1:10]

# Perform the PCA

#pca_depleted <- prcomp(pca_matrix_depleted, retx = TRUE, center = TRUE, scale. = FALSE)

#Plot PCA

#autoplot(pca_depleted, label = TRUE)

#TMM normalization --- ---

# Get log2 counts per million
logcpm <- cpm(y$counts,log=TRUE)
# Check distributions of samples using boxplots
group.col <- rep(c("red","blue"), each = 20)
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,col=group.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(unnormalised, coloured by groups)",cex.main=0.8)

# MDS plots
plotMDS(y,col=group.col)
legend("topright", legend = levels(as.factor(targets$Group)), fill=c("red","blue"))

#Finally, save data and metadata 

#Save metadata

#vroom::vroom_write(sub_metadata, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/RNA_seq_metadata_filtered_DLPFC.txt")

#Save counts

#vroom::vroom_write(sub_counts, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/ROSMAP_RNAseq_filtered_counts_DLPFC.txt")

#END