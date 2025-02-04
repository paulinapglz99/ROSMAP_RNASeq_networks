#Script to build a PCA with RNA-seq data counts 

###Libraries --- ---

pacman::p_load(tidyverse, 
               ggplot2)

#The steps to carry out PCA are:
#1. Center the data
#2. Scale the data
#3. Carry out data reduction
#4.Plot! (Scree, PCA, etc.)

#Get data --- --- 

counts <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/ROSMAP_RNAseq_rawcounts_DLPFC.txt")

colnames(counts)[1] <-"ensembl_gene_id" #change the name for further filtering

#If needed, trim first rows

counts <- counts[-c(1:3, ), ]

#Prepare feature column

counts$ensembl_gene_id <- substr(counts$ensembl_gene_id, 1, nchar(counts$ensembl_gene_id) - 2)

#Obtain factors from metadata

metadata <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/RNA_seq_metadata_DLPFC.txt")
dim(metadata)
#[1] 1141   41

#############PCA###########

#Step 1. Centering the Data

#If your data is not already scaled, scale(scale = T) if you don't need it
#just center with scale(center = T)

#counts <- scale(counts[-1], scale = T, center = T) #omited first column as it contains only names

#This centering and scaling together has the effect of making the column standard deviations equal to 1. 

#apply(counts, 2, sd)

#Step 3. Data Reduction

mat <- as.matrix(counts[-1])
rownames(mat) <- counts$ensembl_gene_id   #Add rownames again
mat <- t(mat) ##we t() to have samples in rows and genes in columns. 
dim(mat)
#[1] 60603   891

#Calculate PCA

pca <- prcomp(mat) #PCA will work on things in rows

#Plot PCA

autoplot(pca, label = TRUE, label.size = 3)

#Scree table 

# Create a data frame with PC number and percentage of variance
variance_table <- data.frame(
  PC = 1:length(pca$sdev),
  Variance_Percentage = pca$sdev^2 / sum(pca$sdev^2) * 100,
  cumulative_percentage = cumsum(pca$sdev^2 / sum(pca$sdev^2) * 100))

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

#Elbow (Scree plot), showing the amount of variance that each PC explains

ggplot(variance_table, aes(x = PC, y = Variance_Percentage)) +
  geom_bar(stat = 'identity', fill = '#B48EAE', color = 'black', position = 'dodge') +
  labs(title = 'Scree plot',
       subtitle = 'filtered data',
       x = 'Principal Components',
       y = 'Variance percentage') +
  scale_x_discrete() +
  theme_minimal()

#Table with the PCs explaining the 95% of data

variance_table_95 <- variance_table %>% 
  filter(cumulative_percentage < 95)

#Save table
#vroom::vroom_write(variance_table_95, file = "name_here.tsv")

#Plotting PCs explaining the 95% of data

ggplot(variance_table_95, aes(x = PC, y = Variance_Percentage)) +
  geom_bar(stat = 'identity', fill = '#95B8D1', color = 'black', position = 'dodge') +
  labs(title = 'Scree plot',
       subtitle = 'for the PCs explaining the 95% of data',
       x = 'Principal Components',
       y = 'Variance percentage') +
  scale_x_discrete() +
  theme_minimal()

#PCA to table

pca_df <- pca$x %>% as.data.frame() %>% rownames_to_column(var = 'specimenID')
pca_df <- pca_df %>%
  as_tibble() %>% 
  left_join(metadata)

#Plot PCA scoreplots --- --- 

#Plot PC1 and PC2

pca_df %>% 
# filter(PC2>-100000 & PC1 < 500000) %>%  #trampeando
  ggplot() +
  aes(x = PC1, y = PC2, colour = as.factor(notes)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA analysis",
       subtitle = "de PC1 vs PC2", x = "PC1 (93.23%)", y = "PC2 (2.65%)")

#Plot PC1 and PC3

pca_df %>% 
 #filter(PC3 < 100000  & PC1 < 500000 ) %>%  #trampeando
  ggplot() +
  aes(x = PC1, y = PC3) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID), colour = "#473144") +
  labs(title = "Gráfico de Dispersión PCA",
       subtitle = "de PC1 vs PC3", x = "PC1(22.98%)", y = "PC3(14.78%)")

#Plot PC2 and PC3

pca_df %>% 
# filter(PC2>-100000 & PC3 < 100000) %>%  #trampeando
  ggplot() +
  aes(x = PC2, y = PC3) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID), colour = '#8a7090')  +
  labs(title = "Gráfico de Dispersión PCA",
       subtitle = "de PC2 vs PC3", x = "PC2 (17.58%)", y = "PC3 (14.78%)")

#Loading Plot

x<- pca$rotation[,1] #extract PC1
x <- sort(x)
x[1]

barplot(pca$rotation[,1], col = '#B48EAE')

# Crear un dataframe con los nombres de las variables y los valores de carga para PC1
loading_data_PC1 <- data.frame(Variable = names(pca$rotation[, 1]),
                               PC1_Loadings = pca$rotation[, 1]) %>% 
                           mutate(loadings_Percentage = PC1_Loadings^2 / sum(PC1_Loadings^2) * 100, 
                                 cumulative_percentages = cumsum(PC1_Loadings^2 / sum(PC1_Loadings^2) * 100))