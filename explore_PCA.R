#Script to build a PCA with RNA-seq data counts 

###Libraries --- ---

pacman::p_load(tidyverse)

#The steps to carry out PCA are:
#1. Center the data
#2. Optionally, scale the data
#3. Carry out data reduction
#4. Optionally, undo any scaling, likely using a limited number of PCs
#5. Optionally, undo the centering, likely using a limited number of PCs

#Get data --- --- 

counts <- vroom::vroom(file = "/datos/rosmap/FPKM_data/filtered_FPKM_matrix_250124.csv")
metadata <- vroom::vroom(file = "/datos/rosmap/metadata/cli_bio_metadata.csv")

#Manage data --- ---

#Extract gene data to further managing

gene_names <- counts$gene_id

#Step 1. Centering the Data

#If your data is not already scaled, scale(scale = T) if you don't need it
#just center with scale(center = T)

counts <- scale(counts[-1], scale = T, center = T) #omited first column as it contains only names
dim(counts)
#[1] 55889   624 genes, samples

#This centering and scaling together has the effect of making the column standard 
#deviations equal to 1. 

apply(counts, 2, sd)

############################PCA 

#Step 3. Data Reduction

mat <- as.matrix(counts)
rownames(mat) <- gene_names  #Add rownames
dim(mat)

#Calculate PCA

pca <- prcomp(t(mat)) #we t() to have samples in rows

#Elbow (Scree plot)

plot(pca, main = 'Scree plot', xlab = 'Principal Components from 1 to 10', 
     col = '#B48EAE',bg = 'lightgray')

#Scree table 

# Calculate the cumulative percentage of variance
cumulative_percentage <- cumsum(pca$sdev^2 / sum(pca$sdev^2) * 100)

# Create a data frame with PC number and percentage of variance
variance_table <- data.frame(
  PC = 1:length(pca$sdev),
  Variance_Percentage = pca$sdev^2 / sum(pca$sdev^2) * 100,
  cumulative_Percentage = cumulative_percentage)

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

variance_table_95 <- variance_table %>% 
  filter(cumulative_percentage <= 95)

#Save table
#vroom::vroom_write(variance_table_95, file = "variance95.tsv")

#PCA to table

pca_df <- pca$x %>% as.data.frame() %>% rownames_to_column(var = 'specimenID')
pca_df <- pca_df %>%
  as_tibble() %>% 
  left_join(metadata)

#Plot PCA scoreplot

#Plot PC1 and PC2

pca_df %>% 
  #filter(PC2>-0.5 & PC1 < 1000000) %>%  #trampeando
  ggplot() +
  aes(x = PC1, y = PC2, colour = as.factor(cogdx)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "Gráfico de Dispersión PCA",
       subtitle = "de PC1 vs PC2", x = "PC1 (17.08%)", y = "PC2 (12.08%)")

#Plot PC1 and PC3

pca_df %>% 
  ggplot() +
  aes(x = PC1, y = PC3) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "Gráfico de Dispersión PCA",
       subtitle = "de PC1 vs PC3", x = "PC1 (17.08%)", y = "PC3 (9.8%)")


#Plot PC2 and PC3

pca_df %>% 
  #filter(PC2>-0.5 & PC1 < 1000000) %>%  #trampeando
  ggplot() +
  aes(x = PC2, y = PC3) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID))  +
  labs(title = "Gráfico de Dispersión PCA",
       subtitle = "de PC2 vs PC3", x = "PC1 (12.08%)", y = "PC3 (9.8%)")

#Loading Plot

barplot(pca$rotation[,1], col = '#B48EAE')

# Crear un dataframe con los nombres de las variables y los valores de carga para PC1
loading_data_PC1 <- data.frame(Variable = names(pca$rotation[, 1]),
                               PC1_Loadings = pca$rotation[, 1]) %>% 
                           mutate(loadings_Percentage = PC1_Loadings^2 / sum(PC1_Loadings^2) * 100, 
                                 cumulative_percentages = cumsum(PC1_Loadings^2 / sum(PC1_Loadings^2) * 100)) 