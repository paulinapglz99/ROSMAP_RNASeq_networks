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
#If your data is not already scaled, then do Steps 1 and 2, if you don't need it
#skip to Step 3

#Step 1. Centering the Data

counts_centered <- scale(counts[-1], scale = F, center = T) 

#Step 2. #2. Scale the data

counts_scaled_centered <- scale(counts_centered, scale = T, center = F) 

#This scaling has the effect of making the column standard deviations equal to one

apply(counts_scaled_centered, 2, sd)

############################PCA 

#Step 3. Data Reduction

mat <- as.matrix(counts_scaled_centered)
rownames(mat) <- counts$gene_id  #Add rownames

#Calculate PCA

pca <- prcomp(t(mat)) #we t() to have samples in rows

#Elbow (Scree plot)

plot(pca)

#Loading Plot

barplot(pca$rotation[,1])

#PCA to table

pca_df <- pca$x %>% as.data.frame() %>% rownames_to_column(var = 'specimenID')
pca_df <- pca_df %>%
  as_tibble() %>% 
  left_join(metadata)

#Plot PCA

pca_df %>% 
  filter(PC2>-0.5 & PC1 < 1000000) %>%  #trampeando
  ggplot() +
  aes(x = PC1, y = PC2) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID))
