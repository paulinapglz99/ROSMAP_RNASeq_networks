#Script to build a PCA with RNA-seq data counts 

###Libraries --- ---

pacman::p_load(dplyr)

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

counts_scaled_centered <- scale(counts[-1], scale = T, center = F) 

#This scaling has the effect of making the column standard deviations equal to one
#(or near to 1)

apply(counts_scaled_centered, 2, sd)

#Add rownames

rownames(counts_scaled_centered) <- counts$gene_id

############################PCA 

#Step 3. Data Reduction

mat <- as.matrix(counts[,-1])
rownames(mat) <- counts$gene_id

pca <- prcomp(t(mat)) #we t() to have samples in rows

#Scree plot

plot(pca, main = "")

#Loading Plot

barplot(pca$rotation[,1], main = "")

#PCA to table

pca_df <- pca$rotation
pca_df %>%
  as_tibble() %>% 
  left_join(pca_counts, metadata, by = )

pca_df %>% 
  filter(PC2>-0.5) %>%  #trampeando
  ggplot() +
  aes(x = PC1, y = PC2, color= as.factor(msex)) +
  geom_point() +
  geom_text(mapping = aes(label = sample))
