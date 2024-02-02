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

#Step 1. Centering the Data

counts_centered <- scale(counts[-1], scale = F, center = T) 

#Step 2. #2. Scale the data

counts_scaled_centered <- scale(counts[-1], scale = T, center = F) 

#This scaling has the effect of making the column standard deviations equal to one
#(or near to 1)

apply(counts_scaled_centered, 2, sd)

#Add rownames

rownames(counts_scaled) <- counts$gene_id

#Step 3. Data Reduction

pca_counts <- prcomp(counts_scaled_centered) #yay

#Scree plot

plot(pca_counts, main = "")

#Loading Plot

barplot(pca_counts$rotation[,1], main = "")