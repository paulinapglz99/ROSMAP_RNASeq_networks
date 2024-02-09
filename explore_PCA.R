#Script to build a PCA with RNA-seq data counts 

###Libraries --- ---

pacman::p_load(tidyverse, 
               ggplot2)

#The steps to carry out PCA are:
#1. Center the data
#2. Optionally, scale the data
#3. Carry out data reduction

#Get data --- --- 

counts <- vroom::vroom(file = "/datos/rosmap/FPKM_data/filtered_FPKM_matrix_250124.csv")

colnames(counts)[1] <-"ensembl_gene_id" #change the name for further filtering

#Extract gene names to further managing

gene_names <- counts$ensembl_gene_id

#Manage data --- ---

#Annotation

### Generate annotation with ensembl ------ ------
#First we generate mart object

mart <- useEnsembl("ensembl",                         
                   dataset="hsapiens_gene_ensembl")

#We create myannot, with GC content, biotype, info for length & names per transcript

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values =  counts$ensembl_gene_id,  #annotate the genes in the count matrix 
                 mart = mart)

#Filtering ------ ------

#left join for further filtering

counts <- myannot %>% left_join(counts, 
                                    by = "ensembl_gene_id")
dim(counts)
#[1] 49399   631

#Filter to obtain only protein coding 

counts <- counts %>% 
  filter(gene_biotype == "protein_coding" & hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
  distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column
dim(expression)
#[1] 18848   631

#Obtain new annotation after filtering

myannot <- counts %>% 
  dplyr::select(1:7)
dim(myannot)
#[1] 18848     7

#Obtain counts 

expression_counts <- counts %>% 
  dplyr::select(ensembl_gene_id, 8:ncol(counts))      
dim(expression_counts)
#[1] 18848   625

#Obtain factors from metadata

metadata <- vroom::vroom(file = "/datos/rosmap/metadata/RNA_seq_metadata_250124.csv")
dim(metadata)
#[1] 624   4

#############PCA###########

#Step 1. Centering the Data

#If your data is not already scaled, scale(scale = T) if you don't need it
#just center with scale(center = T)

expression_counts <- scale(expression_counts[-1], scale = T, center = T) #omited first column as it contains only names
dim(counts)
#[1] 55889   624 genes, samples

#This centering and scaling together has the effect of making the column standard 
#deviations equal to 1. 

apply(expression_counts, 2, sd)

#PCA --- ---

#Step 3. Data Reduction

mat <- as.matrix(expression_counts)
rownames(mat) <- gene_names  #Add rownames
dim(mat)
#[1] 55889   625

#Calculate PCA

pca <- prcomp(t(mat)) #we t() to have samples in rows

#Scree table 

# Create a data frame with PC number and percentage of variance
variance_table <- data.frame(
  PC = 1:length(pca$sdev),
  Variance_Percentage = pca$sdev^2 / sum(pca$sdev^2) * 100,
  cumulative_percentage = cumsum(pca$sdev^2 / sum(pca$sdev^2) * 100))

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

#Elbow (Scree plot)

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
  filter(cumulative_percentage <= 95)

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
  aes(x = PC1, y = PC2) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID), colour = "#2a2e45") +
  labs(title = "Gráfico de Dispersión PCA",
       subtitle = "de PC1 vs PC2", x = "PC1 (22.98%)", y = "PC2 (17.58%)")

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

barplot(pca$rotation[,1], col = '#B48EAE')

# Crear un dataframe con los nombres de las variables y los valores de carga para PC1
loading_data_PC1 <- data.frame(Variable = names(pca$rotation[, 1]),
                               PC1_Loadings = pca$rotation[, 1]) %>% 
                           mutate(loadings_Percentage = PC1_Loadings^2 / sum(PC1_Loadings^2) * 100, 
                                 cumulative_percentages = cumsum(PC1_Loadings^2 / sum(PC1_Loadings^2) * 100))