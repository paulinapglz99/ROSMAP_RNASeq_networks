#
#5.2.matrix_to_edgelist.R
#Script that takes co-expression matrices (adjacency matrices) and converts them 
#into edgelist format, for later visualization and analysis.

#Libraries  --- --- 

pacman::p_load('tidyverse', 
               "Matrix")

#Read adjacency matrix --- ---

#This matrix may be very heavy
matrix <- readRDS("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.rds")
dim(matrix)
#[1] 22070 22070 <- AD

matrix <- as.data.frame(matrix)

#Pivot  ----

#this gives a table of connections between genes (edgelist), needed for the three approaches
#slow

full_edgelist <- matrix %>%
  rownames_to_column(var = "feature") %>%
  pivot_longer(-feature, names_to = "gene_to", 
               values_to = "MI")
dim(full_edgelist)
# [1] 223,532,401         3
#[1] 487,084,900         3

#As matrix is only the upper triangle, we must deplete NAs
class(full_edgelist$MI)
full_edgelist$MI <- unlist(full_edgelist$MI)
class(full_edgelist$MI)

full_edgelist <- full_edgelist %>% filter(!is.na(MI))
dim(full_edgelist)
#[1] 243531415         3

#save edgelist for later

vroom::vroom_write(full_edgelist, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_edgelist.txt')
#This edgelist can be found as /datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.txt.gz

#END