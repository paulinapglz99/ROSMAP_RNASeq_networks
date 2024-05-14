#
#3.rds_to_matrix.R
#Script to convert RDS file to matrix and equal diagonal MI matrix to 0
#Script by Aidee Lashmi and paulinapglz.99@gmail.com

#libraries --- ---

pacman::p_load("tidyverse", 
               "Matrix")

#Get data  --- ---
#Read Mutual information matrix in rds format

#slow
mi_list <- readRDS("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho.rds")
coexpre_mat_unnested<- mi_list
#Unnest RDS --- ---

coexpre_mat_unnested <- as.data.frame(do.call(cbind, mi_list)) %>% as.matrix()
dim(coexpre_mat_unnested) #this must be a square matrix
#[1] 14951 14951
#[1] 22070 22070

#If the diagonal of the matrix does not have the same values, we must verify
#that at least the elements in the diagonal are the highest values of the columns.

IM_test <- coexpre_mat_unnested[, 5] %>% as.data.frame() #place the column you want to test in the brackets [,n].
max(IM_test) #max number must be the diagonal value

#Set the diagonal to 0 --- ---

diag(coexpre_mat_unnested) <- 0

#Keep only the upper triangle of the matrix to avoid double information and double-edged matrices. 

coexpre_mat_unnested[lower.tri(coexpre_mat_unnested, diag = TRUE)] <- NA

#Save matrix --- ---

saveRDS(coexpre_mat_unnested, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho.rds")

#Pivot it to build an edgelist --- ---

#full_edgelist <- as.data.frame(coexpre_mat_unnested) %>%
#  rownames_to_column(var = "ensembl_gene_id") %>%
#  pivot_longer(-ensembl_gene_id, names_to = "gene_to", 
#              values_to = "MI")
#dim(full_edgelist)
# [1] 223532401         3

#full_edgelist$MI <- full_edgelist$MI %>% as.numeric()

#We need to delete NAs, as the MI matrix is only the upper triangle
#full_edgelist <- full_edgelist %>% na.omit()

#vroom::vroom_write(full_edgelist, file = '/datos/rosmap/coexpre_matrix/full_net_ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_edgelist.tsv')


#Next script is 4.my_network_analysis.R