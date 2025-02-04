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

mi_list <- readRDS("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/oldthings/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.rds")
# 
# args <- commandArgs(trailingOnly = TRUE)
# mi_list <- args[1]
# mi_list <- readRDS(mi_list)

#Unnest RDS --- ---

rows <- names(mi_list)
columns <- unique(unlist(lapply(mi_list, names)))

# Create an empty matrix
coexpre_mat <- matrix(0, nrow = length(rows), ncol = length(columns))
rownames(coexpre_mat) <- rows
colnames(coexpre_mat) <- columns

# Fill matrix with unnesting values
for (i in seq_along(mi_list)) {
  coexpre_mat[i, names(mi_list[[i]])] <- unlist(mi_list[[i]])
}

# Convertir a matriz
coexpre_mat <- as.matrix(coexpre_mat)

#If the diagonal of the matrix does not have the same values, we must verify
#that at least the elements in the diagonal are the highest values of the columns.
# 
# IM_test <- coexpre_mat[, 5] %>% as.data.frame() #place the column you want to test in the brackets [,n].
# max(IM_test) #max number must be the diagonal value

#Keep only the upper triangle of the matrix to avoid double information and double-edged matrices. 

coexpre_mat[lower.tri(coexpre_mat, diag = TRUE)] <- NA

#Save matrix --- ---

saveRDS(coexpre_mat,
        file = mi_list)
#END