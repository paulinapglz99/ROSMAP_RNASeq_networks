#
#3.rds_to_matrix.R
#Script to convert RDS file to matrix and equal diagonal MI matrix to 0
#Script by Aidee Lashmi and paulinapglz.99@gmail.com

#libraries --- ---

pacman::p_load(tidyverse, Matrix)

#Get data  --- ---
#Read Mutual information matrix in rds format

mi_list <- read_rds('/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho.rds')

#Unnest RDS --- ---

coexpre_mat_unnested <- as.data.frame(do.call(cbind, mi_list)) %>% as.matrix()
dim(coexpre_mat_unnested) #this must be a square matrix
#[1] 14951 14951

#If the diagonal of the matrix does not have the same values, we must verify
#that at least the elements in the diagonal are the highest values of the columns.

IM_test <- coexpre_mat_unnested[, 5] %>% as.data.frame() #place the column you want to test in the brackets [,n].
max(IM_test) #max number must be the diagonal value

#Set the diagonal to 0 --- ---

diag(coexpre_mat_unnested) <- 0

#Save matrix --- ---

saveRDS(coexpre_mat_unnested, file = '/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_zero.rds')

#Next script is 4.my_network_analysis.R