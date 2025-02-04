#
#Script 4.discretized_matrix.R
#This script takes RNA counts (already QCed and normalized with 2.prepro-mRNA) and calculates a coexpression matrix.
#With this coexpression matrix, the mutual information matrix is performed later on next block

#NOTE: I decided to leave all types of filtering, so use the script as u need to subset by clinical
#or pathology dx, I use only all NIA_reagan dichotomized to subset

#Libraries --- ---

pacman::p_load("tidyverse", 
               "infotheo")

#Read metadata -- ---

metadata <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_filteredQC_DLPFC.txt")
dim(metadata)
#[1] 1105   42

#Subset by diagnosis --- ---

#Filter by NIA-Reagan pathological diagnosis, using a dichotomized approach

specimenID_wNIA_reagan <- metadata %>% 
  filter(!is.na(dicho_NIA_reagan))
dim(specimenID_wNIA_reagan)
#[1] 774  42

AD_pathology <- metadata %>% filter(dicho_NIA_reagan == 1)
dim(AD_pathology)
#[1] 473  42

no_AD_pathology <-metadata %>% filter(dicho_NIA_reagan == 0)
dim(no_AD_pathology)
#[1] 301  42

#Read counts data --- ---

counts <- readRDS(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_filteredQC_counts_DLPFC.rds")
dim(counts) 
#[1] 28263   774

####  Data discretization --- ---

#this generates a discretized expression matrix

mat_dis <- infotheo::discretize(t(counts)) # I t() because we want genes to be columns to calculate MI in next script
dim(mat_dis)
#[1]   774 28263

#Regenerate gene names in cols and specimenIDs in rows

mat_dis <- mat_dis %>% mutate(specimenID = colnames(counts), .before = 1)
mat_dis <- mat_dis  %>% column_to_rownames(var = "specimenID")

# Finally we subset RNAseq counts by  diagnosis --- ---

# Subset RNAseq by pathological var 

#To obtain all NIA_reagan AD counts

AD_pathology_counts <- mat_dis %>% dplyr::filter(rownames(.) %in% AD_pathology$specimenID)
dim(AD_pathology_counts)

#To obtain all NIA_reagan noAD counts

noAD_pathology_counts <- mat_dis %>% dplyr::filter(rownames(.) %in% no_AD_pathology$specimenID)
dim(noAD_pathology_counts)

#Save discretized matrix --- --- 

saveRDS(AD_pathology_counts, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/ROSMAP_DLFPC_AD_NIAReagan_discretizedmatrix.rds")

saveRDS(noAD_pathology_counts, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/ROSMAP_DLFPC_noAD_NIAReagan_discretizedmatrix.rds")

#END