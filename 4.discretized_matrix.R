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

metadata <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_filteredQC_DLPFC.txt')
dim(metadata)
#[1] 497  41

#Subset by diagnosis --- ---

#Filter by NIA-Reagan pathological diagnosis, using a dichotomized approach

specimenID_wNIA_reagan <- metadata %>% 
  filter(!is.na(dicho_NIA_reagan))
dim(specimenID_wNIA_reagan)
#[1] 497  41

AD_pathology <- metadata %>% filter(dicho_NIA_reagan == 1)
dim(AD_pathology)
#[1] 181  41

no_AD_pathology <-metadata %>% filter(dicho_NIA_reagan == 0)
dim(no_AD_pathology)
#[1] 316  41

#Alternatively, for clinical diagnosis --- --- 

#Filter by clinical diagnosis (cogdx variable is blind to pathological post mortem diagnosis)
#In many studies, we can use only antemortem dx variables

#Library for no MCI (no dementia known at the moment of death)

NCI_cogdx <- metadata %>% filter(cogdx == 1)
dim(NCI_cogdx)
# [1] 201  14 

# For MCI (mild cognitive impairment) and NO other cause of CI

MCI_cogdx <- metadata %>% filter(cogdx == 2)
dim(MCI_cogdx)
#[1] 158   14

# Library for AD dementia an no other cause of CI

AD_cogdx <- metadata %>% filter(cogdx == 4)
dim(AD_cogdx)
#[1] 222   14   #254 individuals with probable or possible AD

#Read counts data --- ---

counts <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/ROSMAP_RNAseq_filteredQC_counts_DLPFC.txt') #Aqui es donde van los conteos de RNASeq
dim(counts) 
#[1] 22070   498

####  Data discretization --- ---

#this generates a discretized expression matrix

mat_dis <- infotheo::discretize(t(counts[-1])) # I t() because we want genes to be columns to calculate MI in next script

#Regenerate gene names in cols and specimenIDs in rows

colnames(mat_dis) <- counts$feature

mat_dis <- mat_dis %>% mutate(specimenID = colnames(counts)[-1], .before = 1)

# Finally we subset RNAseq counts by  diagnosis --- ---

# Subset RNAseq by pathological var 

#To obtain all NIA_reagan AD counts

AD_pathology_counts <- mat_dis %>% dplyr::filter(specimenID %in% AD_pathology$specimenID)
dim(AD_pathology_counts)
#[1] 22070   182
#[1]   181 22071 <-

#To obtain all NIA_reagan noAD counts

noAD_pathology_counts <- mat_dis %>% dplyr::filter(specimenID %in% no_AD_pathology$specimenID)
dim(noAD_pathology_counts)
#[1]   316 22071

#Save discretized matrix --- --- 

#vroom::vroom_write(AD_pathology_counts, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/ROSMAP_DLFPC_AD_NIAReagan_discretizedmatrix.txt")

#vroom::vroom_write(noAD_pathology_counts, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/ROSMAP_DLFPC_noAD_NIAReagan_discretizedmatrix.txt")

#Subset RNAseq counts by clinical var --- ---

cogdxNCI_counts <- mat_dis %>% dplyr::filter(specimenID %in% NCI_cogdx$specimenID)
dim(cogdxNCI_counts)
#[1]   157 22071

cogdxMCI_counts <- mat_dis %>% dplyr::filter(specimenID %in% MCI_cogdx$specimenID)
dim(cogdxMCI_counts)
#[1]   124 22071

cogdxAD_counts <- mat_dis%>% dplyr::filter(specimenID %in% AD_cogdx$specimenID)
dim(cogdxAD_counts)
#[1]   170 22071

#Save discretized matrix --- --- 

#vroom::vroom_write(cogdxNCI_counts, file = "/datos/rosmap/discretized_matrix/ROSMAP_NCI_cogdx_discretizedmatrix_10022024.txt")

#vroom::vroom_write(cogdxMCI_counts, file = "/datos/rosmap/discretized_matrix/ROSMAP_MCI_cogdx_discretizedmatrix_10022024.txt")

#vroom::vroom_write(cogdxAD_counts, file = "/datos/rosmap/discretized_matrix/ROSMAP_AD_cogdx_discretizedmatrix_10022024.txt")
 
#END