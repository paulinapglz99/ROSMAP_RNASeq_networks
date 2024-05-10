#
#Script 4.discretized_matrix.R
#This script takes RNA counts (already QCed and normalized with 1.1.prepro-mRNA) and calculates a coexpression matrix.
#With this coexpression matrix, the mutual information matrix is performed later on next block
#Previous script is 1.1.pre-promRNA.R

#NOTE: I decided to leave all types of filtering, so use the script as u need to subset by clinical
#or pathology dx, I use only all NIA_reagan dichotomized to subset

#Libraries --- ---

pacman::p_load("tidyverse", 
               "infotheo")

#Read metadata -- ---

metadata <- vroom::vroom(file = '/datos/rosmap/metadata/RNA_seq_metadata_080224.csv')
dim(metadata)
#[1] 624  14

#Subset by diagnosis --- ---

#Filter by NIA-Reagan pathological diagnosis, using a dichotomized approach

specimenID_wNIA_reagan <- metadata %>% 
  filter(!is.na(dicho_NIA_reagan))
dim(specimenID_wNIA_reagan)
#[1] 434  14

AD_pathology <- metadata %>% 
  filter(dicho_NIA_reagan == 1)
dim(AD_pathology)
#[1] 255  14

no_AD_pathology <-metadata %>% 
  filter(dicho_NIA_reagan == 0)
dim(no_AD_pathology)
#[1] 179  14

#Alternatively, for clinical diagnosis --- --- 

#Filter by clinical diagnosis (cogdx variable is blind to pathological post mortem diagnosis)
#In many studies, we can use only antemortem dx variables

#Library for no MCI (no dementia known at the moment of death)

NCI_cogdx <- metadata %>% 
  filter(cogdx == 1)
dim(NCI_cogdx)
# [1] 201  14 

# For MCI (mild cognitive impairment) and NO other cause of CI

MCI_cogdx <- metadata %>% 
  filter(cogdx == 2)
dim(MCI_cogdx)
#[1] 158   14

# Library for AD dementia an no other cause of CI

AD_cogdx <- metadata %>% 
  filter(cogdx == 4)
dim(AD_cogdx)
#[1] 222   14   #254 individuals with probable or possible AD

#Read counts data --- ---

counts <- vroom::vroom(file = '/datos/rosmap/FPKM_data/ROSMAP_QCed_count_matrixfiltered_090224.tsv') #Aqui es donde van los conteos de RNASeq
dim(counts) 
#[1] 14951   625

####  Data discretization --- ---

#this generates a discretized expression matrix

mat_dis <- infotheo::discretize(t(counts[-1])) # I t() because we want genes to be columns to calculate MI in next script

#Regenerate gene names in cols and specimenIDs in rows

colnames(mat_dis) <- counts$ensembl_gene_id

mat_dis <- mat_dis %>% 
  mutate(specimenID = colnames(counts)[-1], .before = 1)

# Finally we subset RNAseq counts by  diagnosis --- ---

# Subset RNAseq by pathological var 

#To obtain all NIA_reagan_counts
NIA_reagan_counts <- counts %>% 
  dplyr::select(1, all_of(specimenID_wNIA_reagan$specimenID))
dim(NIA_reagan_counts)
#[1] 14951   435

#To obtain all NIA_reagan AD counts

AD_pathology_counts <- counts %>%
  dplyr::select(1, all_of(AD_pathology$specimenID))
dim(AD_pathology_counts)
#[1] 14951   256

#To obtain all NIA_reagan noAD counts

noAD_pathology_counts <- counts %>%
  dplyr::select(1, all_of(no_AD_pathology$specimenID))
dim(noAD_pathology_counts)
#[1] 14951   180

#Save discretized matrix --- --- 

#vroom::vroom_write(AD_pathology_counts, file = "/datos/rosmap/discretized_matrix/ROSMAP_AD_NIAReagan_discretizedmatrix_10022024.tsv")

#vroom::vroom_write(noAD_pathology_counts, file = "/datos/rosmap/discretized_matrix/ROSMAP_noAD_NIAReagan_discretizedmatrix_10022024.tsv")

# Subset RNAseq FPKMS by clinical var --- ---

cogdxNCI_counts <- counts %>% dplyr::select(1, all_of(NCI_cogdx$specimenID))
dim(cogdxNCI_counts)
#[1] 14951   202

cogdxMCI_counts <- counts %>% dplyr::select(1, all_of(MCI_cogdx$specimenID))
dim(cogdxMCI_counts)
# [1] 14951   159

cogdxAD_counts <- counts %>% dplyr::select(1, all_of(AD_cogdx$specimenID))
dim(cogdxAD_counts)
#[1] 14951   223

#Save discretized matrix --- --- 

#vroom::vroom_write(cogdxNCI_counts, file = "/datos/rosmap/discretized_matrix/ROSMAP_NCI_cogdx_discretizedmatrix_10022024.txt")

#vroom::vroom_write(cogdxMCI_counts, file = "/datos/rosmap/discretized_matrix/ROSMAP_MCI_cogdx_discretizedmatrix_10022024.txt")

#vroom::vroom_write(cogdxAD_counts, file = "/datos/rosmap/discretized_matrix/ROSMAP_AD_cogdx_discretizedmatrix_10022024.txt")

###Next script is 3.mutualinformation_matrix.R