#
#Script 2.dis_mat_coexpre.R
#This script takes RNA counts (already QCed and normalized) and calculates a coexpression matrix.
#With this coexpression matrix, the mutual information matrix is performed later on next block
#Previous script is 1.1.pre-promRNA.R

#NOTE: I decided to leave all types of filtering, so use the script as u need to subset by clinical
#or pathology dx, I use only all NIA_reagan dichotomized to subset

#Libraries --- ---

pacman::p_load("tidyverse")

#Subset counts by cognitive diagnosis --- ---

#we have various diagnosis indicators in our metadata: cogdx, ceradsc
#First, we assign libraries for the respective cognitive diagnosis

#In 1.MatchFPKMandClinicalMetadata i already deleted cogdx == 1 and 6

metadata <- vroom::vroom(file = '/datos/rosmap/metadata/RNA_seq_metadata_080224.csv')
dim(metadata)
#[1] 624  14

#Subset by diagnosis --- ---

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

#Read counts data --- ---

counts <- vroom::vroom(file = '/datos/rosmap/FPKM_data/ROSMAP_QCed_count_matrixfiltered_090224.tsv')
dim(counts)
#[1] 14951   625

# Finally we subset RNAseq counts by cogdx ----------

cogdxNCI_counts <- counts %>%
  dplyr::select(1, all_of(NCI_cogdx$specimenID))
dim(cogdxNCI_counts)
#[1] 14951   202

cogdxMCI_counts <- counts %>%
  dplyr::select(1, all_of(MCI_cogdx$specimenID))
dim(cogdxMCI_counts)
# [1] 14951   159

cogdxAD_counts <- counts %>%
  dplyr::select(1, all_of(AD_cogdx$specimenID))
dim(cogdxAD_counts)
#[1] 14951   223

# Subset RNAseq FPKMS by NIA_reagan dichotomous dx----------

NIA_reagan_counts <- counts %>% 
  dplyr::select(1, all_of(specimenID_wNIA_reagan$specimenID))
dim(NIA_reagan_counts)
#[1] 14951   435

AD_pathology_counts <- counts %>%
  dplyr::select(1, all_of(AD_pathology$specimenID))
dim(AD_pathology_counts)
#[1] 14951   256

noAD_pathology_counts <- counts %>%
  dplyr::select(1, all_of(no_AD_pathology$specimenID))
dim(noAD_pathology_counts)
#[1] 14951   180

#Save tables as needed --- --- 

vroom_write(NIA_reagan_counts, 
      file = '/datos/rosmap/FPKM_data/ROSMAP_NIA_reagandicho_counts_notdis.tsv')

####  Data discretization --- ---

#this generates a discretized expression matrix

mat_dis <- infotheo::discretize(t(NIA_reagan_counts[-1])) # I t() because we want genes to be columns to calculate MI in next script
dim(mat_dis)
#[1]   434 14951

#Regenerate gene names in cols and specimenIDs in rows

colnames(mat_dis) <- NIA_reagan_counts$ensembl_gene_id

mat_dis <- mat_dis %>% 
  mutate(specimenID = colnames(NIA_reagan_counts)[-1], .before = 1)

#Save discretized matrix --- --- 

#vroom::vroom_write(mat_dis, file = "/datos/rosmap/discretized_matrix/ROSMAP_allNIAReaganspecimen_discretizedmatrix_10022024.tsv")

###Next script is 3.mutualinformation_matrix.R