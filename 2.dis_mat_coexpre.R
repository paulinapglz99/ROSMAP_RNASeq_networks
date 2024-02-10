#This script takes FPKM normalized data from RNAseq and calculates a coexpression matrix.
#With this coexpression matrix, the mutual information matrix is performed later on next block
#Previous script is 1.1.pre-promRNA.R

#Libraries

pacman::p_load("tidyverse")

#Subset by cognitive diagnosis ------------------

#we have various diagnosis indicators in our metadata: cogdx, ceradsc
#First, we assign libraries for the respective cognitive diagnosis

#In 1.MatchFPKMandClinicalMetadata i already deleted cogdx == 1 and 6

metadata <- vroom::vroom(file = '/datos/rosmap/metadata/RNA_seq_metadata_080224.csv')

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
#[1] 14951   624

# Finally we subset RNAseq FPKMS by cogdx ----------

cogdxNCI_counts <- counts %>%
  dplyr::select(gene_id, all_of(NCI_cogdx$specimenID))
dim(FPKM_noMCI)
#[1] 55889   201  #201 specimenIDs

cogdxMCI_counts <- counts %>%
  dplyr::select(gene_id, all_of(MCI_cogdx$specimenID))
dim(FPKM_MCI)
# [1] 55889   169  #169 specimenIDs

cogdxAD_counts <- counts %>%
  dplyr::select(gene_id, all_of(AD_cogdx$specimenID))
dim(FPKM_AD)
#[1] 55889   255  #255 spcimenIDs

FPKM_confirmedAD <- counts %>% 
  dplyr::select(gene_id, all_of(confirmedAD$specimenID))
dim(FPKM_confirmedAD)
#[1] 55889   193  #193 specimenIDs

#Save tables

#vroom_write(FPKM_AD, 
#        file = 'FPKM_AD.csv', 
#       delim = ',')

FPKM <- FPKM_AD  #For non-hardcoding reasons, change this to FPKM_MCI, FPKM_noMCI, FPKM_confirmedAD, etc.

#Generate annotation with ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", 
                                "gene_biotype"),
                 filters = "ensembl_gene_id", 
                 values = FPKM$gene_id,
                 mart = mart)

##join biomart annotations with data table, to generate the annotated data

expre <-left_join(x = FPKM,
                  y = myannot,
                  by = c('gene_id'='ensembl_gene_id'))                

##Keep only the 'gene coding' type data (since for now they are the only ones we are interested in).

protcod <- filter(expre, 
                  gene_biotype =='protein_coding')
dim(protcod)
#[1] 19064   257

#delete name column, we dont need it for discretization

valores_expre <- protcod %>% 
  dplyr::select(-gene_id, 
                -percentage_gene_gc_content, 
                -gene_biotype)
dim(valores_expre)
#[1] 19064   254  #19064 genes

#pu

gene_names <- pull(protcod, 'gene_id')

#vroom_write(protcod, 
#            file = '/datos/rosmap/discretized_matrix/protcod_AD.txt',  #this for MCI and noMCI also
#            delim = ',')

####  Data discretization
#this generates a discretized expression matrix

mat_dis<-infotheo::discretize(t(valores_expre))
dim(mat_dis)
#[1]   254 19064  #still 19064 genes yay

###Next script is 3.mutualinformation_matrix.R
