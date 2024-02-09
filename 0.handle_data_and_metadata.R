#
#Script 0.Handle_metadata.R
#This pre-pre-process script reads ROSMAP metadata and FPKM counts data and handle both for further analysis
#paulinapglz.99@gmail.com

#libraries  ----- 

pacman::p_load("dplyr")

#Read metadata ---------------
#This data was directly downloaded from https://www.synapse.org/#!Synapse:syn3157322

clinical_metadata <- vroom::vroom(file = '/datos/rosmap/metadata/ROSMAP_clinical.csv')
                          
biospecimen_metadata <- vroom::vroom(file = '/datos/rosmap/metadata/ROSMAP_biospecimen_metadata.csv')

#merge metadata in one dataset

cli_bio_metadata <- left_join(x = clinical_metadata, 
                              y = biospecimen_metadata, 
                              by = "individualID")

#save table for later

#vroom_write(cli_bio_metadata,
#           file = "/datos/rosmap/metadata/cli_bio_metadata.csv", 
#           delim = ",")

#read expression data, there's 2 count archives ---------------
#This data was directly downloaded from https://www.synapse.org/#!Synapse:syn3388564 

#In both archives I delete first column as it is repeated data

FPKM_p1_p6 <- vroom::vroom(file = '/datos/rosmap/FPKM_data/ROSMAP_RNAseq_FPKM_gene_plates_1_to_6_normalized.tsv')[-1] 
FPKM_p7_p8 <- vroom::vroom(file = '/datos/rosmap/FPKM_data/ROSMAP_RNAseq_FPKM_gene_plates_7_to_8_normalized.tsv')[-1] 

#merge expression data into one

FPKM_p1_p8<-left_join(x=FPKM_p1_p6,
                      y=FPKM_p7_p8,
                      by="gene_id")

#Here we handle with inopportune data like the copy of the gene ID col
#and two repeated samples with _6 and _7 are deleted (_0 was left)

FPKM_p1_p8 <- FPKM_p1_p8 %>% 
  dplyr::select(-c("492_120515_6", "492_120515_7")) #remove inopportune data 
          
#note: I lost 2 samples (492_120515_6 and 492_120515_7)

#save table for later

#vroom_write(FPKM_p1_p8,
#           file = "/datos/rosmap/FPKM_data/RNAseq_FPKM_1_to_8_merged.csv", 
#           delim = ",")

#Next script is 1.MatchFPKMandClinicalMetadata.R