#
#Script 0.Handle_metadata.R
#This pre-pre-process script reads ROSMAP metadata and FPKM counts data and handle both for further analysis
#By paulinapglz.99@gmail.com

#libraries  ----- 

pacman::p_load("dplyr")

#Read metadata --- ---

#This data was directly downloaded from https://www.synapse.org/#!Synapse:syn3157322

clinical_metadata <- vroom::vroom(file = '/datos/rosmap/metadata/ROSMAP_clinical.csv')
                          
biospecimen_metadata <- vroom::vroom(file = '/datos/rosmap/metadata/ROSMAP_biospecimen_metadata.csv')

#merge metadata in one dataset

cli_bio_metadata <- left_join(x = clinical_metadata, 
                              y = biospecimen_metadata, 
                              by = "individualID")
dim(cli_bio_metadata)
#[1] 14133    37

#save table for later

#vroom_write(cli_bio_metadata,
#           file = "/datos/rosmap/metadata/cli_bio_metadata.csv", 
#           delim = ",")

#Add other dx variants to metadata
#NIA reagan was constructed as https://www.sciencedirect.com/science/article/pii/S0197458097000572?via%3Dihub and
#https://www.radc.rush.edu/docs/var/detail.htm?category=Pathology&subcategory=Alzheimer%27s+disease&variable=niareagansc

cli_bio_metadata <- cli_bio_metadata %>%
  mutate(NIA_reagan_ADLikelihood = case_when(
    ceradsc == 4 ~ "4",  #No AD (0)
    (ceradsc == 3 & (braaksc == 1 | braaksc == 2)) ~ "3", #Low (0)
    (ceradsc == 2 & (braaksc == 3 | braaksc == 4)) ~ "2", #Intermediate
    (ceradsc == 1 & (braaksc == 5 | braaksc == 6)) ~ "1", #High
    TRUE ~ NA_character_  # Handle no-specified cases
  )) %>% 
  mutate(dicho_NIA_reagan = case_when(
    (NIA_reagan_ADLikelihood == 1 | NIA_reagan_ADLikelihood == 2) ~ "1", 
    (NIA_reagan_ADLikelihood == 3 | NIA_reagan_ADLikelihood == 4) ~ "0"
  )) %>% 
  mutate(final_diagnostic = case_when(
    cogdx == 1 & (NIA_reagan_ADLikelihood == 4) ~ "NCI",
    cogdx == 2 & (NIA_reagan_ADLikelihood == 3) ~ "MCI",
    cogdx == 4  & (NIA_reagan_ADLikelihood == 1 | NIA_reagan_ADLikelihood == 2) ~ "AD",
    TRUE ~ NA_character_  # Handle no-specified cases
  ))
dim(cli_bio_metadata)
#[1] 14133    40

#Questions about metadata --- ---

#N of participants with omic assays

cli_bio_metadata <- cli_bio_metadata %>% 
  filter(!is.na(assay))
dim(cli_bio_metadata)
#[1] 13107    40

N <- unique(cli_bio_metadata$individualID)
length(N)
#[1] 2558

#Demographics by sex

#Males

metadata_males <- cli_bio_metadata %>% filter(msex == 1)
dim(metadata_males)
#[1] 4228   40

N_males <- unique(metadata_males$individualID)
length(N_males)
#[1] 747
  
#Females

metadata_females <- cli_bio_metadata %>% filter(msex == 0)
dim(metadata_females)
#[1] 8873   40

N_females <- unique(metadata_females$individualID)
length(N_females)
#[1] 1810

#How much AD and no AD per sex --- ---

#Males
AD_males <- metadata_males %>% filter(dicho_NIA_reagan == 1)
length(unique(AD_males$individualID))
#[1] 201

#
noAD_males <- metadata_males %>% filter(dicho_NIA_reagan != 1 | is.na(dicho_NIA_reagan))
length(unique(noAD_males$individualID))  
#[1] 322

#Females
  
AD_females <- metadata_females  %>% filter(dicho_NIA_reagan == 1)
length(unique(AD_females$individualID))
#[1] 480




noAD_females <- metadata_females %>% filter(dicho_NIA_reagan != 1 | is.na(dicho_NIA_reagan))
length(unique(noAD_females$individualID))
#[1] 1330

#APOE per sex --- ---

#Males
APOE_males <- metadata_males %>% filter(apoe_genotype == "24" |apoe_genotype == "34" | apoe_genotype == "44" )
length(unique(APOE_males$individualID))
#[1] 187

no_APOE_males <- metadata_males %>% filter(apoe_genotype == "23" | apoe_genotype == "33" | apoe_genotype == "22" )
length(unique(no_APOE_males$individualID))
#[1] 556

#Females

APOE_females <- metadata_females %>% filter(apoe_genotype == "24" |apoe_genotype == "34" | apoe_genotype == "44" )
APOE_females <- unique(APOE_females$individualID)
length(APOE_females)
#[1] 449

no_APOE_females <- metadata_females %>% filter(apoe_genotype == "23" | apoe_genotype == "33" | apoe_genotype == "22" )
no_APOE_females <- unique(no_APOE_females$individualID)
length(no_APOE_females)
#[1] 1345

#Years of education per sex--- ---

#Males

educ_males_AD <- mean(AD_males$educ)
#[1] 17.23493

educ_males_noAD <- mean(noAD_males$educ)
#[1] 16.98173

#Females

educ_females_AD <- mean(AD_females$educ, na.rm = T)
#[1]  15.78124

educ_females_noAD <- mean(noAD_females$educ, na.rm = T)
#[1] 15.77209

#Number of reported omics trials --- ---

#Biocrates Bile Acids --- ---

#Males

biocrates_bile_males_AD <- metadata_males %>% filter(assay == "Biocrates Bile Acids")
length(unique(biocrates_bile_males_AD$individualID))
#[1] 31

#Females

biocrates_bile_females_AD <- metadata_females  %>% filter(assay == "Biocrates Bile Acids")
length(unique(biocrates_bile_females_AD$individualID))
#[1] 80
  
#Biocrates p180 --- ---

#Males

biocrates_p180_males_AD <- metadata_males %>% filter(assay == "Biocrates p180")
length(unique(biocrates_p180_males_AD$individualID))
#[1] 127

biocrates_p180_males_noAD <- metadata_males %>% filter(assay == "Biocrates p180")
length(unique(biocrates_p180_males_AD$individualID))
#[1] 127

#Females

biocrates_p180_females_AD <- metadata_females  %>% filter(assay == "Biocrates p180")
length(unique(biocrates_p180_females_AD$individualID))
#[1] 452

#Females

#ChIPSeq --- ---

#Males

#Females

#Metabolon --- ---

#Males

#Females

#TMT quantitation --- --- 

#Males

#Females

#Label free mass spectrometry --- ---

#Males

#Females

#MethylationArray --- ---
#Males

#Females

#miRNA Array  --- ---

#Males

#Females

#RNA Array  --- ---

#Males

#Females

#RNAseq  --- ---

#Males

#Females

#scRNASeq  --- ---

#Males

#Females

#SNP Array --- ---

#Males

#Females

#Whole Genome Sequencing --- ---

#Males

#Females

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