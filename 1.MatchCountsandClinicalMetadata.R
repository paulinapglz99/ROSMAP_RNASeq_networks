#
#Script 1.MatchCountsandClinicalMetadata.R
#This code  Matchs FPKM files and and ClinicalMetadata to prepare counts for the RNA-Seq QC 
#By paulinapglz.99@gmail.com

#Libraries  ----- 

pacman::p_load('dplyr', 
               'stringr')

#Read data ---------

counts <- vroom::vroom(file = '/datos/rosmap/FPKM_data/originalFPKM_RNASeq/RNAseq_FPKM_1_to_8_merged.csv')
dim(counts)
#[1] 55889   639
metadata <- vroom::vroom(file = '/datos/rosmap/metadata/cli_bio_metadata.csv') #merged metadata for all assays, archive from script 0

#Simplify metadata to have only the RNAseq assays in brain  ------

metadata <- metadata %>% 
  filter(assay == 'rnaSeq', 
         organ == 'brain') %>% 
  dplyr::select(individualID, specimenID, Study, assay,age_death, msex, pmi, 
                cogdx, ceradsc, braaksc) #I also added columns for the factors later

#Take the last two characters (batch info) from specimenID identifiers to match them with the RNASeq metadata

batches <- data.frame(
  specimenID = colnames(counts)[-1]
)

# Creates a new column named 'batch' with the last two characters
batches <- batches %>%
  mutate(batch = substr(specimenID, nchar(specimenID) - 1, nchar(specimenID)),
         specimenID = substr(specimenID, 1, nchar(specimenID) - 2))

#Delete _ from the batch column
batches$batch <- gsub("_", "", batches$batch)

# Add the new columns to metadata with batch info

metadata <- metadata %>%
  left_join(batches, by = "specimenID")

#Add other dx variants to metadata
#NIA reagan was constructed as https://www.sciencedirect.com/science/article/pii/S0197458097000572?via%3Dihub and
#https://www.radc.rush.edu/docs/var/detail.htm?category=Pathology&subcategory=Alzheimer%27s+disease&variable=niareagansc

metadata <- metadata %>%
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

#Fix count matrix column names to match it with metadata

colnames(counts) <- substr(colnames(counts)[],1, nchar(colnames(counts))-2) #delete batch info from count matrix
colnames(counts)[1] <- "gene_id"  #as we remove the last two characters we have to recompose the gene_id

#Create new count matrix only with columns present in our metadata

counts_matrix <- counts[,(colnames(counts) %in% metadata$specimenID)]  %>%  #only specimenIDs in the metadata
  mutate(gene_id = counts$gene_id, .before = 1)
counts_matrix$gene_id <- str_sub(counts_matrix$gene_id,1, 15) #delete version of genes for annotation purposes
dim(counts_matrix)
#[1] 55889   637

#Filter again metadata to maintain only specimenIDs in the count matrix

metadata <- metadata %>%
  filter(specimenID %in% colnames(counts_matrix))%>% 
  filter(cogdx != 6)  #we don't need cogdx == 6, other types of dementia
dim(metadata)
#[1] 636   13

#Count of dx

NIA_reagan_counts <- table(metadata$NIA_reagan_ADLikelihood, useNA = "always") %>% as.data.frame()

fdx_counts <- table(metadata$final_diagnostic, useNA = "always") %>% as.data.frame()

#Finally we filter again RNAseq count matrix to have only specimenIDs with metadata and viceversa  --- ---

counts_matrix <- counts_matrix %>%
  dplyr::select(all_of(metadata$specimenID)) %>% 
  mutate(gene_id = counts_matrix$gene_id, .before = 1) #added gene names to a first column 
dim(counts_matrix)
#[1] 55889   625

#Save data --- ---

#Save filtered metadata

#vroom::vroom_write(metadata, file = '/datos/rosmap/metadata/RNA_seq_metadata_080224.csv')

#Save count table

#vroom::vroom_write(counts_matrix, file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_250124.csv')

#next script 1.1.prepro-mRNA.R