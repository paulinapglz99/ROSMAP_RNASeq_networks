#Script that Matchs FPKM files and and ClinicalMetadata to prepare counts for the RNA-Seq QC 
#paulinapglz.99@gmail.com

#Libraries  ----- 

pacman::p_load('dplyr', 
               'stringr')

#Read data ---------

counts <- vroom::vroom(file = '/datos/rosmap/FPKM_data/RNAseq_FPKM_1_to_8_merged.csv') #ROSMAP all plates RNAseq normalized counts, archive from script 0
metadata <- vroom::vroom(file = '/datos/rosmap/metadata/cli_bio_metadata.csv') #merged metadata for all assays, archive from script 0

#Remove the last two characters from specimenID identifiers to match them with the RNASeq metadata

colnames(counts) <- substr(colnames(counts),1, nchar(colnames(counts))-2)
colnames(counts)[1] <- "gene_id"  #as we remove the last two characters we have to recompose the gene_id

#Simplify metadata to have only the RNAseq assays in brain  ------

metadata <- metadata %>% 
  filter(assay == 'rnaSeq', 
         organ == 'brain') %>% 
  dplyr::select(individualID, specimenID, assay, 
                msex, cogdx, ceradsc, braaksc) #I also added columns for the factors later

# filter specific columns based on the values of "specimenID"
#in RNA_seq_metadata and then add a new column "gene_id" in the first position of 
#the result.

counts_matrix <- counts[,(colnames(counts) %in% metadata$specimenID)] %>%  #only specimenIDs in the metadata
  mutate(gene_id = counts$gene_id, .before = 1)
counts_matrix$gene_id <- str_sub(counts_matrix$gene_id,1, 15)
dim(counts_matrix)
#[1] 55889   637

#Simplify metadata to maintain only specimenIDs in the count matrix

metadata <- metadata %>%
  filter(specimenID %in% colnames(counts_matrix))%>% 
  filter(cogdx != 6)  #we don't need cogdx == 6, other types of dementia
dim(metadata)
#[1] 636   5

#Now we have only specimenIDs from cognitive diagnosis we are interested in
#Finally we subset RNAseq count tables to have only specimenIDs with cogdx wanted  ----------

counts_matrix <- counts_matrix %>%
  dplyr::select(all_of(RNA_seq_metadata$specimenID)) %>% 
  mutate(gene_id = counts_matrix$gene_id, .before = 1)
dim(counts_matrix)
#[1] 55889   625

#Generate new variable based on ante & post mortem data ----------

metadata <- metadata %>%
  mutate(final_diagnostic = case_when(
    cogdx == 1 ~ "NCI",
    cogdx == 2 & (ceradsc == 1 | ceradsc == 2) ~ "AD",
    cogdx == 2 ~ "MCI",
    TRUE ~ NA_character_  # Handle no-specified cases
  ))

#Count of dx

dx_counts <- table(RNA_seq_metadata$final_diagnostic, useNA = "always") %>% as.data.frame()

#Save data --- ---

#Save filtered metadata

vroom::vroom_write(RNA_seq_metadata, file = '/datos/rosmap/metadata/RNA_seq_metadata_250124.csv')

#Save count table

vroom::vroom_write(counts_matrix, file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_250124.csv')

#next script 1.1.prepro-mRNA.R