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

RNA_seq_metadata <- metadata %>% 
  filter(assay == 'rnaSeq', 
         organ == 'brain') %>% 
  dplyr::select(individualID, specimenID, msex, cogdx, ceradsc) #I also added columns for the factors later

# filter specific columns based on the values of "specimenID"
#in RNA_seq_metadata and then add a new column "gene_id" in the first position of 
#the result.

fpkm_matrix <- counts[,(colnames(counts) %in% RNA_seq_metadata$specimenID)] %>%  #only specimenIDs in the metadata
  mutate(gene_id = counts$gene_id, .before = 1)
fpkm_matrix$gene_id <- str_sub(fpkm_matrix$gene_id,1, 15)

#Simplify metadata to maintain only specimenIDs in the count matrix

RNA_seq_metadata <-  RNA_seq_metadata %>%
  filter(specimenID %in% colnames(counts))
dim(RNA_seq_metadata)
#[1] 636   5

#Subset by cognitive diagnosis ------------------

#First, we assign libraries for the respective cognitive diagnosis
#as the metadata dictionary says:

##1 NCI: No cognitive impairment (No impaired domains)
##2 MCI: Mild cognitive impairment (One impaired domain) and NO other cause of CI
##3 MCI: Mild cognitive impairment (One impaired domain) AND another cause of CI
##4 AD: Alzheimer’s dementia and NO other cause of CI (NINCDS PROB AD) 
##5 AD: Alzheimer’s dementia AND another cause of CI (NINCDS POSS AD)
##6 Other dementia: Other primary cause of dementia

#we don't need cogdx == 6, so we filter the metadata

AD_MI_cogdx <- RNA_seq_metadata %>% 
  filter(cogdx != 6)         
dim(AD_MI_cogdx)
# [1] 624   5

#Now we have only specimenIDs from cognitive diagnosis we are interested in
# Finally we subset RNAseq count tables to have only specimenIDs with cogdx wanted  ----------

counts_cogdx <- fpkm_matrix %>%
  dplyr::select(all_of(AD_MI_cogdx$specimenID))

#next script 1.1.prepro-mRNA.R
