#
#Script 0.Handle_data_and_metadata.R
#This pre-pre-process script reads ROSMAP metadata and FPKM counts data and handle both for further analysis
#By paulinapglz.99@gmail.com

#libraries  ----- 

pacman::p_load("dplyr", 
               "ggplot2", 
               "gridExtra", 
               "ggcorrplot")

#Read metadata --- ---

#This data was directly downloaded from https://www.synapse.org/#!Synapse:syn27000096

metadata <- vroom::vroom(file = '/datos/rosmap/data_by_counts/metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv')
dim(metadata)
#[1] 3400   38

#Merge metadata in one dataset
#Add other dx variants to metadata

#NIA reagan was constructed as https://www.sciencedirect.com/science/article/pii/S0197458097000572?via%3Dihub 2.B.Neuropathological Assessment
#and https://www.radc.rush.edu/docs/var/detail.htm?category=Pathology&subcategory=Alzheimer%27s+disease&variable=adnc

metadata <- metadata %>%
  filter(assay == "rnaSeq") %>% 
  mutate(NIA_reagan_ADLikelihood = case_when(         
    ceradsc == 4 ~ "0",  #No AD (0)
    (ceradsc == 1 & (braaksc == 5 | braaksc ==  6)) ~ "3", #High likelihood
    (ceradsc == 2 & (braaksc == 3 | braaksc == 4)) ~ "2", #Intermediate likelihood
    (ceradsc == 3 & (braaksc == 1 | braaksc == 2)) ~ "1", #Low likelihood
    TRUE ~ NA_character_  # Handle no-specified cases
  )) %>% 
  mutate(dicho_NIA_reagan = case_when(
    (NIA_reagan_ADLikelihood == 0 | NIA_reagan_ADLikelihood == 1) ~ "0", #no AD pathology
    (NIA_reagan_ADLikelihood == 2 | NIA_reagan_ADLikelihood == 3) ~ "1"  #AD pathology
  )) %>% 
  mutate(is_resilient = case_when(
    cogdx == 1 & (braaksc != 0 & (ceradsc == 1 | ceradsc ==2)) ~ "resilient", 
    TRUE ~ NA_character_ 
      ))
dim(metadata)
#[1] 2809   41
# 
# #Explore metadata --- ---
# 
# # If i'd like to observe only data per individual I need to delete duplicates, then only for exploration purposes
# 
# metadata_ind <- metadata %>% distinct(individualID, .keep_all = TRUE)
# dim(metadata)
# #[1] 1169   41
# 
# #Types of brain regions
# 
# unique(metadata_ind$tissue)
# 
# #I have data from brain regions
# #[1] "frontal cortex"                 "temporal cortex"                "dorsolateral prefrontal cortex"
# #[4] NA                               "Head of caudate nucleus"        "posterior cingulate cortex"    
# 
# #Exploring metadata by tissue --- ---
# 
# dim(metadata %>% filter(!is.na(tissue)))
# 
# met_tissue <- ggplot(metadata, aes(x = as.factor(tissue), fill = as.factor(tissue))) +
#   geom_bar() +
#   geom_text(stat='count', aes(label=..count..), vjust=-0.5) +  
#   labs(
#     title = "RNA-seq Samples sequenced by tissue",
#     x = "",
#     y = "") + 
#   scale_x_discrete(labels = c("dorsolateral prefrontal cortex" = "Dorsal prefrontal cortex", 
#                               "frontal cortex" = "Frontal cortex", 
#                               "Head of caudate nucleus" = "Head of caudate nucleus", 
#                               "posterior cingulate cortex" = "Posterior cingulate cortex", 
#                               "temporal cortex" = "Temporal cortex")) +
#   guides(fill = "none")  +
#   scale_fill_brewer(palette="Set1") +
#   theme_minimal()
# 
# ggsave("metadata_by_tissue.png", plot = met_tissue, 
#        device = "png",
#        width = 30, 
#        height = 15,
#        units = "cm", 
#        dpi = "print")

#Stratification by diagnosis --- ---

#For AD
AD_metadata <- metadata_ind %>% filter(dicho_NIA_reagan == 1)
dim(AD_metadata)
#[1] 594  41

N_AD <- length(unique(AD_metadata$individualID))
#[1] 594

#For no AD

noAD_metadata <- metadata_ind %>% filter(dicho_NIA_reagan == 0)
dim(noAD_metadata)
#[1] 312  41

N_noAD <- length(unique(noAD_metadata$individualID))
#[1] 312

#Distribution of sex

#For AD

AD_sex <- ggplot(AD_metadata, aes(x = as.factor(msex), fill = as.factor(msex))) +
         geom_bar() +
         geom_text(stat='count', aes(label=..count..), vjust=-0.5) +  
           labs(
         title = "For AD pathology patients",
         x = "Sex",
         y = "Count",
         fill = "Sex") + 
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_AD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
  scale_fill_brewer(palette="Set1") +
  theme_minimal()

#Vis

AD_sex

#For no AD

no_ADsex <-  ggplot(noAD_metadata, aes(x = as.factor(msex), fill = as.factor(msex))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +  
  labs(title = "For no AD pathology patients",
    x = "Sex",
    y = "Count", 
    fill = "Sex") +
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_noAD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
  scale_fill_brewer(palette="Set1") +  
  theme_minimal()

#Vis

no_ADsex

#Grid both graphs

sex_grid <- grid.arrange(AD_sex, no_ADsex, ncol = 2)

#Vis

sex_grid

#Distribution of education

AD_educ <- ggplot(AD_metadata, aes(x = as.factor(educ), fill = as.factor(educ))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +  
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_AD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
  theme(legend.position = "none") +
  labs(
    title = "For AD pathology patients",
    x = "Years of education",
    y = "Count"
  ) +  
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_AD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
  theme_minimal()
  
#For noAD

noAD_educ <-  ggplot(noAD_metadata, aes(x = as.factor(educ), fill = as.factor(educ))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) + 
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_noAD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
  theme(legend.position = "none") +
  labs(
    title = "For no AD pathology patients",
    x = "Years of education",
    y = "Count"
  ) +
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_noAD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
  theme_minimal()

grid.arrange(AD_educ, noAD_educ, ncol = 2)

#Race 

#For AD

AD_race <- ggplot(AD_metadata, aes(x = as.factor(race), fill = as.factor(race))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +  
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_AD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
  theme(legend.position = "none") +
  labs(
    title = "For AD pathology patients",
    x = "Race",
    y = "Count") +
  scale_fill_brewer(palette="Set1") +  
  theme_minimal() +
  scale_x_discrete(labels = c("White", "African American", "American Indian or Alaska Native"))

#For no AD

noAD_race <-  ggplot(noAD_metadata, aes(x = as.factor(race), fill = as.factor(race))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) + 
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_noAD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
  theme(legend.position = "none") +
  labs(
    title = "For no AD pathology patients",
    x = "Years of education",
    y = "Count"
  ) + 
  scale_fill_brewer(palette="Set1") +  
  theme_minimal() +
  scale_x_discrete(labels = c("White", "African American", "American Indian or Alaska Native"))


grid.arrange(AD_race, noAD_race, ncol = 2)

#Tissue 

AD_tissue <- ggplot(AD_metadata, aes(x = as.factor(tissue), fill = as.factor(tissue))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +  
  labs(
    title = "For AD pathology patients",
    x = "Tissue",
    y = "Count", 
    fill = "Tissue") +
  scale_fill_brewer(palette="Set1") +  
  theme_minimal()

noAD_tissue <- ggplot(noAD_metadata, aes(x = as.factor(tissue), fill = as.factor(tissue))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +  
  labs(
    title = "For AD pathology patients",
    x = "Tissue",
    y = "Count", 
    fill = "Tissue") +
  scale_fill_brewer(palette="Set1") +  
  theme_minimal()

grid.arrange(AD_tissue, noAD_tissue, ncol = 1)

#Distribution of Age of death

metadata_ind$age_death <- gsub("90\\+", "90", metadata_ind$age_death)

ade_death <- ggplot(metadata_ind, aes(x = dicho_NIA_reagan, y = as.numeric(age_death), fill = dicho_NIA_reagan)) +
  geom_violin() +
  theme_minimal() +
  labs(title = "Age of death",
    x = "Dichotomized NIA Reagan",
    y = "Age of death") +
  scale_fill_brewer(palette="Set1") +  
  theme_minimal()
  
#APOE genotype

AD_APOE <- ggplot(AD_metadata, aes(x = as.factor(apoe_genotype), fill = as.factor(apoe_genotype))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +  
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_AD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
    labs(
    title = "For AD pathology patients",
    x = "APOE genotype",
    y = "Count") +
    scale_fill_brewer(palette="Set1") +  
    theme_minimal()

noAD_APOE <-  ggplot(noAD_metadata, aes(x = as.factor(apoe_genotype), fill = as.factor(apoe_genotype))) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-0.5) +
  annotate("text", x = Inf, y = Inf, label = paste("Total N =", N_noAD), hjust = 1, vjust = 1, size = 4) +  # Agregar el total en una esquina
    labs(
    title = "For no AD pathology patients",
    x = "APOE genotype",
    y = "Count") +
  scale_fill_brewer(palette="Set1") +  
  theme_minimal()

grid.arrange(AD_APOE, noAD_APOE, ncol = 1)

#Distribution of post-mortem interval

PMI_interval <- ggplot(metadata_ind, aes(x = dicho_NIA_reagan, y = as.numeric(pmi), fill = dicho_NIA_reagan)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "Post mortem Inverval",
    x = "Dichotomized NIA Reagan",
    y = "Post mortem inverval"
  ) +
  scale_fill_brewer(palette="Set1") +  
  theme_minimal()

#RIN

RIN <- ggplot(metadata_ind, aes(x = dicho_NIA_reagan, y = as.numeric(RIN), fill = dicho_NIA_reagan)) +
  geom_boxplot() +
  theme_minimal() +
  labs(
    title = "RIN",
    x = "Dichotomized NIA Reagan",
    y = "RIN" 
  ) + 
  scale_fill_brewer(palette="Set1") +  
  theme_minimal()


RIN

#Define metadata by brain region --- ---

#Metadata for frontal cortex (FC)
metadata_FC <- metadata %>% filter(tissue == "frontal cortex")
dim(metadata_FC)
#[1] 123  41

#Metadata for temporal cortex (TC)
metadata_TC <- metadata %>% filter(tissue == "temporal cortex")
dim(metadata_TC)
#[1] 125  41

#Metadata Dorsoral Prefrontal Cortex (DLPFC)
metadata_DLPFC<- metadata %>% filter(tissue == "dorsolateral prefrontal cortex")
dim(metadata_DLPFC)
#[1] 1141   41

#Metadata for Head of caudate nucleus (HCN)
metadata_HCN <- metadata %>% filter(tissue == "Head of caudate nucleus")
dim(metadata_HCN)
#[1] 749  41

#Metadata for posterior cingulate cortex (PCC)
metadata_PCC <- metadata %>% filter(tissue == "posterior cingulate cortex")
dim(metadata_PCC)
#[1] 671  41

#Save metadata --- --- 

#Metadata for frontal cortex (FC)
# 
# vroom::vroom_write(metadata_FC, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/FC/RNA_seq_metadata_FC.txt")
# 
# #Metadata for temporal cortex (TC)
# 
# vroom::vroom_write(metadata_TC, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/TC/RNA_seq_metadata_TC.txt")
# 
# #Metadata Dorsoral Prefrontal Cortex (DLPFC)
# 
# vroom::vroom_write(metadata_DLPFC, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
# 
# #Metadata for Head of caudate nucleus (HCN)
# 
# vroom::vroom_write(metadata_HCN, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/HCN//RNA_seq_metadata_HCN.txt")
# 
# #Metadata for posterior cingulate cortex (PCC)
# 
# vroom::vroom_write(metadata_PCC, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/PCC/RNA_seq_metadata_PCC.txt")

#Read expression data, there's 4 count archives --- ---

#This data was directly downloaded from https://www.synapse.org/#!Synapse:syn3388564 

counts_one <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/raw_counts/ROSMAP_batch1_gene_all_counts_matrix_clean.txt')
counts_two <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/raw_counts/ROSMAP_batch2_gene_all_counts_matrix_clean.txt') 
counts_three <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/raw_counts/ROSMAP_batch3_gene_all_counts_matrix_clean.txt')
counts_four <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/raw_counts/ROSMAP_batch4_gene_all_counts_matrix_clean.txt')

#Merge expression data into one

counts <- dplyr::left_join(counts_one, counts_two, by = 'feature') %>%
  dplyr::left_join(counts_three, by = 'feature') %>% 
  dplyr::left_join(counts_four, by = 'feature')
dim(counts)
#[1] 60607  2911

#Counts from the frontal cortex

counts_FC <- counts[, (colnames(counts) %in% metadata_FC$specimenID)] %>% mutate(features = counts[1], .before = 1)
dim(counts_FC)
#[1] 60607     124

#Counts from the temporal cortex

counts_TC <- counts[, (colnames(counts) %in% metadata_TC$specimenID)] %>% mutate(counts[1], .before = 1)
dim(counts_TC)
#[1] 60607     126

#Counts for Dorsoral Prefrontal Cortex
counts_DLPFC <- counts[, (colnames(counts) %in% unique(metadata_DLPFC$specimenID))] %>% mutate(counts[1], .before = 1)
dim(counts_DLPFC)
#[1] 60607   1142

#Counts for  Head of caudate nucleus 
counts_HCN <- counts[, (colnames(counts) %in% unique(metadata_HCN$specimenID))] %>% mutate(counts[1], .before = 1)
dim(counts_HCN)
#[1] 60607   750

#Counts for posterior cingulate cortex
counts_PCC <- counts[, c(colnames(counts) %in% metadata_PCC$specimenID)] %>% mutate(counts[1], .before = 1)
dim(counts_PCC)
#[1] 60607   672

#Mapping counts

mapping_counts <- data.frame(
  Region = c("DLPFC", "HCN", "PCC", "TC", "FC"),
  samples = c(dim(counts_DLPFC)[2], dim(counts_HCN)[2], dim(counts_PCC)[2], dim(counts_TC)[2], dim(counts_FC)[2]), 
  features = c(dim(counts_DLPFC)[1], dim(counts_HCN)[1], dim(counts_PCC)[1], dim(counts_TC)[1], dim(counts_FC)[1])
  )


#Plot

mapping_counts.p <- ggplot(mapping_counts, aes(x = Region, y = samples, fill = Region)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = samples), vjust = -0.5, size = 3.5) + # Añadir números encima de las barras
  labs(
    title = "Number of specimens sequenced",
    x = "Brain region",
    y = "SpecimenID sequenced") +
  annotate("text", x = Inf, y = Inf, label = paste("Total features = 60,607"), hjust = 1, vjust = 1, size = 4) +
  scale_fill_brewer(palette="Set1") +  
  theme_minimal()

#Vis

mapping_counts.p

#Save count data ---- ---

#Counts for Dorsoral Prefrontal Cortex

vroom::vroom_write(counts_DLPFC, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.txt")

#Counts for  Head of caudate nucleus 
# 
# vroom::vroom_write(counts_HCN, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/HCN/ROSMAP_RNAseq_rawcounts_HCN.txt")
# 
# #Counts for posterior cingulate cortex
# 
# vroom::vroom_write(counts_PCC, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/PCC/ROSMAP_RNAseq_rawcounts_PCC.txt")
# 
# #Counts for Frontal Cortex
# 
# vroom::vroom_write(counts_FC, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/FC/ROSMAP_RNAseq_rawcounts_FC.txt")
# 
# #Counts for Temporal cortex
# 
# vroom::vroom_write(counts_TC, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/TC/ROSMAP_RNAseq_rawcounts_TC.txt")

#END