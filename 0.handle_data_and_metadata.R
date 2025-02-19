#
#Script 0.Handle_data_and_metadata.R
#This pre-pre-process script reads ROSMAP metadata and counts data and handle both for further analysis
#By paulinapglz.99@gmail.com
#Modified by keilaperezf99@gmail.com

#This data was directly downloaded from https://www.synapse.org/#!Synapse:syn27000096

#libraries  ----- 

pacman::p_load("dplyr",
               "vroom", 
               "purrr",
               "ggplot2", 
               "viridis", 
               "gridExtra")

#Functions --- ---

# Function to filter and summarize data by brain region

summarize_by_tissue <- function(metadata, tissue_name) {
  tissue_data <- metadata %>% filter(tissue == tissue_name)
  cat(paste0("Metadata dim ", tissue_name, ": "), dim(tissue_data), "\n")
  cat("NIA-Reagan diagnosis table:\n")
  print(table(tissue_data$dicho_NIA_reagan, useNA = "ifany"))
  return(tissue_data)
}


summarise_metadata <- function(data, group_var, fill_var, title, fill_labels = NULL) {
  # Resumir los datos
  sum_rosmap <- data %>%
    select(tissue, all_of(fill_var)) %>%
    mutate(across(all_of(fill_var), ~ ifelse(is.na(.), "NA", as.character(.)))) %>%
    group_by(tissue, !!sym(fill_var)) %>%
    summarise(count = n(), .groups = 'drop')
  
  # Si se proporcionan etiquetas de llenado, úsalas para recodificar
  if (!is.null(fill_labels)) {
    sum_rosmap <- sum_rosmap %>%
      mutate(!!sym(fill_var) := recode(!!sym(fill_var), !!!fill_labels))
  }
  
  # Crear el gráfico horizontal
  plot <- ggplot(sum_rosmap, aes(x = count, y = tissue, fill = !!sym(fill_var))) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = count), position = position_stack(vjust = 0.5), hjust = -0.1) + 
    theme_minimal() +
    labs(x = "Count", y = "Tissue", fill = title) +
    ggtitle(title) +
    scale_fill_manual(values = c("red", "blue", "yellow", "green", "purple")) + 
    scale_x_continuous(limits = c(0, 950)) +
    theme(axis.text.y = element_text(hjust = 1), 
          legend.position = "bottom")  
  return(plot)
}

##################################### ROSMAP ##################################### 

#Read metadata --- ---

metadata_ROSMAP <- vroom::vroom(file = '/datos/rosmap/data_by_counts/metadata/RNAseq_Harmonization_ROSMAP_combined_metadata.csv')
dim(metadata_ROSMAP)
#[1] 3400   38

#Add other dx variants to metadata

#NIA reagan was constructed as https://www.sciencedirect.com/science/article/pii/S0197458097000572?via%3Dihub 2.B.Neuropathological Assessment
#and https://www.radc.rush.edu/docs/var/detail.htm?category=Pathology&subcategory=Alzheimer%27s+disease&variable=adnc

metadata_ROSMAP <- metadata_ROSMAP %>%
  filter(assay == "rnaSeq") %>%  
  mutate(NIA_reagan_ADLikelihood = case_when(         
    (ceradsc == 1 & (braaksc == 5 | braaksc ==  6)) ~ "3", #High likelihood
    (ceradsc == 2 & (braaksc == 3 | braaksc == 4)) ~ "2", #Intermediate likelihood
    (ceradsc == 3 & (braaksc == 1 | braaksc == 2)) ~ "1", #Low likelihood
    ceradsc == 4 ~ "0",  #No AD (0)
    TRUE ~ NA_character_  # Handle no-specified cases
  )) %>% 
  mutate(dicho_NIA_reagan = case_when(
    (NIA_reagan_ADLikelihood == 0 | NIA_reagan_ADLikelihood == 1) ~ "0", #no AD pathology
    (NIA_reagan_ADLikelihood == 2 | NIA_reagan_ADLikelihood == 3) ~ "1"  #AD pathology
  )) %>% 
  mutate(is_resilient = case_when(
    cogdx == 1 & (braaksc != 0 & (ceradsc == 1 | ceradsc ==2)) ~ "resilient", 
    TRUE ~ NA_character_ 
  )) %>% 
  mutate(is_AD = case_when(
    cogdx == 1 & (braaksc != 0 & (ceradsc == 1 | ceradsc ==2)) ~ "resilient",
    cogdx == 1 | ceradsc == 4 ~ "noAD",
    (cogdx %in% c(4, 5) & ceradsc == 1) ~ "AD",
    cogdx %in% c(2, 3) ~ "MCI",
    TRUE ~ NA_character_
  ))
dim(metadata_ROSMAP)
#[1] 2809   42

table(metadata_ROSMAP$is_resilient, useNA = "ifany")

table(metadata_ROSMAP$dicho_NIA_reagan, useNA = "ifany")

table(metadata_ROSMAP$is_AD, useNA = "ifany")

#Define metadata by brain region --- ---

tissues_ROSMAP <- unique(metadata_ROSMAP$tissue)
#[1] "frontal cortex"                 "temporal cortex"  "dorsolateral prefrontal cortex" "Head of caudate nucleus"       
#[5] "posterior cingulate cortex"    

metadata_tissue_ROSMAP <- lapply(tissues_ROSMAP, 
                                 function(tissue) summarize_by_tissue(metadata_ROSMAP, tissue))
names(metadata_tissue_ROSMAP) <- tissues_ROSMAP
# 
# Metadata dim frontal cortex:  123 42 
# NIA-Reagan diagnosis table:
#   
#   0    1 <NA> 
#   33   49   41 
# Metadata dim temporal cortex:  125 42 
# NIA-Reagan diagnosis table:
#   
#   0    1 <NA> 
#   33   51   41 
# Metadata dim dorsolateral prefrontal cortex:  1141 42 
# NIA-Reagan diagnosis table:
#   
#   0    1 <NA> 
#   307  486  348 
# Metadata dim Head of caudate nucleus:  749 42 
# NIA-Reagan diagnosis table:
#   
#   0    1 <NA> 
#   206  327  216 
# Metadata dim posterior cingulate cortex:  671 42 
# NIA-Reagan diagnosis table:
#   
#   0    1 <NA> 
#   196  272  203 

#Save metadata --- --- 

#Metadata for frontal cortex (FC)
# 
# vroom::vroom_write(metadata_tissue_ROSMAP[["frontal cortex"]],
#                    file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/FC/RNA_seq_metadata_FC.txt")
# #
# # #Metadata for temporal cortex (TC)
# #
# vroom::vroom_write(metadata_tissue_ROSMAP[["temporal cortex"]],
#                    file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/TC/RNA_seq_metadata_TC.txt")
#
# #Metadata Dorsolateral Prefrontal Cortex (DLPFC)
#
# vroom::vroom_write(metadata_tissue_ROSMAP[["dorsolateral prefrontal cortex"]],
#                    file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
#
# #Metadata for Head of caudate nucleus (HCN)
#
# vroom::vroom_write(metadata_tissue_ROSMAP[["Head of caudate nucleus"]],
#                    file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/HCN/RNA_seq_metadata_HCN.txt")
#
# #Metadata for posterior cingulate cortex (PCC)
#
# vroom::vroom_write(metadata_tissue_ROSMAP[["posterior cingulate cortex"]],
#                    file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/PCC/RNA_seq_metadata_PCC.txt")
#
#Read expression data, there's 4 count archives --- ---

#This data was directly downloaded from https://www.synapse.org/#!Synapse:syn3388564 

# Load all batch files
batch_files <- c('ROSMAP_batch1_gene_all_counts_matrix_clean.txt.gz',
                 'ROSMAP_batch2_gene_all_counts_matrix_clean.txt.gz',
                 'ROSMAP_batch3_gene_all_counts_matrix_clean.txt.gz',
                 'ROSMAP_batch4_gene_all_counts_matrix_clean.txt.gz')

counts_ROSMAP <- reduce(map(batch_files, vroom), left_join, by = 'feature')

#Define names and paths
tissue_names <- c("FC", "TC", "DLPFC", "HCN", "PCC")
tissue_paths <- c(
  "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/FC/ROSMAP_RNAseq_rawcounts_FC.rds",
  "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/TC/ROSMAP_RNAseq_rawcounts_TC.rds",
  "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.rds",
  "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/HCN/ROSMAP_RNAseq_rawcounts_HCN.rds",
  "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/PCC/ROSMAP_RNAseq_rawcounts_PCC.rds"
)

#Extract and save counts per tissue

Map(function(tissue, path, i) {
  tissue_counts <- counts_ROSMAP[, colnames(counts_ROSMAP) %in% metadata_tissue_ROSMAP[[i]]$specimenID] %>%
    mutate(counts_ROSMAP[1], .before = 1)
  
  saveRDS(tissue_counts, file = path)
}, tissue_names, tissue_paths, seq_along(tissue_names))

#Summarize ROSMAP --- ---

#By NIA-Reagan
rosmap_nia_reagan.p <- summarise_metadata(
  data = metadata_ROSMAP,
  group_var = "tissue",
  fill_var = "dicho_NIA_reagan",
  title = "NIA Reagan proportion by tissue - ROSMAP"
)
rosmap_nia_reagan.p

#by CERAD score

rosmap_ceradsc.p <- summarise_metadata(
  data = metadata_ROSMAP,
  group_var = "tissue",
  fill_var = "ceradsc",
  title = "CERAD score proportion by tissue - ROSMAP"
)
rosmap_ceradsc.p

#By diagnostic

rosmap_dx.p <- summarise_metadata(
  data = metadata_ROSMAP,
  group_var = "tissue",
  fill_var = "is_AD",
  title = "AD proportion by tissue - ROSMAP"
)
rosmap_dx.p

#DLFPC --- ---

metadata_DLFPC_ROSMAP <- metadata_ROSMAP %>% 
  filter(tissue == "dorsolateral prefrontal cortex")

#Covariate table --- ---

#Diagnosis

table(metadata_DLFPC_ROSMAP$dicho_NIA_reagan)

#Age by dx

age <- as.numeric(gsub("\\+", "", metadata_DLFPC_ROSMAP$age_death))
age_mean <- mean(age, na.rm = TRUE)
age_sd <- sd(age, na.rm = TRUE)

age_summary_by_dx <- metadata_DLFPC_ROSMAP %>%
  mutate(age_death = as.numeric(gsub("\\+", "", age_death))) %>% 
  group_by(dicho_NIA_reagan) %>%
  summarise(
    mean_age = mean(age_death, na.rm = TRUE),
    sd_age = sd(age_death, na.rm = TRUE)
  )

table(metadata_DLFPC_ROSMAP$dicho_NIA_reagan)

#Schooling

mean_schooling <-  metadata_DLFPC_ROSMAP %>%
  group_by(dicho_NIA_reagan) %>% 
  summarise(
    mean_educ = mean(educ, na.rm = TRUE),
    sd_educ = sd(educ, na.rm = TRUE)
  )

#Individuals by sex

sex_dx_counts <- metadata_DLFPC_ROSMAP %>%
  count(dicho_NIA_reagan, msex)  # Contar el número de individuos por diagnóstico y sexo

#Post-mortem interval

pmi_summary_by_dx <- metadata_DLFPC_ROSMAP %>%
  group_by(dicho_NIA_reagan) %>%
  summarise(
    mean_pmi = mean(pmi, na.rm = TRUE),
    sd_pmi = sd(pmi, na.rm = TRUE)
  )

#APOE genotype

apoe_summary_by_dx <- metadata_DLFPC_ROSMAP %>%
  count(dicho_NIA_reagan, apoe_genotype) 

#Race 

race_counts <- metadata_DLFPC_ROSMAP %>%
  group_by(dicho_NIA_reagan) %>% 
  mutate(race = case_when(
    race == 1 ~ "White",
    race == 2 ~ "Black or African American",
    race == 3 ~ "American Indian or Alaska Native",
    race == 4 ~ "Native Hawaiian or Other Pacific Islander",
    race == 5 ~ "Asian",
    race == 6 ~ "Other",
    race == 7 ~ "Unknown",
    TRUE ~ "Missing"
  )) %>% 
  count(race)

#Study

study_summary_by_dx <- metadata_DLFPC_ROSMAP %>%
  count(dicho_NIA_reagan, Study) 

#END