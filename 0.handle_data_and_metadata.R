#
#Script 0.Handle_data_and_metadata.R
#This pre-pre-process script reads ROSMAP metadata and counts data and handle both for further analysis
#By paulinapglz.99@gmail.com
#Modified by keilaperezf99@gmail.com

#This data was directly downloaded from https://www.synapse.org/#!Synapse:syn27000096

#libraries  ----- 

pacman::p_load("dplyr",
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

counts_one <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/raw_counts/ROSMAP_batch1_gene_all_counts_matrix_clean.txt.gz')
counts_two <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/raw_counts/ROSMAP_batch2_gene_all_counts_matrix_clean.txt.gz') 
counts_three <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/raw_counts/ROSMAP_batch3_gene_all_counts_matrix_clean.txt.gz')
counts_four <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/raw_counts/ROSMAP_batch4_gene_all_counts_matrix_clean.txt.gz')

#Merge expression data into one

counts_ROSMAP <- dplyr::left_join(counts_one, counts_two, by = 'feature') %>%
  dplyr::left_join(counts_three, by = 'feature') %>% 
  dplyr::left_join(counts_four, by = 'feature')
dim(counts_ROSMAP)
#[1] 60607  2911

#Counts from the frontal cortex

counts_FC_ROSMAP <- counts_ROSMAP[, (colnames(counts_ROSMAP) %in% metadata_tissue_ROSMAP[[1]]$specimenID)]%>% 
  mutate(counts_ROSMAP[1], .before = 1)
dim(counts_FC_ROSMAP)
#[1] 60607     124

#Counts from the temporal cortex

counts_TC_ROSMAP <- counts_ROSMAP[, (colnames(counts_ROSMAP) %in% metadata_tissue_ROSMAP[[2]]$specimenID)] %>% 
  mutate(counts_ROSMAP[1], .before = 1)
dim(counts_TC_ROSMAP)
#[1] 60607     126

#Counts for Dorsoral Prefrontal Cortex
counts_DLPFC_ROSMAP <- counts_ROSMAP[, (colnames(counts_ROSMAP) %in% unique(metadata_tissue_ROSMAP[[3]]$specimenID))] %>%
  mutate(counts_ROSMAP[1], .before = 1)
dim(counts_DLPFC_ROSMAP)
#[1] 60607   1142

#Counts for  Head of caudate nucleus 
counts_HCN_ROSMAP <- counts_ROSMAP[, (colnames(counts_ROSMAP) %in% unique(metadata_tissue_ROSMAP[[4]]$specimenID))] %>% 
  mutate(counts_ROSMAP[1], .before = 1)
dim(counts_HCN_ROSMAP)
#[1] 60607   750

#Counts for posterior cingulate cortex
counts_PCC_ROSMAP <- counts_ROSMAP[, c(colnames(counts_ROSMAP) %in% metadata_tissue_ROSMAP[[5]]$specimenID)] %>%
  mutate(counts_ROSMAP[1], .before = 1)
dim(counts_PCC_ROSMAP)
#[1] 60607   672

#Save count data for ROSMAP ---- ---

#Counts for Dorsoral Prefrontal Cortex
# # 
# saveRDS(counts_DLPFC_ROSMAP, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.rds")
# # 
# # #Counts for  Head of caudate nucleus
# # 
# saveRDS(counts_HCN_ROSMAP, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/HCN/ROSMAP_RNAseq_rawcounts_HCN.rds")
# # 
# # #Counts for posterior cingulate cortex
# # 
# saveRDS(counts_PCC_ROSMAP, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/PCC/ROSMAP_RNAseq_rawcounts_PCC.rds")
# # 
# # #Counts for Frontal Cortex
# # 
# saveRDS(counts_FC_ROSMAP, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/FC/ROSMAP_RNAseq_rawcounts_FC.rds")
# # 
# # #Counts for Temporal cortex
# # 
# saveRDS(counts_TC_ROSMAP, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/TC/ROSMAP_RNAseq_rawcounts_TC.rds")

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

##################################### MSBB ##################################### 

metadata_MSBB <- vroom::vroom(file = "/datos/rosmap/data_by_counts/metadata/RNAseq_Harmonization_MSBB_combined_metadata.csv")
metadata_MSBB <- metadata_MSBB %>% rename(ceradsc = CERAD, braaksc = Braak)

#Filter MSSB metadata 
#Based on CDR classification, subjects are grouped as no cognitive deficits (CDR = 0),
#questionable dementia (CDR = 0.5), mild dementia (CDR = 1.0), moderate dementia (CDR = 2.0), and severe to terminal dementia (CDR = 3.0–5.0). 

metadata_MSBB <- metadata_MSBB %>%
  filter(assay == "rnaSeq") %>%  
  mutate(NIA_reagan_ADLikelihood = case_when(         
    (ceradsc == 1 & (braaksc == 5 | braaksc ==  6)) ~ "3", #High likelihood
    (ceradsc == 2 & (braaksc == 3 | braaksc == 4)) ~ "2", #Intermediate likelihood
    (ceradsc == 3 & (braaksc == 1 | braaksc == 2)) ~ "1", #Low likelihood
    ceradsc == 4 ~ "0",  #No AD (0)
    TRUE ~ NA_character_  # Handle no-specified cases
  )) %>% 
  mutate(dicho_NIA_reagan = case_when(
    (NIA_reagan_ADLikelihood == "0" | NIA_reagan_ADLikelihood == "1") ~ "0", #no AD pathology
    (NIA_reagan_ADLikelihood == "2" | NIA_reagan_ADLikelihood == "3") ~ "1"  #AD pathology
  )) %>%
  mutate(CDR_dicho = case_when(
    CDR < 2  ~ "no_AD", 
    CDR >= 0 ~ "dementia", 
    TRUE ~ NA_character_ 
  )) %>%
  mutate(is_resilient = case_when(
    CDR %in% c(3.0, 4.0, 5.0) & (braaksc != 0 & (ceradsc == 1 | ceradsc ==2)) ~ "resilient", 
    TRUE ~ NA_character_ 
  )) %>% 
  mutate(is_AD = case_when(
    CDR %in% c(3.0, 4.0, 5.0) & (braaksc != 0 & (ceradsc == 1 | ceradsc ==2)) ~ "resilient",
    CDR == 0 | ceradsc == 4 ~ "noAD",
    (CDR %in% c(3.0, 4.0, 5.0) &  (ceradsc == 1 | ceradsc ==2)) ~ "AD",
    CDR %in% c(1.0, 2.0) ~ "MCI",
    TRUE ~ NA_character_
  ))
table(metadata_MSBB$CDR_dicho, useNA = "ifany")

table(metadata_MSBB$NIA_reagan_ADLikelihood, useNA = "ifany")

table(metadata_MSBB$is_AD, useNA = "ifany")

#Define metadata by brain region --- ---

tissues_MSBB <- unique(metadata_MSBB$tissue)

metadata_tissue_MSBB <- lapply(tissues_MSBB,
                               function(tissues_MSBB) summarize_by_tissue(metadata_MSBB, tissues_MSBB))

names(metadata_tissue_MSBB) <- tissues_MSBB

# Metadata dim superior temporal gyrus:  334 36 
# NIA-Reagan diagnosis table:
#   
#   0    1 <NA> 
#   49   20  265 
# Metadata dim parahippocampal gyrus:  315 36 
# NIA-Reagan diagnosis table:
#   
#   0    1 <NA> 
#   43   16  256 
# Metadata dim frontal pole:  310 36 
# NIA-Reagan diagnosis table:
#   
#   0    1 <NA> 
#   47   19  244 
# Metadata dim inferior frontal gyrus:  308 36 
# NIA-Reagan diagnosis table:
#   
#   0    1 <NA> 
#   45   19  244 
# Metadata dim prefrontal cortex:  15 36 
# NIA-Reagan diagnosis table:
#   
#   1 <NA> 
#   2   13 

#Save metadata --- --- 

# #STG
# vroom::vroom_write(metadata_tissue_MSBB[["superior temporal gyrus"]],
#                    file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/STG/metadata/MSBB_RNAseq_metadata_STG.txt")
# #PHG
# vroom::vroom_write(metadata_tissue_MSBB[["parahippocampal gyrus"]],
#                    file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/PHG/MSBB_RNAseq_metadata_PHG.txt")
# 
# #FP
# 
# vroom::vroom_write(metadata_tissue_MSBB[["frontal pole"]],
#                    file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/FP/MSBB_RNAseq_metadata_FP.txt")
# 
# #inferior frontal gyrus IFG
# 
# vroom::vroom_write(metadata_tissue_MSBB[["inferior frontal gyrus"]],
#                    file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/IFG/MSBB_RNAseq_metadata_IFG.txt")
# #PFC
# 
# vroom::vroom_write(metadata_tissue_MSBB[["prefrontal cortex"]],
#                    file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/PFC/MSBB_RNAseq_metadata_PFC.txt")

#Read expression data --- --- 

counts_MSSB <-  vroom::vroom(file = "/datos/rosmap/data_by_counts/MSBB_counts/MSBB_gene_all_counts_matrix_clean.txt")

#Stratify data by brain region --- --- 

#Counts "superior temporal gyrus" (STG)
counts_STG_MSBB <- counts_MSSB[, (colnames(counts_MSSB) %in% metadata_tissue_MSBB[[1]]$specimenID)] %>%
  mutate(counts_MSSB[1], .before = 1)
dim(counts_STG_MSBB)

#Counts from parahippocampal gyrus (PHCG)

counts_PHG_MSBB <- counts_MSSB[, (colnames(counts_MSSB) %in% metadata_tissue_MSBB[[2]]$specimenID)] %>%
  mutate(counts_MSSB[1], .before = 1)
dim(counts_PHG_MSBB)

#Counts for frontal pole (FP)
counts_FP_MSBB <- counts_MSSB[, (colnames(counts_MSSB) %in% unique(metadata_tissue_MSBB[[3]]$specimenID))] %>%
  mutate(counts_MSSB[1], .before = 1)
dim(counts_FP_MSBB)

#Counts for inferior frontal gyrus (IFG)
counts_IFG_MSBB <- counts_MSSB[, (colnames(counts_MSSB) %in% unique(metadata_tissue_MSBB[[4]]$specimenID))] %>%
  mutate(counts_MSSB[1], .before = 1)
dim(counts_IFG_MSBB)

#Counts for prefrontal cortex (PFC)
counts_PFC_MSBB <- counts_MSSB[, c(colnames(counts_MSSB) %in% metadata_tissue_MSBB[[5]]$specimenID)] %>% 
  mutate(counts_MSSB[1], .before = 1)
dim(counts_PFC_MSBB)

#Save count data for MSBB ---- ---
# 
# # Counts "superior temporal gyrus"
# saveRDS(counts_STG_MSBB, file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/STG/MSBB_RNAseq_rawcounts_STG.rds")
#  
# # #Counts for parahippocampal gyrus (PHG)
#  
# saveRDS(counts_PHG_MSBB, file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/PHG/MSBB_RNAseq_rawcounts_PHG.rds")
#  
# # #Counts for frontal pole (FP)
# 
# saveRDS(counts_FP_MSBB, file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/FP/MSBB_RNAseq_rawcounts_FP.rds")
# 
# # Counts for inferior frontal gyrus (IFG)
# 
# saveRDS(counts_IFG_MSBB, file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/IFG/MSBB_RNAseq_rawcounts_IFG.rds")
#  
# # #Counts for prefrontal cortex (PFC)
#  
# saveRDS(counts_PFC_MSBB, file = "/datos/rosmap/data_by_counts/MSBB_counts/counts_by_tissue/PFC/MSBB_RNAseq_rawcounts_counts_PFC_MSBB.rds")

#Summarize MSBB --- ---

MSBB_cerad.p <- summarise_metadata(
  data = metadata_MSBB,
  group_var = "tissue",
  fill_var = "ceradsc",
  title = "CERAD proportion by tissue - MSBB"
)

MSBB_cerad.p

#By AD

MSBB_AD.p <- summarise_metadata(
  data = metadata_MSBB,
  group_var = "tissue",
  fill_var = "is_AD",
  title = "AD proportion by tissue - MSBB"
)

MSBB_AD.p

##################################### Mayo Clinic ##################################### 

metadata_Mayo <- vroom::vroom(file = "/datos/rosmap/data_by_counts/metadata/RNAseq_Harmonization_Mayo_combined_metadata.csv")

tissues_Mayo <- unique(metadata_Mayo$tissue)

metadata_tissue_Mayo <- lapply(tissues_Mayo,
                               function(tissues_Mayo) summarize_by_tissue(metadata_Mayo, tissues_Mayo))

names(metadata_tissue_Mayo) <- tissues_Mayo

#Metatada Mayo save --- ---

# vroom::vroom_write(metadata_tissue_Mayo[["cerebellum"]],
#        file = "/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/cerebellum/Mayo_RNAseq_metadata_CRB.txt")
# 
# 
# vroom::vroom_write(metadata_tissue_Mayo[["temporal cortex"]],
#       file = "/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/TC/Mayo_RNAseq_metadata_TC.txt")

#Get expression data --- ---

counts_Mayo <-  vroom::vroom(file = "/datos/rosmap/data_by_counts/Mayo_counts/Mayo_gene_all_counts_matrix_clean.txt")

#Stratify by brain region --- ---

#counts for cerebellum
counts_CRB_Mayo <- counts_Mayo[, (colnames(counts_Mayo) %in% metadata_tissue_Mayo[[1]]$specimenID)] %>% 
  mutate(counts_Mayo[1], .before = 1)
dim(counts_CRB_Mayo)

#Counts from Temporal cortex

counts_TC_Mayo <- counts_Mayo[, (colnames(counts_Mayo) %in% metadata_tissue_Mayo[[2]]$specimenID)] %>% 
  mutate(counts_Mayo[1], .before = 1)
dim(counts_TC_Mayo)

#Save count data for Mayo ---- ---

# Counts for cerebellum
# 
# saveRDS(counts_CRB_Mayo, file = "/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/cerebellum/Mayo_RNAseq_rawcounts_CRB.rds")
# 
# # #Counts for Temporal cortex
# 
# saveRDS(counts_TC_Mayo, file = "/datos/rosmap/data_by_counts/Mayo_counts/counts_by_tissue/TC/Mayo_RNAseq_rawcounts_TC.rds")

#Summarize Mayo --- ---

dx_Mayo.p <- summarise_metadata(
  data = metadata_Mayo,
  group_var = "tissue",
  fill_var = "diagnosis",
  title = "dx proportion by tissue - Mayo"
)

dx_Mayo.p

#Sumarize everything --- 

grid <- grid.arrange(rosmap_nia_reagan.p, rosmap_ceradsc.p, rosmap_dx.p, MSBB_cerad.p, MSBB_AD.p, dx_Mayo.p, ncol = 2)

#Save 

ggsave(filename = "proportions_diagnosis_lo.png", plot = grid, 
       device = "png", width = 17, height = 10, 
       dpi = 300
)

#END