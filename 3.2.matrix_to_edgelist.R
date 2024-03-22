#
#3.2.matrix_to_edgelist.R
#Script that analyzes co-expression matrices (adjacency matrices) and converts them 
#into graph format with igraph, using heuristics for a cut in the MI, for later visualization and analysis.

#Libraries  --- --- 

pacman::p_load('tidyverse', 
               'igraph')

#Read adjacency matrix --- ---

matrix <- read_rds('/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_zero.rds')

#Pivot  ----

#this gives a table of connections between genes (edgelist), needed for the three approaches

full_edgelist <- as.data.frame(matrix) %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  pivot_longer(-ensembl_gene_id, names_to = "gene_to", 
               values_to = "MI")
dim(full_edgelist)
# [1] 223532401         3

#save edgelist for later
#vroom::vroom_write(full_edgelist, file = '/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.tsv')

full_edgelist <- vroom::vroom(file = '/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.tsv')

#END