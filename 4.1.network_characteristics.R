#
#This script takes a matrix of NxN and makes iterative cuts of mutual information by percentiles.
#Then, it calculates characteristics such as number of nodes, number of vertices, clustering coefficient, etc. for each graph.
#Generates scatter plots describing the topology of each of the networks, comparing healthy and diseased subjects.

#Libraries --- ---

pacman::p_load('tidyverse', 
               'igraph',
               'data.table',
               'ggplot2')

#If needed, read adjacency matrix and pivot it to build an edgelist --- ---

#matrix <- read_rds('/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_zero.rds')

#Pivot

#this gives a table of connections between genes (edgelist)

#full_edgelist <- as.data.frame(matrix) %>%
#  rownames_to_column(var = "ensembl_gene_id") %>%
#  pivot_longer(-ensembl_gene_id, names_to = "gene_to", 
#               values_to = "MI")
#dim(full_edgelist)
# [1] 223532401         3

#full_edgelist$MI <- full_edgelist$MI %>% as.numeric()

#vroom::vroom_write(full_edgelist, file = '/datos/rosmap/coexpre_matrix/full_net_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.tsv')

#Get edgelist --- ---

full_edgelist <- vroom::vroom(file = '/datos/rosmap/coexpre_matrix/full_net_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.tsv') %>% as.data.frame()
edgelist_prueba <- full_edgelist[1:100, ]

#Function to obtain the edgelists per percentile

calcular_percentiles <- function(data, percentiles) {
  result_list <- list()
  
  for (p in percentiles) {
    percentile_value <- quantile(as.numeric(data$MI), p)
    table_subset <- subset(data, as.numeric(MI) > percentile_value)
    
    result_list[[paste0("percentile_", gsub("\\.", "", as.character(p)))]] <- list(
      percentile_value = percentile_value,
      table = table_subset
    )
  }
  
  return(result_list)
}

#Calcular las tablas de los percentiles para guardarlos en una lista

percentiles <- c(0.1, 0.2, 0.3)

tablas_percentiles <- calcular_percentiles(edgelist_prueba, percentiles)

#Generar una lista de grafos a partir de las tablas de percentiles en la lista

generar_grafo <- function(tablas_edgelists) {
  grafos <- list()
  
  for (nombre_percentil in names(tablas_edgelists)) {
    tabla <- tablas_edgelists[[nombre_percentil]]$table
    grafo <- graph_from_data_frame(tabla, directed = FALSE)
    
    # Puedes personalizar el nombre del grafo según el percentil
    nombre_grafo <- paste("grafo_", nombre_percentil, sep = "")
    
    grafos[[nombre_grafo]] <- grafo
  }
  
  return(grafos)
}

#

grafos_resultados <- generar_grafo(tablas_edgelists = tablas_percentiles)

#Los saco de la lista

for (nombre_grafo in names(grafos_resultados)) {
  assign(nombre_grafo, grafos_resultados[[nombre_grafo]])
}

###########HASTA AKA TODO CHIDO