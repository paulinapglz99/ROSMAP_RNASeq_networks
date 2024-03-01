#
#This script takes a matrix of NxN and makes iterative cuts of mutual information by percentiles.
#Then, it calculates characteristics such as number of nodes, number of vertices, clustering coefficient, etc. for each graph.
#Generates scatter plots describing the topology of each of the networks, comparing healthy and diseased subjects.

#Libraries --- ---

pacman::p_load('tidyverse', 
               'igraph',
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

calculate_percentiles <- function(data, percentiles) {
  result_list <- list()
  
  for (p in percentiles) {
    percentile_value <- quantile(as.numeric(data$MI), p)
    table_subset <- subset(data, as.numeric(MI) > percentile_value)
    
    result_list[[paste0("percentile_", gsub("\\.", "", as.character(p*100)))]] <- list(
      percentile_value = percentile_value,
      table = table_subset
    )
  }
  
  return(result_list)
}

#Calculate the percentile tables to store them in a list.

#First assign the percentiles, already divided by 10. e.g, if you want percentile 70.5, write 0.705

percentiles <- c(0.99, 0.3, 0.4)

percentile_tables <- calculate_percentiles(edgelist_prueba, percentiles)

#Generate a list of networks from the percentile tables in the list

generate_graph <- function(tables_edgelists) {
  networks <- list()
  
  for (percentile_name in names(tables_edgelists)) {
    table <- tables_edgelists[[percentile_name]]$table
    graph <- graph_from_data_frame(table, directed = FALSE)
    
    # You can customize the network name according to the percentile.
    name_graph <- paste("graph_", percentile_name, sep = "")
    
    networks[[name_graph]] <- graph
  }
  
  return(networks)
}

#Generate networks results

results_networks <- generate_graph(tables_edgelists = percentile_tables)

###########HASTA AKA TODO CHIDO



############################################


  length_v <- length(V(grafo))
  length_E <- length(E(grafo))
  clusters_no <- clusters(grafo)$no
  clustering_coefficient <- transitivity(grafo)
  max_weight <- max(E(grafo)$weight)
  min_weight <- min(E(grafo)$weight)

  length(V(grafito))
  length(E(grafito))
  plot(grafito)
  
  grafete <- results_networks[[2]]
  
  length(V(grafete))
  length(E(grafete))
  plot(grafete)
  
  calcular_caracteristicas <- function(grafo) {
    length_v <- length(V(grafo))
    length_E <- length(E(grafo))
    clusters_no <- clusters(grafo)$no
    clusteringcoefficient <- transitivity(grafo)
            
    data.frame(
      length_v = length_v,
      length_E = length_E,
      clusters = clusters_no,
      clusteringcoefficient = clusteringcoefficient,
              )
  }
  
  # Aplica la función a cada grafo en la lista
  
  resultados_caracteristicas <- lapply(results_networks, calcular_caracteristicas)

  
  # Combina los resultados en una lista de tablas
  tabla_resultados <- do.call(rbind, Map(cbind, nombres_grafos = names(grafos_resultados), resultados_caracteristicas))
  
  # Visualiza la tabla
  print(tabla_resultados)