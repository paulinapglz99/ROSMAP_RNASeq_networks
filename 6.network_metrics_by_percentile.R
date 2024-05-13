#
#6.network_metrics_by_percentile.R
#This script takes an edgelist and makes iterative cuts of mutual information by percentiles.
#Then, it calculates characteristics such as number of nodes, number of vertices, clustering coefficient, etc. for each graph.
#The next script, 4.1.plto_network_characteristics.R Generates scatter plots describing the topology of each of the networks, comparing healthy and diseased subjects.

#Libraries --- ---

pacman::p_load('tidyverse', 
               'igraph')

#Set timer 

tempus <- Sys.time()

#Get edgelist --- ---

full_edgelist <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.txt')

#Declare functions --- ---

#Function to obtain the edgelists per percentile

calculate_percentiles <- function(data, percentiles) {
  result_list <- list()
  
  for (p in percentiles) {
    percentile_value <- quantile(as.numeric(data$MI), p)
    table_subset <- data %>% 
      filter(as.numeric(MI) > percentile_value)    
    result_list[[paste0("percentile_", gsub("\\.", "", as.character(p*100)))]] <- list(
      percentile_value = percentile_value,
      table = table_subset
    )
  }
  
  return(result_list)
}

#Function that generates a list of networks from the percentile tables in the list

generate_graph <- function(tables_edgelists) {
  networks <- list()
  
  for (percentile_name in names(tables_edgelists)) {
    table <- tables_edgelists[[percentile_name]]$table
    graph <- graph_from_data_frame(table, directed = FALSE)
    # Assigning MI values to the edges
    E(graph)$MI <- table$MI
    
    # You can customize the network name according to the percentile.
    name_graph <- paste("graph_", percentile_name, sep = "")
    
    networks[[name_graph]] <- graph
  }
  
  return(networks)
}

#Function that calculates the metrics for each network

calculate_metrics <- function(graph) {
  # Number of vertices
  length_v <- length(V(graph))
  # Number of edges
  length_E <- length(E(graph))
  # Number of clusters
  components_no <- components(graph)$no
  # Clustering coefficient
  clustering_coefficient <- transitivity(graph, type = 'undirected')
  # Max and min MI 
  max_MI <- max(E(graph)$MI)
  min_MI <- min(E(graph)$MI)
  
  #Finding the size of the biggest connected component
  max_comp_size <- max(components(graph)$csize)
  #Percentage of genes in larger components
  percentage_genes_in_larger_component <- max_comp_size / vcount(graph) * 100
  
  #Calculate the number of communities derived from infomap algorithm
  cluster_infomap <- cluster_infomap(graph)
  membership_infomap <- membership(cluster_infomap)
  
  #Extract list of nodes by community
  infomap_nodes_by_community <- split(V(graph)$name, membership_infomap)
  no_cluster_infomap <- length(infomap_nodes_by_community)
  
  # Output with metrics
  data.frame(
    length_v = length_v,
    length_E = length_E,
    components_no = components_no,
    clustering_coefficient = clustering_coefficient,
    max_MI = max_MI,
    min_MI = min_MI,
    no_cluster_infomap = no_cluster_infomap,
    max_comp_size = max_comp_size, 
    percentage_genes_in_larger_component = percentage_genes_in_larger_component
  )
}

#Applying functions --- ---

# First assign the percentiles, already divided by 100. e.g, if you want percentile 70.5, write 0.705

percentiles <- c(0.999999, 0.99999, 0.9999, 0.999, 0.99, 0.98, 0.9, 0.8)

# Calculate percentile tables

percentile_tables <- calculate_percentiles(full_edgelist, percentiles)

# Generate networks results

results_networks <- generate_graph(tables_edgelists = percentile_tables)

# Calculate metrics for each graph

results_metrics <- lapply(results_networks, calculate_metrics)

# Combine the results in a table

metric_table <- do.call(rbind.data.frame, results_metrics) 
metric_table <- metric_table %>% 
  mutate(percentile_no = rownames(metric_table), .before = 1)

# Show the table

print(metric_table)

# Print time elapsed

print(Sys.time() - tempus)

#Save table

#vroom::vroom_write(metric_table, file = "/datos/rosmap/cuts_by_heuristics/AD_graphs/metrics_percentiles_normalizedMI_AD_ROSMAP_RNAseq_MutualInfo_NIA_Reagan_dicho.txt")

#If you want to save graphs

graph_to_save <- results_networks[[3]] #Indicate the graph to save by the index

#Save in graphml format

write_graph(graph_to_save, file = '/datos/rosmap/graphs/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml',
            format = "graphml")

#Save edgelist

edgelist_to_save <- percentile_tables[[3]] %>% as.data.frame() #Indicate the graph to save by the index
edgelist_to_save<- edgelist_to_save[-1]

#Save graph in edgelist format

#vroom::vroom_write(edgelist_to_save,  file = '/datos/rosmap/cuts_by_MI/noAD_graphs/percentile99.99_ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_edgelist.tsv')

#Obtain adjacency matrix

adj_matrix <- as_adjacency_matrix(results_networks[[3]], sparse = FALSE, attr = "MI") %>% as.data.frame()

vroom::vroom_write(adj_matrix, file = '/datos/rosmap/cuts_by_heuristics/noAD_graphs/percentile99.99_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_we_adjacency_matrix.tsv')

#END