#
#This script takes a matrix of NxN and makes iterative cuts of mutual information by percentiles.
#Then, it calculates characteristics such as number of nodes, number of vertices, clustering coefficient, etc. for each graph.
#Generates scatter plots describing the topology of each of the networks, comparing healthy and diseased subjects.

#Libraries --- ---

pacman::p_load('tidyverse', 
               'igraph')

#If needed, read adjacency matrix and pivot it to build an edgelist --- ---

#matrix <- read_rds('/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_zero.rds')

#Pivot

#this gives a table of connections between genes (edgelist)

#full_edgelist <- as.data.frame(matrix) %>%
#  rownames_to_column(var = "ensembl_gene_id") %>%
#  pivot_longer(-ensembl_gene_id, names_to = "gene_to", 
#              values_to = "MI")
#dim(full_edgelist)
# [1] 223532401         3

#full_edgelist$MI <- full_edgelist$MI %>% as.numeric()

#vroom::vroom_write(full_edgelist, file = '/datos/rosmap/coexpre_matrix/full_net_ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_edgelist.tsv')

#Get edgelist --- ---

full_edgelist <- vroom::vroom(file = '/datos/rosmap/coexpre_matrix/full_net_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.tsv') %>% as.data.frame()
full_edgelist <- full_edgelist[1:100,]

#Declare functions --- ---

#Function to calculate p-value

pvalue<-function(MI, n=100){
  alfa = 1.062
  beta = -48.7
  gamma = -0.634
  p = exp(alfa -MI*(-beta + (-gamma * n)))
  return(p)
}

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
  clusters_no <- components(graph)$no
  # Clustering coefficient
  clustering_coefficient <- transitivity(graph, type = 'undirected')
  #Max and min weigth
  
  max_weight <- max(E(graph)$MI)
  min_weight <- min(E(graph)$MI)
  
  #MI p-value
  
  p_value <- pvalue(MI = E(graph)$MI)
  
  # Output with metrics
  data.frame(
    length_v = length_v,
    length_E = length_E,
    clusters_no = clusters_no,
    clustering_coefficient = clustering_coefficient,
    max_weight = max_weight,
    min_weight = min_weight,
    p_value = p_value
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
  
#Save table
  
#vroom::vroom_write(metric_table, file = '/datos/rosmap/cuts_by_MI/AD_graphs/metrics_percentiles_AD_ROSMAP_RNAseq_MutualInfo_NIA_Reagan_dicho.txt')

#END