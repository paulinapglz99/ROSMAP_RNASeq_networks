#
#5.network_coreness_analysis.R
#This script makes cuts by network coreness
#paulinapglz.99@gmail.com

#Libraries --- ---
pacman::p_load('tidyverse', 
               'igraph')

#Read graphs ----
graph_0 <- read.graph(file = "graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml",
                      format = "graphml")

graph_1 <- read.graph(file = "graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml",
                      format = "graphml")

#Graph list ----
graphList <- list(graph_0 = graph_0, 
                  graph_1 = graph_1)

#Function for processing each network  ----
process_graph <- function(graph, min_coreness, output_filename) {
  #Calculate coreness
  coreness_df <- coreness(graph) %>% 
    as.data.frame() %>% 
    rename("." = "core_by_node")
  
  coreness_df$gene <- rownames(coreness_df)
  rownames(coreness_df) <- NULL
  
  #Histogram of coreness distribution 
  hist_coreness <- ggplot(coreness_df, aes(x = core_by_node)) +
    geom_histogram(fill = "skyblue", color = "white") +
    labs(title = "Coreness per node",
         subtitle = "General distribution of coreness",
         x = "Coreness by node", y = "Frequency") +
    scale_x_continuous(breaks = seq(0, max(coreness_df$core_by_node), by = 10)) +
    theme_light()
  
  print(hist_coreness)
  
  #Filter
  coreness_filter <- coreness_df %>% 
    filter(core_by_node >= min_coreness)
  
  #Histogram after filtering
  hist_coreness_filter <- ggplot(coreness_filter, aes(x = core_by_node)) +
    geom_histogram(fill = "skyblue", color = "white") +
    labs(title = "Histograma de coreness por nodo con filtro",
         subtitle = paste0("Nodos con coreness >= ", min_coreness),
         x = "Coreness by node", y = "Frecuencia") +
    scale_x_continuous(breaks = seq(0, max(coreness_filter$core_by_node), by = 10)) +
    theme_light()
  
  print(hist_coreness_filter)
  
  #Subgraph of filtedes core
  coreness_filter_v <- coreness_filter$gene
  subgraph_kcore <- induced_subgraph(graph, vids = coreness_filter_v)
  
  #Save graph
  write_graph(subgraph_kcore, 
              file = output_filename, 
              format = "graphml")
  
  return(subgraph_kcore)
}

#Apply function to graphics list --- ---
output_graphs <- lapply(seq_along(graphList), function(i) {
  graph_name <- names(graphList)[i]
  graph <- graphList[[i]]
  output_filename <- paste0(graph_name, "_filtered_coreness.graphml")
  
  process_graph(graph, min_coreness = 20, output_filename = output_filename)
})

#list of processed subnetworks ----
names(output_graphs) <- names(graphList)

#END
