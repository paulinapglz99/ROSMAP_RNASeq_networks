# Librerías ----
pacman::p_load('tidyverse', 
               'igraph')

# Leer gráficos ----
graph_0 <- read.graph(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml",
                      format = "graphml")

graph_1 <- read.graph(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml",
                      format = "graphml")

# Crear lista de gráficos ----
graphList <- list(graph_0 = graph_0, 
                  graph_1 = graph_1)

# Función para procesar cada grafo ----
process_graph <- function(graph, min_coreness, output_filename) {
  # Calcular coreness
  coreness_df <- coreness(graph) %>% 
    as.data.frame() %>% 
    rename("." = "core_by_node")
  
  coreness_df$gene <- rownames(coreness_df)
  rownames(coreness_df) <- NULL
  
  # Histograma de distribución de coreness
  hist_coreness <- ggplot(coreness_df, aes(x = core_by_node)) +
    geom_histogram(fill = "skyblue", color = "white") +
    labs(title = "Histograma de coreness por nodo",
         subtitle = "Distribución general de coreness",
         x = "Coreness by node", y = "Frecuencia") +
    scale_x_continuous(breaks = seq(0, max(coreness_df$core_by_node), by = 10)) +
    theme_light()
  
  print(hist_coreness)
  
  # Filtrar por coreness
  coreness_filter <- coreness_df %>% 
    filter(core_by_node >= min_coreness)
  
  # Histograma después del filtro
  hist_coreness_filter <- ggplot(coreness_filter, aes(x = core_by_node)) +
    geom_histogram(fill = "skyblue", color = "white") +
    labs(title = "Histograma de coreness por nodo con filtro",
         subtitle = paste0("Nodos con coreness >= ", min_coreness),
         x = "Coreness by node", y = "Frecuencia") +
    scale_x_continuous(breaks = seq(0, max(coreness_filter$core_by_node), by = 10)) +
    theme_light()
  
  print(hist_coreness_filter)
  
  # Subgrafo con nodos filtrados
  coreness_filter_v <- coreness_filter$gene
  subgraph_kcore <- induced_subgraph(graph, vids = coreness_filter_v)
  
  # Guardar subgrafo
  write_graph(subgraph_kcore, 
              file = output_filename, 
              format = "graphml")
  
  return(subgraph_kcore)
}

# Aplicar función a la lista de gráficos ----
output_graphs <- lapply(seq_along(graphList), function(i) {
  graph_name <- names(graphList)[i]
  graph <- graphList[[i]]
  output_filename <- paste0(graph_name, "_filtered_coreness.graphml")
  
  process_graph(graph, min_coreness = 20, output_filename = output_filename)
})

# Resultado: lista de subgrafos procesados ----
names(output_graphs) <- names(graphList)
