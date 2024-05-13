#
#5.network_coreness_analysis.R
#This script makes cuts by network coreness
#paulinapglz.99@gmail.com

#Libraries ----

pacman::p_load('tidyverse', 
               'igraph')

#Read graph ----

graph <- read.graph(file = "noMCI_17112023_graph_MI0.532.graphml",
                    format = "graphml")

#Cut by coreness ----

coreness <- coreness(graph) %>% 
  as.data.frame() %>% 
  rename("core_by_node" = ".") 

coreness$gene <- rownames(coreness)

rownames(coreness) <- NULL

#Histogram of the coreness distribution

hist_coreness <- ggplot(coreness, aes(x = core_by_node)) +
  geom_histogram(fill = "skyblue", color = "white") +
  labs(title = "Histograma de coreness por nodo",
       subtitle = "Para sujetos ROSMAP sin demencia,despues del corte de percolacion de MI",
       x = "Coreness by node", y = "Frecuencia") +
  scale_x_continuous(breaks = seq(0, max(coreness$core_by_node), by = 10)) +
  theme_light()

# Filtering by coreness -----

coreness_filter <- coreness %>% 
  filter(core_by_node >= 20) #filtrar para tener un grafo con nodos con un core max de __

#Histogram of filtered coreness----

hist_coreness_filter <- ggplot(coreness_filter, aes(x = core_by_node)) +
  geom_histogram(fill = "skyblue", color = "white") +
  labs(title = "Histograma de coreness por nodo con filtro",
       subtitle = "Para sujetos ROSMAP sin demencia, despues del corte de percolacion de MI",
       x = "Coreness by node", y = "Frecuencia") +
  scale_x_continuous(breaks = seq(0, max(coreness_filter$core_by_node), by = 10)) +
  theme_light()

#Now the subgraph ----

# Create a vector with the name of nodes I want to subset

coreness_filter_v <- coreness_filter$gene

#Making a subgraph only with the vertices of the core we want

subgraph_kcore <- induced_subgraph(graph,
                                   vids = coreness_filter_v)

#Plot 

subraph_kcore_filter.p <- plot(subgraph_kcore)

#write graphs in graphml ----

#write_graph(subgraph_kcore, 
#           file = 'coexpre_noMCI_18112023_graph_MI0.5_coreness_20.graphml', 
#           'graphml') #'edgelist.txt' / '.graphml'

#Automatize to have other kcore graphs -----

#####################################################################################################

# Function to filter and save the subgraph

coreness_maker <- function(coreness, graph, n) {
  
  # Filter by the value of core_by_node
  coreness_filter <- coreness %>% 
    filter(core_by_node >= n)
  
  # Obtain the vector with the names of nodes
  coreness_filter_v <- coreness_filter$gene
  
  # Create the subgraph
  subgraph_kcore <- induced_subgraph(graph, vids = coreness_filter_v)
  
  # Write the network in graphml format
  filename <- paste0('coexpre_noMCI_18112023_graph_MI0.5_coreness_', n, '.graphml')
  write_graph(subgraph_kcore,
              file = filename, 
              format = 'graphml')
}

# 'coreness' and 'graph' are initial data

# Iterate over core_by_node values from 25 to 145 with increments of 10

for (i in seq(25, 145, by = 10)) {
  coreness_maker(coreness, graph, i)
}
