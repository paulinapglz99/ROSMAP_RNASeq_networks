#
#4.my_network_cut_minimal_single_component_heuristic.R

#paulinapglz.99@gmail.com

#NOTE for this script
#There are several approaches to cut the number of interactions in a network. In this pipeline, we look at three heuristics:

#1. Mutual Information cut-off according to the MI distribution of a network
#2. The single component minimum graph. the minimum one-component graph is obtained by making iterative cuts to find the interaction cut to obtain a single-component graph 
#topological heuristic that looks for a graph that maintains a unique component with as few interactions as possible. 
#3.The percentile-cut approach is explained as a method to obtain the percentiles with edges with the highest
#mutual information (top edges). 

#You can choose an approach depending on what you want to analyze or compare them. I recommend to get your data and only use one 
#heuristic at a time

#Here we code the 

########--- 2nd HEURISTIC: minimal single component --- ########

#Libraries  --- --- 

pacman::p_load('tidyverse', 
               'igraph', 
               'ggplot2')

#Get data --- ---


#Get data --- ---

full_edgelist <- vroom::vroom(file = '/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.tsv')

#Mutual information cuts -------

#Note: The MI cut heuristic for finding a minimum graph of a single component involves making empirical
#iterative cuts until the critical percolation point is identified, i.e., the point at which sufficient 
#connections are lost for the system to cease to be a single module.  

#You can try using slices starting at the 'tail' of the histogram and moving upward.

#MI>0.2

matrix_MI_0.2 <- full_edgelist %>% filter(MI >= 0.2)
dim(matrix_MI_0.2)
# [1] 18318810        3

#Build graph with MI long matrix

graph0.2 <- graph_from_data_frame(matrix_MI_0.2, directed = F)

#Number of components

components(graph0.2)

#$no [1] 5

#MI>0.1

matrix_MI_0.1 <- full_edgelist %>% filter(MI >= 0.1)
dim(matrix_MI_0.1)
# [1] 62423892        3

#Build graph with MI long matrix

graph0.1 <- graph_from_data_frame(matrix_MI_0.1, directed = F)  

#Number of components

components(graph0.1)

####

#Iterative cut --- ---

#This function takes a sequence of slices and applies it to an edgelist, to perform a
#filtering and obtain the edges between these values. Then it calculates the components 
#according to the graphs of each cut. The result is a table with the components according to the MI cuts.

generate_network_and_components <- function(edgelist, cutoff) {
  
  matrix_MI_filtrado <- full_edgelist %>%
    filter(MI >= cutoff)
  
  graph <- graph_from_data_frame(matrix_MI_filtrado, directed = FALSE)
  num_components <- components(graph)$no
  
  result <- tibble(cutoff_value = cutoff, num_components = num_components)
  return(result)
}

# Iterate over cutoff values from 0.1 to 0.8 in 0.01 increments.

cutoff_values <- seq(0.1, 0.2, by = 0.01) #add here the two cut-off points between which the critical percussive point is located

# Create a list to store the results

results_list <- lapply(cutoff_values, function(cutoff) {
  generate_network_and_components(full_edgelist, cutoff)
})

# Table with the cutoff value and number of components of the net

cutoff_table <- bind_rows(results_list)

#cutoff_value num_components
#<dbl>          <dbl>
#  1         0.1               1
#2         0.11              1
#3         0.12              1
#4         0.13              1
#5         0.14              1 <--------
#6         0.15              2 <--------
#7         0.16              2
#8         0.17              2
#9         0.18              5
#10         0.19              4
#11         0.2               5

#The network ceases to be a single component between the MI cut of 0.14 and 0.15. 
#If finer cuts are needed, the focus can be changed in the cutoff_values line e.g
#cutoff_values <- seq(0.14, 0.5, by = 0.001)

#
# Iterate over cutoff values from 0.1 to 0.8 in 0.01 increments.

cutoff_values <- seq(0.14, 0.15, by = 0.001) #add here the two cut-off points between which the critical percussive point is located

# Create a list to store the results

results_list <- lapply(cutoff_values, function(cutoff) {
  generate_network_and_components(full_edgelist, cutoff)
})

# Table with the cutoff value and number of components of the net

cutoff_table <- bind_rows(results_list)

#  cutoff_value num_components
#  1        0.14               1
#2        0.141              1
#3        0.142              1
#4        0.143              1
#5        0.144              1 <--------
#6        0.145              2
#7        0.146              2
#8        0.147              2
#9        0.148              2
#10        0.149              2
#11        0.15               2

#Build graph with MI leading to one component

matrix_MI_sc <- full_edgelist %>% filter(MI >= 0.144)

#Histogram 

png("min_sin_comp_MI_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(matrix_MI_sc$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Mutual information histogram",   # Título principal
     sub = "patients with AD, minimal single component",   # Subtítulo
     xlab = "Mutual Information",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Generate graph

min_sin_comp_graph <- graph_from_data_frame(matrix_MI_sc, directed = FALSE)

#Graph info

summary(min_sin_comp_graph)

#IGRAPH f79c29c UN-- 13350 35620678 -- 
#  + attr: name (v/c), MI (e/x)

#number of components

components(min_sin_comp_graph)

#Degree by node 

degreee_by_node_min_sin_g <- degree(min_sin_comp_graph) %>% as.data.frame()

png("min_sin_comp_MI_degree_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(degreee_by_node_min_sin_g$.), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Degree by node histogram",   # Título principal
     sub = "patients with AD, minimal single component",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()


#Coreness by node

coreness_by_node_min_sin_g <- coreness(min_sin_comp_graph) %>% as.data.frame()


png("min_sin_comp_MI_coreness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_min_sin_g$.), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Coreness by node histogram",   # Título principal
     sub = "patients with AD, minimal single component",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()


#edge betweenneess by node

betweeness_by_node_min_sin_g <-betweenness(min_sin_comp_graph, directed = F) %>% as.data.frame()


png("min_sin_comp_MI_betweeness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_min_sin_g$.), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Betweeness by node histogram",   # Título principal
     sub = "patients with AD, minimal single component",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Plot graph, slow

png('min_sin_comp_graph_MI_0.144_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.png')
plot(min_sin_comp_graph, vertex.color="green")
dev.off()

#write graphs in graphml --- ---

#write_graph(min_sin_comp_graph,
#           file = '/datos/rosmap/cut_by_MI_one_component/min_sin_comp_graph_MI_0.144_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
#          'graphml') #'edgelist.txt' / '.graphml'

#END