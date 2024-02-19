#
#4.my_network_analysis.R
#script that analyzes coexpression matrices (adjacency matrices) and converts them 
#into graph format with igraph, for later visualization and analysis.
#paulinapglz.99@gmail.com

#Libraries  --- --- 

pacman::p_load('tidyverse', 
               'igraph', 
               'ggplot2')

#Read adjacency matrix --- ---

matrix <- read_rds('/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_zero.rds') %>% as.data.frame()

#Pivot  ----
#this gives a table of connections between genes

MI_matriz_long <- matrix %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  pivot_longer(-ensembl_gene_id, names_to = "gene_to", 
               values_to = "MI")
dim(MI_matriz_long)

#Plot MI distribution

#This histogram is necessary to decide on the following mutual information cutoffs 

hist(as.numeric(MI_matriz_long$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Mutual information histogram",   # Título principal
     sub = "patients with AD, no MI cut",   # Subtítulo
     xlab = "Mutual Information",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)

########--- 

#Mutual information cuts -------

#MI>0.5

matrix_MI_0.5 <- MI_matriz_long %>% 
 filter(MI >= 0.5)
dim(matrix_MI_0.5)
#[1] 440646      3

#Plot MI distribution

hist(as.numeric(matrix_MI_0.5$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Mutual information histogram",   # Título principal
     sub = "patients with AD, cut on MI >= 0.5",   # Subtítulo
     xlab = "Mutual Information",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)


#MI>0.8

matrix_MI_0.8 <- matrix_MI_0.5 %>% 
  filter(MI >= 0.8)

#Plot MI distribution

hist(as.numeric(matrix_MI_0.8$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Mutual information histogram",   # Título principal
     sub = "patients with AD, cut on MI >= 0.8",   # Subtítulo
     xlab = "Mutual Information",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)

#Build net from pivot data --- ---

graph <- graph_from_data_frame(MI_matriz_long, 
                                 directed = F)  #full net

graph0.5 <- graph_from_data_frame(matrix_MI_0.5, 
                               directed = F)    #cut on 0.5

graph0.8 <- graph_from_data_frame(matrix_MI_0.8, 
                                     directed = F) #cut on 0.8

#write graphs in graphml --- ---

#write_graph(graph0.8,
#           file = '/datos/rosmap/cut_by_MI_one_component/graph0.8_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
#          'graphml') #'edgelist.txt' / '.graphml'

####

#Network generalities --- ---

#Summary 

summary(graph0.5) #resumen

#Number of components 

components(graph0.5) #saber cuantos componentes tiene

#Degree by node 

degree(graph0.5) #devuelve un vector donde para cada nodo tengo el valor de grado, es decir el numero de L

#edge betweenneess

V(graph0.8)$betweenneess <- betweenness(graph0.8, directed = F)

#Iterative cut --- ---

# Function to build the network and obtain the number of components

generate_network_and_components <- function(matrix_MI, cutoff) {
  matrix_MI_filtrado <- matrix_MI %>%
    filter(mi >= cutoff)
  
  graph <- graph_from_data_frame(matrix_MI_filtrado, directed = FALSE)
  num_components <- components(graph)$no
  
  result <- tibble(cutoff_value = cutoff, num_components = num_components)
  return(result)
}

# Iterate over cutoff values from 0.01 to 0.8 in 0.01 increments.

cutoff_values <- seq(0.5, 0.8, by = 0.01)

# Create a list to store the results

results_list <- lapply(cutoff_values, function(cutoff) {
  generate_network_and_components(matrix_MI, cutoff)
})

# Table with the cutoff value and number of components of the net

cutoff_table <- bind_rows(results_list)

#Now I want to use the same for a fineblanking between 0.53 and 0.54

fineblanking_cutoff_values <- seq(0.53, 0.54, by = 0.001)

# Create a list to store the results

fineblanking_list <- lapply(fineblanking_cutoff_values, function(cutoff) {
  generate_network_and_components(matrix_MI, cutoff)
})

# Table with the cutoff value and number of components of the net

fineblanking_list_table <- bind_rows(fineblanking_list)

#Cut will be at 0.532

#Construct graph --- ---

#MI>0.532

matrix_MI_0.532 <- matrix_MI %>% 
  filter(mi >= 0.532)

graphMI_0.532 <- graph_from_data_frame(matrix_MI_0.532, 
                               directed = F)
#Components of new graph

components(graphMI_0.532) # one component $csize of 18503

#save graph ----

#write_graph(graphMI_0.532, 
#           file = 'noMCI_17112023_graph_MI0.532.graphml', 
#           'graphml') #'edgelist.txt' / '.graphml'

#Network generalities, run this only if you want to know ---

#Summary 

summary(graphMI_0.532) #nodes and edges

#Degree by node 

degree(graphMI_0.532) #devuelve un vector donde para cada nodo tengo el valor de grado, es decir el numero de L

#edge betweenneess

E(graph0.5)$betweenneess <- betweenness(graph0.5,
                                          directed = F)

#Next script is 5.network_coreness_analysis.R 