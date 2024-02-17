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
matrix2<- vroom::vroom('/datos/rosmap/coexpre_matrix/matriz_coexpre_noMCI_11052023_zero.txt')
#Pivot  ----
#this gives a table of connections between genes

MI_matriz_long <- matrix %>%
  rownames_to_column(var = "gene") %>%
  pivot_longer(-gene, names_to = "gene_to", values_to = "MI")

########--- 

#Mutual information cuts -------

#MI>0.5

matrix_MI_0.5 <- MI_matriz_long %>% 
 filter(MI >= 0.5)
dim(matrix_MI_0.5)

#Plotting

ggplot(matrix_MI_0.5, aes(x = MI)) +
  geom_bar(fill = "skyblue", color = "white") +
  labs(title = "mutual information histogram",
       subtitle = "patients with AD",
       x = "Mutual Information", y = "freq") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))


#MI>0.8

matrix_MI_0.8 <- matrix_MI %>% 
  filter(mi >= 0.8)

hist_MI0_8 <- ggplot(matrix_MI_0.8, aes(x = mi)) +
  geom_histogram(fill = "skyblue", color = "white") +
  labs(title = "Histograma de Informacion Mutua",
       subtitle = "Para sujetos ROSMAP con cogdx de AD, corte MI > 0.8",
       x = "Informacion Mutua", y = "Frecuencia") +
  scale_y_continuous(labels = function(x) format(x, scientific = FALSE))

#Build net from pivot data ----

graph <- graph_from_data_frame(matrix_MI, 
                                 directed = F)  #full net

graph0.5 <- graph_from_data_frame(matrix_MI_0.5, 
                               directed = F)    #cut on 0.5

graph0.8 <- graph_from_data_frame(matrix_MI_0.8, 
                                     directed = F) #cut on 0.8

#write graphs in graphml ----

#write_graph(graphAD_0.5, 
#           file = 'coexpre_noMCI_16112023_graph_MI0.5.graphml', 
#           'graphml') #'edgelist.txt' / '.graphml'

#Network generalities ---

#Summary 

summary(graph0.5) #resumen

#Number of components 

components(graphAD0.8) #saber cuantos componentes tiene

#Degree by node 

degree(graph0.5) #devuelve un vector donde para cada nodo tengo el valor de grado, es decir el numero de L

#edge betweenneess

E(graph0.5)$betweenneessAD <- betweenness(graph0.5,
            directed = F)

#Iterative cutt --------

# Función para construir la red y obtener el número de componentes
generate_network_and_components <- function(matrix_MI, cutoff) {
  matrix_MI_filtrado <- matrix_MI %>%
    filter(mi >= cutoff)
  
  graph <- graph_from_data_frame(matrix_MI_filtrado, directed = FALSE)
  num_components <- components(graph)$no
  
  result <- tibble(cutoff_value = cutoff, num_components = num_components)
  return(result)
}

# Iterar sobre los valores de corte de 0.01 a 0.8 en incrementos de 0.01

cutoff_values <- seq(0.5, 0.8, by = 0.01)

# Crear una lista para almacenar los resultados

results_list <- lapply(cutoff_values, function(cutoff) {
  generate_network_and_components(matrix_MI, cutoff)
})

# Table with the cutoff value and number of components of the net

cutoff_table <- bind_rows(results_list)

#Now I want to use the same for a fineblanking between 0.53 and 0.54

fineblanking_cutoff_values <- seq(0.53, 0.54, by = 0.001)

# Crear una lista para almacenar los resultados

fineblanking_list <- lapply(fineblanking_cutoff_values, function(cutoff) {
  generate_network_and_components(matrix_MI, cutoff)
})

# Table with the cutoff value and number of components of the net

fineblanking_list_table <- bind_rows(fineblanking_list)

#Cut will be at 0.532

#Construct graph ----

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