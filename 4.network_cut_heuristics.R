#
#4.my_network_analysis.R
#Script that analyzes co-expression matrices (adjacency matrices) and converts them 
#into graph format with igraph, using heuristics for a cut in the MI, for later visualization and analysis.

#paulinapglz.99@gmail.com

#NOTE for this script

#This script is designed to filter the interactions of a complete graph using three heuristics:
#1. Mutual Information cut-off according to the MI distribution of a network
#2. The single component minimum graph. The minimal graph approach is a 
#topological heuristic that looks for a graph that maintains a unique component with as few interactions as possible. 
#3.The percentile-cut approach is explained as a method to obtain the percentiles with edges with the highest
#mutual information (top edges). 

#You can choose an approach depending on what you want to analyze or compare them. I recommend to get your data and only use one 
#heuristic at a time

#Libraries  --- --- 

pacman::p_load('tidyverse', 
               'igraph', 
               'ggplot2')

#Read adjacency matrix --- ---

matrix <- read_rds('/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_zero.rds') %>% as.data.frame()

#Pivot  ----

#this gives a table of connections between genes (edgelist), needed for the three approaches

full_edgelist <- matrix %>%
  rownames_to_column(var = "ensembl_gene_id") %>%
  pivot_longer(-ensembl_gene_id, names_to = "gene_to", 
               values_to = "MI")
dim(full_edgelist)
# [1] 223532401         3

########--- 1st HEURISTIC: cutoff according to MI distribution --- ########

#Plot MI distribution

#This histogram is necessary to decide on the following mutual information cutoffs 
png("full_net_MI_histogram_ROSMAP_RNA_seq.png")
hist(as.numeric(full_edgelist$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Mutual information histogram",   # Título principal
     sub = "patients with AD, no MI cut",   # Subtítulo
     xlab = "Mutual Information",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#If we look at the histogram, we will notice that there is a point in the MI distribution where we see a dip, 
#which is called the 'tail of the histogram'. The first approach is to make a cut in this tail. 

matrix_MI_tailcut<- full_edgelist %>% 
  filter(MI >= 0.5)   
dim(matrix_MI_tailcut)
#[1] 440646      3

#Build graph

tailcut_graph <-  graph_from_data_frame(matrix_MI_tailcut, directed = F)

#Components

components(tailcut_graph)

#Save graph

#write graphs in graphml --- ---

#write_graph(tailcut_graph,
#            file = '/datos/rosmap/cut_by_MI_one_component/tailcut_graph_MI_0.5_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
#            'graphml') #'edgelist.txt' / '.graphml'


########--- 2nd HEURISTIC: minimal single component --- ########

#Mutual information cuts -------

#Note: The MI cut heuristic for finding a minimum graph of a single component involves making empirical
#iterative cuts until the critical percolation point is identified, i.e., the point at which sufficient 
#connections are lost for the system to cease to be a single module.  

#You can try using slices starting at the 'tail' of the histogram and moving upward.
#Here I start at MI>0.3

#MI>0.3

matrix_MI_0.3<- full_edgelist %>% 
  filter(MI >= 0.3)
dim(matrix_MI_0.3)
#[1] 5618420       3

#Build graph with MI long matrix

graph0.3 <- graph_from_data_frame(matrix_MI_0.3, directed = F)    #cut on 0.3

#Number of components

components(graph0.3)

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

min_sin_comp_graph <- graph_from_data_frame(matrix_MI_sc, directed = FALSE)

#Add graph characteristics

#Number of components

min_sin_comp_graph$comp <- components(min_sin_comp_graph)

#Number of edges

min_sin_comp_graph$E <- E(min_sin_comp_graph)

#Name of vertices

min_sin_comp_graph$V <- V(min_sin_comp_graph)

#Degree by node 

min_sin_comp_graph$degree <- degree(min_sin_comp_graph) 

#Coreness by node

min_sin_comp_graph$coreness <- coreness(min_sin_comp_graph)

#edge betweenneess by node

V(min_sin_comp_graph)$betweenneess <- betweenness(min_sin_comp_graph, directed = F)

#write graphs in graphml --- ---

write_graph(min_sin_comp_graph,
           file = '/datos/rosmap/cut_by_MI_one_component/min_sin_comp_graph_MI_0.144_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
          'graphml') #'edgelist.txt' / '.graphml'

########--- 3rd HEURISTIC: percentile cut-off --- ########

#The P percentile indicates the value below which a specific percentage of the data falls.

# Percentile 90 --- --- 

percentile_90 <- quantile(as.numeric(full_edgelist$MI), 0.9)
#     90%  of observations are above 0.1830432 

percentile_90_table <- subset(full_edgelist, as.numeric(MI) > percentile_90)
dim(percentile_90_table)
#[1] 22353240        3

#write graph

graph_percentile90 <- graph_from_data_frame(percentile_90_table, directed = F)

#analyze

components(graph_percentile90)

#write graphs --- ---

write_graph(graph_percentile90,
            file = '/datos/rosmap/cut_by_MI_one_component/graph_percentile90_MI_0.18_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
            'graphml') #'edgelist.txt' / '.graphml'

# Percentile 80 --- --- 

percentile_80 <- quantile(as.numeric(full_edgelist$MI), 0.8)
#     80%  0.1256503 

percentile_80_table <- subset(full_edgelist, as.numeric(MI) > percentile_80)
dim(percentile_80_table)
#[1] 44706480        3

# graph

graph_percentile80 <- graph_from_data_frame(percentile_80_table, directed = F)

#write graphs --- ---

write_graph(graph_percentile80,
            file = '/datos/rosmap/cut_by_MI_one_component/graph_percentile80_MI_0.126_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
            'graphml') #'edgelist.txt' / '.graphml'


#You can performm further analysis to graphs in cytoscape or with igraph

#Next script is 5.network_coreness_analysis.R 