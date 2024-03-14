#
#5.1 modularity analysis
#This script takes a list of genes and performsa variety of analysis with igraph and INFOMAP

#Libraries --- ---

pacman::p_load('igraph',
               'igraphdata',
               'dplyr')

#Get data --- ---

graph <- read_graph(file = '~/redesROSMAP/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', 
                    format = 'graphml')

#Analyze modules --- --- 

#Calculate clustering coefficient of whole network

clustering_coefficient <- transitivity(graph, type = 'undirected')

#Calculations on components

components <- components(graph)

#Number of components

no_components <- components(graph)$no

#Calculate community membership

membership <- membership(components)

#List of nodes by community

nodes_by_community <- split(V(graph)$name, membership)

#Extract bigger component

bigger_comp <- nodes_by_community[[1]]

#Save list of genes

vroom::vroom_write(bigger_comp, file = '~/redesROSMAP/graphs/A')

#END