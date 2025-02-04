#
#7.topological_comparison.R
#This script makes the topological comparison of coexpression graphs 
#constructed for people with AD and people without pathological AD.

#paulinapglz.99@gmail.com
#Part of code adapted from https://github.com/guillermodeandajauregui/BiologicalModuleComparison/blob/master/comparisonParameters.R

#Libraries --- ---
#install.packages("svglite")
pacman::p_load("igraph", 
               "ggraph",
               "tidyverse", 
               "gridExtra", 
               "svglite", 
               "tidygraph")

#Set seed --- --- 

set.seed(10)

#Declare functions --- --- 

#Declare function that compares edges

jaccard_edges <- function(g1, g2){
  return(length(E(igraph::intersection(g1, g2)))/length(E(igraph::union(g1, g2))))
}

#Declare function that compares nodes

jaccard_nodes <- function(g1,g2){
  a = sort(vertex.attributes(graph = g1)[["name"]])
  b = sort(vertex.attributes(graph = g2)[["name"]])

  deLaCueva = length(intersect(a,b))/length(union(a,b))
  return(deLaCueva)
}

#Get data --- --- 

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

#Save graphs in a list

graphLists <- list(graphAD = graphAD,
                  graphnoAD = graphnoAD)

#Network topological comparison --- --- 

# Extraer los nombres de los nodos de cada red
nodes_AD <- V(graphAD)$name
nodes_noAD <- V(graphnoAD)$name

# Encontrar los nodos comunes
common_nodes <- intersect(nodes_AD, nodes_noAD)
num_common_nodes <- length(common_nodes)

# Calcular el porcentaje de coincidencia respecto a cada red
percent_common_AD <- (num_common_nodes / length(nodes_AD)) * 100
percent_common_noAD <- (num_common_nodes / length(nodes_noAD)) * 100

# Imprimir resultados
cat("Número de nodos en común:", num_common_nodes, "\n")
cat("Porcentaje de nodos en común en la red g1:", percent_common_AD, "%\n")
cat("Porcentaje de nodos en común en la red g2:", percent_common_noAD, "%\n")

jaccard_nodes(graphAD, graphnoAD)

#Calculate diameter of both graphs --- --- 

diameter <- sapply(X = graphLists, FUN = diameter)
diameter
# graphAD graphnoAD 
# 13       12 

#Size of larger component --- ---

max(components(graphAD)$csize)

max(components(graphnoAD)$csize)


#Eigenvector centrality of the network --- ---

eigen <- sapply(graphLists, FUN = eigen_centrality)
eigen

# graphAD      graphnoAD   
# vector  numeric,1113 numeric,1074
# value   178.6451     187.1376    
# options list,20      list,20 

#Calculate clustering coefficient of both networks ---- ---

clustering_coefficient <- sapply(X = graphLists, FUN = transitivity)
clustering_coefficient

# graphAD graphnoAD 
# 0.5956005 0.6252799 

infomap_modularity <- sapply(X = graphLists, FUN = cluster_infomap)
infomap_modularity

# Extract information from the modules

membership_modularity <- sapply(X = infomap_modularity, FUN = membership)

# Assign the modules as attributes of the vertices

V(graphAD)$modules <- membership_modularity[[1]]

V(graphnoAD)$modules <- membership_modularity[[2]]

#Comparison of modular structures between networks --- --- 

#To compare two networks at the modular level, it would be optimal to keep the same set of nodes. 
#Because the heuristic cut was made on the edges, we have an unequal number of nodes in the smaller
#network, so we will add the missing nodes to the smaller network. 

#graphs with equal node set

ADnodes <- V(graphAD)

noADnodes <- V(graphnoAD)

# Find elements in noADnodes but not in ADnodes
missing_elements_in_ADnodes <- setdiff(names(noADnodes), names(ADnodes))
missing_elements_in_ADnodes <- as.list(missing_elements_in_ADnodes)
names(missing_elements_in_ADnodes) <- missing_elements_in_ADnodes
length(missing_elements_in_ADnodes)

# Find elements in ADnodes but not in noADnodes
missing_elements_in_noADnodes <- setdiff(names(ADnodes), names(noADnodes))
missing_elements_in_noADnodes <- as.list(missing_elements_in_noADnodes)
names(missing_elements_in_noADnodes) <- missing_elements_in_noADnodes
length(missing_elements_in_noADnodes)

#Add nodes missing to AD graph

graphAD_plus <- add_vertices(graphAD, nv = length(missing_elements_in_ADnodes), attr = missing_elements_in_ADnodes)

#Add nodes missing to noAD graph

graphnoAD_plus <- add_vertices(graphnoAD, nv = length(missing_elements_in_noADnodes), attr = missing_elements_in_noADnodes)

#Set list of graphs

graphLists_plus <- list(graphAD_plus = graphAD_plus, 
                        graphnoAD_plus = graphnoAD_plus)

#Do they have similar nodes? --- --- 

#To make the topological comparison for nodes and edges we use the Jaccard index

#Apply function to compare nodes

NodesJaccard <- sapply(X = graphLists_plus, FUN = jaccard_nodes, g1 = graphLists_plus[[1]])
NodesJaccard

#graphAD graphnoAD  <-----
#0.6465551 1.0000000 

#Do they have similar edges? --- ---

#Apply function to compare edges

EdgesJaccard <- sapply(X = graphLists, FUN = jaccard_edges, g1 = graphLists[[1]])
EdgesJaccard

# graphAD graphnoAD 
# 1.0000000 0.6839084 

#Apply modularity algorithm --- ---

#Mutual information must be taken into account as weights of the edges

infomap_modularity <- list()

for (i in 1:length(graphLists)) {
  # Aplicar cluster_infomap a cada grafo y almacenar el resultado en results_list
  infomap_modularity[[i]] <- cluster_infomap(graph = graphLists[[i]], e.weights = graphLists[[i]]$mut_info_norm)
}

# Extract information from the modules

graphAD_plus_modu <- cluster_infomap(graphAD_plus, e.weights = graphAD_plus$mut_info_norm)

graphnoAD_plus_modu <- cluster_infomap(graphnoAD_plus, e.weights = graphnoAD_plus$mut_info_norm)

#Modules from plus networks

graphAD_plus_modules <- membership(graphAD_plus_modu)

graphnoAD_plus_modules <-  membership(graphnoAD_plus_modu)

# Assign the modules as attributes of the vertices

V(graphAD_plus)$modules <- graphAD_plus_modules

V(graphnoAD_plus)$modules <- graphnoAD_plus_modules

#Compare modularity between graphs, applying 
#variation of information "vi"
#normalized mutual information "nmi"
#split-join distance "split-join distance"
#Rand index "Rand index"
#adjusted Rand index "adjusted Rand index"

possible_algos <- c("vi", "nmi", "split.join", "rand", "adjusted.rand")

comparison_methods <- sapply(X = possible_algos, FUN = function(i){
  igraph::compare(comm1 = graphAD_plus_modules,
                  comm2 = graphnoAD_plus_modules,
                  method = i
  )
})

comparison_methods

#END