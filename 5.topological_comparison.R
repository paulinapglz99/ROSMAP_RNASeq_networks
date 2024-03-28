#
#topological_comparison.R
#This script makes the topological and modular comparison of coexpression graphs 
#constructed for people with AD and people without pathological AD.

#paulinapglz.99@gmail.com
#Part of code adapted from https://github.com/guillermodeandajauregui/BiologicalModuleComparison/blob/master/comparisonParameters.R

#Libraries --- ---

pacman::p_load("igraph", 
               "ggraph",
               "tidyverse")

#Get data --- --- 

graphAD <- read_graph(file = '/datos/rosmap/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

#Save graphs in a list

graphLists = list(graphAD = graphAD,
                  graphnoAD = graphnoAD)

#Network topological comparison --- --- 

#Calculate clustering coefficient of both networks

clustering_coefficient <- sapply(X = graphLists, FUN = transitivity)

#graphAD graphnoAD 
#0.4522337 0.4454295 

#Degree distributions look the same? --- --- 

#Table of degree distribution

#Degree distributions of all graphs in my list
degree_distributions <- sapply(X = graphLists, FUN = degree)

#Build degree distribution dataframe

ADdegree <- degree_distributions[["graphAD"]]

noADdegree <- degree_distributions[["graphnoAD"]]

degree_distributions <- rbind(data_frame(gene = names(ADdegree), degree = ADdegree, dx = "AD"), 
                              data_frame(gene = names(noADdegree), degree = noADdegree, dx = "noAD"))

#Plot both degree distributions

ggplot(degree_distributions, aes(x = degree, fill = dx)) +
  geom_histogram(binwidth = 1, color = "black", alpha =0.5, position = "identity") +
  labs(title = "Node degree distributions",
       subtitle = "for AD and noAD coexpression networks", 
       x = "Degree",
       y = "Freq") +
  theme_minimal() +
  scale_fill_manual(values = c("AD" = "hotpink", "noAD" = "blue4")) +
  guides(fill = guide_legend(title = "Diagnosis"))

#Do they have similar nodes? --- --- 

#To make the topological comparison for nodes and edges we use the Jaccard index

#Declare function that compares nodes

jaccard_nodes <- function(g1,g2){
  a = sort(vertex.attributes(graph = g1)[["name"]])
  b = sort(vertex.attributes(graph = g2)[["name"]])
  
  deLaCueva = length(intersect(a,b))/length(union(a,b))
  return(deLaCueva)
}

#Apply function to compare nodes

NodesJaccard = sapply(X = graphLists, FUN = jaccard_nodes, g1 = graphnoAD)
NodesJaccard

#graphAD graphnoAD 
#0.603139  1.000000 

#Do they have similar edges? --- ---

#Declare function that compares edges

jaccard_edges <- function(g1, g2){
  return(length(E(igraph::intersection(g1, g2)))/length(E(igraph::union(g1, g2))))
}

#Apply function to compare edges

EdgesJaccard = sapply(X = graphLists, FUN = jaccard_edges, g1 = graphnoAD)
EdgesJaccard

# graphAD graphnoAD 
#0.3561419 1.0000000 

#Apply modularity algorithm --- ---

infomap_modularity <- sapply(X = graphLists, FUN = cluster_infomap)

# Extract information from the modules

membership_modularity <- sapply(X = infomap_modularity, FUN = membership)

# Assign the modules as attributes of the vertices

V(graphAD)$modules <- membership_modularity[["graphAD"]]

V(graphnoAD)$modules <- membership_modularity[["graphnoAD"]]

# Plot the graph with ggraph and color by module.

#For AD graph
ggraph(graphAD, layout = 'kk') +
  geom_edge_link(alpha = 0.6, size = 0.3) +
  geom_node_point(aes(color = factor(modules))) +
  scale_color_manual(values = rainbow(max(membership_modularity[["graphAD"]]))) +
  theme_void() +
  ggtitle("Nodes colored by module for AD graph")

#For noADgraph

ggraph(graphnoAD, layout = 'kk') +
  geom_edge_link(alpha = 0.6, size = 0.3) +
  geom_node_point(aes(color = factor(modules))) +
  scale_color_manual(values = rainbow(max(membership_modularity[["graphnoAD"]]))) +
  theme_void() +
  ggtitle("Nodes colored by module for noAD graph")

#Comparison of modular structures --- ---

#The modularity index Q modularity(), is a measure of the proportion of edges that occur within communities,
#relative to the expected proportion if all edges were placed randomly.

modularity_scoreAD <- modularity(graphAD, membership_modularity[["graphAD"]])
#[1] 0.3491475

modularity_scorenoAD <- modularity(graphnoAD, membership_modularity[["graphnoAD"]])
#[1] 0.3790371

#Comparison of modular structures between networks

#To compare two networks at the modular level, it would be optimal to keep the same set of nodes. 
#Because the heuristic cut was made on the edges, we have an unequal number of nodes in the smaller
#network, so we will add the missing nodes to the smaller network. 

#graphs with equal node set

ADnodes <- V(graphAD)

noADnodes <- V(graphnoAD)

# Find elements in noADnodes but not in ADnodes
missing_elements_in_ADnodes <- setdiff(names(noADnodes), names(ADnodes))
length(missing_elements_in_ADnodes)
#[1] 296

# Find elements in ADnodes but not in noADnodes
missing_elements_in_noADnodes <- setdiff(names(ADnodes), names(noADnodes))
length(missing_elements_in_noADnodes)
#[1] 235

#Add nodes missing to AD graph

graphAD_plus <- add_vertices(graphAD, nv = length(missing_elements_in_ADnodes))

#Add nodes missing to noAD graph

graphnoAD_plus <- add_vertices(graphnoAD, nv = length(missing_elements_in_noADnodes))

# Extract information from the modules

graphAD_plus_modu <- cluster_infomap(graphAD_plus)

graphnoAD_plus_modu <- cluster_infomap(graphnoAD_plus)

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

#vi           nmi    split.join          rand adjusted.rand 
#3.534727e+00  5.633515e-01  1.339000e+03  8.429085e-01  6.274289e-02 

#END
#Next script is the 5.1functional_comparison_of_modules.R