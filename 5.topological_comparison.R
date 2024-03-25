#
#topological_comparison.R
#This script makes the topological and modular comparison of coexpression graphs 
#constructed for people with AD and people without pathological AD.

#paulinapglz.99@gmail.com
#Functions from https://github.com/guillermodeandajauregui/BiologicalModuleComparison/blob/master/comparisonParameters.R

#Libraries --- ---

pacman::p_load("igraph", 
               "tidyverse")

#Get data --- --- 

graphAD <- read_graph(file = '/datos/rosmap/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

#Save graphs in a list

graphLists = list(graphAD = graphAD,
                  graphnoAD = graphnoAD)

#Network topological comparison --- --- 

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
  guides(fill = guide_legend(title = "Diagnóstico"))

#To make the topological comparison for nodes and edges we use the Jaccard index

#Do they have similar nodes? --- --- 

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




