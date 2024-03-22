#
#jaccard_index_probe.R

#Libraries --- ---


pacman::p_load("igraph", 
               "tidyverse")

#Get data --- --- 


graphAD <- read_graph(file = '~/redesROSMAP/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '~/redesROSMAP/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

#do they have similar nodes?

jaccard_nodes <- function(g1,g2){
  a = sort(vertex.attributes(graph = g1)[["name"]])
  b = sort(vertex.attributes(graph = g2)[["name"]])
  
  deLaCueva = length(intersect(a,b))/length(union(a,b))
  return(deLaCueva)
}

# NodesJaccard = c(
#           jaccard_nodes(g, g_alt),
#           jaccard_nodes(g, g_mix),
#           jaccard_nodes(g, g_rew)
# )

NodesJaccard = sapply(X = graphLists, FUN = jaccard_nodes, g1 = g)

igraph::compare(comm1 = graphAD, comm2 = graphnoAD, method = 'nmi')

