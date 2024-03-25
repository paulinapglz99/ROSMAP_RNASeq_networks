#
#randomize_data.R
#When we build a network partition, we must make sure that the modularization we
#observe is not due to chance. To do this, we can build a null model of partitions, 
#and construct a distribution of such graphs generated from randomizations of our
#original matrix. We will perform a hypothesis test to indicate whether there is a 
#difference or no difference between the randomized networks and our problem network.

#Libraries --- ---

pacman::p_load('dplyr')

#Get data --- ---

graphAD <- read_graph(file = '~/redesROSMAP/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '~/redesROSMAP/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

#We need to get the adjacency matrix from the graphs

#As I need the weighted matrix to have a more accurrate model, I'd like to have the full matrix and then filter by the 
#node names I'd like to keep, in this case, the node names in our graphs

node_listAD <- V(graphAD)$name #all the genes in the AD network

node_listnoAD <-  V(graphnoAD)$name #all the genes in the noAD network


