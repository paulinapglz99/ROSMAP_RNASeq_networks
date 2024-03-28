#
#5.network_modularity.R
#This script applies the INFOMAP modularity algorithm to detect hierarchical communities in networks.

#Libraries --- --- 

#Get data --- --- 

graphAD <- read_graph(file = '/datos/rosmap/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

infomap_AD <- cluster_infomap(graphAD)

sizes(infomap_AD)

is_hierarchical(infomap_AD)
