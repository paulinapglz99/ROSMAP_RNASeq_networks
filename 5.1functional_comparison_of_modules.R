#
#5.1functional_comparison_of_modules.R

#This script makes the functional modular comparison of coexpression graphs 
#constructed for people with AD and people without pathological AD.

#Modules will be compared in terms of the associated biological functions identified by an enrichment analysis.
#This comparison involves answering three complementary questions: 

#1.How similar are the sets of biological functions that are associated to the whole network, through the enrichment of individual modules?
#2.How similar are the modules found in each network, in terms of the sets of associated biological functions?
#3.In how many modules is represented each biological process?


#paulinapglz.99@gmail.com
#Part of code adapted from https://github.com/guillermodeandajauregui/BiologicalModuleComparison/blob/master/comparisonParameters.R

#Libraries --- ---

pacman::p_load("igraph", 
               "tidyverse",
               "clusterProfiler")

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- --- 

graphAD <- read_graph(file = '/datos/rosmap/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

#Save graphs in a list

graphLists <- list(graphAD = graphAD,
                  graphnoAD = graphnoAD)

#Modularity algorithm 

modularity <- sapply(X = graphLists, FUN = cluster_infomap)

nodes_membership <- sapply(infomap_modularity, FUN = membership)

# Extract list of nodes by community for each graph

nodes_by_community_list <- lapply(seq_along(infomap_membership), function(i) {
  split(V(graphLists[[i]])$name, infomap_membership[[i]])
})

#I'd like to add names to lists of graphs again

names(nodes_by_community_list)[1] <- "graphAD"
names(nodes_by_community_list)[2] <- "graphnoAD"

#Extraer lista de vectores 

nodes_by_community_AD <- nodes_by_community_list["graphAD"]

##########################################

#cluster_informap find community structure that minimizes the expected description length of a random walker trajectory 
cluster_infomap <- cluster_infomap(graph)

#See membership of nodes
membership_infomap <- membership(cluster_infomap)

#Extract list of nodes by community
infomap_nodes_by_community <- split(V(graph)$name, membership_infomap)

##########################################

#Enrichment --- ---

#Define functions to use

#Similarity of Enriched Processes, Jaccard Index

jaccard_simplex <- function(a,b){
  length(intersect(a,b))/length(union(a,b))
}

#Enrichment of full graphs

#For AD graph

enrichment_fullnet_AD <- enrichGO(gene = V(graphAD)$name,
         OrgDb = "org.Hs.eg.db", 
         keyType = 'ENSEMBL',
         readable = TRUE,
         ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
         pvalueCutoff = 0.05, 
         qvalueCutoff = 0.10)

f <- names(enrichment_fullnet_AD@geneSets)

#For noAD graph

enrichment_fullnet_noAD <- enrichGO(gene = V(graphnoAD)$name,
                                  OrgDb = "org.Hs.eg.db", 
                                  keyType = 'ENSEMBL',
                                  readable = TRUE,
                                  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.10)

f1 <- names(enrichment_fullnet_noAD@geneSets)


full_net_enrichment_list <- list(enrichment_fullnet_AD = enrichment_fullnet_AD, 
                                 enrichment_fullnet_noAD = enrichment_fullnet_noAD)

#Comparison of enrichments

EnrichedProcessJ = jaccard_simplex(names(enrichment_fullnet_AD@geneSets), names(enrichment_fullnet_noAD@geneSets))

####################################

for (i in seq_along(nodes_by_community_AD$graphAD)) {
  community <- nodes_by_community_AD$graphAD[[i]]
}

results_list <- lapply(nodes_by_community_AD$graphAD, function(community) {
  enrichment_results <- enrichGO(gene = community,
                                 OrgDb = "org.Hs.eg.db",
                                 keyType = 'ENSEMBL',
                                 readable = TRUE,
                                 ont = "BP",
                                 pvalueCutoff = 0.05,
                                 qvalueCutoff = 0.10)
  
  return(enrichment_results)
})




enrichment_results <- enrichGO(gene = community,
                               OrgDb = "org.Hs.eg.db",
                               keyType = 'ENSEMBL',
                               readable = TRUE,
                               ont = "BP",
                               pvalueCutoff = 0.05,
                               qvalueCutoff = 0.10)

#Functional module comparison
#
#By analyzing the sets of identified biological functions for each
#network, we can compare and contrast which functions are
#shared by both networks or uniquely found in each


