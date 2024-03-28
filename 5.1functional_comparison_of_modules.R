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

#Comparison of biological function sets associated to the overall network --- ---

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

#For noAD graph

enrichment_fullnet_noAD <- enrichGO(gene = V(graphnoAD)$name,
                                  OrgDb = "org.Hs.eg.db", 
                                  keyType = 'ENSEMBL',
                                  readable = TRUE,
                                  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.10)

#Comparison of enrichments

OverallEnrichedProcessJ <- jaccard_simplex(names(enrichment_fullnet_AD@geneSets), names(enrichment_fullnet_noAD@geneSets))
#[1] 0.8042686

#Similarity of modules in terms of associated biological functions --- ---

#Modularity algorithm 

modularity <- sapply(X = graphLists, FUN = cluster_infomap)

#Split lists of nodes by module

nodes_membership <- sapply(infomap_modularity, FUN = membership)

# Extract list of nodes by community for each graph

nodes_by_community_list <- lapply(seq_along(infomap_membership), function(i) {
  split(V(graphLists[[i]])$name, infomap_membership[[i]])
})

#I'd like to add names to lists of graphs again

names(nodes_by_community_list)[1] <- "graphAD"
names(nodes_by_community_list)[2] <- "graphnoAD"

#Extract list of nodes by community 

nodes_by_community_AD <- nodes_by_community_list[["graphAD"]]

nodes_by_community_noAD <- nodes_by_community_list[["graphnoAD"]]

#Enrichment to all modules of a network

#Define function that enrichs 

enrichGO_for_vector1 <- function(gene_vector) {
  if(length(gene_vector) >= 10) {
    enrichGO_result <- enrichGO(gene = gene_vector,
                                OrgDb = org.Hs.eg.db, 
                                keyType = 'ENSEMBL',
                                readable = TRUE,
                                ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                                pvalueCutoff = 0.05, 
                                qvalueCutoff = 0.10)
    return(enrichGO_result)
  } else {
    return(NULL)  # Return NULL if the gene vector is too small
  }
}

#Enrichment of graphAD modules
enriched_results_AD <- lapply(nodes_by_community_AD, enrichGO_for_vector)

#Enrichment of graphAD modules
enriched_results_noAD <- lapply(nodes_by_community_noAD, enrichGO_for_vector)

#It is possible to compare each module of the network of interest to the modules
#of another network, in terms of their associated biological function





#Functional module comparison
#
#By analyzing the sets of identified biological functions for each
#network, we can compare and contrast which functions are
#shared by both networks or uniquely found in each


