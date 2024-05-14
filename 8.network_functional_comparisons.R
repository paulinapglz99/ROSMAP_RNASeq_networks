#
#8.network_functional_comparisons.R

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
               "clusterProfiler", 
               "gridExtra")

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- --- 

graphAD <- read_graph(file = '/datos/rosmap/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

#Save graphs in a list

graphLists <- list(graphAD = graphAD,
                   graphnoAD = graphnoAD)

#Set seed for modularity algorithm --- ---

set.seed(10)

#Comparison of biological function sets associated to the overall network --- ---

#Enrichment of full graphs

#For AD graph

enrichment_fullnet_AD <- enrichGO(gene = V(graphAD)$name,
                                  OrgDb = "org.Hs.eg.db", 
                                  keyType = 'ENSEMBL',
                                  readable = TRUE,
                                  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.10)

#Gene concept network enrichment

enrichment_fullnet_AD_cnet <- cnetplot(enrichment_fullnet_AD, circular = TRUE, colorEdge = TRUE, showCategory= 10)

#Dotplot enrichment

enrichment_fullnet_AD_dot <- dotplot(enrichment_fullnet_AD)

#Save plots

ggsave("enrichment_fullnet_AD_cnet.png", enrichment_fullnet_AD_cnet, width = 15, height = 8)

#For noAD graph

enrichment_fullnet_noAD <- enrichGO(gene = V(graphnoAD)$name,
                                    OrgDb = "org.Hs.eg.db", 
                                    keyType = 'ENSEMBL',
                                    readable = TRUE,
                                    ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

#Barplot enrichment

barplot(enrichment_fullnet_noAD)

#Dotplot enrichment

dotplot(enrichment_fullnet_noAD)


#Define imilarity of Enriched Processes, Jaccard Index function --- ---

jaccard_simplex <- function(a,b){
  length(intersect(a,b))/length(union(a,b))
}

#Comparison of enrichments

OverallEnrichedProcessJ <- jaccard_simplex(names(enrichment_fullnet_AD@geneSets), names(enrichment_fullnet_noAD@geneSets))
#[1] 0.8042686

#Similarity of modules in terms of associated biological functions --- ---

#Modularity algorithm 

modularity <- sapply(X = graphLists, FUN = cluster_infomap)
#This can't be seen if it's not called from console

#Split lists of nodes by module

nodes_membership <- sapply(modularity, FUN = membership)

# Extract list of nodes by community for each graph

nodes_by_community_list <- lapply(seq_along(nodes_membership), function(i) {
  split(V(graphLists[[i]])$name, nodes_membership[[i]])
})

#I'd like to add names to lists of graphs again

names(nodes_by_community_list)[1] <- "graphAD"
names(nodes_by_community_list)[2] <- "graphnoAD"

#Extract list of nodes by community 

nodes_by_community_AD <- nodes_by_community_list[["graphAD"]]

nodes_by_community_noAD <- nodes_by_community_list[["graphnoAD"]]

#Enrichment to all modules of a network --- ---

#Define function that performs enrichment

enricher <- function(nodes_by_community) {
  enrichGO_result <- enrichGO(gene = nodes_by_community,
                              OrgDb = org.Hs.eg.db, 
                              keyType = 'ENSEMBL',
                              readable = TRUE,
                              ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)
  return(enrichGO_result)
}

#These functions may print a " --> No gene can be mapped...." on the console. 
#This happens when there are lists of genes that cannot be enriched. It will still generate a result

# Define a function to replace NULL by an empty S4 object "enrichResult" to handle NULLs in the enrichment lists

replace_null <- function(x) {
  if (is.null(x)) {
    return(new("enrichResult"))
  } else {
    return(x)
  }
}

#Enrichment of graphAD modules

enriched_results_AD <- lapply(nodes_by_community_AD, enricher) #slow
enriched_results_AD <- lapply(enriched_results_AD, replace_null) # Exchange NULLs for empty S4 objects

#Enrichment of graphAD modules
enriched_results_noAD <- lapply(nodes_by_community_noAD, enricher) #slow
enriched_results_noAD <- lapply(enriched_results_noAD, replace_null) # Exchange NULLs for empty S4 objects

#Functional module comparison between two graphs --- ---

#Similarity of modules in terms of associated biological functions

#It is possible to compare each module of the network of interest to the modules
#of another network, in terms of their associated biological function

# Create an empty array to store the results of the Jaccard index.
num_modules_AD <- length(enriched_results_AD)
num_modules_noAD <- length(enriched_results_noAD)

similarity_matrix <- matrix(NA, nrow = num_modules_AD, ncol = num_modules_noAD)

# Asignar nombres de filas y columnas
rownames(similarity_matrix) <- paste("AD", 1:num_modules_AD, sep = "_")
colnames(similarity_matrix) <- paste("noAD", 1:num_modules_noAD, sep = "_")

#Iterate to obtain similarity matrix 

for (i in 1:num_modules_AD) {
  for (j in 1:num_modules_noAD) {
    # Acceder a los geneSets de las redes AD y no AD
    geneSets_AD <- names(enriched_results_AD[[i]]@geneSets)
    geneSets_noAD <- names(enriched_results_noAD[[j]]@geneSets)
    similarity_matrix[i, j] <- jaccard_simplex(geneSets_AD, geneSets_noAD)
  }
}

head(similarity_matrix)

#Find the highest similarity of the matrix  
max(similarity_matrix, na.rm = TRUE)

#Number of modules in networks that have Jaccard index J=1 with a module of the Main network --- ---

# Function to count the number of columns with value equal to 1
count_ones <- function(row) {
  sum(row == 1)
}

# Apply the function to each row of the matrix
ones_count <- apply(similarity_matrix, 1, count_ones)

# Crear una tabla con los resultados
num_equal_nodes <- data.frame(module_number = rownames(similarity_matrix), modules_with_jaccardindex1 = ones_count)

#Number of modules associated to a given biological function --- ---

# Obtener la longitud del vector más corto
min_length <- min(length(v1), length(v2))

# Calcular la distancia euclidiana solo para las dimensiones compartidas
euclidean_distance <- sqrt(sum((v1[1:min_length] - v2[1:min_length])^2))

print(euclidean_distance)

#

euclidean_distance <- matrix(NA, nrow = num_modules_AD, ncol = num_modules_noAD)

# Asignar nombres de filas y columnas
rownames(euclidean_distance) <- paste("AD", 1:num_modules_AD, sep = "_")
colnames(euclidean_distance) <- paste("noAD", 1:num_modules_noAD, sep = "_")

for (i in 1:num_modules_AD) {
  for (j in 1:num_modules_noAD) {
    # Acceder a los geneSets de las redes AD y no AD
    geneSets_AD <- length(enriched_results_AD[[i]]@geneSets)
    geneSets_noAD <- length(enriched_results_noAD[[j]]@geneSets)
    euclidean_distance[i, j] <- sqrt(sum((v1[1:num_modules_AD] - v2[1:num_modules_noAD])^2))
  }
}
