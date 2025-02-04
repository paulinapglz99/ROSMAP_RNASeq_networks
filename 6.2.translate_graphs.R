#
#6.2.translate_graphs.R
#This function allows you to change the vertex names of a graph. It translates vertices in ensembl to symbol. 
#paulinapglz.99@gmail.com

pacman::p_load('igraph',
               'dplyr', 
               'biomaRt')

# Define function to change names of vertex for symbol names --- ---

# Connecting to the Ensembl database through biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart)
}

translate_vertex_names <- function(graph) {
  #Extract vertex names
  graph_vnames <- V(graph)$name
  
  #Translate names
  graph_vnames_trad <- convert_ens_to_symbol(graph_vnames)
  
  #Replace the missing values in the column 'external_gene_name' with the values of 'ensembl_gene_id'.
  graph_vnames_trad$external_gene_name <- ifelse(graph_vnames_trad$external_gene_name == "", graph_vnames_trad$ensembl_gene_id, graph_vnames_trad$external_gene_name)
  
  # Create a vector of translated names using the dictionary
  # We need to ensure that the actual names of the network are in the dictionary
  graph_vnames_trad <- setNames(graph_vnames_trad$external_gene_name, graph_vnames_trad$ensembl_gene_id)
  
  # Sort graph_vnames_trad according to the order of graph_vnames
  sorted_graph_vnames_trad <- graph_vnames_trad[match(graph_vnames, names(graph_vnames_trad))]
    #Assign the new names to the network vertices.
  V(graph)$name_trad <- sorted_graph_vnames_trad
  
  return(graph)
}

#Translate AD graphs

#Get data --- ---

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml',
                        format = 'graphml')

graphs <- list(graphAD = graphAD, 
               graphnoAD = graphnoAD)

#Translate graphs --- ---

graphs_trad <- sapply(graphs, translate_vertex_names)

#Save graphs --- ---
 
write_graph(graphs_trad[["graphAD"]], file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml',
            format = "graphml")
# 
write_graph(graphs_trad[["graphnoAD"]], file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml',
             format = "graphml")

#END