#non sequitur
#
#7.2.hub_and_betweeness_genes.R
#This script makes an analysis of hub genes and high betweeness genes

#Libraries --- ---

pacman::p_load('igraph',
               'ggplot2', 
               'dplyr', 
               'gridExtra', 
               'biomaRt', 
               "ggVennDiagram", 
               "ClusterProfiler")

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- ---

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml',
                        format = 'graphml')

graphs <- list(graphAD = graphAD, 
               graphnoAD = graphnoAD)

#Define functions --- --- 

#1. Function to convert gene names

# Connecting to the Ensembl database through biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart)
}

# Define function to change names of vertex for symbol names --- ---

translate_vertex_names <- function(graph) {
  graph_vnames <- V(graph)$name   # Extract vertex names
  graph_vnames_trad <- convert_ens_to_symbol(graph_vnames)   # Translate names
  # Replace the missing values in the column 'external_gene_name' with the values of 'ensembl_gene_id'.
  graph_vnames_trad$external_gene_name <- ifelse(graph_vnames_trad$external_gene_name == "", graph_vnames_trad$ensembl_gene_id, graph_vnames_trad$external_gene_name)
  # Create a vector of translated names using the dictionary
  # We need to ensure that the actual names of the network are in the dictionary
  graph_vnames_trad <- setNames(graph_vnames_trad$external_gene_name, graph_vnames_trad$ensembl_gene_id)
  # Sort graph_vnames_trad according to the order of graph_vnames
  sorted_graph_vnames_trad <- graph_vnames_trad[match(graph_vnames, names(graph_vnames_trad))]
  # Assign the new names to the network vertices.
  V(graph)$name <- sorted_graph_vnames_trad
  return(graph)
}

#2. Jaccard Index function

jaccard_simplex <- function(a,b){
  length(intersect(a,b))/length(union(a,b))
}

#3. Create degree data frames

process_distribution <- function(degree_distribution) {
  # Convert to data frame and count frequencies
  distribution.df <- as.data.frame(table(degree_distribution$degree))
  # Rename the column
  distribution.df <- distribution.df %>% rename(Var1 = "degree")
  # Order by degree
  distribution.df <- distribution.df[order(distribution.df$degree, decreasing = F), ]
  # Calculate the cumulative sum
  distribution.df$CumulativeDegree <- cumsum(distribution.df$Freq)
  # Logarithm of the cumulative sum
  distribution.df$logCumulativeDegree <- log10(distribution.df$CumulativeDegree)
  # Add threshold
  CumulativeDegree_threshold <- quantile(distribution.df$logCumulativeDegree, probs = 0.95)
  # Puedes devolver el umbral si es necesario, por ejemplo:
  return(list(degree_distribution = distribution.df, threshold = CumulativeDegree_threshold))
  return(distribution.df)
}

########## HUB GENES ########## 

#Calculate degree of nodes
nodes_degree <- sapply(X = graphs, FUN = degree)

#Table of degree distribution

degree_distribution <- list()

for (i in 1:length(nodes_degree)) {
  degree_distribution[[i]] <- data.frame(gene = names(nodes_degree[[i]]), degree = nodes_degree[[i]])
}

#Calculate percentile 95 of genes with higher degree ---- ---

q_threshold.de <- list()

for (i in 1:length(nodes_degree)) {
  # Aplicar cluster_infomap a cada grafo y almacenar el resultado en results_list
  q_threshold.de[[i]] <- quantile(nodes_degree[[i]], probs = 0.95)
}

#Plot degree distribution --- ---
# Aplicar la función a cada elemento de la lista
processed_distribution <- lapply(degree_distribution, FUN = process_distribution)

# Processed distribution of both 
processed_distribution_AD <- processed_distribution[[1]]$degree_distribution
processed_distribution_AD$dx <- "AD"
processed_distribution_noAD <-  processed_distribution[[2]]$degree_distribution
processed_distribution_noAD$dx<- "no AD"

#Plot both distributions

degree_distr.p <- ggplot(d_distribution, aes(x = degree, y = logCumulativeDegree, group = dx, color = dx)) +
  geom_point() +
  geom_line() +
  labs(title = "Degree Distribution",
       x = "Degree",
       y = "logCumulativeDegree") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Rotar los nombres del eje X
    theme_light()
    
#Save plot 

ggsave("d_distribution.png", plot = degree_distr.p, width = 10, height = 4, units = "in", dpi = 300)



