#
#randomize_data.R
#When we build a network partition, we must make sure that the modularization we
#observe is not due to chance. To do this, we can build a null model, 
#and construct a distribution of such graphs generated from randomizations of our
#original matrix. We will perform a hypothesis test to indicate whether there is a 
#difference or no difference between the randomized networks and our problem network.

#Libraries --- ---

pacman::p_load('tidyverse', 
               "igraph")

#Get data --- ---

#Adjacency matrices from both conditions 

adjm_AD <- vroom::vroom(file = '/datos/rosmap/cuts_by_heuristics/AD_graphs/percentile99.99_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_we_adjacency_matrix.tsv')

adjm_noAD <- vroom::vroom(file = '/datos/rosmap/cuts_by_heuristics/noAD_graphs/percentile99.99_ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_we_adjacency_matrix.tsv')

# Function to randomize the rows of the adjacency matrix

randomize_rows <- function(adj_matrix) {
  num_rows <- nrow(adj_matrix)
  shuffled_indices <- sample(1:num_rows)
  adj_matrix_shuffled <- adj_matrix[shuffled_indices, ]
  return(adj_matrix_shuffled)
}

# Matriz de adyacencia original (sustitúyela por tu propia matriz)
adj_matrix_original <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), nrow = 3)

#If you want to iterate in only one adjacecncy matrix --- ---

# Número de iteraciones
num_iterations <- 5

# Iterar para crear grafos randomizados
for (i in 1:num_iterations) {
  # Randomizar las filas de la matriz de adyacencia
  adj_matrix_randomized <- randomize_rows(adjm_AD)
  
  # Crear un objeto de grafo a partir de la matriz de adyacencia randomizada
  graph <- graph_from_adjacency_matrix(adj_matrix_randomized, mode = "undirected", weighted = TRUE)
  
  # Imprimir el grafo
  print(graph)
}

#If you want to iterate in a list of adjacecncy matrices --- ---

#List of matrices
adjacency_matrices <- list(adjm_AD = adjm_AD, 
                           adjm_noAD = adjm_noAD)

# Function to randomize the rows of a list of adjacency matrices.

randomize_rows_multiple <- function(adj_matrices) {
  randomized_matrices <- lapply(adj_matrices, function(adj_matrix) {
    num_rows <- nrow(adj_matrix)
    shuffled_indices <- sample(1:num_rows)
    adj_matrix_shuffled <- adj_matrix[shuffled_indices, ]
    return(adj_matrix_shuffled)
  })
  return(randomized_matrices)
}

# Número de iteraciones
num_iterations <- 5

# Iterar para crear grafos randomizados
for (i in 1:num_iterations) {
  # Randomizar las filas de cada matriz de adyacencia en la lista
  adj_matrices_randomized <- randomize_rows_multiple(adjacency_matrices)
  
  # Crear y mostrar los grafos correspondientes a las matrices randomizadas
  for (j in seq_along(adj_matrices_randomized)) {
    adj_matrix_randomized <- adj_matrices_randomized[[j]]
    graphs_randomized_matrices <- graph.adjacency(adj_matrix_randomized, mode = "undirected", weighted = TRUE)
    print(graph)
  }
}

