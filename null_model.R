#
#randomize_data.R
#When we build a network partition, we must make sure that the modularization we
#observe is not due to chance. To do this, we can build a null model, 
#and construct a distribution of such graphs generated from randomizations of our
#original matrix. We will perform a hypothesis test to indicate whether there is a 
#difference or no difference between the randomized networks and our problem network.

#paulinapglz.99@gmail.com

#Libraries --- ---

pacman::p_load('tidyverse', 
               "igraph")

#Get data --- ---

#Adjacency matrices from both conditions 

adjm_AD <- vroom::vroom(file = '/datos/rosmap/cuts_by_heuristics/AD_graphs/percentile99.99_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_we_adjacency_matrix.tsv')

adjm_noAD <- vroom::vroom(file = '/datos/rosmap/cuts_by_heuristics/noAD_graphs/percentile99.99_ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_we_adjacency_matrix.tsv')

#If you want to iterate in a list of adjacecncy matrices --- ---

# Función para randomizar las filas de una matriz de adyacencia
randomize_rows <- function(adj_matrix) {
  num_nodes <- nrow(adj_matrix)
  random_index <- sample(num_nodes)
  adj_matrix_randomized <- adj_matrix[random_index, random_index]
  return(adj_matrix_randomized)
}

# Lista de matrices de adyacencia
lista_matrices <- list(
  matrix(c(0, 1, 1, 0,
           1, 0, 0, 1,
           1, 0, 0, 1,
           0, 1, 1, 0), nrow = 4, byrow = TRUE),
  matrix(c(0, 1, 0, 1,
           1, 0, 1, 0,
           0, 1, 0, 1,
           1, 0, 1, 0), nrow = 4, byrow = TRUE)
)

# Lista para almacenar las matrices randomizadas
lista_matrices_randomizadas <- list()
# Lista para almacenar los grafos randomizados
lista_grafos_randomizados <- list()

#Number of iterations 

num_iter <- 5

# Iterar sobre cada matriz en la lista original
for (i in seq_along(lista_matrices)) {
  matriz_original <- lista_matrices[[i]]
  # Crear y almacenar grafos para la matriz original
  grafo_original <- graph_from_adjacency_matrix(matriz_original, mode = "undirected")
  lista_grafos_randomizados[[paste0("Grafo_", i, "_original")]] <- grafo_original
  
  # Randomizar filas y crear grafos para cada matriz randomizada
  for (j in 1:num_iter) { 
    matriz_randomizada <- randomize_rows(matriz_original)
    lista_matrices_randomizadas[[paste0("Matriz_", i, "_randomizada_", j)]] <- matriz_randomizada
    grafo_randomizado <- graph_from_adjacency_matrix(matriz_randomizada, mode = "undirected")
    lista_grafos_randomizados[[paste0("Grafo_", i, "_randomizado_", j)]] <- grafo_randomizado
  }
}

# Ahora, lista_matrices_randomizadas contiene las matrices randomizadas
# y lista_grafos_randomizados contiene los grafos correspondientes

# Verificar los nombres de las listas
print(names(lista_matrices_randomizadas))
print(names(lista_grafos_randomizados))
