#
#null_model.R
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
rownames(adjm_AD) <- colnames(adjm_AD)

adjm_noAD <- vroom::vroom(file = '/datos/rosmap/cuts_by_heuristics/noAD_graphs/percentile99.99_ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_we_adjacency_matrix.tsv')
rownames(adjm_noAD) <- colnames(adjm_noAD)

########################

# Función para randomizar las filas de una matriz de adyacencia
randomizar_filas <- function(matriz) {
  return(matriz[sample(nrow(matriz)), ])
}

# Crear una matriz de adyacencia de ejemplo
matriz_adyacencia <- matrix(c(0, 1, 1, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0), nrow = 4, byrow = TRUE)

# Número de iteraciones
num_iteraciones <- 5

# Lista para almacenar las matrices con filas randomizadas
matrices_randomizadas <- list()

# Iterar para randomizar las filas y almacenar las matrices resultantes
for (i in 1:num_iteraciones) {
  matriz_randomizada <- randomizar_filas(adjm_AD)
  rownames(matriz_randomizada) <- colnames(matriz_randomizada) # Asignar nombres de columnas como nombres de filas
  matrices_randomizadas[[i]] <- matriz_randomizada
}



######################

random_adjm_AD <- adjm_AD[sample(nrow(adjm_AD)), ]
rownames(random_adjm_AD) <- colnames(random_adjm_AD)

random_adjm_AD.g <- graph_from_adjacency_matrix(random_adjm_AD)

#If you want to iterate in a list of adjacecncy matrices --- ---

# Define la función para randomizar las filas de una matriz de adyacencia
randomize_adj_matrix <- function(adj_matrix) {
  randomized_indices <- sample(nrow(adj_matrix)) # Permuta los índices de las filas
  randomized_matrix <- adj_matrix[randomized_indices, randomized_indices] # Aplica la permutación a las filas y columnas
  return(randomized_matrix)
}

# Define la función para crear un grafo a partir de una matriz de adyacencia
create_graph <- function(adj_matrix) {
  graph_from_adjacency_matrix(adj_matrix, mode = "undirected")
}

# Crea una lista de listas para almacenar las matrices randomizadas
randomized_matrices <- list()

# Crea una lista para almacenar los grafos randomizados
randomized_graphs <- list()

# Suponiendo que tienes una matriz de adyacencia llamada adj_matrix
# Puedes reemplazar esto con tu propia matriz

# Define el número de iteraciones para generar matrices randomizadas
num_iterations <- 5

# Itera para generar matrices randomizadas y grafos
for (i in 1:num_iterations) {
  randomized_matrix <- randomize_matrix(adjm_AD)
  randomized_graph <- graph_from_adjacency_matrix(randomized_matrix, directed = F)
  
  # Agrega la matriz randomizada a la lista de matrices randomizadas
  randomized_matrices[[i]] <- randomized_matrix
  
  # Agrega el grafo randomizado a la lista de grafos randomizados
  randomized_graphs[[i]] <- randomized_graph
}

# Accede a las matrices y grafos randomizados en la lista
# Por ejemplo, para acceder a la tercera matriz randomizada:
# third_randomized_matrix <- randomized_matrices[[3]]

# Y para acceder al segundo grafo randomizado:
# second_randomized_graph <- randomized_graphs[[2]]

###########################################################################

#Second approach, el rewire

#Libraries --- ---

pacman::p_load("igraph", 
               "ggraph",
               "tidyverse")

#Get data --- --- 

graphAD <- read_graph(file = '/datos/rosmap/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml', format = 'graphml')

# 1. Cargar el grafo existente
# Supongamos que tienes tus datos y los cargas en un grafo existente llamado 'grafo_existente'
# Puedes cargar tu grafo de diferentes maneras, dependiendo del formato de tus datos

# Definir el número de iteraciones
num_iteraciones <- 20

# Lista para almacenar los grafos randomizados
grafos_random <- list()

# Iterar para generar grafos aleatorios
for (i in 1:num_iteraciones) {
  # Generar una copia del grafo existente
  grafo_random <- graphAD
  
  # Aplicar alguna técnica de randomización
  # Por ejemplo, puedes utilizar rewire() para realizar un rewiring aleatorio de las aristas
  grafo_random <- rewire(grafo_random, each_edge(prob = 0.5))
  
  # Almacenar el grafo aleatorio en la lista
  grafos_random[[i]] <- grafo_random
}

#Apply modularization to all randomized graphs

infomap_modularity <- sapply(X = grafos_random, FUN = cluster_infomap)

#Esto genera grafos que contienen un solo modulo, por lo que no sirve para hacer el grafico de distribucion de modulos que queria hacer

