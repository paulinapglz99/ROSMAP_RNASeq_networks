#3.mutual_information_matrix.R
#Script that calculates mutual information in discretized matrices from RNA-seq expression data
#paulinapglz.99@gmail.com

pacman::p_load('future', 
               'tidyverse', 
               'furrr', 
               'infotheo')

# Function to calculate mutual information matrix for a given file
mutual_information <- function(file_path, output_path) {
  # Read data
  mat_dis <- readRDS(file = file_path)
  
  # Indexing data
  my_index <- colnames(mat_dis) %>% as.character()
  
  # Create vector with named indexes
  my_index_i <- seq_along(my_index) # Index number of my_index
  names(my_index_i) <- my_index  # Giving names to index
  
  # Set multicore plan
  plan(multicore, workers = 40)
  
  # Set timer
  tempus <- Sys.time()
  
  # Calculate mutual information
  MI_MI <- future_map(
    .x = my_index_i,      # Named vector where we are going to apply a function .f
    .f = function(k){     # Create function that calculates MI between a vector k...
      kk = mat_dis[k]
      map(
        .x = my_index_i,     # Same named vector
        .f = function(m){
          mm = mat_dis[m]       # ...and a vector m
          mutinformation(kk, mm)
        })
    }, .progress = TRUE)
  
  # Print time
  print(Sys.time() - tempus)
  
  # Write matrix
  saveRDS(MI_MI, output_path)
}

# Paths for the matrices
file_path_AD <- "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/ROSMAP_DLFPC_AD_NIAReagan_discretizedmatrix.rds"
file_path_noAD <- "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/ROSMAP_DLFPC_noAD_NIAReagan_discretizedmatrix.rds"

# Output paths
output_path_AD <- "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.rds"
output_path_noAD <- "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho.rds"

# Calculate mutual information for both matrices
mutual_information(file_path_AD, output_path_AD)
mutual_information(file_path_noAD, output_path_noAD)

#END