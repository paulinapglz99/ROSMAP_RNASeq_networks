# 3.mutual_information_matrix.R
# Script that calculates mutual information in a discretized matrix from RNA-seq expression data
# paulinapglz.99@gmail.com

pacman::p_load('future', 
               'tidyverse', 
               'furrr', 
               'infotheo')

# Read data ----------
args <- commandArgs(trailingOnly = TRUE)
mat_dis <- args[1]
output_file <- args[2] # Segundo argumento para el nombre del archivo de salida

mat_dis <- readRDS(mat_dis)

# Indexing data --------------

my_index <- colnames(mat_dis) %>% as.character()

# Create vector with named indexes

my_index_i <- seq_along(my_index) # index number of my_index

names(my_index_i) <- my_index  # giving names to index

# The dimensions of the number of variables in mat_dis and the vector size must be the same

# Set multicore plan --------

plan(multicore, workers = 20)

# Set timer ---- 

tempus <- Sys.time()

# Calculate mutual information ------

# Parallel MI with future_map() and map()

MI_MI <- future_map(
  .x = my_index_i,      # named vector where we are going to apply a function .f
  .f = function(k){     # create function that calculates MI between a vector k...
    kk = mat_dis[k]
    map(
      .x = my_index_i,     # same named vector
      .f = function(m){
        mm = mat_dis[m]       # ...and a vector m
        mutinformation(kk, mm)
      })
  }, .progress = TRUE)

# Print time ----
print(Sys.time() - tempus)

# Write matrix ----

saveRDS(MI_MI, output_file) # Usar el argumento para el nombre del archivo de salida

# Next script is the 3.1.rds_to_matrix.R