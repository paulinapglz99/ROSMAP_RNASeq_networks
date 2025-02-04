# Libraries  --- ---
pacman::p_load('tidyverse', 
               'Matrix', 
               'furrr')

# Read adjacency matrix --- ---
args <- commandArgs(trailingOnly = TRUE)

# This matrix may be very heavy
matrix <- args[1]
matrix <- readRDS(matrix)

output <- args[2]

# matrix <- readRDS("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.rds")
dim(matrix)

matrix <- as.data.frame(matrix)

# Extract universe of genes
universe <- colnames(matrix)
length(universe)

# Save genes
write(universe,
      file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/universe.txt')

# Parallel setup --- ---
# Set up the number of cores to use
plan(multicore, workers = parallel::detectCores() - 1)

# Pivot and parallel processing  ----
# Split matrix rows into chunks for parallel processing
matrix_chunks <- split(matrix, seq(nrow(matrix)))

# Use future_map_dfr to process each chunk in parallel and bind the results into a dataframe
full_edgelist <- future_map_dfr(matrix_chunks, function(chunk) {
  chunk %>%
    rownames_to_column(var = "feature") %>%
    pivot_longer(-feature, names_to = "gene_to", values_to = "MI") %>%
    filter(!is.na(MI))
}, .progress = TRUE)

dim(full_edgelist)
# Expected output size reduction after filtering NAs
#[1] 53369946        3

# Save edgelist for later
saveRDS(full_edgelist, file = output)

# END