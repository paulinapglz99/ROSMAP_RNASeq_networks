#
#4.1.1.normalize_edgelist.R
#This script takes an edgelist and normalizes it 

#Libraries --- ---

pacman::p_load( "dplyr",
                "microbenchmark")

# Record the start time

start_time <- Sys.time()

#Get the data --- ---

edgelist <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.txt')

# get the min
minx <- min( edgelist$MI )
maxx <- max( edgelist$MI )
diffx <- maxx - minx

# now normalize all the data
edgelist_normalized <- edgelist %>% 
  mutate( mut_info_norm = ( MI - minx ) / diffx )

# Record the end time
end_time <- Sys.time()

# Calculate the elapsed time
elapsed_time <- end_time - start_time

# Print the elapsed time
print(paste("Elapsed Time: ", elapsed_time))

#Write data
vroom::vroom_write(edgelist_normalized, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.txt')

#END