#
#5.3.normalize_edgelist.R
#This script takes an edgelist and normalizes it 
#By APPG

#Libraries --- ---

pacman::p_load( "dplyr",
                "microbenchmark")

# Record the start time

start_time <- Sys.time()

#Get the data --- ---

#edgelist <- readRDS(file = 'ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.rds')

edgelist <- args[1]
edgelist <- readRDS(edgelist)

output <- args[2]

#get the min
minx <- min( edgelist$MI )
maxx <- max( edgelist$MI )
diffx <- maxx - minx

#now normalize all the data
edgelist_normalized <- edgelist %>% 
  mutate( mut_info_norm = ( MI - minx ) / diffx )

#Record the end time
end_time <- Sys.time()

#Calculate the elapsed time
elapsed_time <- end_time - start_time

#Print the elapsed time
print(paste("Elapsed Time: ", elapsed_time))

#Write data
#saveRDS(edgelist_normalized, file = 'ROSMAP_DLFPC_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.rds')


saveRDS(edgelist_normalized, file = output)

#END