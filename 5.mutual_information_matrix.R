#
#3.mutual_information_matrix.R
#Script that calculates mutual information in a discretized matrix from RNA-seq expression data
#paulinapglz.99@gmail.com

pacman::p_load('future', 
               'tidyverse', 
               'furrr', 
               'infotheo')

#read data ----------

mat_dis <- readRDS(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/ROSMAP_DLFPC_AD_NIAReagan_discretizedmatrix.rds")

# indexing data --------------

my_index <- colnames(mat_dis) %>% as.character()

#Create vector with named indexes

my_index_i <- seq_along(my_index) #index number of my_index

names(my_index_i) <- my_index  #giving names to index

#the dimensions of the number of variables in mat_dis and the vector size must be the same

#set multicore plan --------

plan(multicore, workers = 40)

#set timer ---- 

tempus <- Sys.time()

#Calculate mutual information ------

#Parallel MI with future_map() and map()

MI_MI <- future_map(
  .x = my_index_i,      #named vector where we are going to apply a function .f
  .f = function(k){     #create function that calculates MI between a vector k...
    kk = mat_dis[k]
    map(
      .x = my_index_i,     #same named vector
      .f = function(m){
        mm = mat_dis[m]       #...and a vector m
        mutinformation(kk,mm)
      })
  }, .progress = TRUE)

#print time ----
print(Sys.time() - tempus)

#write matrix ----

saveRDS(MI_MI, "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/MI_matrices_NIA_Reagan/ROSMAP_DLFPC_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho.rds")

#Next script is the 3.1.rds_to_matrix.R