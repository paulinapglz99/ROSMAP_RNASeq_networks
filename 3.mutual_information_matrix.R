#
#3.mutual_information_matrix.R
#Script que hace mutual information para datos de expresion RNAseq de ROSMAP
#paulinapglz.99@gmail.com

pacman::p_load('future', 
             'tidyverse', 
             'furrr')
#read data ----------

mat_dis <- vroom::vroom(file = '/datos/rosmap/discretized_matrix/ROSMAP_allNIAReaganspecimen_discretizedmatrix_10022024.tsv')

# indexing data --------------

my_index <- pull(mat_dis, 1) #pull vector from 1st col of mat_dis

#Create vector with named indexes

my_index_i <- seq_along(my_index) #index number of my_index

names(my_index_i) <- my_index  #giving names to index

#set multicore plan --------

plan(multicore, workers = 40)

#set timer ---- 

tempus <- Sys.time()

#Calculate mutual information ------

MI_MI <- future_map(
  .x = my_index_i, 
  .f = function(k){
kk = mat_dis[k]
 map(
   .x = my_index_i,
     .f = function(m){
  mm = mat_dis[m]
   mutinformation(kk,mm)
 })
})

#write matrix ----

#saveRDS(MI_MI, "ROSMAP_RNAseq_MutualInfo_allNIA_Reagan_dicho.rds")

#print time ----

print(Sys.time() - tempus)

#Next script is the 4.rds_to_matrix.R and 5.zero_diag.R