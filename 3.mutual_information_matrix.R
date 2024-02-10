#
#3.mutual_information_matrix.R
#Script que hace mutual information para datos de expresion RNAseq de ROSMAP
#paulinapglz.99@gmail.com

pacman::p_load('future', 
             'tidyverse', 
             'furrr', 
             'infotheo')
#read data ----------

mat_dis <- vroom::vroom(file = '/datos/rosmap/discretized_matrix/ROSMAP_allNIAReaganspecimen_discretizedmatrix_10022024.tsv')
dim(mat_dis)
#[1] 14951   435

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

#Mutual information can't calculate for 1st column, gene names, so 
mat_dis <- mat_dis[-1]

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
})

#print time ----

print(Sys.time() - tempus)

#write matrix ----

#saveRDS(MI_MI, "ROSMAP_RNAseq_MutualInfo_allNIA_Reagan_dicho.rds")

#Next script is the 4.rds_to_matrix.R and 5.zero_diag.R