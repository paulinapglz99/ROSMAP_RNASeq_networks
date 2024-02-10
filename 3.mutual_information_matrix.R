#Script que hace mutual information para datos de expresion RNAseq de ROSMAP
#paulinapglz.99@gmail.com

pacman::p_load('future', 
             'tidyverse', 
             'infotheo',
             'furrr')

#set timer ---- 

tempus <- Sys.time()

#read data ----------

protcod <- vroom(file = '/datos/rosmap/discretized_matrix/protcod_AD.txt')  #change for every desired diagnosis 

mat_dis <- vroom(file = '/datos/rosmap/discretized_matrix/mat_dis_AD.txt')

# indexing data --------------

my_index <- pull(protcod, 'gene_id') #to rename matrix

#Quiero un vector que tenga los conteos de cada elemento de otro vector

my_index_i <- seq_along(my_index) #indice numero de my_index

names(my_index_i) <- my_index  #etiquetar con indices, le doy nombre a los elementos del vector

#set multicore plan --------

plan(multicore, workers = 40)

#Calculate mutual information ------

MI_MI <- future_map(.x = my_index_i, .f = function(k){
kk = mat_dis[k]
 map(.x = my_index_i, .f = function(m){
  mm = mat_dis[m]
   mutinformation(kk,mm)
 })
})

#write matrix ----

#write_rds(x = MI_MI,
#         file = "ROSMAP_RNAseq_MutualInfo_allMCI_matrix.rds") #this one did not work

saveRDS(MI_MI, "ROSMAP_RNAseq_MutualInfo_allMCI_matrix.rds")

#print time ----

print(Sys.time() - tempus)

#Next script is the 4.rds_to_matrix.R and 5.zero_diag.R
