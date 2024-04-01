#pruebas para mutual information
#iterar sobre el path my_index

my_index1_10 <- my_index[1:10]

#map pide un vector .x 

#se usan funciones anonimas o estandar

map(.x = my_index1_10, .f = print) #iterando

#en lugar de print

map(.x = my_index1_10, .f = function(i){
  
  paste0('yo soy ', i)
  
}

) #iterando la funcion()sobre un vector my_index1_10

#iterando sobre renglones

map(.x = my_index1_10, .f = function(i){
  
  paste0('yo soy ', i)
}
) 




my_index1_10[7]

mat_dis[,"V725"]

#si quisiera extraer el gene_id de la posicion 7 de mi my_index1_10 pero de la matriz mat_dis

mat_dis[,my_index1_10[7]]

#Quiero un vector que tenga los conteos de cada elemento de otro vector

my_index1_10_num <- seq_along(my_index1_10) #indice numero de my_index

names(my_index1_10_num) <- my_index1_10  #etiquetar con indices

my_index1_10_num


#Recapitulando

map(.x = my_index1_10_num, .f = function(i){
  
  paste0('yo soy ', i)
  
}

) #me devuelve una lista nombrada, el nombre de la lista es el nombre del vector sobre el que itere



#####

map(.x = my_index1_10_num, .f = function(i){
  
  mat_dis[,i]   #iterando sobre cada indice numerico (los indices numeros estan asociados a un valor)
  
}
)


############################################################################

### quiero iterar sobre mi indice numero y devolver los 10 primeros valores de mi 
###matriz discretizada

#hacer mis indices

my_index <- pull(protcod, 'gene_id') #de otro dataset not here pero que tenia en mi ambiente, ese esta en recap.R

#Quiero un vector que tenga los conteos de cada elemento de otro vector

my_index1_10_num <- seq_along(my_index1_10) #indice numero de my_index

names(my_index1_10_num) <- my_index1_10  #etiquetar con indices, le doy nombre a los elementos del vector

#mision: sacar la info mutua de los genes del 1 al 10, cuando lo tenga listo

MI_1_10 <- map(.x = my_index1_10_num, .f = function(k){
    
 kk = mat_dis[k]
 
 map(.x = my_index1_10_num, .f = function(m){
   
   mm = mat_dis[m]
   
   mutinformation(kk,mm)
   
 })
  })

#write the rds file

write_rds(x =MI_1_10, file = "MI_AD_matrix_prueba1_10.rds")
