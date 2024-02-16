#
#3.rds_to_matrix.R
#Script to convert RDS file to matrix
#Script by Aidee Lashmi and paulinapglz.99@gmail.com

#libraries --- ---

library(tidyverse)
library(Matrix)

#Get data  --- ---
#Read Mutual information matrix

mi_list <- read_rds('/datos/rosmap/coexpre_matrix/')

#Unneest RDS --- ---

coexpre_mat_unnested <- as.data.frame(do.call(cbind, mi_list)) %>% as.matrix() 

mm2 <- mm

mm2[upper.tri(mm2)]  <- mm2[lower.tri(mm2)]

#mm2   %>% as_tibble() %>% pull(1) %>% unlist() %>% as.numeric() %>% length() #vroom::vroom_write(x = "/datos/rosmap/matriz_coexpre_20231011.txt")

mm3 <- mm2   %>% as_tibble() %>% mutate(V1 = colnames(.)) %>% pivot_longer(cols = -V1, names_to = "V2", values_to = "value") 
mm3 <- mm3 %>% mutate(value = mm3$value %>% sapply(FUN = function(i){ifelse(is.null(i), 0, i)}))

#Write mew matrix

mm3 %>% vroom::vroom_write("matriz_coexpre_allMCI_11052023.txt") #change name file every time

mm3 %>% 
  rename(gene = V1) %>% 
  pivot_wider(id_cols = gene, names_from = V2, values_from = value) %>% 
  vroom::vroom_write("matriz_coexpre_allMCI_11052023.txt")
