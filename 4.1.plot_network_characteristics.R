#
#Plot 4.1.network_characteristics.R

#Libraries --- ---

pacman::p_load('dplyr',
               'ggplot2')

#Get ddata --- ---

metrics_AD <- vroom::vroom(file = '/datos/rosmap/cuts_by_MI/noAD_graphs/metrics_percentiles_noAD_ROSMAP_RNAseq_MutualInfo_NIA_Reagan_dicho.txt')

metrics_noAD <- vroom::vroom(file = '/datos/rosmap/cuts_by_MI/noAD_graphs/metrics_percentiles_AD_ROSMAP_RNAseq_MutualInfo_NIA_Reagan_dicho.txt')

#Manage data --- ---

#For AD
metrics_AD <- metrics_AD %>% 
    mutate(dx = 'AD', .before = 1) 

#For noAD
metrics_noAD <- metrics_noAD %>% 
  mutate(dx = 'noAD', .before = 1) 

#Merge tables

metrics <- rbind(metrics_AD, metrics_noAD)

#Plot --- ---

#Plot number of vertices 

ggplot(datos, aes(x = percentile_no, y = length_v)) +
  geom_point() +
  labs(title = "Gráfico de Dispersión: length_v vs. percentile_no",
       x = "Percentile number",
       y = "Number of vertices") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

