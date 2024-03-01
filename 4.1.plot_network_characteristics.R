#
#Plot 4.1.network_characteristics.R

#Libraries --- ---

pacman::p_load('dplyr',
               'ggplot2', 
               'gridExtra')

#Get ddata --- ---

metrics_AD <- vroom::vroom(file = '/datos/rosmap/cuts_by_MI/AD_graphs/metrics_percentiles_AD_ROSMAP_RNAseq_MutualInfo_NIA_Reagan_dicho.txt')

metrics_noAD <- vroom::vroom(file = '/datos/rosmap/cuts_by_MI/noAD_graphs/metrics_percentiles_noAD_ROSMAP_RNAseq_MutualInfo_NIA_Reagan_dicho.txt')

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

v_plot <- ggplot(metrics, aes(x = percentile_no, y = length_v, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # Puedes elegir diferentes paletas
  labs(title = "Number of vertices per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Number of vertices") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot number of edges 
E_plot <- ggplot(metrics, aes(x = percentile_no, y = length_E, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # Puedes elegir diferentes paletas
  labs(title = "Number of edges per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Number of edges") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot number of clusters

cluster_plot <- ggplot(metrics, aes(x = percentile_no, y = clusters_no, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # Puedes elegir diferentes paletas
  labs(title = "Number of clusters per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Number of clusters") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot clustering_coefficient

clusteringcoe_plot <- ggplot(metrics, aes(x = percentile_no, y = clustering_coefficient, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # Puedes elegir diferentes paletas
  labs(title = "Clustering coefficient per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Clustering coefficient ") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot max weight

max_weight_plot <- ggplot(metrics, aes(x = percentile_no, y = max_weight, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # Puedes elegir diferentes paletas
  labs(title = "Max weight per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Max weight") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot max weight

min_weight_plot <- ggplot(metrics, aes(x = percentile_no, y = min_weight, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # Puedes elegir diferentes paletas
  labs(title = "Minimal weight per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Min weight") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Arrange in a grid for presentation

grid.arrange(v_plot, E_plot, cluster_plot, clusteringcoe_plot, max_weight, min_weight,ncol = 2)

#END