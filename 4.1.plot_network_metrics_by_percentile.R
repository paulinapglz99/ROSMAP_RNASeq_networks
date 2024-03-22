#
#Plot 4.1.network_characteristics.R

#Libraries --- ---

pacman::p_load('dplyr',
               'ggplot2', 
               'gridExtra')

#Get ddata --- ---

metrics_AD <- vroom::vroom(file = '/datos/rosmap/cuts_by_MI/AD_graphs/metrics_percentiles_normalizedMI_AD_ROSMAP_RNAseq_MutualInfo_NIA_Reagan_dicho.txt')

metrics_noAD <- vroom::vroom(file = '/datos/rosmap/cuts_by_MI/noAD_graphs/metrics_percentiles_normalizedMI_noAD_ROSMAP_RNAseq_MutualInfo_NIA_Reagan_dicho.txt')

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
  geom_line(aes(group = dx), linewidth = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
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
E_plot <- ggplot(metrics, aes(x = percentile_no, y = log(length_E), color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Number of edges per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Number of edges, in log scale") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot number of components

components_plot <- ggplot(metrics, aes(x = percentile_no, y = components_no, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Number of components per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Number of components") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot number of components with INFOMAP


components_INFOMAP_plot <- ggplot(metrics, aes(x = percentile_no, y = no_cluster_infomap, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Number of components calculated with INFOMAP per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Number of components calculated with INFOMAP") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot clustering coefficient

clusteringcoe_plot <- ggplot(metrics, aes(x = percentile_no, y = clustering_coefficient, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") + # You can choose different palettes
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

max_weight_plot <- ggplot(metrics, aes(x = percentile_no, y = max_MI, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Max weight per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Max weight") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot min weight

min_weight_plot <- ggplot(metrics, aes(x = percentile_no, y = min_MI, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +   # You can choose different palettes
  labs(title = "Minimal weight per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Min weight") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot max comp size

max_comp_size_plot <- ggplot(metrics, aes(x = percentile_no, y = max_comp_size, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Maximal component size per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Larger component size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot percentage of genes in the max size comp

per_genes_max_size_plot <- ggplot(metrics, aes(x = percentile_no, y = percentage_genes_in_larger_component, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Percentage of genes in larger component per cut percentile",
       subtitle = "by pathological diagnosis of Alzheimer's disease",
       x = "Percentile number",
       y = "Percentage of genes in larger component ") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Arrange in a grid for presentation

png('network_normalized_characteristics_bypercentile.png',
    width     = 35,
    height    = 15,
    units     = "in",
    res       = 600,
    pointsize = 4 )
grid.arrange(v_plot, E_plot, 
             components_plot, clusteringcoe_plot,
             max_weight_plot, min_weight_plot,
             max_comp_size_plot, per_genes_max_size_plot, 
             components_INFOMAP_plot,
             ncol = 4)
dev.off()

#END