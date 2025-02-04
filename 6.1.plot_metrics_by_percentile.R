#
#Plot 4.1.network_characteristics.R

#Libraries --- ---

pacman::p_load('dplyr',
               'ggplot2', 
               'gridExtra', 
               'stringr')

#Get ddata --- ---

metrics_AD <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/metrics_percentiles_normalizedMI_AD_ROSMAP_RNAseq_DLFPC_MutualInfo_NIA_Reagan_dicho.txt')

metrics_noAD <- vroom::vroom(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/metrics_percentiles_normalizedMI_noAD_ROSMAP_RNAseq_DLFPC_MutualInfo_NIA_Reagan_dicho.txt')

#Manage data --- ---

#For AD
metrics_AD <- metrics_AD %>% 
  mutate(dx = 'AD', .before = 1) 

#For noAD
metrics_noAD <- metrics_noAD %>% 
  mutate(dx = 'noAD', .before = 1) 

#Merge tables

metrics <- rbind(metrics_AD, metrics_noAD)

metrics <- metrics %>%
  mutate(
    # Extrae los números después del guion bajo
    raw_percentile = str_extract(percentile_no, "\\d+"),
    # Determina cuántos dígitos hay y dónde agregar el punto decimal
    percentile = case_when(
      nchar(raw_percentile) == 2 ~ paste0(raw_percentile, "%"),
      nchar(raw_percentile) == 3 ~ paste0(substr(raw_percentile, 1, 2), ".", substr(raw_percentile, 3, 3), "%"),
      nchar(raw_percentile) == 4 ~ paste0(substr(raw_percentile, 1, 2), ".", substr(raw_percentile, 3, 4), "%"),
      nchar(raw_percentile) == 5 ~ paste0(substr(raw_percentile, 1, 2), ".", substr(raw_percentile, 3, 5), "%"),
      nchar(raw_percentile) == 6 ~ paste0(substr(raw_percentile, 1, 2), ".", substr(raw_percentile, 3, 6), "%")
    )
  ) %>%
  select(-raw_percentile) %>%   
  mutate(
    # Convertir percentile en un factor ordenado
    percentile = factor(percentile, levels = c("85%", "90%", "98%", "99%", "99.9%", "99.99%", "99.999%", "99.9999%"), ordered = TRUE)
  )


#Plot --- ---

#Plot number of vertices 

v_plot <- ggplot(metrics, aes(x = percentile, y = length_v, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), linewidth = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Number of vertices (genes)",
       x = "Percentile",
       y = "Number of vertices") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot number of edges 
E_plot <- ggplot(metrics, aes(x = percentile, y = log(length_E), color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Number of edges",
       x = "Percentile",
       y = "Number of edges (log)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot number of components

components_plot <- ggplot(metrics, aes(x = percentile, y = components_no, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Number of components",
       x = "Percentile",
       y = "Number of components") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot number of clusters with INFOMAP

components_INFOMAP_plot <- ggplot(metrics, aes(x = percentile, y = no_cluster_infomap, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Number of clusters (Infomap)",
       x = "Percentile",
       y = "Number of clusters") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot clustering coefficient

clusteringcoe_plot <- ggplot(metrics, aes(x = percentile, y = clustering_coefficient, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") + # You can choose different palettes
  labs(title = "Clustering coefficient",
       x = "Percentile",
       y = "Clustering coefficient ") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot max weight

max_weight_plot <- ggplot(metrics, aes(x = percentile, y = max_MI, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Maximal weight",
       x = "Percentile",
       y = "Max weight") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot min weight

min_weight_plot <- ggplot(metrics, aes(x = percentile, y = min_MI, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +   # You can choose different palettes
  labs(title = "Minimal weight",
       x = "Percentile",
       y = "Min weight") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot max comp size

max_comp_size_plot <- ggplot(metrics, aes(x = percentile, y = max_comp_size, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "Maximal component size",
       x = "Percentile",
       y = "Larger component size") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Plot percentage of genes in the max size comp

per_genes_max_size_plot <- ggplot(metrics, aes(x = percentile, y = percentage_genes_in_larger_component, color = dx)) +
  geom_point(size = 3) +
  geom_line(aes(group = dx), size = 1) +
  scale_y_continuous(limits = c(min(metrics$percentage_genes_in_larger_component), 100), breaks = seq(50, 100, by = 20)) +
  scale_color_brewer(palette = "Set1") +  # You can choose different palettes
  labs(title = "% of genes in larger component",
       x = "Percentile",
       y = "% of genes in larger component") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1))

#Arrange in a grid for presentation

grid <- grid.arrange(v_plot, E_plot, 
             components_plot, clusteringcoe_plot,
             max_weight_plot, min_weight_plot,
             max_comp_size_plot, per_genes_max_size_plot, 
             components_INFOMAP_plot,
             ncol = 3)

ggsave(filename = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/network_characteristics_bypercentile.jpg',
       plot = grid,  # Esto asume que el último gráfico generado es el que quieres guardar
       width = 30,
       height = 15,
       units = "in",
       dpi = 600,
       pointsize = 4)


#END