
#Libraries  --- --- 

pacman::p_load('tidyverse', 
               'igraph', 
               'ggplot2')

#edge betweenneess by node

betweeness_by_node_min_sin_g <-betweenness(min_sin_comp_graph, directed = F) %>% as.data.frame()

png("min_sin_comp_MI_betweeness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_min_sin_g$.), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Betweeness by node histogram",   # Título principal
     sub = "patients with AD, minimal single component",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()