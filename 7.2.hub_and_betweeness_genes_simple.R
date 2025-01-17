#
#7.2.hub_and_betweeness_genes.R
#This script makes an analysis of hub genes and high betweeness genes

#Libraries --- ---

pacman::p_load('igraph',
               'ggplot2', 
               'dplyr', 
               'gridExtra', 
               "ggVennDiagram", 
               "clusterProfiler", 
               "ggraph", 
               "biomaRt")

library("org.Hs.eg.db", character.only = TRUE)

#Function to translate genes --- ---
# Connecting to the Ensembl database through biomaRt
#mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://asia.ensembl.org")

convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart)
}

convert_symbol_to_ens <- function(gene_ids) {
  # Separar los IDs que ya son Ensembl (asumen formato típico como "ENSG000001...")
  ensembl_ids <- gene_ids[grepl("^ENSG[0-9]+", gene_ids)]
  
  # Identificar los IDs que no son Ensembl
  symbols <- gene_ids[!gene_ids %in% ensembl_ids]
  
  # Realizar la conversión solo para los símbolos
  if (length(symbols) > 0) {
    converted <- getBM(
      attributes = c("external_gene_name", "ensembl_gene_id"),
      filters = "external_gene_name",
      values = symbols,
      mart = mart
    )
    
    # Crear un diccionario para la conversión
    conversion_dict <- setNames(converted$ensembl_gene_id, converted$external_gene_name)
    
    # Reemplazar los símbolos con sus Ensembl IDs
    symbols_converted <- conversion_dict[symbols]
  } else {
    symbols_converted <- character(0)  # Si no hay símbolos, devolver vacío
  }
  
  # Combinar los Ensembl IDs originales y los convertidos
  result <- c(ensembl_ids, symbols_converted)
  
  # Mantener el orden original
  result <- result[match(gene_ids, c(ensembl_ids, symbols))]
  
  # Devolver el resultado
  return(result)
}


translate_vertex_names <- function(graph) {
  #Extract vertex names
  graph_vnames <- V(graph)$name
  
  #Translate names
  graph_vnames_trad <- convert_ens_to_symbol(graph_vnames)
  
  #Replace the missing values in the column 'external_gene_name' with the values of 'ensembl_gene_id'.
  graph_vnames_trad$external_gene_name <- ifelse(graph_vnames_trad$external_gene_name == "", graph_vnames_trad$ensembl_gene_id, graph_vnames_trad$external_gene_name)
  
  # Create a vector of translated names using the dictionary
  # We need to ensure that the actual names of the network are in the dictionary
  graph_vnames_trad <- setNames(graph_vnames_trad$external_gene_name, graph_vnames_trad$ensembl_gene_id)
  
  # Sort graph_vnames_trad according to the order of graph_vnames
  sorted_graph_vnames_trad <- graph_vnames_trad[match(graph_vnames, names(graph_vnames_trad))]
  #Assign the new names to the network vertices.
  V(graph)$name_trad <- sorted_graph_vnames_trad
  
  return(graph)
}

# Función para generar subgrafos inducidos
generate_induced_subgraphs <- function(graph, nodes_of_interest, exclusive_nodes) {
  subgraphs <- list()
  
  for (node in nodes_of_interest) {
    # Crear subgrafo inducido
    neighbors <- neighbors(graph, node)
    induced_nodes <- c(node, names(neighbors))
    subgraph <- induced_subgraph(graph, vids = induced_nodes)
    
    # Etiquetar nodos exclusivos
    V(subgraph)$label <- ifelse(V(subgraph)$name %in% exclusive_nodes, "AD-exclusive", "Shared")
    
    # Guardar subgrafo
    subgraphs[[node]] <- subgraph
  }
  
  return(subgraphs)
}


#Function to find new coexpressed genes
find_new_coexpressed_genes <- function(subgraph, control_graph, main_node) {
  # Vecinos en el subgrafo
  neighbors_AD <- neighbors(subgraph, main_node)
  neighbors_AD_names <- V(subgraph)[neighbors_AD]$name_trad
  
  # Vecinos en la red de control
  neighbors_noAD <- neighbors(control_graph, main_node)
  neighbors_noAD_names <- V(control_graph)[neighbors_noAD]$name
  
  # Nuevos genes coexpresados
  new_genes <- setdiff(neighbors_AD_names, neighbors_noAD_names)
  
  # Devolver resultados
  return(new_genes)
}

#Get data --- ---

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml',
                        format = 'graphml')

graphs <- list(graphAD = graphAD, 
               graphnoAD = graphnoAD)

universe <- scan(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/universe.txt", what = character())

##### HUB GENES ####

#Calculate degree of nodes
nodes_degree <- sapply(X = graphs, FUN = degree)

#Threshold
cutoff_degrees <- lapply(nodes_degree, function(deg) {
  quantile(deg, probs = 0.90)
})

#Hub nodes

hub_nodes <- list()
for (name in names(graphs)) {
  deg <- nodes_degree[[name]]
  cutoff <- cutoff_degrees[[name]]
  hub_nodes[[name]] <- V(graphs[[name]])[deg >= cutoff]$name
}

#What genes do both networks share?

shared_hub_genes <- intersect(hub_nodes[[1]], hub_nodes[[2]])

#What genes are in the AD network but not in the no AD network?
control_hubs_notAD <- setdiff(hub_nodes[[2]], hub_nodes[[1]])

convert_ens_to_symbol(control_hubs_notAD)


AD_hubs_notcontrol <- setdiff(hub_nodes[[1]], hub_nodes[[2]])
length(AD_hubs_notcontrol)

AD_hubs_notcontrol.df1 <- data.frame(gene = names(nodes_degree$graphAD[names(nodes_degree$graphAD) %in% 
                                                                        AD_hubs_notcontrol]),
                                    degree = nodes_degree$graphAD[names(nodes_degree$graphAD) %in% 
                                                                    AD_hubs_notcontrol]) 

AD_hubs_notcontrol.df <- convert_ens_to_symbol(AD_hubs_notcontrol.df1$gene)


AD_hubs_notcontrol.df <- AD_hubs_notcontrol.df %>%
  left_join(AD_hubs_notcontrol.df1, by = c("ensembl_gene_id" = "gene")) %>% arrange(desc(degree))

#I want to know what is the degree of this same nodes in the control network

degree_in_control <- degree(graph = graphnoAD, v = AD_hubs_notcontrol)
degree_in_AD <- degree(graph = graphAD, v = AD_hubs_notcontrol)

#Calculating the empirical cumulative distribution function (ECDF)
ecdf_degree <- ecdf(degree(graphnoAD))  # Basado en el grado de todos los nodos en la red

# Calcular el percentil para cada nodo
percentile_ecdf <- ecdf_degree(as.numeric(degree_in_control)) * 100  # Asegurarse de que sean números

#df
degree_in_control <- data.frame(
  ensembl_gene_id = names(degree_in_control), 
  degree_in_control = degree_in_control, 
  degree_percentile = percentile_ecdf, 
  degree_percentile1 = (100 - percentile_ecdf), 
  difference = degree
  
)

degree_in_control[2,5] <- "PDE4DIP"

degree_in_control.x <- convert_ens_to_symbol(degree_in_control)

degree_in_control <- degree_in_control %>% left_join(degree_in_control.x)

degree_in_AD_df <- data.frame(
  ensembl_gene_id = names(degree_in_AD), 
  degree_in_AD = degree_in_AD
)

# Combinar degree_in_AD con la tabla degree_in_control
degree_changes <- degree_in_control %>%
  left_join(degree_in_AD_df, by = "ensembl_gene_id")  %>%
  mutate(degree_difference = degree_in_AD - degree_in_control) %>% 
  arrange(by = desc(degree_difference))

# Convertir el dataframe a formato largo
degree_changes.l <- degree_changes %>%
  pivot_longer(
    cols = c(degree_in_AD, degree_in_control), 
    names_to = "network", 
    values_to = "degree"
  )


#Plot
degree_changes.p <- ggplot(degree_changes.l, aes(x = factor(external_gene_name, levels = degree_changes$external_gene_name), y = degree,
                             color = network )) +
  geom_point() + 
  geom_segment(aes(y=0, yend=degree), size = 1) +
  geom_text(data = subset(degree_changes.l, network == "degree_in_AD"), 
            aes(label = degree_difference), vjust = -0.5, size = 3, color = "black") +# Agregar texto encima de los puntos
  #coord_polar(theta = "x") +                                    # Coordenadas polares para gráfico circular
  theme_minimal() +                                             # Tema limpio
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),          # Rotar etiquetas para mejor lectura
    legend.position = "bottom"                                  # Leyenda abajo
  ) +
  labs(
    title = "", 
    x = NULL, 
    y = "Degree",
    fill = "Red"
  )  +
  scale_color_manual(values = c("degree_in_AD" = "red", "degree_in_control" = "blue"))

#save graph
ggsave("degree_changes.jpg",
       plot = degree_changes.p, 
       device = "jpg",
       width = 8, #55
       height = 10, #30
       units = "in", 
       dpi = 300
)

#I want to know something about CROCC

# # Especificar el nodo de interés
# CROCC <- "ENSG00000226321"
# 
# #Neighbours
# CROCC_neigh <- neighbors(graphAD, CROCC)
# CROCC_neigh_names <- V(graphAD)$name[CROOC_neigh]
# 
# #Create subgraph
# CROCC_neigh.s <- induced_subgraph(graphAD, c("ENSG00000226321", CROCC_neigh_names))
# 
# # Dibujar el subgrafo
# plot(CROCC_neigh.s, vertex.size = 30, vertex.label.cex = 1.5)

#Generate induced subgraphs of the hub genes and their neighbours
AD_hubs_notcontrol.x <- hub_nodes$graphAD[!(hub_nodes$graphAD %in% hub_nodes$graphnoAD)]
subgraphs_AD <- generate_induced_subgraphs(
  graph = graphAD,
  nodes_of_interest = AD_hubs_notcontrol,
  exclusive_nodes = AD_hubs_notcontrol
)
# Aplicar la traducción de nombres a todos los subgrafos
subgraphs_AD <- lapply(subgraphs_AD, translate_vertex_names)

#Save subgraphs
for (node in names(subgraphs_AD)) {
  write_graph(subgraphs_AD[[node]], file = paste0("subgraph_", node, ".graphml"), format = "graphml")
}

#Find new coexpressed subgraphs --- ---

#Create a list to store results
new_coexpressed_genes.l <- list()

#Iterate for each subgraph
for (main_node in names(subgraphs_AD)) {
  # Obtener el subgrafo
  subgraph <- subgraphs_AD[[main_node]]
  
  # Encontrar nuevos genes coexpresados
  new_genes <- find_new_coexpressed_genes(
    subgraph = subgraph,
    control_graph = graphnoAD,
    main_node = main_node
  )
  
  # Concatenar los genes separados por comas
  new_genes_combined <- paste(new_genes, collapse = ", ")
  
  # Almacenar los resultados
  new_coexpressed_genes.l[[main_node]] <- data.frame(
    Main_Node = main_node,
    New_Genes = new_genes_combined
  )
}

# Combinar todos los resultados en una sola tabla
new_coexpressed_genes.df <- do.call(rbind, new_coexpressed_genes.l)

# Separar los genes y contar su frecuencia
new_coexpressed_genes.fr <- results_table %>%
  # Dividir la columna New_Genes en genes individuales
  dplyr::mutate(New_Genes = strsplit(New_Genes, ", ")) %>%
  # Expandir la lista en filas
  tidyr::unnest(New_Genes) %>%
  # Contar las apariciones de cada gen
  dplyr::count(New_Genes, sort = TRUE)

(new_coexpressed_genes.fr)

# Gráfico de barras para los 10 genes más frecuentes
new_coexpressed_genes.fr %>%
 # head(50) %>%
  ggplot(aes(x = reorder(New_Genes, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +
  labs(
    title = "Top 10 Genes Más Frecuentes",
    x = "Genes",
    y = "Frecuencia"
  ) +
  theme_minimal()

#Enrichment of common newly expressed genes

new_coexpressed_genes.fr$New_Genes_ense <- convert_symbol_to_ens(new_coexpressed_genes.fr$New_Genes)

new_coexpressed_genes.enr <-  enrichGO(
  gene = gene_frequencies$New_Genes_ense,
  OrgDb = "org.Hs.eg.db", 
  keyType = 'ENSEMBL',
  readable = TRUE,
  universe = universe,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

new_coexpressed_genes.cnet <- cnetplot(gene_frequencies_enr, showCategory = 20)



#Venn Diagram

hubsvenn <- ggVennDiagram(list(AD = hub_nodes[[1]], C =  hub_nodes[[2]])) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "a) Shared Hub genes") +
  theme( plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
hubsvenn

#Enrichment of these genes

AD_hubs_notcontrol_enrichment <- enrichGO(
  gene = AD_hubs_notcontrol,
  OrgDb = "org.Hs.eg.db", 
  keyType = 'ENSEMBL',
  readable = TRUE,
  universe = universe,
  ont = "MF",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

AD_hubs_notcontrol_enrichment.df <- as.data.frame(AD_hubs_notcontrol_enrichment)

#Plot enrichment of hubs

cnetplot(AD_hubs_notcontrol_enrichment, circular = F, colorEdge = TRUE) 

#Induced subgraph of hub genes AD_hubs_notcontrol

V(graphAD)$degree <- nodes_degree[["graphAD"]]

hub_subgraph_AD <- induced_subgraph(graph = graphAD, 
                                 vids = V(graphAD)[name %in% hub_nodes[["graphAD"]]])

plot(hub_subgraph_AD)

hub_subgraph_AD.p <- ggraph(hub_subgraph_AD, layout = 'kk') + 
  geom_edge_link(color = "grey89", show.legend = FALSE) +  # Dibujar las aristas
  geom_node_point(aes(color = degree), size = 5) +  # Dibujar y colorear los nodos
  scale_color_gradient("Degree",low = "lightblue", high = "red") +  # Escala de color para betweenness
  geom_node_text(aes(label = name_trad), repel = TRUE, size = 3, color = "black") +  # Etiquetar los nodos con sus nombres
  ggtitle("a)") +
  theme_void() 
hub_subgraph_AD.p 

#Induced subgraph of hub genes AD_hubs_notcontrol

hub_subgraph <- induced_subgraph(graph = graphAD, 
                                 vids = V(graphAD)[name %in% AD_hubs_notcontrol])

hub_subgraph.p <- ggraph(hub_subgraph, layout = 'kk') + 
  geom_edge_link(color = "grey89", show.legend = FALSE) +  # Dibujar las aristas
  geom_node_point(aes(color = degree), size = 5) +  # Dibujar y colorear los nodos
  scale_color_gradient("Degree", low = "lightblue", high = "red") +  # Escala de color para betweenness
  geom_node_text(aes(label = name_trad), repel = TRUE, size = 3, color = "black") +  # Etiquetar los nodos con sus nombres
  ggtitle("c)") +
  theme_void() 
hub_subgraph.p

#Save graph 

write_graph(hub_subgraph, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_hubs_induced_subgraph.graphml',
             format = "graphml")

#Enrichment with KEGG

kk <- enrichKEGG(gene = AD_hubs_notcontrol,
                 organism = 'hsa',
                 keyType = 'ENSEMBL',
                 pvalueCutoff = 0.05, 
                 universe = universe)

##### HIGH BETWEENESS GENES ####

#Calculate degree of nodes
nodes_betwe <- sapply(X = graphs, FUN = betweenness)

#Threshold
cutoff_betwe <- lapply(nodes_betwe, function(x) {
  quantile(x, probs = 0.88)
})

highbe_nodes <- list()
for (name in names(graphs)) {
  betwe <- nodes_betwe[[name]]
  cutoff <- cutoff_betwe[[name]]
  highbe_nodes[[name]] <- V(graphs[[name]])[betwe >= cutoff]$name
}

#What genes do both networks share?

shared_highbe_genes <- intersect(highbe_nodes[[1]], highbe_nodes[[2]])
shared_highbe_genes
#What genes are in the AD network but not in the no AD network?

AD_highbe_notcontrol <- setdiff(highbe_nodes[[1]], highbe_nodes[[2]])
length(AD_highbe_notcontrol)

#Venn Diagram

highbevenn <- ggVennDiagram(list(AD = highbe_nodes[[1]], C =  highbe_nodes[[2]])) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "b) Shared high betweeness genes") +
  theme( plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

#Enrichment of these genes

AD_highbe_notcontrol_enrichment <- enrichGO(
  gene = AD_highbe_notcontrol,
  OrgDb = "org.Hs.eg.db", 
  keyType = 'ENSEMBL',
  readable = TRUE,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

AD_highbe_notcontrol_enrichment.df <- as.data.frame(AD_highbe_notcontrol_enrichment)

vroom::vroom_write(AD_highbe_notcontrol_enrichment.df,
                   "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_highbe_notcontrol_enrichment.csv")

#Plot enrichment of hubs

AD_highbe_notcontrol_enrichment.p <- cnetplot(AD_highbe_notcontrol_enrichment, circular = TRUE, colorEdge = TRUE) +
  guides(edge_color = "none") +
  ggtitle("e)") 

ggsave("AD_highbe_notcontrol_enrichment.jpg",
       plot = AD_highbe_notcontrol_enrichment.p, 
       device = "jpg",
       width = 15, #55
       height = 10, #30
       units = "in", 
       dpi = 300
)

#Induced subgraph of hub genes

V(graphAD)$betweenness <- nodes_betwe[["graphAD"]]

highbe_subgraph_AD <- induced_subgraph(graph = graphAD, 
                                    vids = V(graphAD)[name %in% highbe_nodes[["graphAD"]]])

highbe_subgraph_AD.p <- ggraph(highbe_subgraph_AD, layout = 'kk', maxiter = 1000) + 
  geom_edge_link(color = "grey89", show.legend = FALSE) +  # Dibujar las aristas
  geom_node_point(aes(color = betweenness), size = 5) +  # Dibujar y colorear los nodos
  scale_color_gradient("Betweenness centrality", low = "lightblue", high = "red") +  # Escala de color para betweenness
  geom_node_text(aes(label = name_trad), repel = TRUE, size = 3, color = "black") +  # Etiquetar los nodos con sus nombres
  ggtitle("c)") +
  theme_void() 

highbe_subgraph_AD.p
#Genes only in AD graph

highbeb_subgraph <- induced_subgraph(graph = graphAD, 
                                 vids = V(graphAD)[name %in% AD_highbe_notcontrol])

plot(highbeb_subgraph)

#Save graph

write_graph(highbeb_subgraph, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_highbe_induced_subgraph.graphml',
            format = "graphml")

#Only in AD

highbeb_subgraph.p <- ggraph(highbeb_subgraph, layout = 'kk') + 
  geom_edge_link(aes(alpha = 0.8), color = "grey89", show.legend = FALSE) +  # Dibujar las aristas
  geom_node_point(aes(color = betweenness), size = 5) +  # Dibujar nodos sin leyenda
  scale_color_gradient("Betweenness", low = "lightblue", high = "red") +  # Eliminar la leyenda del color
  geom_node_text(aes(label = name_trad), repel = TRUE, size = 3, color = "black") +  # Etiquetar los nodos con sus nombres
  ggtitle("d)") +
  theme_void() +
  theme(legend.title = element_text(hjust = 0.9))

#Plot graphs and Venn diagrams in a grid

lay <- rbind(c(3,1,1,2,2),  # Primera fila con tres gráficos
             c(5,4,4,2,2))  # Segunda fila con dos gráficos más, debajo del gráfico central

 grid <- grid.arrange(hub_subgraph.p,
                      highbeb_subgraph.p, 
                      hubsvenn,
                      AD_highbe_notcontrol_enrichment.p, 
                      highbevenn,
                      layout_matrix = lay)

 ggsave("grid_hubs_and_highbe1.jpg",
        plot = grid, 
        device = "jpg",
        width = 18, #55
        height = 12, #30
        units = "in", 
        dpi = 300)
 
 #NEXT QUESTION IS What modules do these genes belong to? --- ---

#Save hub genes to explore them in the partitions

write(AD_hubs_notcontrol,
      file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_hubs_notcontrolens.txt')

#Save high betweeness genes to explore them in the partitions

write(AD_highbe_notcontrol,
      file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_highbe_notcontrol_ens.txt')

#END