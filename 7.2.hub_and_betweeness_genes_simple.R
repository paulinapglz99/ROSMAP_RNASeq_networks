#
#7.2.hub_and_betweeness_genes.R
#This script makes an analysis of hub genes and high betweeness genes

#Libraries --- ---

pacman::p_load('igraph',
               'ggplot2', 
               'tidyverse', 
               'gridExtra', 
               "ggVennDiagram", 
               "clusterProfiler", 
               "ggraph", 
               "biomaRt")

library("org.Hs.eg.db", character.only = TRUE)

#Function to translate genes --- ---
# Connecting to the Ensembl database through biomaRt
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
#If this does not work, try mirrors
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", 
#                 host = "https://useast.ensembl.org")

#Function to convert ensembl IDs to symbols
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = mart)
}

#Function to convert symbols to ensembl IDs

convert_symbol_to_ens <- function(gene_ids) {
  # Separate IDs that are already Ensembl (assume typical format such as “ENSG000001...”)
  ensembl_ids <- gene_ids[grepl("^ENSG[0-9]+", gene_ids)]
  
  #Identify the IDs that are not Ensembl
  symbols <- gene_ids[!gene_ids %in% ensembl_ids]
  
  #Perform conversion for symbols only
  if (length(symbols) > 0) {
    converted <- getBM(
      attributes = c("external_gene_name", "ensembl_gene_id"),
      filters = "external_gene_name",
      values = symbols,
      mart = mart
    )
    #Create a dictionary for conversion
    conversion_dict <- setNames(converted$ensembl_gene_id, converted$external_gene_name)
    
    #Replace symbols with their Ensembl IDs
    symbols_converted <- conversion_dict[symbols]
  } else {
    symbols_converted <- character(0)  # Si no hay símbolos, devolver vacío
  }
  
  #Combining the original and converted Ensembl IDs
  result <- c(ensembl_ids, symbols_converted)
  
  #Maintain original order
  result <- result[match(gene_ids, c(ensembl_ids, symbols))]
  
  return(result)
}

#Function to translate graph vertex names

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

#Function to create induced subgraphs
generate_induced_subgraphs <- function(graph, nodes_of_interest, exclusive_nodes) {
  subgraphs <- list()
  
  for (node in nodes_of_interest) {
    #Create induced subgrapgh
    neighbors <- neighbors(graph, node)
    induced_nodes <- c(node, names(neighbors))
    subgraph <- induced_subgraph(graph, vids = induced_nodes)
    #Label exclusive nodes
    V(subgraph)$label <- ifelse(V(subgraph)$name %in% exclusive_nodes, "AD-exclusive", "Shared")
    #Save subgraphs
    subgraphs[[node]] <- subgraph
  }
  
  return(subgraphs)
}

#Function to find new coexpressed genes
find_new_coexpressed_genes <- function(subgraph, control_graph, main_node) {
  #Neighbors in subgraph
  neighbors_AD <- neighbors(subgraph, main_node)
  neighbors_AD_names <- V(subgraph)[neighbors_AD]$name_trad
  #Neighbors in the control net
  neighbors_noAD <- neighbors(control_graph, main_node)
  neighbors_noAD_names <- V(control_graph)[neighbors_noAD]$name
  #New coexpressed genes
  new_genes <- setdiff(neighbors_AD_names, neighbors_noAD_names)
  
  return(new_genes)
}

#Get data --- ---

graphAD <- read_graph(file = 'ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = 'ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml',
                        format = 'graphml')

graphs <- list(graphAD = graphAD,  #graph1
               graphnoAD = graphnoAD) #graph2

universe <- scan(file = "universe.txt", what = character())

##### HUB GENES ####

#Calculate degree of nodes
nodes_degree <- sapply(X = graphs, FUN = degree)

#Threshold
cutoff_degrees <- lapply(nodes_degree, function(deg) {
  quantile(deg, probs = 0.90)  #to get the 10% of higher degree nodes
})

#Hub nodes

hub_nodes <- list()
for (name in names(graphs)) {
  deg <- nodes_degree[[name]]
  cutoff <- cutoff_degrees[[name]]
  hub_nodes[[name]] <- V(graphs[[name]])[deg >= cutoff]$name
}

#What genes do both networks share?

shared_hub_genes <- intersect(hub_nodes[[1]],  
                              hub_nodes[[2]])  #but not in control
shared_hub_genes

#What genes are in the control network but not in the AD network?

control_hubs_notAD <- setdiff(hub_nodes[[2]], #elements in control
                              hub_nodes[[1]]) #But not in AD

#What genes are in the AD network but not in the no AD network?

AD_hubs_notcontrol <- setdiff(hub_nodes[[1]],  #elements in AD
                              hub_nodes[[2]])  #but not in control
length(AD_hubs_notcontrol)

#As data frames
AD_hubs_notcontrol.df1 <- data.frame(gene = names(nodes_degree$graphAD[names(nodes_degree$graphAD) %in% 
                                                                        AD_hubs_notcontrol]),
                                    degree = nodes_degree$graphAD[names(nodes_degree$graphAD) %in% 
                                                                    AD_hubs_notcontrol]) 

AD_hubs_notcontrol.df <- convert_ens_to_symbol(AD_hubs_notcontrol.df1$gene)


AD_hubs_notcontrol.df <- AD_hubs_notcontrol.df %>%
  left_join(AD_hubs_notcontrol.df1, by = c("ensembl_gene_id" = "gene")) %>% 
  arrange(desc(degree))

#Venn Diagram --- --- 

hubsvenn <- ggVennDiagram(list(AD = hub_nodes[[1]], C =  hub_nodes[[2]])) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "a) Shared Hub genes") +
  theme( plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")
hubsvenn

#Enrichment of genes AD-only --- ---

AD_hubs_notcontrol_enrichment <- enrichGO(
  gene = AD_hubs_notcontrol,
  OrgDb = "org.Hs.eg.db", 
  keyType = 'ENSEMBL',
  readable = TRUE,
  universe = universe,
  ont = "MF",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

#As data frame
AD_hubs_notcontrol_enrichment.df <- as.data.frame(AD_hubs_notcontrol_enrichment)

#Plot enrichment of hubs

cnetplot(AD_hubs_notcontrol_enrichment, circular = F, colorEdge = TRUE) 

#Induced subgraph of hub genes AD_hubs_notcontrol

V(graphAD)$degree <- nodes_degree[["graphAD"]]

hub_subgraph_AD <- induced_subgraph(graph = graphAD, 
                                    vids = V(graphAD)[name %in% AD_hubs_notcontrol])
 
# plot(hub_subgraph_AD)
# 
# hub_subgraph_AD.p <- ggraph(hub_subgraph_AD, layout = 'kk') + 
#   geom_edge_link(color = "grey89", show.legend = FALSE) +  # Dibujar las aristas
#   geom_node_point(aes(color = degree), size = 5) +  # Dibujar y colorear los nodos
#   scale_color_gradient("Degree",low = "lightblue", high = "red") +  # Escala de color para betweenness
#   geom_node_text(aes(label = name_trad), repel = TRUE, size = 3, color = "black") +  # Etiquetar los nodos con sus nombres
#   ggtitle("a)") +
#   theme_void() 
# hub_subgraph_AD.p 

#Save graph 

# write_graph(hub_subgraph_AD, file = 'ROSMAP_RNAseq_DLPFC_AD_hubs_induced_subgraph.graphml',
#             format = "graphml")

#I want to know what is the degree of AD hub genes nodes in the control network

#Degree of only-AD hubs in the control network
degree_in_control <- degree(graph = graphnoAD, v = AD_hubs_notcontrol)
#Degree of only-AD hubs in AD network
degree_in_AD <- degree(graph = graphAD, v = AD_hubs_notcontrol)

#Calculating the empirical cumulative distribution function (ECDF)
ecdf_degree <- ecdf(degree(graphnoAD))  #Based in degree of all nodes in the network

#Percentile for each node
percentile_ecdf <- ecdf_degree(as.numeric(degree_in_control)) * 100  # Asegurarse de que sean números

#df
degree_dif <- data.frame(
  ensembl_gene_id = names(degree_in_control), 
  degree_in_control = degree_in_control, 
  degree_in_AD = degree_in_AD,
  degree_percentile = percentile_ecdf, 
  degree_percentile1 = (100 - percentile_ecdf), 
  degree_diff = degree_in_AD- degree_in_control
)

degree_diff.x <- convert_ens_to_symbol(degree_dif)

degree_dif <- degree_dif %>% 
  left_join(degree_diff.x) %>% 
  arrange(by = desc(degree_diff))
degree_dif$external_gene_name[13] <- "PDE4DIP"

#Long format for plot
degree_changes.l <- degree_dif %>%
  pivot_longer(
    cols = c(degree_in_AD, degree_in_control), 
    names_to = "network", 
    values_to = "degree"
  )

#Plot
degree_changes.p <- ggplot(degree_changes.l, aes(x = factor(external_gene_name, 
                                                            levels = degree_dif$external_gene_name), y = degree,
                             color = network )) +
  geom_point() + 
  geom_segment(aes(y=0, yend=degree), size = 1.5) +
  geom_text(data = subset(degree_changes.l, network == "degree_in_AD"), 
            aes(label = degree_diff), vjust = -0.5, size = 3, color = "black") +# Agregar texto encima de los puntos
  #coord_polar(theta = "x") +                                    # Coordenadas polares para gráfico circular
  theme_minimal() +                                             # Tema limpio
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),          # Rotar etiquetas para mejor lectura
    legend.position = "none"                                  # Leyenda abajo
  ) +
  labs(
    title = "", 
    x = NULL, 
    y = "Degree",
    fill = "Red"
  )  +
  scale_color_manual(values = c("degree_in_AD" = "red", 
                                "degree_in_control" = "blue"))
degree_changes.p

#save plot
# ggsave("degree_changes.jpg",
#        plot = degree_changes.p, 
#        device = "jpg",
#        width = 8, #55
#        height = 10, #30
#        units = "in", 
#        dpi = 300
# )

################################################

#I want to know what is the degree of control hub genes nodes in the AD network --- ---

#Degree of control-only hubs in control network
degree_in_control_control <- degree(graph = graphnoAD, v = control_hubs_notAD)
#Degree of control-only hubs in AD network
degree_in_AD_control <- degree(graph = graphAD, v = control_hubs_notAD)

#Calculating the empirical cumulative distribution function (ECDF)
#ecdf_degree <- ecdf(degree(graphnoAD))  #Based in degree of all nodes in the network

#Percentile for each node
percentile_ecdf <- ecdf_degree(as.numeric(degree_in_control_control)) * 100

#df
degree_dif_control <- data.frame(
  ensembl_gene_id = names(degree_in_control_control), 
  degree_in_control = degree_in_control_control, 
  degree_in_AD = degree_in_AD_control,
  degree_percentile = percentile_ecdf, 
  degree_percentile1 = (100 - percentile_ecdf), 
  degree_diff = degree_in_AD_control - degree_in_control_control
)

degree_diff_control.x <- convert_ens_to_symbol(degree_dif_control$ensembl_gene_id)

degree_dif_control <- degree_dif_control %>% 
  left_join(degree_diff_control.x) %>% 
  arrange(by = desc(degree_diff))
degree_dif_control$external_gene_name[7] <- "SNORD56"
degree_dif_control$external_gene_name[8] <- "HINT1P"
degree_dif_control$external_gene_name[9] <- "LOC100506405"

#Long format for plot
degree_changes_control.l <- degree_dif_control %>%
  pivot_longer(
    cols = c(degree_in_AD, degree_in_control), 
    names_to = "network", 
    values_to = "degree"
  )

#Plot
degree_changes_control.p <- ggplot(degree_changes_control.l, aes(x = factor(external_gene_name),
                                                            y = degree,
                                                 color = network)) +
  geom_point() + 
  geom_segment(aes(y = 0, yend = degree, color = network), size = 1.5) +
  geom_point(aes(color = network), size = 2) + 
  geom_segment(data = subset(degree_changes_control.l, network == "degree_in_AD"), 
               aes(y = 0, yend = degree), size = 1.5, color = "red") +
  geom_text(data = subset(degree_changes_control.l, network == "degree_in_control"), 
            aes(label = degree_diff), vjust = -0.5, size = 3, color = "black") + 
  theme_minimal() +                                             
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),        
    legend.position = "none"                                  
  ) +
  labs(
    title = "", 
    x = NULL, 
    y = "Degree",
    fill = "Red"
  )  +
  scale_color_manual(values = c("degree_in_AD" = "red", 
                                "degree_in_control" = "blue"))
degree_changes_control.p

################################################

#Generate induced subgraphs of the hub genes and their neighbors ---- ---

subgraphs_AD <- generate_induced_subgraphs(
  graph = graphAD,
  nodes_of_interest = AD_hubs_notcontrol,
  exclusive_nodes = AD_hubs_notcontrol)

#Translate
subgraphs_AD <- lapply(subgraphs_AD, translate_vertex_names)

#Save subgraphs
# for (node in names(subgraphs_AD)) {
#   write_graph(subgraphs_AD[[node]], file = paste0("subgraph_", node, ".graphml"), format = "graphml")
# }

#Find new coexpressed genes --- ---

#Create a list to store results
new_coexpressed_genes.l <- list()

#Iterate for each subgraph
for (main_node in names(subgraphs_AD)) {
  #Obtain subraph Obtener el subgrafo
  subgraph <- subgraphs_AD[[main_node]]
  
  #Find new coexpressed genes
  new_genes <- find_new_coexpressed_genes(
    subgraph = subgraph,
    control_graph = graphnoAD,
    main_node = main_node
  )
  
  #Concatenate genes
  new_genes_combined <- paste(new_genes, collapse = ", ")
  
  #df
  new_coexpressed_genes.l[[main_node]] <- data.frame(
    Main_Node = main_node,
    New_Genes = new_genes_combined
  )
}

#Combine results in a df
new_coexpressed_genes.df <- do.call(rbind, new_coexpressed_genes.l)

#Gene frequence in the table
new_coexpressed_genes.fr <- new_coexpressed_genes.df %>%
  dplyr::mutate(New_Genes = strsplit(New_Genes, ", ")) %>%
  tidyr::unnest(New_Genes) %>%
  dplyr::count(New_Genes, sort = TRUE)

#Bar chart for the 10 most frequent genes
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
  gene = new_coexpressed_genes.fr$New_Genes_ense,
  OrgDb = "org.Hs.eg.db", 
  keyType = 'ENSEMBL',
  readable = TRUE,
  universe = universe,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

new_coexpressed_genes.cnet <- cnetplot(new_coexpressed_genes.enr,
                                       showCategory = 20)
# 
# #I want to know something about CROCC
# 
# CROCC <- "ENSG00000226321"
# 
# # Obtener los vecinos de CROCC en ambas redes
# CROCC_neigh_AD <- neighbors(graphAD, CROCC)
# CROCC_neigh_noAD <- neighbors(graphnoAD, CROCC)
# 
# # Convertir los índices de los vecinos a nombres
# CROCC_neigh_AD_names <- V(graphAD)$name[CROCC_neigh_AD]
# CROCC_neigh_noAD_names <- V(graphnoAD)$name[CROCC_neigh_noAD]
# 
# # Identificar los vecinos diferentes entre las dos redes
# CROCC_diff_neigh_AD <- setdiff(unique(CROCC_neigh_AD_names), 
#                          unique(CROCC_neigh_noAD_names)) # Vecinos de CROCC en la red AD pero no en la red noAD
# #enrich
# CROCC.er <- enrichGO(
#   gene = CROCC_diff_neigh_AD,
#   OrgDb = "org.Hs.eg.db", 
#   keyType = 'ENSEMBL',
#   readable = TRUE,
#   universe = universe,
#   ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
#   pvalueCutoff = 0.05, 
#   qvalueCutoff = 0.10)
# 
# CROCC.er.he <- heatplot(CROCC.er, showCategory=30)

#Differences grid

#diff.grid <- grid.arrange(degree_changes.p, new_coexpressed_genes.cnet, ncol = 2)

diff.grid <-cowplot::plot_grid(degree_changes.p,new_coexpressed_genes.cnet,
                               ncol=2, labels=letters[1:2])
diff.grid

#save graph
ggsave("hub_neig_changes_grid.jpg",
       plot = diff.grid, 
       device = "jpg",
       width = 20, #55
       height = 10, #30
       units = "in", 
       dpi = 300
)

################################################

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

#Extract betweenness values from vector
AD_highbe_betweenness <- data.frame(
  ensembl_gene_id = AD_highbe_notcontrol,
  betweenness = nodes_betwe[[1]][AD_highbe_notcontrol]
)

AD_highbe_betweenness.x <- convert_ens_to_symbol(AD_highbe_betweenness)

AD_highbe_betweenness <- AD_highbe_betweenness %>%
  left_join(AD_highbe_betweenness.x, by = "ensembl_gene_id") %>% 
  arrange(desc(betweenness))

AD_highbe_betweenness$external_gene_name[1:10]

#Venn Diagram

highbevenn <- ggVennDiagram(list(AD = highbe_nodes[[1]], C =  highbe_nodes[[2]])) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = "b) Shared high betweeness genes") +
  theme( plot.title = element_text(size=15, hjust = 0.5), legend.position = "none")

#Enrichment of only AD genes

AD_highbe_notcontrol_enrichment <- enrichGO(
  gene = AD_highbe_notcontrol,
  OrgDb = "org.Hs.eg.db", 
  keyType = 'ENSEMBL',
  readable = TRUE,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

AD_highbe_notcontrol_enrichment.df <- as.data.frame(AD_highbe_notcontrol_enrichment)

# vroom::vroom_write(AD_highbe_notcontrol_enrichment.df,
#                    "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_highbe_notcontrol_enrichment.csv")

#Plot enrichment of hubs

AD_highbe_notcontrol_enrichment.p <- cnetplot(AD_highbe_notcontrol_enrichment, circular = TRUE, colorEdge = TRUE) +
  guides(edge_color = "none") +
  ggtitle("e)") 

# ggsave("AD_highbe_notcontrol_enrichment.jpg",
#        plot = AD_highbe_notcontrol_enrichment.p, 
#        device = "jpg",
#        width = 15, #55
#        height = 10, #30
#        units = "in", 
#        dpi = 300
# )

#Induced subgraph of hub genes only in AD

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
# 
# write_graph(highbeb_subgraph, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_highbe_induced_subgraph.graphml',
#             format = "graphml")

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

 # ggsave("grid_hubs_and_highbe1.jpg",
 #        plot = grid, 
 #        device = "jpg",
 #        width = 18, #55
 #        height = 12, #30
 #        units = "in", 
 #        dpi = 300)
 # 
 
 #NOW ONLY GENES IN CONTROL --- ---
 
 #What genes are in the CONTROL network but not in the AD network?
 
 control_highbe_notAD <- setdiff(highbe_nodes[[2]], highbe_nodes[[1]])
 length(control_highbe_notAD)

 control_highbe_notAD.enr <- enrichGO(
   gene = control_highbe_notAD,
   OrgDb = "org.Hs.eg.db", 
   keyType = 'ENSEMBL',
   readable = TRUE,
   ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
   pvalueCutoff = 0.05, 
   qvalueCutoff = 0.05)
 
control_highbe_notAD.cnet <- cnetplot(control_highbe_notAD.enr)

#Extract betweenness values from vector
control_highbe_betweenness <- data.frame(
  ensembl_gene_id = control_highbe_notAD,
  betweenness = nodes_betwe[["graphnoAD"]][control_highbe_notAD]
)

control_highbe_notAD.x <- convert_ens_to_symbol(control_highbe_betweenness$ensembl_gene_id)

control_highbe_betweenness <- control_highbe_betweenness %>%
  left_join(control_highbe_notAD.x, by = "ensembl_gene_id") %>% 
  arrange(desc(betweenness))

control_highbe_betweenness$external_gene_name[1:10]

#I want to know what is the betweenness of AD hub genes nodes in the control network --- ---

#Calculate betweenness for control
betweenness_control <- betweenness(graphnoAD, v = V(graphnoAD), directed = FALSE)
betweenness_control_df <- data.frame(
  ensembl_gene_id = names(betweenness_control),
  betweenness_in_control = betweenness_control
)

#Calculate betweenness for AD
betweenness_AD <- betweenness(graphAD, v = V(graphAD), directed = FALSE)
betweenness_AD_df <- data.frame(
  ensembl_gene_id = names(betweenness_AD),
  betweenness_in_AD = betweenness_AD
)

#Join dataframes by 'ensembl_gene_id', allowing for missing nodes in some network
merged_bet <- merge(
  betweenness_control_df, 
  betweenness_AD_df, 
  by = "ensembl_gene_id", 
  all = TRUE
)

merged_bet.x <- convert_ens_to_symbol(merged_bet$ensembl_gene_id)

merged_bet<- merged_bet %>% left_join(merged_bet.x, by = "ensembl_gene_id")

#Replace NA values with 0 (assuming missing nodes have betweenness 0).
merged_bet[is.na(merged_bet)] <- 0
 
#Calculate percentiles of betweenness for each network
merged_bet$betweenness_percentile_control <- rank(merged_bet$betweenness_in_control) / nrow(merged_bet)
merged_bet$betweenness_percentile_AD <- rank(merged_bet$betweenness_in_AD) / nrow(merged_bet)

#Calculate the difference between the values of betweenness
merged_bet$betweenness_diff <-  merged_bet$betweenness_in_control - merged_bet$betweenness_in_AD

#only in control

only_bet_control <- merged_bet %>% 
  filter(ensembl_gene_id %in% control_highbe_betweenness$ensembl_gene_id)

only_bet_control <- only_bet_control %>%
  arrange(desc(betweenness_diff))

only_bet_control$external_gene_name[34]<- "ENSG00000249877"
only_bet_control$external_gene_name[28]<- "ENSG00000281383"

only_bet_control$betweenness_diff <- round(only_bet_control$betweenness_diff,
                                           digits = 2)

#Long format for gpglot
only_bet_control.l <- only_bet_control %>%
  tidyr::pivot_longer(
    cols = c(betweenness_in_control, betweenness_in_AD), 
    names_to = "network", 
    values_to = "betweenness")

only_bet_control.l <- only_bet_control.l %>%
  arrange(desc(betweenness_diff))

only_bet_control.l$external_gene_name <- factor(
  only_bet_control.l$external_gene_name, 
  levels = only_bet_control$external_gene_name
)

#Change the names of the levels for better visualisation
only_bet_control.l$network <- factor(
  only_bet_control.l$network, 
  levels = c("betweenness_in_control", "betweenness_in_AD"),
  labels = c("Control", "AD")
)

#Plot
diff_be_only_control_.p <- ggplot(only_bet_control.l, aes(x = external_gene_name, y = betweenness, 
                             color = network)) +
  geom_point(size = 3) + 
  geom_segment(
    aes(x = external_gene_name, xend = external_gene_name, 
        y = 0, yend = betweenness), 
    size = 1.5) + 
  geom_text_repel(data = subset(only_bet_control.l, network == "Control"), 
                  aes(label = betweenness_diff), 
                  size = 3, color = "black", max.overlaps = 10) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1), 
    legend.position = "top"
  ) +
  labs(
    title = "Betweenness of nodes control-network-specific",
    x = NULL,
    y = "Betweenness",
    color = "Network"
  ) +
  scale_color_manual(values = c("Control" = "blue", "AD" = "red"))
diff_be_only_control_.p

#only in AD

only_bet_AD <- merged_bet %>% 
  filter(ensembl_gene_id %in% AD_highbe_betweenness$ensembl_gene_id) 
only_bet_AD$betweenness_diff <-abs(only_bet_AD$betweenness_diff)
only_bet_AD$betweenness_diff <- round(only_bet_AD$betweenness_diff, digits = 2)

only_bet_AD <- only_bet_AD %>% arrange(desc(betweenness_diff))

only_bet_AD$external_gene_name[37]<- "AURKAP"
only_bet_AD$external_gene_name[1]<- "LOC124903002"
only_bet_AD$external_gene_name[24]<- "ENSG00000288049"

#Long format for ggplot 
only_bet_AD.l <- only_bet_AD %>%
  tidyr::pivot_longer(
    cols = c(betweenness_in_control, betweenness_in_AD), 
    names_to = "network", 
    values_to = "betweenness"
  )

#Sort external_gene_name based on betweenness_diff
only_bet_AD.l$external_gene_name <- factor(
  only_bet_AD.l$external_gene_name,
  levels = only_bet_AD$external_gene_name[order(only_bet_AD$betweenness_diff, decreasing = TRUE)]
 )

# Cambiar los nombres de los niveles para mejor visualización
only_bet_AD.l$network <- factor(
  only_bet_AD.l$network, 
  levels = c("betweenness_in_control", "betweenness_in_AD"),
  labels = c("Control", "AD")
)

#Plot
diff__be_only_AD_.p <- ggplot() +
geom_segment(data = subset(only_bet_AD.l, network == "AD"),
             aes(x = external_gene_name, xend = external_gene_name,
                 y = 0, yend = betweenness, color = network), 
             size = 1.5) +
  geom_segment(data = subset(only_bet_AD.l, network == "Control"),
               aes(x = external_gene_name, xend = external_gene_name,
                   y = 0, yend = betweenness, color = network), 
               size = 1.5) +
  geom_point(data = subset(only_bet_AD.l, network == "AD"),
             aes(x = external_gene_name, y = betweenness, color = network), 
             size = 3) +
  geom_point(data = subset(only_bet_AD.l, network == "Control"),
             aes(x = external_gene_name, y = betweenness, color = network), 
             size = 3) +
  geom_text_repel(data = subset(only_bet_AD.l, network == "AD"),
                  aes(x = external_gene_name, y = betweenness, label = betweenness_diff), 
                  size = 3, color = "black", max.overlaps = 10) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1),  
    legend.position = "top") +
  labs(title = "Betweenness of nodes AD-network-specific",
    x = NULL,
    y = "Betweenness",
    color = "Network") +
  scale_color_manual(values = c("AD" = "red", 
                                "Control" = "blue"))

diff__be_only_AD_.p

#Grid

diff.grid_be <- cowplot::plot_grid(diff__be_only_AD_.p,
                                   diff_be_only_control_.p,
                                   ncol = 2, labels = letters[1:2])
diff.grid_be

#save plot
# ggsave("betweenness_changes.jpg",
#        plot = diff.grid_be,
#        device = "jpg",
#        width = 15, #55
#        height = 10, #30
#        units = "in",
#        dpi = 300)

#Save hub genes to explore them in the partitions
# 
# write(AD_hubs_notcontrol,
#       file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_hubs_notcontrolens.txt')
# 
# #Save high betweeness genes to explore them in the partitions
# 
# write(AD_highbe_notcontrol,
#       file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/AD_highbe_notcontrol_ens.txt')
# 
# #END