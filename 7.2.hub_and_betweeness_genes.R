#
#5.1 modularity analysis
#This script takes a list of genes and performs network communities analysis 
#with igraph and INFOMAP

#Libraries --- ---

pacman::p_load('igraph',
               'ggplot2', 
               'dplyr', 
               'gridExtra', 
               'biomaRt')

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- ---

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml',
                        format = 'graphml')

graphs <- list(graphAD = graphAD, 
               graphnoAD = graphnoAD)

#Define function to convert gene names --- ---

# Connecting to the Ensembl database through biomaRt
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Obtener la conversión de ENSMBL a SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = ensembl)
}

##Change names of vertex from both graphs --- ---

# Extract vertex names

graphAD_vnames <- V(graphAD)$name
length(graphAD_vnames)

#Translate names

graphAD_vnames_trad <-convert_ens_to_symbol(graphAD_vnames)

# Reemplazar los valores faltantes en la columna 'external_gene_name' con los valores de 'ensembl_gene_id'
graphAD_vnames_trad <- graphAD_vnames_trad %>%
  mutate(external_gene_name = ifelse(external_gene_name == "", ensembl_gene_id, external_gene_name))


# Create a vector of translated names using the dictionary
# We need to ensure that the actual names of the network are in the dictionary

graphAD_vnames_trad <- setNames(graphAD_vnames_trad$external_gene_name, graphAD_vnames_trad$ensembl_gene_id)

# Ordena graphAD_vnames_trad según el orden de graphAD_vnames
sorted_graphAD_vnames_trad <- graphAD_vnames_trad[match(graphAD_vnames, names(graphAD_vnames_trad))]

# Assign the new names to the network vertices.

V(graphAD)$name <- sorted_graphAD_vnames_trad

###TRYING HARD ######

# Extraer los nombres de los vértices
graphAD_vnames <- V(graphAD)$name

# Traduce los nombres
graphAD_vnames_trad <- convert_ens_to_symbol(graphAD_vnames)

# Crea un vector de nombres traducidos utilizando el diccionario
graphAD_vnames_trad <- setNames(graphAD_vnames_trad$external_gene_name, graphAD_vnames_trad$ensembl_gene_id)

# Ordena graphAD_vnames_trad según el orden de graphAD_vnames y conserva los nombres originales cuando no se encuentra una traducción
sorted_graphAD_vnames_trad <- sapply(graphAD_vnames, function(x) {
  if (x %in% names(graphAD_vnames_trad)) {
    return(graphAD_vnames_trad[[x]])
  } else {
    return(x)
  }
})

# Asigna los nuevos nombres a los vértices del grafo
V(graphAD)$name <- sorted_graphAD_vnames_trad

# Verifica los nombres de los vértices
print(V(graphAD)$name)

###TRYING HARD ######








# Extraer los nombres actuales de los vértices
graphAD_vnames <- V(graphAD)$name

# Crear un vector de nombres traducidos utilizando el diccionario
# Necesitamos asegurar que los nombres actuales del grafo estén en el diccionario
name_map <- setNames(diccionario$external_gene_name, diccionario$ensembl_gene_id)

# Ahora mapeamos los nombres actuales a los nuevos nombres
new_names <- sapply(graphAD_vnames, function(x) {
  if (x %in% names(name_map)) {
    return(name_map[[x]])
  } else {
    return(x)
  }
})

# Asignar los nuevos nombres a los vértices del grafo
V(graphAD)$name <- new_names


#####################

#Identify hub genes --- ---

#Calculate degree of nodes
nodes_degree <- sapply(X = graphs, FUN = degree)

#Table of degree distribution

degree_distributionAD <- data.frame(gene = names(nodes_degree[[1]]), degree = nodes_degree[[1]])

degree_distributionnoAD <- data.frame(gene = names(nodes_degree[[2]]), degree = nodes_degree[[2]])

degree_distribution <- list()

for (i in 1:length(nodes_degree)) {
  # Aplicar cluster_infomap a cada grafo y almacenar el resultado en results_list
  degree_distribution[[i]] <- data.frame(gene = names(nodes_degree[[i]]), degree = nodes_degree[[i]])
}

#Calculate percentile 95 of genes with higher degree ---- ---

q_threshold.de <- list()

for (i in 1:length(nodes_degree)) {
  # Aplicar cluster_infomap a cada grafo y almacenar el resultado en results_list
  q_threshold.de[[i]] <- quantile(nodes_degree[[i]], probs = 0.95)
}

#Plot degree distribution

AD_distribution <- ggplot(degree_distribution[[1]], aes(x = degree, y = ..density..)) +
  geom_histogram(binwidth = 1, fill = "#00688B", color = "black") +
  labs(title = "Degree distribution histogram",
       subtitle = "for patients with AD",
       x = "Grado", y = "Freq") +
  geom_vline(xintercept = q_threshold.de[[1]], color = "red", linetype = "dashed") +
  geom_text(aes(x = q_threshold.de[[1]], y = 0.06, label = "95th percentile"), color = "red", hjust = -0.1) +
  theme_minimal()

noAD_distribution <- ggplot(degree_distribution[[2]], aes(x = degree, y = ..density..)) +
  geom_histogram(binwidth = 1, fill = "#00688B", color = "black") +
  labs(title = "Degree distribution histogram",
       subtitle = "for patients with no AD",
       x = "Grado", y = "Freq") +
  geom_vline(xintercept = q_threshold.de[[2]], color = "red", linetype = "dashed") +
  geom_text(aes(x = q_threshold.de[[2]], y = 0.06, label = "95th percentile"), color = "red", hjust = -0.1) +
  theme_minimal()

grid.arrange(AD_distribution, noAD_distribution)

#Hub genes --- ---

#Hub genes for AD patients

hub_genes_AD <- V(graphAD)$name[nodes_degree[[1]] > q_threshold.de[[1]]]

hub_genes_AD <- nodes_degree[[1]][hub_genes_AD]

hub_gene_AD <- data.frame(ensembl_gene_id = names(hub_genes_AD), degree = hub_genes_AD)

#Order by degree

hub_gene_AD <- hub_gene_AD[order(hub_gene_AD$degree, decreasing = TRUE), ]

#Plot hub genes for AD

#Df
hub_genes_AD.df <- data.frame(ensembl_gene_id = names(nodes_degree[[1]]), degree = nodes_degree[[1]])
hub_genes_AD.df$hub_gene <- ifelse(hub_genes_AD.df$ensembl_gene_id %in% hub_gene_AD$ensembl_gene_id, "Hub Gene", "Not Hub Gene")

#Translate again
hub_genes_AD.df_trads <- convert_ens_to_symbol(hub_genes_AD.df$ensembl_gene_id)

hub_genes_AD.df <- hub_genes_AD.df %>% left_join(hub_genes_AD.df_trads, by= "ensembl_gene_id")

#Order factors for plotting
chrom_order <- c(1:22, "X", "MT")

hub_genes_AD.df$chromosome_name <- factor(hub_genes_AD.df$chromosome_name, levels = chrom_order)

#p
hub_genes_AD.p <- ggplot(hub_genes_AD.df,
                         aes(x = chromosome_name, y = degree, color = hub_gene)) +
  geom_text(data = subset(hub_genes_AD.df, hub_gene == "Hub Gene"), aes(label = external_gene_name), 
            color = 'black', hjust = 0.5, vjust = -0.8, size = 3, check_overlap = TRUE) +  # Añadir etiquetas solo para los "hub genes"
  geom_point(size = 3) +
  theme(axis.text.x = element_blank(), legend.position = "none") +  
  labs(x = "Genes", y = "Degree", title = "Degree of Genes", subtitle = "with hub genes (95%) indicated in patients with pathological AD") +
  scale_color_manual(values = c("Hub Gene" = "#8B3A3A", "Not Hub Gene" = "#1874CD")) +
  theme_light()

#For no AD patients --- --- 

#Hub genes --- ---

#Hub genes for AD patients

hub_genes_noAD <- V(graphnoAD)$name[nodes_degree[[2]] > q_threshold.de[[2]]]

hub_genes_noAD <- nodes_degree[[2]][hub_genes_noAD]

hub_genes_noAD <- data.frame(ensembl_gene_id = names(hub_genes_noAD), degree = hub_genes_noAD)

#Plot hub genes for AD

#Df
hub_genes_noAD.df <- data.frame(ensembl_gene_id = names(nodes_degree[[2]]), degree = nodes_degree[[2]])
hub_genes_noAD.df$hub_gene <- ifelse(hub_genes_noAD.df$ensembl_gene_id %in% hub_genes_noAD$ensembl_gene_id, "Hub Gene", "Not Hub Gene")

#Translate again
hub_genes_AD.df_trads <- convert_ens_to_symbol(hub_genes_noAD.df$ensembl_gene_id)

hub_genes_noAD.df <- hub_genes_noAD.df %>% left_join(hub_genes_AD.df_trads, by= "ensembl_gene_id")

#p

#Order factors for plotting

hub_genes_noAD.df$chromosome_name <- factor(hub_genes_noAD.df$chromosome_name, levels = chrom_order)

hub_genes_noAD.p <- ggplot(hub_genes_noAD.df,
                         aes(x = chromosome_name, y = degree, color = hub_gene)) +
  geom_text(data = subset(hub_genes_noAD.df, hub_gene == "Hub Gene"), aes(label = external_gene_name), 
            color = 'black', hjust = 0.5, vjust = -0.8, size = 3, check_overlap = TRUE) +  # Añadir etiquetas solo para los "hub genes"
  geom_point(size = 3) +
  theme(axis.text.x = element_blank(), legend.position = "none") +  
  labs(x = "Chromosoe", y = "Degree", title = "Degree of Genes", subtitle = "with hub genes (95%) indicated in patients without pathological AD") +
  scale_color_manual(values = c("Hub Gene" = "#8B3A3A", "Not Hub Gene" = "#1874CD")) +
  theme_light()

#Grid both plots --- ---

grid_hub_genes <- grid.arrange(hub_genes_AD.p, hub_genes_noAD.p)

#Save plot

ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNA-seq_DLFPC_NIA_reagan_hub_genes.png", 
grid_hub_genes, width = 20, height = 10, dpi = 300)

#Analysis of hub genes --- --- 

#What genes do both networks share?

shared_hub_genes <- intersect(hub_genes_AD.df$external_gene_name, hub_genes_noAD.df$external_gene_name)

# Genes in AD but not in noAD

hub_gene_AD <- hub_gene_AD %>% left_join(hub_genes_AD.df, by = "ensembl_gene_id")
hub_genes_noAD <- hub_genes_noAD <- hub_genes_noAD %>% left_join(hub_genes_noAD.df, by = "ensembl_gene_id")

#
genes_en_AD_no_noAD_sym <- setdiff(hub_gene_AD$external_gene_name, hub_genes_noAD$external_gene_name)
genes_en_AD_no_noAD_ens<- setdiff(hub_gene_AD$ensembl_gene_id, hub_genes_noAD$ensembl_gene_id)

#Save hub genes to explore them in the partitions

vroom::vroom_write(genes_en_AD_no_noAD_ens)

#Comparison of hub genes not found in healthy people

genes_en_AD_no_noAD.df <- hub_gene_AD %>% filter(external_gene_name %in% genes_en_AD_no_noAD_sym)
dim(genes_en_AD_no_noAD.df)
#[1] 25  6

#Network of hub genes in the AD graph but not in the noAD graph

indexed_ver <- which(V(graphAD)$name %in% genes_en_AD_no_noAD_ens)

genes_en_AD_no_noAD.g <- induced_subgraph(graphAD, indexed_ver)

plot(genes_en_AD_no_noAD.g)

#Give alternative names to the nodes --- --- 

# Convertir el grafo a un data frame
df_grafo <- as.data.frame(get.data.frame(genes_en_AD_no_noAD.g, what = "vertices"))

df_grafo <- df_grafo %>% left_join(genes_en_AD_no_noAD.df, by = c("name" = "ensembl_gene_id"))
df_grafo<- df_grafo[-4]
df_grafo <- df_grafo %>% rename(degree.x = "degree")
genes_en_AD_no_noAD.g$alt_name <- df_grafo[,c(1, 5)]
df_grafo <- df_grafo[order(df_grafo$degree, decreasing = TRUE), ]

ay <- set_vertex_attr(genes_en_AD_no_noAD.g, "alt_name", value=df_grafo$external_gene_name)

#
ggraph(genes_en_AD_no_noAD.g, layout = "auto") +
  geom_node_text(aes(label = genes_en_AD_no_noAD.g$alt_name))

#Save graph

write_graph(genes_en_AD_no_noAD.g, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_hub_genes_AD_noAD.graphml',
            format = "graphml")

################################################################################

#Identify high betweenness nodes --- ---

betweenness_values <- sapply(X = graphs, FUN = betweenness)

q_threshold.be <- list()

for (i in 1:length(betweenness_values)) {
  # Aplicar cluster_infomap a cada grafo y almacenar el resultado en results_list
  q_threshold.be[[i]] <- quantile(betweenness_values[[i]], probs = 0.95)
}

#High betweeness nodes for AD patients --- ---

high_betweenness_nodes_AD <- V(graphAD)$name[betweenness_values[[1]] > q_threshold.be[[1]]]

high_betweenness_nodes_AD <- betweenness_values[[1]][high_betweenness_nodes_AD]

high_betweenness_nodes_AD <- data.frame(ensembl_gene_id = names(high_betweenness_nodes_AD), betweenness_values = high_betweenness_nodes_AD)

# Convert ENSMBL a SYMBOL

traductions_be_AD <- convert_ens_to_symbol(high_betweenness_nodes_AD)

# Merge gene names traduction

high_betweenness_nodes_AD <- high_betweenness_nodes_AD %>% left_join(traductions_be_AD, by = "ensembl_gene_id")

high_betweenness_nodes_AD <- high_betweenness_nodes_AD[order(high_betweenness_nodes_AD$betweenness_values, decreasing = TRUE), ]

#Plot --- ---

#df

#Create data frame 
high_betweenness_nodes_AD.df <- data.frame(ensembl_gene_id = names(betweenness_values[[1]]), betweeness = betweenness_values[[1]])
#Add check column
high_betweenness_nodes_AD.df$high_betweenness <- ifelse(high_betweenness_nodes_AD.df$ensembl_gene_id %in% high_betweenness_nodes_AD$ensembl_gene_id, "High betweeness Gene", "Not High betweeness Gene")
#Translate again
high_betweenness_nodes_AD.df_trads <- convert_ens_to_symbol(hub_genes_noAD.df)
#Paste both dfs
high_betweenness_nodes_AD.df <- high_betweenness_nodes_AD.df %>% left_join(high_betweenness_nodes_AD.df_trads, by = "ensembl_gene_id")
#Arrange chromosome names in order
high_betweenness_nodes_AD.df$chromosome_name <- factor(high_betweenness_nodes_AD.df$chromosome_name, levels = chrom_order)

#p

high_betweenness_nodes_AD.p <- ggplot(high_betweenness_nodes_AD.df,
                                      aes(x = chromosome_name, y = betweeness, color = high_betweenness)) +
  geom_text(data = subset(high_betweenness_nodes_AD.df, high_betweenness == "High betweeness Gene"), aes(label = external_gene_name), 
            color = 'black', hjust = 1.1, vjust = 0.2, size = 3) +  # Añadir etiquetas solo para los "hub genes"
  geom_point(size = 3) +
  scale_fill_manual(values = c("High betweeness Gene" = "red", "Other" = "blue")) +  # Definir los colores para las categorías
  theme(axis.text.x = element_blank(),legend.position = "none") +  # Rotar las etiquetas del eje x
  labs(x = "Chromosome", y = "Betweeness", title = "Betweeness of Genes", subtitle = "with central genes indicated in patients with pathological AD") +
  scale_color_manual(values = c("High betweeness Gene" = "#8B3A3A", "Not High betweeness Gene" = "#1874CD")) +
  theme_light()

high_betweenness_nodes_AD.p

#High betweeness nodes for noAD patients --- --- 

high_betweenness_nodes_noAD <- V(graphnoAD)$name[betweenness_values[[2]] > q_threshold.be[[2]]]

high_betweenness_nodes_noAD <- betweenness_values[[2]][high_betweenness_nodes_noAD]

high_betweenness_nodes_noAD <- data.frame(ensembl_gene_id = names(high_betweenness_nodes_noAD), betweenness_values = high_betweenness_nodes_noAD)

# Convert ENSMBL a SYMBOL

traductions_be_noAD <- convert_ens_to_symbol(high_betweenness_nodes_noAD$ensembl_gene_id)

# Merge gene names traduction

high_betweenness_nodes_noAD <- high_betweenness_nodes_noAD %>% left_join(traductions_be_noAD, by = "ensembl_gene_id")

high_betweenness_nodes_noAD <- high_betweenness_nodes_noAD[order(high_betweenness_nodes_noAD$betweenness_values, decreasing = TRUE), ]

#Plot --- ---

#df

#Create data frame 
high_betweenness_nodes_noAD.df <- data.frame(ensembl_gene_id = names(betweenness_values[[2]]), betweeness = betweenness_values[[2]])
#Add check column
high_betweenness_nodes_noAD.df$high_betweenness <- ifelse(high_betweenness_nodes_noAD.df$ensembl_gene_id %in% high_betweenness_nodes_noAD$ensembl_gene_id, "High betweeness Gene", "Not High betweeness Gene")
#Translate again
high_betweenness_nodes_noAD.df_trads <- convert_ens_to_symbol(high_betweenness_nodes_noAD.df)
#Paste both dfs
high_betweenness_nodes_noAD.df <- high_betweenness_nodes_noAD.df %>% left_join(high_betweenness_nodes_noAD.df_trads, by = "ensembl_gene_id")
#Arrange chromosome names in order
high_betweenness_nodes_noAD.df$chromosome_name <- factor(high_betweenness_nodes_noAD.df$chromosome_name, levels = chrom_order)

#p

high_betweenness_nodes_noAD.p <- ggplot(high_betweenness_nodes_noAD.df,
                                      aes(x = chromosome_name, y = betweeness, color = high_betweenness)) +
  geom_text(data = subset(high_betweenness_nodes_noAD.df, high_betweenness == "High betweeness Gene"), aes(label = external_gene_name), 
            color = 'black', hjust = 1.1, vjust = 0.2, size = 3) +  # Añadir etiquetas solo para los "hub genes"
  geom_point(size = 3) +
  scale_fill_manual(values = c("High betweeness Gene" = "red", "Other" = "blue")) +  # Definir los colores para las categorías
  theme(axis.text.x = element_blank(),legend.position = "none") +  # Rotar las etiquetas del eje x
  labs(x = "Chromosome", y = "Betweeness", title = "Betweeness of Genes", subtitle = "with central genes indicated in patients without pathological AD") +
  scale_color_manual(values = c("High betweeness Gene" = "#8B3A3A", "Not High betweeness Gene" = "#1874CD")) +
  theme_light()

high_betweenness_nodes_noAD.p

#Grid both plots

grid_high_be_genes <- grid.arrange(high_betweenness_nodes_AD.p, high_betweenness_nodes_noAD.p)

#Save plot
ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNA-seq_DLFPC_NIA_reagan_high_be_genes.png", 
       grid_high_be_genes, width = 20, height = 10, dpi = 300)

#What genes do both networks share?

shared_high_be_genes <- intersect(high_betweenness_nodes_AD$external_gene_name, high_betweenness_nodes_noAD$external_gene_name)
#

high_be_genes_en_AD_no_noAD_sym <- setdiff(high_betweenness_nodes_AD$external_gene_name, high_betweenness_nodes_noAD$external_gene_name)
high_be_genes_en_AD_no_noAD_ens<- setdiff(high_betweenness_nodes_AD$ensembl_gene_id, high_betweenness_nodes_noAD$ensembl_gene_id)

#Comparison of hub genes not found in healthy people

high_be_genes_en_AD_no_noAD.df <- high_betweenness_nodes_AD %>% filter(external_gene_name %in% high_be_genes_en_AD_no_noAD_sym)
dim(high_be_genes_en_AD_no_noAD.df)
#[1] 66  4

#Network of hub genes in the AD graph but not in the noAD graph --- ---

indexed_ver_be <- which(V(graphAD)$name %in% high_be_genes_en_AD_no_noAD_ens)

high_be_genes_en_AD_no_noAD.g <- induced_subgraph(graphAD, indexed_ver_be)

plot(high_be_genes_en_AD_no_noAD.g)

#

name <- setNames(high_be_genes_en_AD_no_noAD.df$ensembl_gene_id, high_be_genes_en_AD_no_noAD.df$external_gene_name )
#name <- paste(names(name), "=", unname(name), collapse = ", ")
name <- paste0('"', names(name), '" = "', unname(name), '"', collapse = ", ")
name <- c(name)

V(high_be_genes_en_AD_no_noAD.g)$name <- name[V(high_be_genes_en_AD_no_noAD.g)$name]

plot(high_be_genes_en_AD_no_noAD.g)

sym_names <- high_be_genes_en_AD_no_noAD.df$ensembl_gene_id
names(sym_names) <- high_be_genes_en_AD_no_noAD.df$external_gene_name

sym_names <- c("114" = "Lisa", "115" = "Joe", "113" = "Peter", "112" = "Barbara", "113" = "Bart")
V(high_be_genes_en_AD_no_noAD.g)$sym_names <- sym_names[V(high_be_genes_en_AD_no_noAD.g)$sym_names]

texto <- paste(names(sym_names), "=", unname(sym_names), collapse = ", ")


names <- c(sym_names)
V(high_be_genes_en_AD_no_noAD.g)$label <- names[V(high_be_genes_en_AD_no_noAD.g)$name]
plot(high_be_genes_en_AD_no_noAD.g)

V(high_be_genes_en_AD_no_noAD.g)$additional_name <- high_be_genes_en_AD_no_noAD.df$external_gene_name

#Save graph

write_graph(genes_en_AD_no_noAD.g, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_hub_genes_AD_noAD.graphml',
            format = "graphml")

#NEXT QUESTION IS What modules do these genes belong to? --- ---



#END