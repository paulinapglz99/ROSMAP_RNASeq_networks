#
#5.1 modularity analysis
#This script takes a list of genes and performs network communities analysis 
#with igraph and INFOMAP

#Libraries --- ---

pacman::p_load('igraph',
               'ggplot2', 
               'dplyr')

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- ---

graph <- read_graph(file = '~/redesROSMAP/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml',
                    format = 'graphml')

#Identify hub genes --- ---

nodes_degree <- degree(graph)

#Table of degree distribution

degree_distribution <- data.frame(gene = names(nodes_degree), degree = nodes_degree)

#Percentile 90 of genes with higher degree

q_threshold.de <- quantile(nodes_degree, probs = 0.90)

#Plot degree distribution

ggplot(degree_distribution, aes(x = degree, y = ..density..)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  labs(title = "Degree distribution histogram",
       subtitle = "for patients with no AD",
       x = "Grado", y = "Freq") +
  geom_vline(xintercept = q_threshold.de, color = "red", linetype = "dashed") +
  geom_text(aes(x = q_threshold.de, y = 0.06, label = "90th percentile"), color = "red", hjust = -0.1) +
   theme_minimal()

#Hub genes are 

hub_genes <- V(graph)$name[nodes_degree > q_threshold.de]

hub_genes <- nodes_degree[hub_genes]

hub_gene_distribution <- data.frame(gene = names(hub_genes), degree = hub_genes)

ggplot(hub_gene_distribution, aes(x = degree, y = ..density..)) +
  geom_histogram(binwidth = 1, fill = "skyblue", color = "black") +
  scale_x_continuous(breaks = seq(0, max(hub_gene_distribution$degree), by = 2)) +
  labs(title = "Degree distribution histogram", x = "Grado", y = "Freq") +
  theme_minimal()

#Plot the hub genes

hub_genes.df <- data.frame(genes = names(nodes_degree), degree = nodes_degree)
hub_genes.df$hub_gene <- ifelse(hub_genes.df$genes %in% names(hub_genes), "Hub Gene", "Not Hub Gene")

#Plot hub genes

ggplot(hub_genes.df,
                      aes(x = genes, y = degree, color = hub_gene)) +
  geom_text(data = subset(hub_genes.df, hub_gene == "Hub Gene"), aes(label = genes), 
            color = 'black', hjust = -0.1, vjust = 0.2, size = 3) +  # Añadir etiquetas solo para los "hub genes"
  geom_point(size = 3) +
  theme(axis.text.x = element_blank()) +  # Rotar las etiquetas del eje x
  labs(x = "Genes", y = "Degree", title = "Degree of Genes", subtitle = "with hub genes indicated in patients with pathological AD") 

#Identify high betweenness nodes --- ---

betweenness_values <- betweenness(graph)

q_threshold.be <- quantile(betweenness_values, probs = 0.90)

#High betweenness nodes are 

high_betweenness_nodes <- V(graph)$name[betweenness_values > q_threshold.be]

high_betweenness_nodes <- betweenness_values[high_betweenness_nodes]

#Plot the hub genes

high_betweenness_nodes.df <- data.frame(genes = names(betweenness_values), betweenness_values = betweenness_values)
high_betweenness_nodes.df$high_betweenness <- ifelse(high_betweenness_nodes.df$genes %in% names(high_betweenness_nodes), "High betweeness Gene", "Not High betweeness Gene")

# Crear el gráfico de puntos
high_betweenness_nodes.p <- ggplot(high_betweenness_nodes.df,
                      aes(x = genes, y = betweenness_values, color = high_betweenness)) +
  geom_text(data = subset(high_betweenness_nodes.df, high_betweenness == "High betweeness Gene"), aes(label = genes), 
            color = 'black', hjust = -0.1, vjust = 0.2, size = 3) +  # Añadir etiquetas solo para los "hub genes"
  geom_point(size = 3) +
  theme(axis.text.x = element_blank()) +  # Rotar las etiquetas del eje x
  labs(x = "Genes", y = "Betweeness", title = "Betweeness of Genes", subtitle = "with central brokerage genes indicated in patients with pathological AD") 

high_betweenness_nodes.p

#END
