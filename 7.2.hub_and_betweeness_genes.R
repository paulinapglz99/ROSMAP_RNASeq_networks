#
#7.2.hub_and_betweeness_genes.R
#This script makes an analysis of hub genes and high betweeness genes

#Libraries --- ---

pacman::p_load('igraph',
               'ggplot2', 
               'dplyr', 
               'gridExtra', 
               'biomaRt', 
               "ggVennDiagram", 
               "ClusterProfiler")

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

# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = ensembl)
}

# Define function to change names of vertex for symbol names --- ---

translate_vertex_names <- function(graph) {
  # Extract vertex names
  graph_vnames <- V(graph)$name
  # Translate names
  graph_vnames_trad <- convert_ens_to_symbol(graph_vnames)
  # Reemplazar los valores faltantes en la columna 'external_gene_name' con los valores de 'ensembl_gene_id'
  graph_vnames_trad$external_gene_name <- ifelse(graph_vnames_trad$external_gene_name == "", graph_vnames_trad$ensembl_gene_id, graph_vnames_trad$external_gene_name)
  # Create a vector of translated names using the dictionary
  # We need to ensure that the actual names of the network are in the dictionary
  graph_vnames_trad <- setNames(graph_vnames_trad$external_gene_name, graph_vnames_trad$ensembl_gene_id)
  # Ordena graph_vnames_trad según el orden de graph_vnames
  sorted_graph_vnames_trad <- graph_vnames_trad[match(graph_vnames, names(graph_vnames_trad))]
    # Assign the new names to the network vertices.
  V(graph)$name <- sorted_graph_vnames_trad
  return(graph)
}

#Define similarity of Enriched Processes, Jaccard Index function --- ---

jaccard_simplex <- function(a,b){
  length(intersect(a,b))/length(union(a,b))
}

##Change names of vertex from both graphs --- ---
# 
# # Extract vertex names
# 
# graphAD_vnames <- V(graphAD)$name
# length(graphAD_vnames)
# 
# #Translate names
# 
# graphAD_vnames_trad <-convert_ens_to_symbol(graphAD_vnames)
# 
# # Reemplazar los valores faltantes en la columna 'external_gene_name' con los valores de 'ensembl_gene_id'
# graphAD_vnames_trad <- graphAD_vnames_trad %>%
#   mutate(external_gene_name = ifelse(external_gene_name == "", ensembl_gene_id, external_gene_name))
# 
# # Create a vector of translated names using the dictionary
# # We need to ensure that the actual names of the network are in the dictionary
# 
# graphAD_vnames_trad <- setNames(graphAD_vnames_trad$external_gene_name, graphAD_vnames_trad$ensembl_gene_id)
# 
# # Ordena graphAD_vnames_trad según el orden de graphAD_vnames
# sorted_graphAD_vnames_trad <- graphAD_vnames_trad[match(graphAD_vnames, names(graphAD_vnames_trad))]
# identical(names(sorted_graphAD_vnames_trad), graphAD_vnames)#must say true
# 
# # Assign the new names to the network vertices.
# 
# V(graphAD)$name <- sorted_graphAD_vnames_trad

#Identify hub genes --- ---

#Calculate degree of nodes
nodes_degree <- sapply(X = graphs, FUN = degree)

#Table of degree distribution

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

#For AD --- --- 
#Plot degree distribution
AD_distribution.df <- degree_distribution[[1]]

#Order data by degree column
AD_distribution.df <- AD_distribution.df[order(AD_distribution.df$degree, decreasing = TRUE), ]

#Cummulative sum
AD_distribution.df$CumulativeDegree <- cumsum(AD_distribution.df$degree)


AD_distribution.df$logCumulativeDegree <- log10(AD_distribution.df$CumulativeDegree)

#Add threshold
AD_CumulativeDegree_threshold <- quantile(AD_distribution.df$CumulativeDegree, probs = 0.95)

#
AD_distribution_accum <- ggplot(AD_distribution.df, aes(x = reorder(gene, degree), y = CumulativeDegree, group = 1)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = AD_CumulativeDegree_threshold, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Cumulative Distribution",
       x = "Gene",
       y = "CumulativeDegree") +
  theme(axis.text.x = element_blank()) +
  log10()+
  annotate("text", x = length(AD_distribution.df$gene)/2, y = AD_CumulativeDegree_threshold,
           label = "95% Umbral", vjust = -1, color = "red")

grado_freq <- as.data.frame(table(AD_distribution.df$degree))

#For no AD --- --- 

#Plot degree distribution
noAD_distribution.df <- degree_distribution[[2]]

#Order data by degree column
noAD_distribution.df <- noAD_distribution.df[order(noAD_distribution.df$degree), ]

#Cummulative sum
noAD_distribution.df$CumulativeDegree <- cumsum(noAD_distribution.df$degree)

#
noAD_distribution.df$logCumulativeDegree <- log10(noAD_distribution.df$CumulativeDegree)

#Add threshold
noAD_CumulativeDegree_threshold <- quantile(noAD_distribution.df$logCumulativeDegree, probs = 0.95)

#
noAD_distribution_accum <- ggplot(noAD_distribution.df, aes(x = reorder(gene, degree), y = logCumulativeDegree, group = 1)) +
  geom_line(color = "blue") +
  geom_hline(yintercept = noAD_CumulativeDegree_threshold, linetype = "dashed", color = "red", size = 1) +
  labs(title = "Cumulative Distribution",
       x = "Gene",
       y = "Accumulated Value") +
  theme(axis.text.x = element_blank()) +
  annotate("text", x = length(noAD_distribution.df$gene)/2, y = noAD_CumulativeDegree_threshold,
           label = "95% Umbral", vjust = -1, color = "red")

#Left join both

degree_distribution.df <- AD_distribution.df %>% left_join(noAD_distribution.df, by = "gene")

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

hub_genes_AD <- data.frame(ensembl_gene_id = names(hub_genes_AD), degree = hub_genes_AD)

hub_genes_AD_trad <- convert_ens_to_symbol(hub_genes_AD)

hub_genes_AD <- hub_genes_AD %>% left_join(hub_genes_AD_trad, by= "ensembl_gene_id")

#Order by degree

hub_genes_AD <- hub_genes_AD[order(hub_genes_AD$degree, decreasing = TRUE), ]

#Plot hub genes for AD

#Df
hub_genes_AD.df <- data.frame(ensembl_gene_id = names(nodes_degree[[1]]), degree = nodes_degree[[1]])
hub_genes_AD.df$hub_gene <- ifelse(hub_genes_AD.df$ensembl_gene_id %in% hub_genes_AD$ensembl_gene_id, "Hub Gene", "Not Hub Gene")

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

#Enrichment of the hub genes

hub_genes_AD_enrichment <- enrichGO(
  gene = hub_genes_AD$external_gene_name,
  OrgDb = org.Hs.eg.db, 
  keyType = 'SYMBOL',
  readable = TRUE,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

barplot(hub_genes_AD_enrichment)

#For no AD patients --- --- 

#Hub genes --- ---

#Hub genes for AD patients

hub_genes_noAD <- V(graphnoAD)$name[nodes_degree[[2]] > q_threshold.de[[2]]]

hub_genes_noAD <- nodes_degree[[2]][hub_genes_noAD]

hub_genes_noAD <- data.frame(ensembl_gene_id = names(hub_genes_noAD), degree = hub_genes_noAD)

hub_genes_noAD_trad <- convert_ens_to_symbol(hub_genes_noAD)

hub_genes_noAD <- hub_genes_noAD %>% left_join(hub_genes_noAD_trad, by= "ensembl_gene_id")

hub_genes_noAD <- hub_genes_noAD[order(hub_genes_noAD$degree, decreasing = TRUE), ]

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

#Enrichment

hub_genes_noAD_enrichment <- enrichGO(
  gene = hub_genes_noAD$external_gene_name,
  OrgDb = org.Hs.eg.db, 
  keyType = 'SYMBOL',
  readable = TRUE,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

#Enrichment table

hub_genes_noAD_enrichment.df <- as.data.frame(hub_genes_noAD_enrichment@result)

#Enrichment plot

barplot(hub_genes_noAD_enrichment)



#Grid both plots --- ---

grid_hub_genes <- grid.arrange(hub_genes_AD.p, hub_genes_noAD.p)

#Save plot
# 
# ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNA-seq_DLFPC_NIA_reagan_hub_genes.png", 
# grid_hub_genes, width = 20, height = 10, dpi = 300)

#Analysis of hub genes --- --- 

#What genes do both networks share?

shared_hub_genes <- intersect(hub_genes_AD$external_gene_name, hub_genes_noAD$external_gene_name)

venn.hub <- venn.diagram(
  x = list(AD = hub_genes_AD$external_gene_name, noAD = hub_genes_noAD$external_gene_name),
filename = NULL,
fill = c("red", "blue"),
alpha = 0.5,
cex = 2,
cat.cex = 2,
cat.pos = 0
)

ggVennDiagram(list(AD = hub_genes_AD$external_gene_name, noAD = hub_genes_noAD$external_gene_name)) +
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  theme(legend.position = "none")

# Genes in AD but not in noAD

hub_gene_AD <- hub_gene_AD %>% left_join(hub_genes_AD.df, by = "ensembl_gene_id")
hub_genes_noAD <- hub_genes_noAD <- hub_genes_noAD %>% left_join(hub_genes_noAD.df, by = "ensembl_gene_id")

#
genes_en_AD_no_noAD_sym <- setdiff(hub_genes_AD$external_gene_name, hub_genes_noAD$external_gene_name)

genes_en_AD_no_noAD_ens<- setdiff(hub_genes_AD$ensembl_gene_id, hub_genes_noAD$ensembl_gene_id)

#Enrichment of genes_en_AD_no_noAD_sym

genes_en_AD_no_noAD_sym

genes_en_AD_no_noAD_sym_enrichment <- enrichGO(
  gene = genes_en_AD_no_noAD_sym,
  OrgDb = org.Hs.eg.db, 
  keyType = 'SYMBOL',
  readable = TRUE,
  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.10)

#table of enrichment

genes_en_AD_no_noAD_sym_enrichment.df <- as.data.frame(genes_en_AD_no_noAD_sym_enrichment@result)

#Enrichment plots

barplot(genes_en_AD_no_noAD_sym_enrichment, showCategory=10)

dotplot(genes_en_AD_no_noAD_sym_enrichment, showCategory=10)

cnetplot(genes_en_AD_no_noAD_sym_enrichment, circular = TRUE, colorEdge = TRUE) 

#

pathview(gene.data = genes_en_AD_no_noAD_ens, pathway.id="hsa05014", #AD
         species = "hsa", gene.idtype=gene.idtype.list[3])

#Save hub genes to explore them in the partitions

#write(genes_en_AD_no_noAD_ens, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/genes_en_AD_no_noAD_ens.txt')

#Comparison of hub genes not found in healthy people

genes_en_AD_no_noAD.df <- hub_gene_AD %>% filter(external_gene_name %in% genes_en_AD_no_noAD_sym)
dim(genes_en_AD_no_noAD.df)
#[1] 25  6

#Network of hub genes in the AD graph but not in the noAD graph

indexed_ver <- which(V(graphAD)$name %in% genes_en_AD_no_noAD_ens)

genes_en_AD_no_noAD.g <- induced_subgraph(graphAD, indexed_ver)

plot(genes_en_AD_no_noAD.g)

#Translate vertex names

xgenes_en_AD_no_noAD.g <- translate_vertex_names(genes_en_AD_no_noAD.g)

plot(xgenes_en_AD_no_noAD.g)

#Save graph
# 
# write_graph(xgenes_en_AD_no_noAD.g, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_hub_genes_AD_noAD_trad.graphml',
#             format = "graphml")

#Similitud de procesos biologicos de cada hub gene

hub_genes_AD_enrichment

hub_genes_noAD_enrichment

#Comparison of enrichments

ComparisonEnrichedProcessJ_hubs <- jaccard_simplex(names(hub_genes_AD_enrichment@geneSets), names(hub_genes_noAD_enrichment@geneSets))
#[1] 0.8768494


################################################################################

#Identify high betweenness nodes --- ---

betweenness_values <- sapply(X = graphs, FUN = betweenness)

#Tables of betweenness distribution

betweeness_distribution <- list()

for (i in 1:length(betweenness_values)) {
  # Aplicar cluster_infomap a cada grafo y almacenar el resultado en results_list
  betweeness_distribution[[i]] <- data.frame(gene = names(betweenness_values[[i]]), degree = betweenness_values[[i]])
}

#List to sabe quantile thresholds

q_threshold.be <- list()

for (i in 1:length(betweenness_values)) {
  q_threshold.be[[i]] <- quantile(betweenness_values[[i]], probs = 0.95)
}

#Plot betweeness distribution

AD_be_distribution <- ggplot(betweeness_distribution[[1]], aes(x = degree, y = ..density..)) +
  geom_histogram(binwidth = 1, fill = "#00688B", color = "black") +
  labs(title = "Betweeness distribution histogram",
       subtitle = "for patients with AD",
       x = "Betweeness", y = "Freq") +
  geom_vline(xintercept = q_threshold.be[[1]], color = "red", linetype = "dashed") +
  geom_text(aes(x = q_threshold.be[[1]], y = 0.06, label = "95th percentile"), color = "red", hjust = -0.1) +
  theme_minimal()

noAD_be_distribution <- ggplot(betweeness_distribution[[2]], aes(x = degree, y = ..density..)) +
  geom_histogram(binwidth = 1, fill = "#00688B", color = "black") +
  labs(title = "Betweeness distribution histogram",
       subtitle = "for patients with no AD",
       x = "Betweeness", y = "Freq") +
  geom_vline(xintercept = q_threshold.be[[2]], color = "red", linetype = "dashed") +
  geom_text(aes(x = q_threshold.be[[2]], y = 0.06, label = "95th percentile"), color = "red", hjust = -0.1) +
  theme_minimal()

grid.arrange(AD_be_distribution, noAD_be_distribution)

#High betweeness nodes for AD patients --- ---

high_betweenness_nodes_AD <- V(graphAD)$name[betweenness_values[[1]] > q_threshold.be[[1]]]

high_betweenness_nodes_AD <- betweenness_values[[1]][high_betweenness_nodes_AD]

high_betweenness_nodes_AD <- data.frame(ensembl_gene_id = names(high_betweenness_nodes_AD), betweenness_values = high_betweenness_nodes_AD)

# Convert ENSMBL a SYMBOL

high_betweenness_nodes_AD_trad <- convert_ens_to_symbol(high_betweenness_nodes_AD$ensembl_gene_id)

# Merge gene names traduction

high_betweenness_nodes_AD <- high_betweenness_nodes_AD %>% left_join(high_betweenness_nodes_AD_trad, by = "ensembl_gene_id")

high_betweenness_nodes_AD <- high_betweenness_nodes_AD[order(high_betweenness_nodes_AD$betweenness_values, decreasing = TRUE), ]

#Plot --- ---

#df

#Create data frame 
high_betweenness_nodes_AD.df <- data.frame(ensembl_gene_id = names(betweenness_values[[1]]), betweeness = betweenness_values[[1]])
#Add check column
high_betweenness_nodes_AD.df$high_betweenness <- ifelse(high_betweenness_nodes_AD.df$ensembl_gene_id %in% high_betweenness_nodes_AD$ensembl_gene_id, "High betweeness Gene", "Not High betweeness Gene")
dim(high_betweenness_nodes_AD.df)

#Translate again
high_betweenness_nodes_AD.df_trads <- convert_ens_to_symbol(high_betweenness_nodes_AD.df$ensembl_gene_id)
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

high_betweenness_nodes_noAD <- high_betweenness_nodes_noAD[order(high_betweenness_nodes_noAD$betweenness_values, decreasing = TRUE), ]

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

# #Save plot
# ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNA-seq_DLFPC_NIA_reagan_high_be_genes.png", 
#        grid_high_be_genes, width = 20, height = 10, dpi = 300)

#What genes do both networks share?

shared_high_be_genes <- intersect(high_betweenness_nodes_AD$external_gene_name, high_betweenness_nodes_noAD$external_gene_name)

#

high_be_genes_en_AD_no_noAD_sym <- setdiff(high_betweenness_nodes_AD$external_gene_name, high_betweenness_nodes_noAD$external_gene_name)
high_be_genes_en_AD_no_noAD_ens<- setdiff(high_betweenness_nodes_AD$ensembl_gene_id, high_betweenness_nodes_noAD$ensembl_gene_id)

#Comparison of hub genes not found in healthy people

high_be_genes_en_AD_no_noAD.df <- high_betweenness_nodes_AD %>% filter(external_gene_name %in% high_be_genes_en_AD_no_noAD_sym)
high_be_genes_en_AD_no_noAD.df <- high_be_genes_en_AD_no_noAD.df[order(high_be_genes_en_AD_no_noAD.df$betweenness_values, decreasing = TRUE), ]

dim(high_be_genes_en_AD_no_noAD.df)
#[1] 66  5

#Network of hub genes in the AD graph but not in the noAD graph --- ---

indexed_ver_be <- which(V(graphAD)$name %in% high_be_genes_en_AD_no_noAD_ens)

high_be_genes_en_AD_no_noAD.g <- induced_subgraph(graphAD, indexed_ver_be)

plot(high_be_genes_en_AD_no_noAD.g)

#Translate gene names in graph --- ---

xhigh_be_genes_en_AD_no_noAD.g <- translate_vertex_names(high_be_genes_en_AD_no_noAD.g)

#Save graph

# write_graph(xhigh_be_genes_en_AD_no_noAD.g, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_high_be_genes_AD_noAD.graphml',
#             format = "graphml")

#NEXT QUESTION IS What modules do these genes belong to? --- ---

#Save hub genes to explore them in the partitions

#write(genes_en_AD_no_noAD_ens, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/genes_en_AD_no_noAD_ens.txt')

#Save high betweeness genes to explore them in the partitions

#write(high_be_genes_en_AD_no_noAD_ens, file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/high_be_genes_en_AD_no_noAD_ens.txt')

#END