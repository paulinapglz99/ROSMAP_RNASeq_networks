#
#6.enrichment.R
#Script for a gene-list enrichment
#paulinapglz.99@gmail.com

#Installing ClusterProfiler

#Libraries

pacman::p_load("clusterProfiler", 
               "tidyverse", 
               "org.Hs.eg.db")

#BiocManager::install("org.Hs.eg.db", character.only = TRUE)

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- ---

graph <- read_graph(file = '~/redesROSMAP/graphs/AD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml',
                    format = 'graphml')

#Analyze modules with igraph --- --- 

#Number of components

no_components <- components(graph)$no

#Calculate clustering coefficient of whole network

clustering_coefficient <- transitivity(graph, type = 'undirected')

#Extract all genes from whole network

universe_of_genes <- V(graph)$name #all the genes in the network

#See membership of nodes
components <- components(graph)

membership <- membership(components)

#Extract list of nodes by community
nodes_by_community <- split(V(graph)$name, membership)

#Modularity with infomap --- --- 

#cluster_informap find community structure that minimizes the expected description length of a random walker trajectory 
cluster_infomap <- cluster_infomap(graph)

#See membership of nodes
membership_infomap <- membership(cluster_infomap)

#Extract list of nodes by community
infomap_nodes_by_community <- split(V(graph)$name, membership_infomap)

#Calculate the number of communities derived from infomap algorithm
no_communities_infomap <- length(infomap_nodes_by_community)

#Plot

grafo_tbl <- as_tbl_graph(graph)

ggraph(grafo_tbl) +
  geom_edge_link() +  # Color edges in gray
  geom_node_point(color = membership_infomap) # Set node size and color

#Enrichment --- ---

#Create enrichResult object

gene_o_enrich <- enrichGO(gene = universe_of_genes,
                          OrgDb = org.Hs.eg.db, 
                          keyType = 'ENSEMBL',
                          readable = TRUE,
                          ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)

#Observe GO object

head(gene_o_enrich)

#Plot 

upsetplot(gene_o_enrich)

#Barplot

barplot(gene_o_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)

#Dotplot

dotplot(gene_o_enrich)

#heatplot

heatplot(gene_o_enrich)

#Category Netplot

#The cnetplot depicts the linkages of genes and biological concepts (e.g. GO terms or KEGG pathways) 
#as a network (helpful to see which genes are involved in enriched pathways and genes that may belong to multiple annotation categories).

cnetplot(gene_o_enrich, categorySize="pvalue", showCategory = 10)

#Enriched GO induced graph visual representation of the relationships between Gene Ontology (GO) 
#terms and a set of genes of interest.
#directed acyclic graph (DAG) that depicts the relationships between significant GO terms.

#Nodes: Represent the  GO terms identified in the enrichment analysis. Gray nodes represent GO terms that are not significant in the enrichment analysis.
#Edges: Represent the relationships (parent-child) between GO terms in the Gene Ontology hierarchy. The direction of the arrow indicates the "is-a" relationship
#Coloring: Nodes can be colored according to different factors, such as the enrichment p-value (more significant terms are colored darker) or the category of the GO term (biological process, molecular function, cellular component).

goplot(gene_o_enrich)

#Get the similarity matrix, by default using Jaccard’s similarity index (JC)

term_similarities <- pairwise_termsim(gene_o_enrich)

#emaplot, Enrichment Map for enrichment result of over-representation test or gene set enrichment analysis

emapplot(term_similarities)

#Treeplot, Functional grouping tree diagram for enrichment result of over-representation test or gene set en-

treeplot(term_similarities)

#KEGG Over-representation analysis --- ---

# We convert the gene IDs for the function enrichKEGG
# It could be that we lose some genes here since some conversions are not compatible.

ids <- bitr(universe_of_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

# Removemos IDS duplicados (aquí se usa "ENSEMBL", pero debería de ser lo que hayamos usado como keyType)

#
kk <- enrichKEGG(gene = ids$ENTREZID, 
                 organism="hsa",
                 pvalueCutoff = 0.05, 
                 keyType = "ncbi-geneid")
head(kk)

barplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)

dotplot(kk, 
        showCategory = 10, 
        title = "Enriched Pathways",
        font.size = 8)

