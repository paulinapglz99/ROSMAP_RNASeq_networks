#
#5.1 modularity analysis
#This script takes a list of genes and performsa variety of analysis with igraph and INFOMAP

#Libraries --- ---

pacman::p_load('igraph',
               'igraphdata',
               'clusterProfiler', 
               'enrichplot',
               'dplyr')

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- ---

graph <- read_graph(file = '~/redesROSMAP/graphs/graphical_graphs/',
                    format = 'graphml')

#Analyze modules --- --- 

#Calculate clustering coefficient of whole network

clustering_coefficient <- transitivity(graph, type = 'undirected')

#Calculations on components

components <- components(graph)

#Number of components

no_components <- components(graph)$no
no_components

#Calculate community membership

unverse_of_genes <- V(graph)$name #all the genes in the network

#List of nodes by community

nodes_by_community <- split(V(graph)$name, membership)

#Enrichment by community --- ---

#Create enrichResult object

gene_o_enrich <- enrichGO(gene = universe_of_genes,
                      OrgDb = org.Hs.eg.db, 
                      keyType = 'ENSEMBL',
                      readable = TRUE,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)

#Observe GO object

head(gene_o_enrich)

#Plot 

upsetplot(gene_o_enrich)

#Barplot

barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)

#Dotplot

dotplot(go_enrich)

#END