#
#5.1 modularity analysis
#This script takes a list of genes and performs network communities analysis 
#with igraph and INFOMAP

#Libraries --- ---

pacman::p_load('igraph',
               'igraphdata',
               'clusterProfiler', 
               'enrichplot',
               'dplyr', 
               'pathview')

library("org.Hs.eg.db", character.only = TRUE)

#Get data --- ---

graph <- read_graph(file = '~/redesROSMAP/graphs/noAD_ROSMAP_RNAseq_MutualInfograph_percentile99.99.graphml',
                    format = 'graphml')

#Identify hub genes --- ---

nodes_degree <- degree(graph)

#Percentile 90 of genes with higher degree

q_threshold <- quantile(nodes_degree, probs = 0.99)

#Hub genes are 

hub_genes <- V(graph)$name[nodes_degree > q_threshold]

hub_genes <- nodes_degree[hub_genes]

#Analyze modules --- --- 

#Calculate clustering coefficient of whole network

clustering_coefficient <- transitivity(graph, type = 'undirected')

#Extract all genes from whole network

universe_of_genes <- V(graph)$name #all the genes in the network

#Calculations on components

components <- components(graph)

#Number of components

no_components <- components(graph)$no

#See membership of nodes

membership <- membership(components)

#Extract list of nodes by community

nodes_by_community <- split(V(graph)$name, membership)

#Enrichment by community --- ---

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

#Table

head(gene_o_enrich@result)

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

#KEGG Over-representation analysis --- ---

# We convert the gene IDs for the function enrichKEGG
# It could be that we lose some genes here since some conversions are not compatible.

ids <- bitr(universe_of_genes, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")

# Removemos IDS duplicados (aquí se usa "ENSEMBL", pero debería de ser lo que hayamos usado como keyType)

dedup_ids <- ids[!duplicated(ids[c("ENSEMBL")]),]

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

# Pathview analysis --- ---

#BiocManager::install("pathview")
library(pathview)

# Produce una gráfica de KEGG (PNG)
hsa <- pathview(gene.data=universe_of_genes, pathway.id="hsa04740", species = "hsa", gene.idtype=gene.idtype.list[3])

# Produce una gráfica diferente (PDF)
hsa <- pathview(gene.data=gene_list, pathway.id="hsa04740", species = "hsa", gene.idtype=gene.idtype.list[3], kegg.native = FALSE)

#END

#Alternative clustering with INFOMAP

cluster_infomap <- cluster_infomap(graph)

membership_infomap <- membership(cluster_infomap)

comm_infomap <- communities(cluster_infomap)

plot(membership_infomap)
