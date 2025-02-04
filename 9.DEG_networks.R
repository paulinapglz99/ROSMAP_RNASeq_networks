#
#9.1 DEG_induced_subgraph.R
#This script filters a full edgelist to obtain nodes of the DEGS previously obtained

#Libraries --- ---
pacman::p_load("dplyr", 
               "vroom", 
               "igraph")

#Get data --- ---

#Differentially expressed genes
DEGS <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_DEGS_dichoNIAReagan.txt")

#Graphs

graph_0 <- read_graph(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml", format = "graphml")  
graph_1 <-  read_graph(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml", format = "graphml")  

#Convert into graph
graph_0.degree <- data.frame(ensembl_gene_id = names(degree(graph_0)), degree_g_0 = degree(graph_0))

graph_1.degree <- data.frame(ensembl_gene_id = names(degree(graph_1)), degree_g_1 = degree(graph_1))

#See if there's hub genes in the graphs

table(DEGS$ensembl_gene_id %in% graph_0.degree$ensembl_gene_id)
table(DEGS$ensembl_gene_id %in% graph_1.degree$ensembl_gene_id)

#If false, there will not be induced subgraph, sorry

# Left join

DEGS <- DEGS %>% left_join(graph_0.degree, by = "ensembl_gene_id")

DEGS <-DEGS %>% left_join(graph_1.degree, by = "ensembl_gene_id")

#Plot scatter --- ---

#BaseMean

DEGS_BM_0 <- ggplot(DEGS, aes(x=degree_g_0, y=baseMean_NIA_R_0)) + 
  geom_point() +
  geom_vline(xintercept = max(DEGS$degree_g_0)/2, col="red") + 
  geom_hline(yintercept = max(DEGS$baseMean_NIA_R_0)/2, col="red") +
  xlab("Degree") +
  ylab("baseMean") +
  labs(title = "", 
       subtitle = "NIA-Reagan 0")+
  theme_classic()
DEGS_BM_0

DEGS_BM_1 <- ggplot(DEGS, aes(x=degree_g_1, y=baseMean_NIA_R_1)) + 
  geom_point() +
  geom_vline(xintercept = max(DEGS$degree_g_1)/2, col="red") + 
  geom_hline(yintercept = max(DEGS$baseMean_NIA_R_1)/2, col="red") +
  xlab("Degree") +
  ylab("baseMean") +
  labs(title = "", 
       subtitle = "NIA-Reagan 1") +
  theme_classic()
DEGS_BM_1

#logfoldChange

DEGS_lfc_0 <- ggplot(DEGS, aes(x=degree_g_0, y=log2FoldChange)) + 
  geom_point() +
  geom_vline(xintercept = max(DEGS$degree_g_0)/2, col="red") + 
  geom_hline(yintercept = max(DEGS$log2FoldChange)/2, col="red") +
  xlab("Degree") +
  ylab("log2FoldChange") +
  labs(title = "", 
       subtitle = "NIA-Reagan 0") +
  theme_classic()
DEGS_lfc_0

DEGS_lfc_1 <- ggplot(DEGS, aes(x=degree_g_1, y=log2FoldChange)) + 
  geom_point() +
  geom_vline(xintercept = max(DEGS$degree_g_1)/2, col="red") + 
  geom_hline(yintercept = max(DEGS$log2FoldChange)/2, col="red") +
  xlab("Degree") +
  ylab("log2FoldChange") +
  labs(title = "", 
       subtitle = "NIA-Reagan 1") +
  theme_classic()
DEGS_lfc_0
#Save plots



#END
