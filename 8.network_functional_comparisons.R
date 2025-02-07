#
#8.network_functional_comparisons.R

#This script makes the functional modular comparison of coexpression graphs 
#constructed for people with AD and people without pathological AD.

#Modules will be compared in terms of the associated biological functions identified by an enrichment analysis.
#This comparison involves answering three complementary questions: 

#1.How similar are the sets of biological functions that are associated to the whole network, through the enrichment of individual modules?
#2.How similar are the modules found in each network, in terms of the sets of associated biological functions?
#3.In how many modules is represented each biological process?

#paulinapglz.99@gmail.com
#Part of code adapted from https://github.com/guillermodeandajauregui/BiologicalModuleComparison/blob/master/comparisonParameters.R

#Libraries --- ---

pacman::p_load("igraph", 
               "tidyverse",
               "clusterProfiler", 
               "gridExtra", 
               "biomaRt", 
               "ggraph",
               "tidyheatmaps")

library("org.Hs.eg.db", character.only = TRUE)

#Set seed for modularity algorithm --- ---

set.seed(10)

# Define functions --- ---

#Define similarity of Enriched Processes, Jaccard Index function 

jaccard_simplex <- function(a,b){
  length(intersect(a,b))/length(union(a,b))
}

# Calculate Jaccard Index 
jaccard_index <- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(union(set1, set2))
  return(intersection / union)
}

#Create a function to calculate the assortativity per community
calculate_assortativity <- function(graph, nodes) {
  # Subgrafo para el módulo
  subgraph <- induced_subgraph(graph, vids = nodes)
  
  # Calcular la asortatividad por grado
  assortativity_degree <- assortativity(subgraph, 
                                        types1 = degree(subgraph, mode = "all"))
  
  return(assortativity_degree)
}

#Get data --- --- 

graphAD <- read_graph(file =  '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml',
                        format = 'graphml')

universe <- scan(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/universe.txt", what = character())

#Save graphs in a list

graphLists <- list(graphAD = graphAD,
                   graphnoAD = graphnoAD)

#Comparison of biological function sets associated to the overall network --- ---

#Enrichment of full graphs

#For AD graph

enrichment_fullnet_AD <- enrichGO(gene = V(graphAD)$name,
                                  OrgDb = "org.Hs.eg.db", 
                                  keyType = 'ENSEMBL',
                                  readable = TRUE,
                                  universe = universe, 
                                  ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                                  pvalueCutoff = 0.05, 
                                  qvalueCutoff = 0.10)

#Gene concept network enrichment

enrichment_fullnet_AD_cnet <- cnetplot(enrichment_fullnet_AD, showCategory= 5)

#Dotplot enrichment

enrichment_fullnet_AD_dot <- dotplot(enrichment_fullnet_AD)

#Save plots

#ggsave("enrichment_fullnet_AD_cnet.png", enrichment_fullnet_AD_cnet, width = 15, height = 8)

#For noAD graph

enrichment_fullnet_noAD <- enrichGO(gene = V(graphnoAD)$name,
                                    OrgDb = "org.Hs.eg.db", 
                                    keyType = 'ENSEMBL',
                                    readable = TRUE,
                                    universe = universe, 
                                    ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                                    pvalueCutoff = 0.05, 
                                    qvalueCutoff = 0.10)

#Gene concept network enrichment

enrichment_fullnet_noAD_cnet <- cnetplot(enrichment_fullnet_noAD, showCategory= 5)

#Dotplot enrichment

enrichment_fullnet_noAD_dot <- dotplot(enrichment_fullnet_noAD)

#Comparison of enrichments --- ---

OverallEnrichedProcessJ <- jaccard_simplex(names(enrichment_fullnet_AD@geneSets), names(enrichment_fullnet_noAD@geneSets))
#[1] 0.8042686

#This answers the question 2. How similar are the modules found in each network, in terms of the sets of associated biological functions?
  
#Similarity of modules in terms of associated biological functions --- ---

#Modularity algorithm 

modularities <- sapply(X = graphLists, FUN = cluster_infomap)

#Count modules

len_mod <- sapply(X = modularities, FUN = length)
# graphAD graphnoAD 
# 67        72

#Calculate Q score 

Qscore <- sapply(modularities, FUN = modularity)
Qscore
# graphAD graphnoAD 
# 0.2825230 0.2027528 

#Calculate clustering coefficient

clus_coe <- sapply(graphLists, FUN = transitivity)
clus_coe
# graphAD graphnoAD 
# 0.5956005 0.6252799 

#Split lists of nodes by module

nodes_membership <- sapply(modularities, FUN = membership)

#Create df
nodes_membership_AD.df <- data.frame(ensembl_gene_id = names(nodes_membership$graphAD),  membership = nodes_membership$graphAD)

nodes_membership_noAD.df <- data.frame(ensembl_gene_id = names(nodes_membership$graphnoAD),  membership = nodes_membership$graphnoAD)

# Assign membership to the nodes of each network
for (i in seq_along(graphLists)) {
  V(graphLists[[i]])$community <- nodes_membership[[i]]
}

#Save graphs  --- ----

#Loop through each graph and export it with community membership as GraphML
for (i in seq_along(graphLists)) {
  # Generate filename for each graph (e.g., graphAD.graphml, graphnoAD.graphml)
  graph_name <- ifelse(i == 1, "graphAD", "graphnoAD")
  filename <- paste0(graph_name, "nodes_membership.graphml")

  # Export graph to GraphML format
  write_graph(graphLists[[i]], file = filename, format = "graphml")
}

#If you already have the networks

graphAD <- read_graph(file =  'graphADnodes_membership.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = 'graphnoADnodes_membership.graphml',
                        format = 'graphml')

universe <- scan(file = "universe.txt", what = character())

#Save graphs in a list

graphLists <- list(graphAD = graphAD,
                   graphnoAD = graphnoAD)

#Extract list of nodes by community for each graph
nodes_by_community_list <- lapply(graphLists, function(graph) {
  split(V(graph)$name, V(graph)$community)
})

#Extract modules for both graphs
modules_AD <- nodes_by_community_list[[1]]
modules_noAD <- nodes_by_community_list[[2]]

#Create empty matrix
similarity_matrix <- matrix(0, nrow = length(modules_AD), ncol = length(modules_noAD))

#Fill matrix with jaccard inndex
for (i in seq_along(modules_AD)) {
  for (j in seq_along(modules_noAD)) {
    similarity_matrix[i, j] <- jaccard_index(modules_AD[[i]], modules_noAD[[j]])
  }
}

#Name rows and columns
rownames(similarity_matrix) <- paste0("AD_", seq_along(modules_AD))
colnames(similarity_matrix) <- paste0("Control_", seq_along(modules_noAD))

similarity_matrix
dim(similarity_matrix)
#[1] 68 71

#Is there any modules with perfect similitude?

length(similarity_matrix[similarity_matrix == 1])
#[1] 10

#Is there any modules with moderate similitude?

length(similarity_matrix[similarity_matrix > 0.5])
#[1] 24
 
#Positions equal to 1
pairs_with_one <- which(similarity_matrix == 1, arr.ind = TRUE)
dim(pairs_with_one)
#[1] 10  2

#Dataframe with correspondent modules
module_pairs <- data.frame(
  module_AD = rownames(similarity_matrix)[pairs_with_one[, "row"]],
  module_noAD = colnames(similarity_matrix)[pairs_with_one[, "col"]]
)

#Create heatmap

#Long format for ggplot

sim_heatmap.df <- as.data.frame(as.table(similarity_matrix)) %>%
  rename(Var1 = "module_AD", Var2 = "module_noAD", Freq = "similarity")

#Plot heatmap
sim_heatmap.p <- tidyheatmap(
  df = sim_heatmap.df,
  rows = module_AD,
  columns = module_noAD,
  values = similarity,
  scale = "none", # No queremos escalado adicional
  clustering_method = "average", # Opcional: Método de clustering
  annotation_col = NULL,    # Opcional: Anotaciones en columnas si es necesario
  annotation_row = NULL,     # Opcional: Anotaciones en filas si es necesario
  colors =  c("navy", "white", "firebrick"), 
  main = "Gene module correspondence"
)
sim_heatmap.p

# ggsave(
#   "sim_genes_heatmap.pdf",
#   plot = sim_heatmap.p,
#   device = "pdf",
#   width = 15,
#   height = 10,
#   units = "in",
#   dpi = 300
# )

#Calculate assortativity by module --- ---

assortativity_modules <- lapply(graphLists, function(graph) {
  nodes_by_community <- split(V(graph)$name, V(graph)$community)
  sapply(nodes_by_community, function(nodes) {
    calculate_assortativity(graph, nodes)
  })
})

assortativity__AD <- data.frame(
  name = as.factor(names(assortativity_modules$graphAD)),
  value = assortativity_modules$graphAD
) %>% filter(complete.cases(.)) %>% mutate(group = "AD")

assortativity__control <- data.frame(
  name = as.factor(names(assortativity_modules$graphnoAD)),
  value = assortativity_modules$graphnoAD
) %>% filter(complete.cases(.)) %>% mutate(group = "Control")

#normality test

shapiro.test(assortativity_modules$graphAD)

# Shapiro-Wilk normality test
# 
# data:  assortativity_modules$graphAD
# W = 0.86702, p-value = 0.0002845

shapiro.test(assortativity_modules$graphnoAD)

# Shapiro-Wilk normality test
# 
# data:  assortativity_modules$graphnoAD
# W = 0.88078, p-value = 0.00148

sd(assortativity_modules$graphAD, na.rm = TRUE)

sd(assortativity_modules$graphnoAD, na.rm = TRUE)

#T test

t_test <- t.test(assortativity_modules$graphAD, assortativity_modules$graphnoAD,
                 var.equal = T)  # var.equal = FALSE for unequal variances
t_test


#plot assortativity by module

library(ggpubr)

assortativity_AD.p <- ggdotchart(assortativity__AD, 
           x = "name", 
           y = "value",
          # color = "name",                         
          # palette = "blue",                       
           sorting = "descending",                  
           rotate = TRUE,                           
           dot.size = 3,                            
           y.text.col = FALSE,                      
           ggtheme = theme_pubr()                   
) +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(floor(min(assortativity__AD$value)), 
                                  ceiling(max(assortativity__AD$value)), 
                                  by = 0.1)) +
  ggtitle("Assortativity per module", "in AD network") +
  theme_cleveland()                                     
  
assortativity_control.p <- ggdotchart(assortativity__control, 
                                 x = "name", 
                                 y = "value",
                                # color = "name",                                
                                 # palette = "blue",                             
                                 sorting = "descending",                        
                                 rotate = TRUE,                               
                                 dot.size = 3,                                
                                 y.text.col = FALSE,                          
                                 ggtheme = theme_pubr()) +
  ggtitle("Assortativity per module", "in control network") +
  theme(legend.position = "none") +
  scale_y_continuous(breaks = seq(floor(min(assortativity__control$value)), 
                                  ceiling(max(assortativity__control$value)), 
                                  by = 0.1)) +
  theme_cleveland() 

#Grid both plots

assortativity_grid <- cowplot::plot_grid(assortativity_AD.p, assortativity_control.p,
                   align = "h", axis = "tb", labels = "auto",
          ncol = 1)

# ggsave(
#     "assortativity_grid.jpg",
#     plot = assortativity_grid,
#     device = "jpg",
#     width = 15,
#     height = 15,
#     units = "in",
#     dpi = 300
#   )

#Combine the two datasets
assortativity_df <- bind_rows(assortativity__AD, assortativity__control)

assortativity_viol <- ggplot(assortativity_df, aes(x = group, y = value, fill = group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +  # Violin shape
  geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +  # Optional: Add a boxplot inside
  scale_fill_manual(values = c("AD" = "firebrick", 
                               "Control" = "navy")) +  # Custom colors
  labs(title = " ",
       x = "Network",
       y = "Assortativity") +
  theme(legend.position = "none") +
  guides(fill = "none") +
  theme_minimal()

# ggsave(
#   "assortativity_viol.jpg",
#   plot = assortativity_viol,
#   device = "jpg",
#   width = 8,
#   height = 10,
#   units = "in",
#   dpi = 300
# )
  
#Compare modularity between graphs --- ---
#, applying 
#variation of information "vi"
#normalized mutual information "nmi"
#split-join distance "split-join distance"
#Rand index "Rand index"
#adjusted Rand index "adjusted Rand index"

possible_algos <- c("vi", "nmi", "split.join", "rand", "adjusted.rand")

comparison_methods <- sapply(X = possible_algos, FUN = function(i){
  igraph::compare(comm1 = graphAD_plus_modules,
                  comm2 = graphnoAD_plus_modules,
                  method = i
  )
})

comparison_methods

#Enrichment to all modules of a network --- ----

#I'd like to add names to lists of graphs again

names(nodes_by_community_list)[1] <- "graphAD"
names(nodes_by_community_list)[2] <- "graphnoAD"

#Extract list of nodes by community 

nodes_by_community_AD <- nodes_by_community_list[["graphAD"]]

nodes_by_community_noAD <- nodes_by_community_list[["graphnoAD"]]

#Define function that performs enrichment

enricher <- function(nodes_by_community) {
  enrichGO_result <- enrichGO(gene = nodes_by_community,
                              OrgDb = org.Hs.eg.db, 
                              universe = universe , 
                              keyType = 'ENSEMBL',
                              readable = TRUE,
                              ont = "BP",          #type of GO(Biological Process (BP), cellular Component (CC), Molecular Function (MF)
                              pvalueCutoff = 0.05, 
                              qvalueCutoff = 0.10)
  return(enrichGO_result)
}

#These functions may print a " --> No gene can be mapped...." on the console. 
#This happens when there are lists of genes that cannot be enriched. It will still generate an empty S4 object.

# Define a function to replace NULL by an empty S4 object "enrichResult" to handle NULLs in the enrichment lists

replace_null <- function(x) {
  if (is.null(x)) {
    return(new("enrichResult"))
  } else {
    return(x)
  }
}

#Enrichment of graphAD modules

enriched_results_AD <- lapply(nodes_by_community_AD, enricher) #slow
enriched_results_AD <- lapply(enriched_results_AD, replace_null) # Exchange NULLs for empty S4 objects

#Enrichment of graphAD modules
enriched_results_noAD <- lapply(nodes_by_community_noAD, enricher) #slow
enriched_results_noAD <- lapply(enriched_results_noAD, replace_null) # Exchange NULLs for empty S4 objects

#Functional module comparison between two graphs --- ---

#1.How similar are the sets of biological functions that are associated to the whole network, through the enrichment of individual modules?

#Similarity of modules in terms of associated biological functions

#It is possible to compare each module of the network of interest to the modules
#of another network, in terms of their associated biological function

# Create an empty array to store the results of the Jaccard index.
num_modules_AD <- length(enriched_results_AD)
num_modules_noAD <- length(enriched_results_noAD)

similarity_matrix_enri <- matrix(NA, nrow = num_modules_AD, ncol = num_modules_noAD)

#Assign row and column names
rownames(similarity_matrix_enri) <- paste("AD", 1:num_modules_AD, sep = "_")
colnames(similarity_matrix_enri) <- paste("control", 1:num_modules_noAD, sep = "_")

#Iterate to obtain similarity matrix 

for (i in 1:num_modules_AD) {
  for (j in 1:num_modules_noAD) {
    # Acceder a los geneSets de las redes AD y no AD
    geneSets_AD <- names(enriched_results_AD[[i]]@geneSets)
    geneSets_noAD <- names(enriched_results_noAD[[j]]@geneSets)
    similarity_matrix_enri[i, j] <- jaccard_simplex(geneSets_AD, geneSets_noAD)
  }
}

head(similarity_matrix_enri)

#Find the highest similarity of the matrix  
max(similarity_matrix_enri, na.rm = TRUE)

length(similarity_matrix_enri[similarity_matrix_enri == 1])
#[1] 25

length(similarity_matrix_enri[similarity_matrix_enri > 0.5])
#[1] 52

# Create data frame with module names
module_pairs_enri <- data.frame(
  module_AD = rownames(similarity_matrix_enri)[pairs_with_one[, "row"]],
  module_noAD = colnames(similarity_matrix_enri)[pairs_with_one[, "col"]]
)
print(module_pairs_enri)

sim_enri_heatmap.df <- as.data.frame(as.table(similarity_matrix_enri)) %>%
  rename(Var1 = "module_AD", Var2 = "module_control", Freq = "similarity")

# Heatmap
sim_enri_heatmap.p <- tidyheatmap(
  df = sim_enri_heatmap.df,
  rows = module_AD,
  columns = module_control,
  values = similarity,
  scale = "none", #No aditional scaling
  clustering_method = "average", # Clustering method
  annotation_col = NULL,    #Column annotations
  annotation_row = NULL,     #Rows annotations
  colors =  c("navy", "white", "firebrick"), 
  main = "Biological Processes module correspondence"
)

# ggsave(
#   "sim_enri_genes_heatmap.pdf",
#   plot = sim_enri_heatmap.p,
#   device = "pdf",
#   width = 15,
#   height = 10,
#   units = "in",
#   dpi = 300
# )

#Number of modules in networks that have Jaccard index J=1 with a module of the Main network --- ---

# Function to count the number of columns with value equal to 1
count_ones <- function(row) {
  sum(row == 1)
}

# Apply the function to each row of the matrix
ones_count <- apply(similarity_matrix_enri, 1, count_ones)

# Crear una tabla con los resultados
num_equal_nodes <- data.frame(module_number = rownames(similarity_matrix_enri), modules_with_jaccardindex1 = ones_count)

num_equal_nodes.x <- num_equal_nodes %>% filter(modules_with_jaccardindex1 ==1)
dim(num_equal_nodes.x)
#[1] 15  2

#This answers the question of 1. How similar are the sets of biological functions that are associated to the whole network, through the enrichment of individual modules?
  
# Extraer descripciones de los módulos de AD
modulos_ad <- names(enriched_results_AD)
funciones_ad <- sapply(modulos_ad, function(mod) {
  enriched_results_AD[[mod]]@result$Description[1] # Obtener la función principal
})

funciones_ad <- sapply(modulos_ad, function(mod) {
  if (!is.null(enriched_results_AD[[mod]]@result$Description[1])) {
    enriched_results_AD[[mod]]@result$Description[1]  # Obtener la función principal
  } else {
    "Not enriched" 
  }
})

#Build dataframe
tabla_modulos <- data.frame(
  Modulo = paste0("AD_", modulos_ad),
  Module_main_function = funciones_ad,
  stringsAsFactors = FALSE
)
tabla_modulos$membership <- gsub("AD_", "", tabla_modulos$Modulo)

#Circular plot

library(circlize)

module_descriptions <- sapply(enriched_results_AD, function(res) res@result$Description[1])

#Create similarity matrix using jaccard_simplex
num_modules <- length(geneSets_AD)
jaccard_matrix <- matrix(0, nrow = num_modules, ncol = num_modules,
                         dimnames = list(module_descriptions, module_descriptions))

for (i in seq_along(geneSets_AD)) {
  for (j in seq_along(geneSets_AD)) {
    jaccard_matrix[i, j] <- jaccard_simplex(geneSets_AD[[i]], geneSets_AD[[j]])
  }
}
jaccard_matrix[is.na(jaccard_matrix)] <- 0

normalized_matrix <- jaccard_matrix / max(jaccard_matrix)

hc <- as.dendrogram(hclust(normalized_matrix))
circlize_dendrogram(hc,
                    labels_track_height = NA,
                    dend_track_height = 0.5)

dist_matrix <- as.dist(1 - similarity_matrix)  #Convert similarity to distance
clustering <- hclust(dist_matrix, method = "ward.D2")

################################################################################

#Where are hub and high betweeness genes?
#In which modules are found our hub genes? --- ---

library(ggsankey)
library(ggsankeyfier)

#Get hub genes 

hub_genes <- scan("AD_hubs_notcontrolens.txt",
                  what = character())

hub_genes.x <- nodes_membership_AD.df %>% filter(ensembl_gene_id %in% hub_genes)

hub_genes.xy <- convert_ens_to_symbol(hub_genes.x)

hub_genes.x <- hub_genes.x %>% left_join(hub_genes.xy, by ="ensembl_gene_id")
hub_genes.x <- merge(hub_genes.x, tabla_modulos[, c("membership", "Module_main_function")], by = "membership", all.x = TRUE)
hub_genes.x[2, 3] <- "PDE4DIP Pseudogene"

#hub_genes.x$chromosome_name <- as.factor(hub_genes.x$chromosome_name)
unique(hub_genes.x$chromosome_name)

hub_genes.x$membership <- factor(
  hub_genes.x$membership,
  levels = c("1", "2", "3", "6", "12", "17", "19"))


hub_genes.x$chromosome_name <- factor(hub_genes.x$chromosome_name, levels = myOrder)
hub_genes.x$external_gene_name <- as.factor(hub_genes.x$external_gene_name)
hub_genes.x$Module_main_function <- as.factor(hub_genes.x$Module_main_function)
hub_genes.x$membership <- as.factor(hub_genes.x$membership)

#Convert the table to the long format needed for ggsankey

hub_genes_long <- hub_genes.x %>%
  make_long(membership, external_gene_name, Module_main_function)

sankey_hubs <- ggplot(hub_genes_long, aes(x = x, 
                 next_x = next_x, 
                 node = node, 
                 next_node = next_node,
                 fill = factor(node),
                 label = node)) +
  scale_x_discrete(labels = c("membership" = "Module membership",
                              "external_gene_name" = "Gene", 
                              "Module_main_function" = "Module main function"))+
    geom_sankey(flow.alpha = 0.5, node.color = 1) +
    geom_sankey_label(size = 3.5, color = 1, fill = "white") +
  labs(x = " ")+
    scale_fill_viridis_d(option = "A", alpha = 0.95) +
    theme_sankey(base_size = 16) +
    theme(legend.position = 'none')

sankey_hubs

##
library(ggalluvial)

allu_hubs <- ggplot(data = hub_genes.x,
       aes(axis1 = chromosome_name, axis2 = external_gene_name, axis3 = Module_main_function)) +
  geom_alluvium(aes(fill = membership)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_fill_viridis_d() +
  theme_void()

#In which modules are found our high betweeness genes?--- ---

high_be_genes <-  scan(file = 'AD_highbe_notcontrol_ens.txt',   what = character())

high_be_genes.x <- nodes_membership_AD.df %>% filter(ensembl_gene_id %in% high_be_genes)

high_be_genes.xy <- convert_ens_to_symbol(high_be_genes.x)

high_be_genes.x <- high_be_genes.x %>% left_join(high_be_genes.xy, by ="ensembl_gene_id")
high_be_genes.x <- merge(high_be_genes.x, tabla_modulos[, c("membership", "Module_main_function")], by = "membership", all.x = TRUE)

high_be_genes.x [4,3]<- "(AURKA) Pseudogene"
high_be_genes.x [7,3]<- "TXNRD1"
high_be_genes.x [44,3]<-  "ENSG00000288049-001"
high_be_genes.x$Module_main_function <- gsub(
  pattern = "^detection of mechanical stimulus involved in sensory perception of pain",
  replacement = "detection of mechanical stimulus of pain",
  x = high_be_genes.x$Module_main_function
)

unique(high_be_genes.x$chromosome_name)

high_be_genes.x$chromosome_name <- factor(
  high_be_genes.x$chromosome_name,
  levels = c(1, 2, 3, 4, 5, 6, 7, 10, 11, 12, 14, 15, 18, 19, 20, 22, "X", "Y"))

#Parallel Coordinates chart

high_be_genes.x$chromosome_name <- as.factor(high_be_genes.x$chromosome_name)
high_be_genes.x$external_gene_name <- as.factor(high_be_genes.x$external_gene_name)
high_be_genes.x$Module_main_function <- as.factor(high_be_genes.x$Module_main_function)
high_be_genes.x$membership <- as.factor(high_be_genes.x$membership)

#Convert the table to the long format needed for ggsankey
highbe_genes_long <- high_be_genes.x %>%
  make_long(membership, external_gene_name, Module_main_function)

sankey_highbe <- ggplot(highbe_genes_long, aes(x = x, 
                                          next_x = next_x, 
                                          node = node, 
                                          next_node = next_node,
                                          fill = factor(node),
                                          label = node)) +
  scale_x_discrete(labels = c("membership" = "Module membership", 
                              "external_gene_name" = "Gene", 
                              "Module_main_function" = "Module main function"))+
  geom_sankey(flow.alpha = 0.5, node.color = 1, node.width = 0.8) +
  labs(x = "")+
  geom_sankey_label(size = 3.5, fill = "white") +
  scale_fill_viridis_d(option = "A", alpha = 0.95) +
  theme_sankey(base_size = 16) +
  theme(legend.position = 'none')

sankey_highbe

allu_highbe <- ggplot(data = high_be_genes.x,
       aes(axis1 = chromosome_name, axis2 = external_gene_name, axis3 = Module_main_function)) +
  geom_alluvium(aes(fill = membership)) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  scale_fill_viridis_d() +
  theme_void()


#Grids

sankeys <- grid.arrange(sankey_hubs, sankey_highbe, ncol = 2)

# 
# ggsave("~/DLPFC_paper_plots/sankeys_hubs_highs.jpg", 
#        sankeys, device = "jpg", 
#        height = 10, 
#        width = 15, 
#        units = "in", 
#        dpi = 300
#         )


################################################################################

#Function for counting genes and enriched processes
count_genes_and_processes <- function(enrichment_results) {
 
  results_summary <- list()
  
  for (i in seq_along(enrichment_results)) {
      num_genes <- length(enrichment_results[[i]]@gene)
    
      if (nrow(enrichment_results[[i]]) > 0) {
      num_processes <- nrow(enrichment_results[[i]]@result[enrichment_results[[i]]@result$pvalue < 0.05, ])
    } else {
      num_processes <- 0
    }
    
     results_summary[[i]] <- list(
      "num_genes" = num_genes,
      "num_processes" = num_processes
    )
  }
  return(results_summary)
}

# Enrichment summary for AD and non-AD
summary_AD <- count_genes_and_processes(enriched_results_AD)
summary_noAD <- count_genes_and_processes(enriched_results_noAD)

communities_AD <- names(nodes_by_community_AD)
communities_noAD <- names(nodes_by_community_noAD)

#As data frames
df_AD <- data.frame(
  Community = communities_AD,
  Category = "AD",
  Number_of_Genes = sapply(summary_AD, function(x) x$num_genes),
  Number_of_Enriched_Processes = sapply(summary_AD, function(x) x$num_processes)
)


df_noAD <- data.frame(
  Community = communities_noAD,
  Category = "noAD",
  Number_of_Genes = sapply(summary_noAD, function(x) x$num_genes),
  Number_of_Enriched_Processes = sapply(summary_noAD, function(x) x$num_processes)
)

#merge tables
df_summary <- bind_rows(df_AD, df_noAD)

#Number of modules associated to a given biological function --- ---

#3.In how many modules is represented each biological process?

#Extract GO terms for each module in AD
BP_terms_AD <- lapply(enriched_results_AD, function(result) {
  if (length(result@result) > 0) {
    return(result@result$Description)  # Extraer nombres de procesos biológicos
  } else {
    return(NA)  # Manejar módulos sin resultados
  }
})

#Extract GO terms for each module in noAD
BP_terms_noAD <- lapply(enriched_results_noAD, function(result) {
  if (length(result@result) > 0) {
    return(result@result$Description)
  } else {
    return(NA)
  }
})

#Count how many modules are associated with each biological process in AD.
BP_count_AD <- table(unlist(BP_terms_AD)) %>% as.data.frame()
colnames(BP_count_AD) <- c("term", "in_modules_AD")

#Count how many modules are associated with each biological process in control
BP_count_noAD <- table(unlist(BP_terms_noAD)) %>% as.data.frame()
colnames(BP_count_noAD) <- c("term", "in_modules_noAD")

#Unify the terms of biological processes of both networks.
all_BP_terms <- BP_count_AD %>% left_join(BP_count_noAD, by = "term")
all_BP_terms[is.na(all_BP_terms)] <- 0  # Reemplazar NA por 0
#absolute difference
all_BP_terms <- all_BP_terms %>%
  mutate(difference = abs(in_modules_AD - in_modules_noAD))

#Plot

ggplot(all_BP_terms, aes(x = difference)) +
  geom_histogram(binwidth = 1, fill = "cadetblue", color = "black", alpha = 0.6) +  # Bins más pequeños
  scale_x_continuous(breaks = seq(0, max(all_BP_terms$difference), by = 1)) +  # Eje x de 1 en 1
    labs(
    title = "Distribution of absolute differences",
    subtitle = "in the distribution of biological processes between the two models", 
    x = "Absolute difference (number of modules)",
    y = "Freq"
  ) +
  theme_minimal()

#
all_diff_BP_terms.p <- ggplot(all_BP_terms, aes(x = reorder(term, difference))) +  
  geom_point(aes(y = in_modules_AD, color = "AD"), size = 1, alpha = 0.8) +
  geom_point(aes(y = in_modules_noAD, color = "control"), size = 1, alpha = 0.8) +
  coord_flip() +
  geom_segment(aes(xend = term, y = in_modules_AD, yend = 0, color = "AD"), size = 0.5, alpha = 0.8) +
  geom_segment(aes(xend = term, y = in_modules_noAD, yend = 0, color = "control"), size = 0.5, alpha = 0.8) +
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  labs(x = "Biological Process (GO:BP)", y = "Number of Modules in which they are represented", color = "Network") +
  ggtitle("Absolute differences in representation of Biological Processes") +
  scale_color_manual(values = c("AD" = "red", "control" = "blue"))  # Colores personalizados
all_diff_BP_terms.p

all_diff_BP_terms.p <- ggplot(all_BP_terms, aes(x = reorder(term, difference))) +  # Reordenar según la diferencia en orden descendente
  geom_col(aes(y = in_modules_AD, fill = "AD"), stat = "identity", position = "dodge") +
  geom_col(aes(y = in_modules_noAD, fill = "control", alpha = 0.8), stat = "identity", position = "dodge") +
  coord_flip() + 
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  labs(x = "Biological Process (GO:BP)", y = "Number of Modules", fill = "Network") +
  ggtitle("Top 20 Biological Processes with Largest Differences (AD vs control)")+
  scale_fill_manual(values = c("AD" = "red", "control" = "blue"))  # Colores personalizados
all_diff_BP_terms.p

#Select the 50 terms with the biggest difference 

top_diff_BP_terms <- all_BP_terms %>%
  arrange(desc(difference)) %>%
  head(50)

diff_BP_terms.p <- ggplot(top_diff_BP_terms, aes(x = reorder(term, difference))) +  
  geom_point(aes(y = in_modules_AD, color = "AD"), size = 4, alpha = 0.8) +
  geom_point(aes(y = in_modules_noAD, color = "control"), size = 4, alpha = 0.8) +
coord_flip() +
geom_segment(aes(xend = term, y = in_modules_AD, yend = 0, color = "AD"), size = 1, alpha = 0.8) +
  geom_segment(aes(xend = term, y = in_modules_noAD, yend = 0, color = "control"), size = 1, alpha = 0.8) +
  geom_text(aes(y = max(in_modules_AD, in_modules_noAD) + 1, label = difference), 
            hjust = 0, size = 3) +  # Etiquetas al final de las barras
  theme_minimal() +
  labs(x = "Biological Process (GO:BP)", y = "Number of Modules in which they are represented", color = "Network") +
  ggtitle("Top 50 Biological Processes with Largest Differences (AD vs control)") +
  scale_color_manual(values = c("AD" = "red", "control" = "blue"))  # Colores personalizados
diff_BP_terms.p

#Save plot

# ggsave(filename = "diff_BP_terms.jpg",
#        plot = diff_BP_terms.p,
#        device = "jpg", width = 25,
#        height = 20, units = "cm",dpi = 300)

#Generate a list of connections between modules and GO processes
network_edges <- data.frame()

for (i in seq_along(enriched_results_AD)) {
  module_name <- names(nodes_by_community_AD)[i]  #Module name
  go_terms <- enriched_results_AD[[i]]@result$ID  #Go terms
  
  if (length(go_terms) > 0) {
    for (go in go_terms) {
      network_edges <- rbind(network_edges, data.frame(Module = module_name, GO = go))
    }
  }
}

#Create graph object from our network
g <- graph_from_data_frame(network_edges, directed = FALSE)

#Adding properties to nodes: whether it is a module or a GO term
V(g)$type <- ifelse(V(g)$name %in% names(nodes_by_community_AD), "Module", "GO Term")

#Plot
ggraph(g, layout = "fr") +  # Fruchterman-Reingold layout para redes
  geom_edge_link(aes(edge_alpha = 0.5), color = "gray", show.legend = FALSE) +  # Enlaces
  geom_node_point(aes(color = V(g)$type), size = 5) +  # Nodos
  geom_node_text(aes(label = V(g)$name), repel = TRUE, size = 3) +  # Etiquetas
  scale_color_manual(values = c("Module" = module_color, "GO Term" = go_color)) +  # Asignar colores
  theme_void() +  # Sin fondo
  theme(legend.position = "none")  # Sin leyenda

#Create a data frame containing information about the modules and their genes/processes.
df_AD <- data.frame(
  Module = names(nodes_by_community_AD),
  Number_of_Genes = sapply(summary_AD, function(x) x$num_genes),
  Number_of_Enriched_Processes = sapply(summary_AD, function(x) x$num_processes)
)

df_noAD <- data.frame(
  Module = names(nodes_by_community_noAD),
  Number_of_Genes = sapply(summary_noAD, function(x) x$num_genes),
  Number_of_Enriched_Processes = sapply(summary_noAD, function(x) x$num_processes)
)

#Keeping the first 20 modules of each network
df_AD <- df_AD[order(-df_AD$Number_of_Enriched_Processes), ]
df_noAD <- df_noAD[order(-df_noAD$Number_of_Enriched_Processes), ]

#Only top 20
top_modules_AD <- df_AD[1:20, ]
top_modules_noAD <- df_noAD[1:20, ]

#Generate a list of connections between filtered modules and GO processes
filtered_enrichment_AD <- enriched_results_AD[names(enriched_results_AD) %in% top_modules_AD$Module]
filtered_enrichment_noAD <- enriched_results_noAD[names(enriched_results_noAD) %in% top_modules_noAD$Module]
#Generate a list of connections between filtered modules and GO processes
network_edges_filtered <- data.frame()

#For AD
for (i in seq_along(filtered_enrichment_AD)) {
  module_name <- names(filtered_enrichment_AD)[i]  # Nombre de la comunidad o módulo
  go_terms <- filtered_enrichment_AD[[i]]@result$ID  # Procesos GO significativos para esa comunidad
  
  if (length(go_terms) > 0) {
    for (go in go_terms) {
      network_edges_filtered <- rbind(network_edges_filtered, data.frame(Module = module_name, GO = go))
    }
  }
}

#For control
for (i in seq_along(filtered_enrichment_noAD)) {
  module_name <- names(filtered_enrichment_noAD)[i]  # Nombre de la comunidad o módulo
  go_terms <- filtered_enrichment_noAD[[i]]@result$ID  # Procesos GO significativos para esa comunidad
  
  if (length(go_terms) > 0) {
    for (go in go_terms) {
      network_edges_filtered <- rbind(network_edges_filtered, data.frame(Module = module_name, GO = go))
    }
  }
}

#Filtered graph
g_filtered <- graph_from_data_frame(network_edges_filtered, directed = FALSE)

#Add 
V(g_filtered)$type <- ifelse(V(g_filtered)$name %in% c(top_modules_AD$Module, top_modules_noAD$Module), "Module", "GO Term")

#plot

ggraph(g_filtered, layout = "fr") +  
  geom_edge_link(aes(edge_alpha = 0.5), color = "gray", show.legend = FALSE) +  
  geom_node_point(aes(color = V(g_filtered)$type), size = 5) + 
  geom_node_text(aes(label = V(g_filtered)$name), repel = TRUE, size = 3) +  #
  scale_color_manual(values = c("Module" = module_color, "GO Term" = go_color)) +  
  theme_void() +  
  theme(legend.position = "none") 
#END