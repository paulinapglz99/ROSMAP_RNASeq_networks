#
#4.my_network_cut_by_percentile_heuristic.R

#paulinapglz.99@gmail.com

#NOTE for this script
#There are several approaches to cut the number of interactions in a network. In this pipeline, we look at three heuristics:

#1. Mutual Information cut-off according to the MI distribution of a network
#2. The single component minimum graph. the minimum one-component graph is obtained by making iterative cuts to find the interaction cut to obtain a single-component graph 
#topological heuristic that looks for a graph that maintains a unique component with as few interactions as possible. 
#3.The percentile-cut approach is explained as a method to obtain the percentiles with edges with the highest
#mutual information (top edges). 

#You can choose an approach depending on what you want to analyze or compare them. I recommend to get your data and only use one 
#heuristic at a time

#Here we code the 

########--- 3rd HEURISTIC: percentile cut-off --- ########

#Libraries  --- --- 

pacman::p_load('tidyverse', 
               'igraph', 
               'ggplot2')

#Get data --- ---

full_edgelist <- vroom::vroom(file = '/datos/rosmap/coexpre_matrix/full_net_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_normalizedMI_edgelist.tsv.gz')

#The P percentile indicates the value below which a specific percentage of the data falls.

# Percentile 99.9999 --- --- 

percentile_99.9999 <- quantile(as.numeric(full_edgelist$MI), 0.999999)
#   99.9999%  0.9647999 

percentile_99.9999_table <- subset(full_edgelist, as.numeric(MI) > percentile_99.9999)
dim(percentile_99.9999_table)
#[1] 224   3

#Histogram of MI distribution 

png("MI_distribution_percentile99.9999_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(percentile_99.9999_table$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "MI distribution histogram",   # Título principal
     sub = "patients with AD, 99.9999 percentile",   # Subtítulo
     xlab = "MI",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#write graph

graph_percentile99.9999 <- graph_from_data_frame(percentile_99.9999_table, directed = F)

#Graph info

summary(graph_percentile99.9999)

#IGRAPH 04418d4 UN-- 83 224 -- 
#  + attr: name (v/c), MI (e/x)

#number of components

components(graph_percentile99.9999)

#Degree by node 

degreee_by_node_graph_percentile99.9999<- degree(graph_percentile99.9999) %>% as.data.frame()

png("degreee_by_node_graph_percentile99.9999_MI_degree_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(degreee_by_node_graph_percentile99.9999$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Degree by node histogram",   # Título principal
     sub = "patients with AD, graph_percentile 99.9999%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Coreness by node

coreness_by_node_graph_percentile99.9999 <- coreness(graph_percentile99.9999) %>% as.data.frame()

png("graph_percentile99.9999_MI_coreness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_graph_percentile99.9999$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Coreness by node histogram",   # Título principal
     sub = "patients with AD, graph percentile 99.9999%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()


#edge betweenneess by node

betweeness_by_node_graph_percentile99.9999 <-betweenness(min_sin_comp_graph, directed = F) %>% as.data.frame()

png("min_sin_comp_MI_betweeness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_min_sin_g$.), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Betweeness by node histogram",   # Título principal
     sub = "patients with AD, minimal single component",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Plot graph, slow

png('graph_percentile99.9999_graph_MI_0.96_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.png')
plot(graph_percentile99.9999, edge.width = 2, edge.color = "black")
dev.off()

#write graphs in graphml

write_graph(graph_percentile99.9999,
            file = 'graph_percentile99.9999_graph_MI_0.96_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
            'graphml') #'edgelist.txt' / '.graphml'

# Percentile 99.999 --- --- 
percentile_99.999 <- quantile(as.numeric(full_edgelist$MI), 0.99999)
#  99.999% 0.853441

percentile_99.999_table <- subset(full_edgelist, as.numeric(MI) > percentile_99.999)
dim(percentile_99.999_table)
#[1] 22354     3

#Histogram of MI distribution 

png("MI_distribution_percentile9.999_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(percentile_99.999_table$MI), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "MI distribution histogram",   # Título principal
     sub = "patients with AD, 99.999 percentile",   # Subtítulo
     xlab = "MI",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#write graph

graph_percentile99.999 <- graph_from_data_frame(percentile_99.999_table, directed = F)

#Graph info

summary(graph_percentile99.999)

#IGRAPH a66e9f4 UN-- 1042 22354 -- 
#+ attr: name (v/c), MI (e/x)

#number of components

components(graph_percentile99.999)

#Degree by node 

degreee_by_node_graph_percentile99.999<- degree(graph_percentile99.999) %>% as.data.frame()

png("degreee_by_node_graph_percentile99.999_MI_degree_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(degreee_by_node_graph_percentile99.999$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Degree by node histogram",   # Título principal
     sub = "patients with AD, graph_percentile 99.9999%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Coreness by node

coreness_by_node_graph_percentile99.999 <- coreness(graph_percentile99.999) %>% as.data.frame()


png("graph_percentile99.999_MI_coreness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_graph_percentile99.999$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Coreness by node histogram",   # Título principal
     sub = "patients with AD, graph percentile 99.999%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Plot graph, slow

png('graph_percentile99.999_graph_MI_0.853_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.png')
plot(graph_percentile99.999, edge.width = 2, edge.color = "black")
dev.off()

#write graphs in graphml

write_graph(graph_percentile99.999,
            file = 'graph_percentile99.999_graph_MI_0.853_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
            'graphml') #'edgelist.txt' / '.graphml'

# Percentile 99.99 --- --- 

percentile_99.99 <- quantile(as.numeric(full_edgelist$MI), 0.9999)
#   99.99%  0.7117147 

percentile_99.99_table <- subset(full_edgelist, as.numeric(MI) > percentile_99.99)
dim(percentile_99.99_table)
#[1] 22354     3

#Histogram of MI distribution 

png("MI_distribution_percentile99.99_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(percentile_99.99_table$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "MI distribution histogram",   # Título principal
     sub = "patients with AD, 99.99 percentile",   # Subtítulo
     xlab = "MI",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#write graph

graph_percentile99.99 <- graph_from_data_frame(percentile_99.99_table, directed = F)

#Graph info

summary(graph_percentile99.99)

#IGRAPH 02b3333 UN-- 3378 223532 --
#  + attr: name (v/c), MI (e/x)

#number of components

components(graph_percentile99.99)

#Degree by node 

degreee_by_node_graph_percentile99.99<- degree(graph_percentile99.99) %>% as.data.frame()

png("degreee_by_node_graph_percentile99.99_MI_degree_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(degreee_by_node_graph_percentile99.99$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Degree by node histogram",   # Título principal
     sub = "patients with AD, graph_percentile 99.99%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Coreness by node

coreness_by_node_graph_percentile99.99 <- coreness(graph_percentile99.99) %>% as.data.frame()

png("graph_percentile99.99_MI_coreness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_graph_percentile99.99$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Coreness by node histogram",   # Título principal
     sub = "patients with AD, graph percentile 99.99%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Plot graph, slow

png('graph_percentile99.99_graph_MI_0.853_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.png')
plot(graph_percentile99.99, edge.width = 2, edge.color = "black")
dev.off()

#write graphs in graphml

write_graph(graph_percentile99.99,
            file = 'graph_percentile99.99_graph_MI_0.711_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
            'graphml') #'edgelist.txt' / '.graphml'

# Percentile 99.9 --- --- 

percentile_99.9 <- quantile(as.numeric(full_edgelist$MI), 0.999)
#  99.9% 0.5497759 

percentile_99.9_table <- subset(full_edgelist, as.numeric(MI) > percentile_99.9)
dim(percentile_99.9_table)
# [1] 223532      3

#Histogram of MI distribution 

png("MI_distribution_percentile99.9_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(percentile_99.9_table$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "MI distribution histogram",   # Título principal
     sub = "patients with AD, 99.9 percentile",   # Subtítulo
     xlab = "MI",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#write graph

graph_percentile99.9 <- graph_from_data_frame(percentile_99.9_table, directed = F)

#Graph info

summary(graph_percentile99.9)

#IGRAPH ce2b296 UN-- 3378 223532 -- 
#+ attr: name (v/c), MI (e/x)

#number of components

components(graph_percentile99.9)

#Degree by node 

degreee_by_node_graph_percentile99.9<- degree(graph_percentile99.9) %>% as.data.frame()

png("degreee_by_node_graph_percentile99.9_MI_degree_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(degreee_by_node_graph_percentile99.9$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Degree by node histogram",   # Título principal
     sub = "patients with AD, graph_percentile 99.9%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Coreness by node

coreness_by_node_graph_percentile99.9 <- coreness(graph_percentile99.9) %>% as.data.frame()

png("graph_percentile99.9_MI_coreness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_graph_percentile99.9$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Coreness by node histogram",   # Título principal
     sub = "patients with AD, graph percentile 99.9%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Plot graph, slow

png('graph_percentile99.9_graph_MI_0.853_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.png')
plot(graph_percentile99.9, edge.width = 2, edge.color = "black")
dev.off()

#Write graph  

write_graph(graph_percentile99.9,
            file = 'graph_percentile99.9_graph_MI_0.549_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
            'graphml') #'edgelist.txt' / '.graphml'

# Percentile 99 --- --- 

percentile_99 <- quantile(as.numeric(full_edgelist$MI), 0.99)
# 99% 0.3748686 

percentile_99_table <- subset(full_edgelist, as.numeric(MI) > percentile_99)
dim(percentile_99_table)
#[1] 2235324       3

#Histogram of MI distribution 

png("MI_distribution_percentile99_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(percentile_99_table$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "MI distribution histogram",   # Título principal
     sub = "patients with AD, 99 percentile",   # Subtítulo
     xlab = "MI",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#write graph

graph_percentile99<- graph_from_data_frame(percentile_99_table, directed = F)

#Graph info

summary(graph_percentile99)

#IGRAPH 4db644c UN-- 7558 2235324 -- 
#  + attr: name (v/c), MI (e/x)

#number of components

components(graph_percentile99)

#Degree by node 

degreee_by_node_graph_percentile99<- degree(graph_percentile99) %>% as.data.frame()

png("degreee_by_node_graph_percentile99_MI_degree_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(degreee_by_node_graph_percentile99$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Degree by node histogram",   # Título principal
     sub = "patients with AD, graph_percentile 99%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Coreness by node

coreness_by_node_graph_percentile99 <- coreness(graph_percentile99) %>% as.data.frame()

png("graph_percentile99_MI_coreness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_graph_percentile99$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Coreness by node histogram",   # Título principal
     sub = "patients with AD, graph percentile 99%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Plot graph, slow

png('graph_percentile99graph_MI_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.png')
plot(graph_percentile99, edge.width = 2, edge.color = "black")
dev.off()

# Percentile 90 --- --- 

percentile_90 <- quantile(as.numeric(full_edgelist$MI), 0.9)
#     90%  of observations are above 0.1830432 

percentile_90_table <- subset(full_edgelist, as.numeric(MI) > percentile_90)
dim(percentile_90_table)
#[1] 22353240        3

#write graph

graph_percentile90 <- graph_from_data_frame(percentile_90_table, directed = F)

#analyze

components(graph_percentile90)

#write graphs --- ---

write_graph(graph_percentile90,
            file = '/datos/rosmap/cut_by_MI_one_component/graph_percentile90_MI_0.18_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
            'graphml') #'edgelist.txt' / '.graphml'

# Percentile 90 --- --- 

percentile_90 <- quantile(as.numeric(full_edgelist$MI), 0.9)
#     90%  of observations are above 0.1830432 

percentile_90_table <- subset(full_edgelist, as.numeric(MI) > percentile_90)
dim(percentile_90_table)
#[1] 22353240        3

#Histogram of MI distribution 

png("MI_distribution_percentile90_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(percentile_90_table$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "MI distribution histogram",   # Título principal
     sub = "patients with AD, 90 percentile",   # Subtítulo
     xlab = "MI",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#write graph

graph_percentile90 <- graph_from_data_frame(percentile_90_table, directed = F)

#Graph info

summary(graph_percentile90)

#IGRAPH 35dc62f UN-- 12435 22353240 -- 
#  + attr: name (v/c), MI (e/x)

#number of components

components(graph_percentile90)

#Degree by node 

degreee_by_node_graph_percentile90<- degree(graph_percentile90) %>% as.data.frame()

png("degreee_by_node_graph_percentile90_MI_degree_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(degreee_by_node_graph_percentile90$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Degree by node histogram",   # Título principal
     sub = "patients with AD, graph_percentile 90%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Coreness by node

coreness_by_node_graph_percentile90 <- coreness(graph_percentile90) %>% as.data.frame()

png("graph_percentile90_MI_coreness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_graph_percentile90$.), 
     col = "blue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Coreness by node histogram",   # Título principal
     sub = "patients with AD, graph percentile 90%",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Plot graph, slow

png('graph_percentile90graph_MI_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.png')
plot(graph_percentile90, edge.width = 2, edge.color = "black")
dev.off()

# Percentile 80 --- --- 

percentile_80 <- quantile(as.numeric(full_edgelist$MI), 0.8)
#     80%  0.1256503 

percentile_80_table <- subset(full_edgelist, as.numeric(MI) > percentile_80)
dim(percentile_80_table)
#[1] 44706480        3

# graph

graph_percentile80 <- graph_from_data_frame(percentile_80_table, directed = F)

#write graphs --- ---

write_graph(graph_percentile80,
            file = '/datos/rosmap/cut_by_MI_one_component/graph_percentile80_MI_0.126_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
            'graphml') #'edgelist.txt' / '.graphml'

#You can performm further analysis to graphs in cytoscape or with igraph