#
#4.my_network_cut_MI_distribution_heuristic.R

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

########--- 1st HEURISTIC: cutoff according to MI distribution --- ########

#Libraries  --- --- 

pacman::p_load('tidyverse', 
               'igraph', 
               'ggplot2')

#Get data --- ---

full_edgelist <- vroom::vroom(file = '/datos/rosmap/coexpre_matrix/ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho_edgelist.tsv')

#Plot MI distribution

#This histogram is necessary to decide on the following mutual information cutoffs 
png("full_net_MI_histogram_ROSMAP_RNA_seq.png")
hist(as.numeric(full_edgelist$MI), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Mutual information histogram",   # Título principal
     sub = "patients with AD, no MI cut",   # Subtítulo
     xlab = "Mutual Information",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#If we look at the histogram, we will notice that there is a point in the MI distribution where we see a dip, 
#which is called the 'tail of the histogram'. The first approach is to make a cut in this tail. 

matrix_MI_tailcut<- full_edgelist %>% 
  filter(MI >= 0.5)   
dim(matrix_MI_tailcut)
#[1] 440646      3

#Build graph

tailcut_graph <-  graph_from_data_frame(matrix_MI_tailcut, directed = F)

#Graph info

summary(tailcut_graph)

#Components

components(tailcut_graph)

#Degree by node 

degreee_by_node_tailcut_g <- degree(min_sin_comp_graph) %>% as.data.frame()

png("tailcut_graph_MI_0.5_degree_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(degreee_by_node_tailcut_g$.), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Degree by node histogram",   # Título principal
     sub = "patients with AD, tailcut graph",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()

#Coreness by node

coreness_by_node_min_sin_g <- coreness(min_sin_comp_graph) %>% as.data.frame()


png("min_sin_comp_MI_coreness_by_node_histogram_ROSMAP_AD_RNA_seq.png")
hist(as.numeric(coreness_by_node_min_sin_g$.), 
     col = "skyblue",      # Color de las barras
     border = "white",     # Color del borde de las barras
     main = "Coreness by node histogram",   # Título principal
     sub = "patients with AD, minimal single component",   # Subtítulo
     xlab = "nodes",   # Etiqueta del eje x
     ylab = "freq",   # Etiqueta del eje y
     labels = FALSE   # Deshabilitar la notación científica en el eje y
)
dev.off()


#edge betweenneess by node

betweeness_by_node_min_sin_g <-betweenness(min_sin_comp_graph, directed = F) %>% as.data.frame()


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

#Save graph

#write graphs in graphml --- ---

#write_graph(tailcut_graph,
#            file = '/datos/rosmap/cut_by_MI_one_component/tailcut_graph_MI_0.5_ROSMAP_RNAseq_MutualInfo_AD_NIA_Reagan_dicho.graphml', 
#            'graphml') #'edgelist.txt' / '.graphml'

#END