#transitivity_by_module.R


# The AD network exhibits higher modularity (Q value), but there is no significant difference in the global degree distribution (K-S test).
# The study should explore whether the modularity differences arise from local connection reorganization (e.g., increased hub genes within specific modules) rather than global topological changes.

#Libraries --- ---

pacman::p_load("igraph", 
               "tidyverse", 
               "cowplot") 
            
#Functions --- ---

#Calculate modularity Q by module
compute_modularity_per_module <- function(graph, nodes_by_community) {
  total_edges <- gsize(graph) #Total edges in the network
  module_Q_values <- list()
  
  for (module_name in names(nodes_by_community)) {
    module_nodes <- nodes_by_community[[module_name]]  #Nodes in the module
    subgraph <- induced_subgraph(graph, vids = module_nodes)  #Extract subgraphs
    subgraph_edges <- gsize(subgraph)  #Edges in the module
    
    #Calculare modularities by module
    Q_s <- modularity(graph, membership = V(graph)$community, weights = E(graph)$weight) *
      (subgraph_edges / total_edges)
    
    module_Q_values[[module_name]] <- Q_s
  }
  
  return(module_Q_values)
}

#Count hubs by community
count_hubs_in_communities <- function(nodes_by_community, hubs) {
  sapply(nodes_by_community, 
         function(nodes) sum(nodes %in% hubs))
}

#Get networks --- ---

graphAD <- read_graph(file =  'graphADnodes_membership.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = 'graphnoADnodes_membership.graphml',
                        format = 'graphml')

#Save graphs in a list

graphLists <- list(AD = graphAD,
                   Control = graphnoAD)

#Extract list of nodes by community for each graph
nodes_by_community_list <- lapply(graphLists, function(graph) {
  split(V(graph)$name, V(graph)$community)
})

#Count how many modules each network has --- ---

sapply(nodes_by_community_list, length)
# AD Control 
# 68      71 

#Global modularity of networks --- ---

sapply(graphLists, function(g) modularity(g, membership = V(g)$community))
# 
# AD   Control 
# 0.2825230 0.2027528 

#Modularity Q by module --- ---

modularity_per_module <- mapply(function(graph, name) {
  compute_modularity_per_module(graph, nodes_by_community_list[[name]])
}, graphLists, names(graphLists), SIMPLIFY = FALSE)

#As dataframe --- ---

mod_list <- lapply(names(modularity_per_module), function(graph_name) {
  df <- enframe(modularity_per_module[[graph_name]], name = "Module", value = "Modularity")
  df$Graph <- graph_name  # Agregar columna de red
  return(df)
})

#as data frame
mod_df <- do.call(bind_rows, mod_list) %>%
  mutate(Module = factor(Module,levels = sort(unique(as.numeric(Module)))),  #Order modules
         Modularity = as.numeric(Modularity)) %>%
  mutate(modularity_log = log10(Modularity))

mod_df_AD <- mod_df %>% filter(Graph == "AD")
mod_df_control <- mod_df %>% filter(Graph == "Control")

#Hubs by module --- ---

hubs_nodes <- list()

#Load hubs
for (name in c("graphAD", "graphnoAD")) {
  filename_txt <- paste0("hubs_", name, ".txt")
  hubs_nodes[[name]] <- readLines(filename_txt)
}


#Count how many hubs there are by community --- ---
hubs_by_community <- mapply(count_hubs_in_communities, 
                                  nodes_by_community_list, 
                                  hubs_nodes, 
                                  SIMPLIFY = FALSE)

hubs_by_community_AD.df <- data.frame(
  Module = names(hubs_by_community$AD),  # Nombres de las comunidades
  Hubs_Count = unlist(hubs_by_community$AD)  # Conteo de hubs en cada comunidad
)

hubs_by_community_control.df <- data.frame(
  Module = names(hubs_by_community$AD),  # Nombres de las comunidades
  Hubs_Count = unlist(hubs_by_community$AD)  # Conteo de hubs en cada comunidad
)

#Merge with modularity by module

mod_df_AD_hubs<- mod_df_AD %>% 
  left_join(hubs_by_community_AD.df, by = "Module") %>% 
  mutate(Module = factor(Module,levels = sort(unique(as.numeric(Module)))),  #Order modules
          Modularity = as.numeric(Modularity))


mod_df_contr_hubs<- mod_df_control %>% 
  left_join(hubs_by_community_control.df, by = "Module") %>% 
  mutate(Module = factor(Module,levels = sort(unique(as.numeric(Module))), ),  #Order modules
         Modularity = as.numeric(Modularity))

mod_df_contr_hubs[is.na(mod_df_contr_hubs)] <- 0

#PLOT --- ---
# 
# #plot modularity by module 
# 
# ggplot(mod_df, aes(x = Module, y = modularity_log, fill = Graph)) +
#   geom_col(position = "dodge") +
#   theme_minimal() +
#   labs(title = "Modularity Q by module",
#        x = "Module",
#        y = "Modularity Q",
#        fill = "Red") +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Plot heatmap
#Q_module.ht <-
# ggplot(mod_df, aes(x = Module, y = Graph, fill = modularity_log)) +
# geom_tile() +
# scale_fill_gradient(low = "blue", high = "red") +
# theme_minimal() +
# labs(title = "Modularity Q by module",
#      x = "Module",
#      y = "Graph",
#      fill = "logQ")

#Control

Q_module_control.ht <- ggplot(mod_df_contr_hubs, 
                             aes(x = Module, y = Graph, fill = modularity_log)) +
  geom_tile() +
 # geom_text(aes(label = Hubs_Count), vjust = -14, size = 3) + # Añade texto debajo
  #scale_fill_gradient(low = "blue", high = "red") +
  geom_label(aes(label = Hubs_Count), 
            vjust = 1.5, 
            size = 3, 
            fill = "white", 
            color = "black", 
            label.size = 0.2) +
  scale_fill_viridis(option="plasma")+
  theme_minimal() +
  labs(x = " ",
       y = " ",
       fill = "logQ")

#only AD

Q_module_AD.ht <- ggplot(mod_df_AD_hubs, 
                         aes(x = Module, y = Graph, fill = modularity_log)) +
  geom_tile() +
  #geom_text(aes(label = Hubs_Count), vjust = -14, size = 3) + # Añade texto debajo
  # scale_fill_gradient(low = "blue", high = "red") +
  geom_label(aes(label = Hubs_Count), 
             vjust = 1.5, 
             size = 3, 
             fill = "white", 
             color = "black", 
             label.size = 0.2) +
  scale_fill_viridis(option="plasma")+
  theme_minimal() +
  labs(x = "Module",
       y = " ",
       fill = "logQ")


Q_grid <- cowplot::plot_grid(Q_module_control.ht, Q_module_AD.ht , labels = c('a', 'b'), ncol = 2)

# ggsave("Q_grid.jpg", plot = Q_grid,
#        device = "jpg",
#        width = 6.5,
#        height = 15.5,
#        units = "in",
#        dpi = 300)

#I want a vertical version

Q_module_control.ht <- ggplot(mod_df_contr_hubs, 
                              aes(x = Graph , y = Module, fill = modularity_log)) +
  geom_tile() +
  geom_label(aes(label = Hubs_Count), 
             #vjust = -1.5, 
             size = 3, 
             fill = "white", 
             color = "black", 
             label.size = 0.2) +
  scale_fill_viridis(option="plasma")+
  scale_y_discrete(limits = rev(unique(mod_df_contr_hubs$Module))) +  #Invert plot
  theme_minimal() +
  labs(x = " ",
       y = " ",
       fill = "logQ")

#only AD

Q_module_AD.ht <- ggplot(mod_df_AD_hubs, 
                         aes(x =Graph , y = Module, fill = modularity_log)) +
  geom_tile() +
  geom_label(aes(label = Hubs_Count), 
             #vjust = -1.5, 
             size = 3, 
             fill = "white", 
             color = "black", 
             label.size = 0.2) +
  scale_fill_viridis(option="plasma")+
  scale_y_discrete(limits = rev(unique(mod_df_AD_hubs$Module))) +  #Invert plot
  theme_minimal() +
  labs(x = "",
       y = " ",
       fill = "logQ")

Q_grid <- cowplot::plot_grid(Q_module_control.ht, Q_module_AD.ht , labels = c('a', 'b'), ncol = 2)
Q_grid

# ggsave("Q_grid.jpg", plot = Q_grid,
#        device = "jpg",
#        width = 6.5,
#        height = 15.5,
#        units = "in",
#        dpi = 300)

#END