#Comparisons of modularity algorithms 

#Libraries --- ---

pacman::p_load("igraph", 
               "tidyverse", 
               "cowplot", 
               "leiden") 

#Functions --- ---

#Function to apply community detection algorithms
detect_communities <- function(graph) {
  list(
    Louvain = cluster_louvain(graph),
    Infomap = cluster_infomap(graph),
    Leiden = cluster_leiden(graph)
  )
}

#Calculate Newmann Modularity and count communities
get_metrics <- function(graph, community_obj) {
  tibble(
    Algorithm = names(community_obj),
    Modularity = sapply(community_obj, function(x) modularity(graph, membership(x))),  # Corrected!
    Communities = sapply(community_obj, function(x) length(unique(membership(x))))
  )
}

#Get networks --- ---

graphAD <- read_graph(file =  'graphADnodes_membership.graphml',
                      format = 'graphml')

graphnoAD <- read_graph(file = 'graphnoADnodes_membership.graphml',
                        format = 'graphml')

#Save graphs in a list

graphLists <- list(AD = graphAD,
                   Control = graphnoAD)

#Apply community detection algorithms

communities <- lapply(graphLists, detect_communities)

#Get metrics of comparison

results_AD <- get_metrics(communities$AD) %>% mutate(Network = "AD")
results_Control <- get_metrics(communities$Control) %>% mutate(Network = "Control")

# Combine results
final_results <- bind_rows(results_AD, results_Control)

modularity()

