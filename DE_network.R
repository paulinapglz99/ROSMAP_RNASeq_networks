#
#DE_network.R
#This script takes differentially expressed genes and a co-expression network and makes a sub_network with DE genes

#Libraries --- ---

pacman::p_load("igraph", 
               "dplyr",
               "ggplot")

#Get data --- ---

DE_genes <- vroom::vroom(file = "")

graph <- read_graph(file = "", format = "graphml")

#Create sub-graph --- ---

#Save data --- --- 