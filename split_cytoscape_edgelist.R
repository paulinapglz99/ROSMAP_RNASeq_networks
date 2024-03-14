#
#split_cytoscape_edgelist.R
#this script helps to condition the edgelists coming from cytoscape

#Libraries --- --- 

library(dplyr)
library(stringr)

#Get data --- ---

edgelist <- vroom::vroom(file = '~/redesROSMAP/percentile99.99_ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_edgelist.csv')

# Separate the first column into twoedgelist <- edgelist %>%
  
separate(1, into = c("v1", "v2"), sep = " \\(-\\) ")

# Remove the letters 'n' from columns 'v1' and 'v2'.

edgelist$v1 <- gsub("n", "", edgelist$v1)
edgelist$v2 <- gsub("n", "", edgelist$v2)

#Save table

vroom::vroom_write(edgelist, file = '~/redesROSMAP/percentile99.99_ROSMAP_RNAseq_MutualInfo_noAD_NIA_Reagan_dicho_edgelist.csv')

#END