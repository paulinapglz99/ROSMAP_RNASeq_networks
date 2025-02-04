#Libraries --- ---
#install.packages("svglite")
pacman::p_load("igraph", 
               "ggraph",
               "tidyverse", 
               "gridExtra", 
               "svglite", 
               "tidygraph")

#Set seed --- --- 

set.seed(10)

#Get data --- --- 

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml', 
                      format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml', 
                        format = 'graphml')

#Save graphs in a list

graphLists <- list(graphAD = graphAD,
                   graphnoAD = graphnoAD)

#Degree distributions of all graphs in my list
degree_distributions <- sapply(X = graphLists, FUN = degree)

#Build degree distribution dataframe

ADdegree <- degree_distributions[["graphAD"]]
ADdegree.df <- data.frame(gene = names(ADdegree), degree = ADdegree)
ADdegree_freq <- table(ADdegree.df$degree) %>% as.data.frame()
colnames(ADdegree_freq) <- c("degree", "Freq")
ADdegree_freq$degree <- as.numeric(as.character(ADdegree_freq$degree))
ADdegree_freq$Prob <- ADdegree_freq$Freq / sum(ADdegree_freq$Freq) 
ADdegree_freq <- ADdegree_freq[order(ADdegree_freq$degree, decreasing = F), ]
ADdegree_freq$CumulativeDegree <- cumsum(ADdegree_freq$Freq)
# Logarithm of the cumulative sum
ADdegree_freq$logCumulativeDegree <- log10(ADdegree_freq$CumulativeDegree)
ADdegree_freq$dx <- "AD"
ADdegree_freq$degree <- as.numeric(ADdegree_freq$degree)

#noAD

noADdegree <- degree_distributions[["graphnoAD"]]
noADdegree.df <- data.frame(gene = names(noADdegree), degree = noADdegree)
noADdegree_freq <- table(noADdegree.df$degree) %>% as.data.frame()
colnames(noADdegree_freq) <- c("degree", "Freq")
noADdegree_freq$degree <- as.numeric(as.character(noADdegree_freq$degree))
noADdegree_freq$Prob <- noADdegree_freq$Freq / sum(noADdegree_freq$Freq)
noADdegree_freq <- noADdegree_freq[order(noADdegree_freq$degree, decreasing = F), ]
noADdegree_freq$CumulativeDegree <- cumsum(noADdegree_freq$Freq)
# Logarithm of the cumulative sum
noADdegree_freq$logCumulativeDegree <- log10(noADdegree_freq$CumulativeDegree)
noADdegree_freq$dx<- "no AD"
noADdegree_freq$degree <- as.numeric(noADdegree_freq$degree)

# Find the max value of "Freq" in both networks
max_freq_AD <- max(ADdegree_freq$Freq)
max_freq_noAD <- max(noADdegree_freq$Freq)

# Find max degree value in both graphs 
max_degree_AD <- max(as.numeric(as.character(ADdegree_freq$degree)))
max_degree_noAD <- max(as.numeric(as.character(noADdegree_freq$degree)))

# Define limits for make histograms comparable
max_y <- max(max_freq_AD, max_freq_noAD)
max_x <- max(max_degree_AD, max_degree_noAD)

# Bind rows of both degre dis 
dis <- rbind(ADdegree_freq, noADdegree_freq)

#Plot --- ---

degree_distr.p <- ggplot(dis, aes(x = degree, y = logCumulativeDegree, color = dx)) +
  geom_point(size = 1.5) +
  geom_line(size = 1.5) +
  scale_x_continuous(limits = c(min(dis$degree), max(dis$degree)), breaks = seq(min(dis$degree), max(dis$degree), by = 5)) +  # Ajusta el eje X
  labs(title = "Degree Distribution",
       x = "Degree",
       y = "log CumulativeDegree") +
  scale_color_manual(values = c("AD" = "#A3333D", "no AD" = "#477998")) + # Custom colors for AD and no AD
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # Rotate X-axis labels at 45 degrees
  theme_minimal()

degree_distr.p

#Save plot

ggsave(filename = "bothdx_cumulative_degree_distributions_coexpression_NIAReagan_histogram.png",
       plot = degree_distr.p,
       width = 20,
       height = 10,
       units = "in",
       dpi = 300)

#

degree_distr_x.p <- ggplot(dis, aes(x = as.numeric(degree), y = as.numeric(CumulativeDegree), color = dx)) +
  geom_point(size = 1.5) +
  geom_line(size = 1.5) +
  scale_x_continuous(limits = c(min(dis$degree), max(dis$degree)), breaks = seq(min(dis$degree), max(dis$degree), by = 5)) +  # Ajusta el eje X
  labs(title = "Degree Distribution",
       x = "Degree",
       y = "Cumulative Degree") +
  #scale_x_continuous(breaks = seq(1, max(dis$degree), by =4)) + # Adjust the x-axis to show only integer degrees
  scale_color_manual(values = c("AD" = "#A3333D", "no AD" = "#477998")) + # Custom colors for AD and no AD
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # Rotate X-axis labels at 45 degrees
  theme_minimal()

degree_distr_x.p

ggsave(filename = "bothdx_nolog_cumulative_degree_distributions_coexpression_NIAReagan_histogram.jpg",
       plot = degree_distr_x.p,
       width = 20,
       height = 10,
       units = "in",
       dpi = 300)

#END