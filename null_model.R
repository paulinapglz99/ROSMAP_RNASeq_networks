#null_model.R
#This script builds two null models for each graph, by randomizing the networks while preserving the degree distribution. 
#Specifically, we rewired the edges of the original networks 1,000 times using a degree-preserving randomization algorithm. 
pacman::p_load("igraph", 
        "ggplot2")

#Declare functions ----

# Function to calculate various metrics
calculate_metrics <- function(graph) {
  list(
    clustering = transitivity(graph, type = "global"),  # Global clustering coefficient
    assortativity = assortativity_degree(graph),
    avg_path_length = mean_distance(graph, directed = FALSE),  # Average path length
    modularity = modularity(cluster_infomap(graph))  # Modularity of the network, using infomap
  )
}

#Function to randomize networks ----

randomize_graph <- function(graph, num_randomizations, niter = 1000) {
  random_graphs <- vector("list", num_randomizations)
  for (i in 1:num_randomizations) {
    random_graphs[[i]] <- rewire(graph, with = keeping_degseq(niter = niter))
  }
  return(random_graphs)
}

#Function to calculate metrics from ranzomized networks

calculate_metrics_random_graphs <- function(graph_list) {
  lapply(graph_list, calculate_metrics)
}

#Function to plot any metric -----

plot_metric_distribution <- function(random_metrics, observed_metric, metric_name, graph_type, ci_lower, ci_upper) {
  random_df <- data.frame(value = random_metrics, type = paste0("Random (", graph_type, ")")) # DataFrame de valores randomizados
  observed_df <- data.frame(value = observed_metric, type = paste0("Observed (", graph_type, ")"))   # Valor observado
  combined_df <- rbind(random_df, observed_df)   # Combina los DataFrames
  
  # Plot
  ggplot(random_df, aes(x = value)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.001, fill = "skyblue4", color = "black", alpha = 0.7) +
    geom_density(color = "blue", lwd = 1) +
    geom_vline(aes(xintercept = observed_metric), color = "red", linetype = "dashed", size = 0.8) +
    geom_vline(aes(xintercept = ci_lower), color = "darkgreen", linetype = "dotted", size = 0.8) +  # Línea para el límite inferior del IC
    geom_vline(aes(xintercept = ci_upper), color = "darkgreen", linetype = "dotted", size = 0.8) +  # Línea para el límite superior del IC
    labs(title = paste(metric_name, "(", graph_type, ")"), 
         x = metric_name, 
         y = "Density") + 
    theme_minimal()
}

#Read data ----

graphAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

graphnoAD <- read_graph(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/counts_by_NIA_Reagan/graphs_NIA_Reagan/ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99.graphml', format = 'graphml')

# Calculate metrics for your original networks
metrics_AD <- calculate_metrics(graphAD)
metrics_noAD <- calculate_metrics(graphnoAD)

# Number of random networks to generate
num_randomizations <- 1000

#Randomize

random_graphs_AD <- randomize_graph(graphAD, num_randomizations)
random_graphs_noAD <- randomize_graph(graphnoAD, num_randomizations)

# Calculate metrics for randomized graphs
metrics_AD_random <- calculate_metrics_random_graphs(random_graphs_AD)
metrics_noAD_random <- calculate_metrics_random_graphs(random_graphs_noAD)

# Extract metrics -----

# 1. Extract clustering coefficient
clustering_AD_random <- sapply(metrics_AD_random, function(x) x$clustering)
clustering_noAD_random <- sapply(metrics_noAD_random, function(x) x$clustering)

# Extract average path length 

asort_AD_random <- sapply(metrics_AD_random, function(x) x$assortativity)
asort_noAD_random <- sapply(metrics_noAD_random, function(x) x$assortativity)

# Extract modularity for the random networks (AD and noAD)
modularity_AD_random <- sapply(metrics_AD_random, function(x) x$modularity)
modularity_noAD_random <- sapply(metrics_noAD_random, function(x) x$modularity)

# Confidence invtervals

#If calculate the 2.5 and 97.5 percentile for the randomized network metrics, we generate a 95% confidence interval. 
#If the observed metrics fall outside this interval, this will indicate that the observed metrics are not likely to be the product of chance.

# Intervalo de confianza al 95% para clustering coefficient
clustering_ci_AD <- quantile(clustering_AD_random, probs = c(0.025, 0.975))
clustering_ci_noAD <- quantile(clustering_noAD_random, probs = c(0.025, 0.975))

# Verificar si el valor observado cae fuera del intervalo de confianza
print(clustering_ci_AD)
print(metrics_AD$clustering)
print(clustering_ci_noAD)
print(metrics_noAD$clustering)

# Intervalo de confianza al 95% para asortatividad

aso_ci_AD <- quantile(asort_AD_random, probs = c(0.025, 0.975))
aso_ci_noAD <- quantile(asort_noAD_random, probs = c(0.025, 0.975))

# Verificar si el valor observado cae fuera del intervalo de confianza
print(aso_ci_AD)
print(metrics_AD$assortativity)
print(aso_ci_noAD)
print(metrics_noAD$assortativity)

# Intervalo de confianza al 95% para modularidad

modu_ci_AD <- quantile(modularity_AD_random, probs = c(0.025, 0.975))
modu_ci_noAD <- quantile(modularity_noAD_random, probs = c(0.025, 0.975))

# Verificar si el valor observado cae fuera del intervalo de confianza
print(modu_ci_AD)
print(metrics_AD$modularity)
print(modu_ci_noAD)
print(metrics_noAD$modularity)

# Plot clustering, modularity and assortativity

#Modularity

modu_AD <- plot_metric_distribution(
  random_metrics = modularity_AD_random,
  observed_metric = metrics_AD$modularity,
  metric_name = "Modularity",
  graph_type = "AD",
  ci_lower = modu_ci_AD[1],
  ci_upper = modu_ci_AD[2]
)
  
modu_noAD <-  plot_metric_distribution(
  random_metrics = modularity_noAD_random,
  observed_metric = metrics_noAD$modularity,
  metric_name = "Modularity",
  graph_type = "control",
  ci_lower = modu_ci_noAD[1],
  ci_upper = modu_ci_noAD[2]
)

# Clustering
#AD
clust_AD <- plot_metric_distribution(
  random_metrics = clustering_AD_random,
  observed_metric = metrics_AD$clustering,
  metric_name = "Clustering Coefficient",
  graph_type = "AD",
  ci_lower = clustering_ci_AD[1],
  ci_upper = clustering_ci_AD[2]
)
#noAD
clust_noAD <- plot_metric_distribution(
  random_metrics = clustering_noAD_random,
  observed_metric = metrics_noAD$clustering,
  metric_name = "Clustering Coefficient",
  graph_type = "control",
  ci_lower = clustering_ci_noAD[1],
  ci_upper = clustering_ci_noAD[2]
)

# Asortativity
#AD
asort_AD <- plot_metric_distribution(
  random_metrics = asort_AD_random,
  observed_metric = metrics_AD$assortativity,
  metric_name = "Assortativity",
  graph_type = "AD",
  ci_lower = aso_ci_AD[1],
  ci_upper = aso_ci_AD[2]
)
#noAD
asort_noAD <- plot_metric_distribution(
  random_metrics = asort_noAD_random,
  observed_metric = metrics_AD$assortativity,
  metric_name = "Assortativity",
  graph_type = "control",
  ci_lower = aso_ci_noAD[1],
  ci_upper = aso_ci_noAD[2]
)

# Función para calcular gamma a partir de la distribución de grados ----

randomize_graph_nodis <- function(graph, num_randomizations, niter = 1000) {
  random_graphs <- vector("list", num_randomizations)
  for (i in 1:num_randomizations) {
    random_graphs[[i]] <- rewire(graph, with = each_edge(
      0.8,
      loops = FALSE,
      multiple = TRUE
    ) )
  }
  return(random_graphs)
}

#Randomize sin mantener la degree dis

random_graphs_AD <- randomize_graph_nodis(graphAD, num_randomizations)
random_graphs_noAD <- randomize_graph_nodis(graphnoAD, num_randomizations)

# Function to calculate gamma (scaling exponent) ----
calculate_gamma <- function(graph) {
  degree_distribution <- degree(graph)  # Get degree of each node
  degree_freq <- as.data.frame(table(degree_distribution))  # Frequency table of degrees
  colnames(degree_freq) <- c("degree", "Freq")
  degree_freq$degree <- as.numeric(as.character(degree_freq$degree))
  
  # Calculate probability and log values for log-log regression
  degree_freq$Prob <- degree_freq$Freq / sum(degree_freq$Freq)
  degree_freq <- degree_freq[degree_freq$degree > 0, ]  # Exclude zero degrees
  degree_freq$log_degree <- log(degree_freq$degree)
  degree_freq$log_Prob <- log(degree_freq$Prob)
  
  # Fit linear model to log-log degree distribution to get gamma
  fit <- lm(log_Prob ~ log_degree, data = degree_freq)
  gamma <- -coef(fit)["log_degree"]
  return(gamma)
}

# Calculate observed gamma values for AD and noAD networks ----
gamma_AD_observed <- calculate_gamma(graphAD)
gamma_noAD_observed <- calculate_gamma(graphnoAD)

# Randomize networks and calculate gamma for each random network ----
gamma_random_AD <- sapply(random_graphs_AD, calculate_gamma)
gamma_random_noAD <- sapply(random_graphs_noAD, calculate_gamma)


# Calculate 95% confidence intervals for gamma ----
gamma_ci_AD <- quantile(gamma_random_AD, probs = c(0.025, 0.975))
gamma_ci_noAD <- quantile(gamma_random_noAD, probs = c(0.025, 0.975))

# Check if observed gamma falls outside the confidence interval
cat("95% CI for Gamma (AD):", gamma_ci_AD, "\n")
cat("Observed Gamma (AD):", gamma_AD_observed, "\n")
cat("95% CI for Gamma (control):", gamma_ci_noAD, "\n")
cat("Observed Gamma (control):", gamma_noAD_observed, "\n")

# Plot gamma distribution ----
plot_gamma_distribution <- function(random_gammas, observed_gamma, graph_type, gamma_ci) {
  random_df <- data.frame(value = random_gammas, type = paste0("Random (", graph_type, ")"))
  observed_df <- data.frame(value = observed_gamma, type = paste0("Observed (", graph_type, ")"))
  combined_df <- rbind(random_df, observed_df)
  
  # Plot
  ggplot(random_df, aes(x = value)) +
    geom_histogram(aes(y = ..density..), binwidth = 0.01, fill = "skyblue4", color = "black", alpha = 0.7) +
    geom_density(color = "blue", lwd = 1) +
    geom_vline(aes(xintercept = observed_gamma), color = "red", linetype = "dashed", size = 0.8) +
    # Add the confidence interval lines
    geom_vline(aes(xintercept = gamma_ci[1]), color = "green", linetype = "dotted", size = 0.8) +
    geom_vline(aes(xintercept = gamma_ci[2]), color = "green", linetype = "dotted", size = 0.8) +
    labs(title = paste("Gamma (", graph_type, ")"), 
         x = "Gamma", 
         y = "Density") + 
    theme_minimal()
}

# Plot gamma distribution for AD and noAD
gamma_AD_plot <- plot_gamma_distribution(random_gammas = gamma_random_AD, observed_gamma = gamma_AD_observed, 
                        graph_type = "AD", gamma_ci = gamma_ci_AD)
gamma_AD_plot
gamma_noAD_plot <- plot_gamma_distribution(random_gammas = gamma_random_noAD, observed_gamma = gamma_noAD_observed, 
                                           graph_type = "control", gamma_ci = gamma_ci_noAD)
gamma_noAD_plot

#Grid null model

grid_nullm <- grid.arrange(clust_AD, clust_noAD, 
                           asort_AD, asort_noAD, 
                           modu_AD, modu_noAD,
                           gamma_AD_plot, gamma_noAD_plot,
                           ncol = 2)

ggsave(
  "null_model_grid.jpg",
  plot = grid_nullm,
  device = "jpg",
  width = 15,
  height = 10,
  units = "in",
  dpi = 300
)

END