#Libraries --- ---
#install.packages("svglite")
pacman::p_load("igraph", 
               "ggraph",
               "tidyverse", 
               "gridExtra", 
               "svglite", 
               "tidygraph")

#Functions --- ---

#Calculate degree distribution and freq table
calculate_degree_distribution <- function(graph, label) {
  degree_data <- degree(graph)
  degree_df <- data.frame(gene = names(degree_data), degree = degree_data)
  degree_freq <- table(degree_df$degree) %>% as.data.frame()
  colnames(degree_freq) <- c("degree", "Freq")
  degree_freq$degree <- as.numeric(as.character(degree_freq$degree))
  degree_freq$Prob <- degree_freq$Freq / sum(degree_freq$Freq)
  degree_freq <- degree_freq[order(degree_freq$degree, decreasing = FALSE), ]
  degree_freq$CumulativeDegree <- cumsum(degree_freq$Freq)
  degree_freq$logCumulativeDegree <- log10(degree_freq$CumulativeDegree)
  degree_freq$dx <- label
  degree_freq$CDF <- degree_freq$CumulativeDegree / sum(degree_freq$Freq)
  return(degree_freq)
}

#Set seed --- --- 

set.seed(10)

#Get data --- --- 

graphAD <- read_graph(file = 'ROSMAP_RNAseq_DLPFC_AD_MutualInfograph_percentile99.99_trad.graphml', format = 'graphml')

graphnoAD <- read_graph(file = 'ROSMAP_RNAseq_DLPFC_noAD_MutualInfograph_percentile99.99_trad.graphml', format = 'graphml')


#Save graphs in a list

graphLists <- list(graphAD = graphAD,
                   graphnoAD = graphnoAD)

#Centrality of networks --- ---

degree_centralizations <- lapply(graphLists, 
                                 function(g) centr_degree(g, mode = "all", normalized = TRUE)$centralization)
degree_centralizations

#Network topological comparison --- --- 

#Table of degree distribution

#Degree distributions of all graphs in my list
degree_distributions <- sapply(X = graphLists, FUN = degree)

#Build degree distribution dataframe

ADdegree_freq <- calculate_degree_distribution(graphAD, "AD")

#noAD
noADdegree_freq <- calculate_degree_distribution(graphnoAD, "control")

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

###### log-log graph ######

log_log_AD <- ggplot(ADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
  scale_x_log10() +  # Escala logarítmica en el eje x
  scale_y_log10() +  # Escala logarítmica en el eje y
  labs(x = expression("log<k>"), y = expression(logp(k)), 
       title = "", 
       subtitle = "AD network \ndegree distribution") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", 
               size = 1.2, formula = y ~ x) +  # Ajuste lineal con estilo modificado
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )


log_log_noAD <-  ggplot(noADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
 scale_x_log10() + 
 scale_y_log10() + 
  labs(x = expression("log<k>"), y = expression(logp(k)), 
       title = "", 
       subtitle = "Control network \ndegree distribution") +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", 
               size = 1.2, formula = y ~ x) +  # Ajuste lineal con estilo modificado
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  
    axis.title = element_text(size = 14),  
    axis.text = element_text(size = 12),   
    panel.grid = element_line(size = 0.5, color = "grey80") 
  )

log_grid <- grid.arrange(log_log_AD, log_log_noAD, ncol =2)

#Find gammas

#Logarithms of the values of k and P(k) are taken.
ADdegree_freq$log_degree <- log(ADdegree_freq$degree)
ADdegree_freq$log_Prob <- log(ADdegree_freq$Prob)

#Log-log linear regression fitting
ajuste_AD <- lm(log_Prob ~ log_degree, data = ADdegree_freq)

#Results
summary(ajuste_AD)

# Calcular el valor de gamma (la pendiente negativa de la regresión)
gamma_AD <- -coef(ajuste_AD)["log_degree"]
cat("El valor de gamma para AD es:", gamma_AD, "\n")

#control

noADdegree_freq$log_degree <- log(noADdegree_freq$degree)
noADdegree_freq$log_Prob <- log(noADdegree_freq$Prob)

# Ajuste de regresión lineal en escala log-log
ajuste_noAD <- lm(log_Prob ~ log_degree, data = noADdegree_freq)

# Mostrar los resultados del ajuste
summary(ajuste_noAD)

# Calcular el valor de gamma (la pendiente negativa de la regresión)
gamma_noAD <- -coef(ajuste_noAD)["log_degree"]
cat("El valor de gamma noAD es:", gamma_noAD, "\n")

# ggsave(filename = "degree_dis_loglog.jpg",
#        plot = log_grid,
#        width = 25,
#        height = 25,
#        units = "cm",
#        dpi = 300,
# )

#Superponed

# Crear el gráfico combinado
log_log_combined <- ggplot(dis, aes(x = degree, y = Prob, color = dx)) +
  geom_point(size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
  scale_x_log10() +  # Escala logarítmica en el eje x
  scale_y_log10() +  # Escala logarítmica en el eje y
  labs(x = expression("log<k>"), y = expression(logp(k)), 
       title = " ") +
  geom_smooth(method = "lm", se = FALSE, linetype = "longdash", 
              size = 1.2, aes(color = dx), formula = y ~ x) +  # Ajuste lineal con colores diferenciados
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  annotate("text", 
           x = min(dis$degree) * 1.1,           # Ajusta cerca del límite inferior del eje x
           y = min(dis$Prob) * 1.1,             # Ajusta cerca del límite inferior del eje y
           label = "γcontrol = 0.7584 \n γAD = 0.7908", 
           color = "black", 
           size = 5, 
           vjust = 1,                # Ajusta la posición vertical
           hjust = 0)  +              # Ajusta la posición horizontal
  theme(
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  ) +
  scale_color_manual(values = c("AD" = "blue", "Control" = "red")) +
  theme(legend.position = "none")
log_log_combined

# ggsave(filename = "combined_degree_dis_loglog.jpg",
#        plot = log_log_combined,
#        width = 25,
#        height = 25,
#        units = "cm",
#        dpi = 300,
# )

#Statistical significance test
#Are my distributions significantly different?

#Create cumulative distribution functions
AD_CDF <- approxfun(ADdegree_freq$degree, ADdegree_freq$CDF, method = "linear", yleft = 0, yright = 1)
noAD_CDF <- approxfun(noADdegree_freq$degree, noADdegree_freq$CDF, method = "linear", yleft = 0, yright = 1)
ks_result <- ks.test(AD_CDF(ADdegree_freq$degree), noAD_CDF(noADdegree_freq$degree))

# Asymptotic two-sample Kolmogorov-Smirnov test
# 
# data:  AD_CDF(ADdegree_freq$degree) and noAD_CDF(noADdegree_freq$degree)
# D = 0.037916, p-value = 0.9971
# alternative hypothesis: two-sided

#. Con este valor p, no podemos concluir que las distribuciones son significativamente diferentes.

#Cummulative degree

degree_distr_x.p <- ggplot(dis, aes(x = as.numeric(degree), y = as.numeric(CumulativeDegree), color = dx)) +
  #geom_point(size = 1.5) +
  geom_line(size = 2) +
  #geom_smooth(method = "loess", se = FALSE, size = 1.5) +  # Suavizado LOESS para suavizar la línea
  scale_x_continuous(limits = c(min(dis$degree), max(dis$degree)), breaks = seq(min(dis$degree), max(dis$degree), by = 15)) +  # Ajusta el eje X
  labs(x = "Degree",
       y = "Cumulative frequency", 
       color = NULL) +  # Elimina el título de la leyenda
  #scale_x_continuous(breaks = seq(1, max(dis$degree), by =4)) + # Adjust the x-axis to show only integer degrees
  scale_color_manual(values = c("AD" = "#A3333D", "control" = "#477998")) + # Custom colors for AD and no AD
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  # Rotate X-axis labels at 45 degrees
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8),  # Tamaño reducido para etiquetas en el eje X
    panel.grid.major = element_blank(),  # Elimina las líneas de la cuadrícula principal
  )
degree_distr_x.p

#Save
# ggsave(filename = "accum_distr.jpg",
#        plot = degree_distr_x.p,
#        width = 50,
#        height = 25,
#        units = "cm",
#        dpi = 300,
# )

#Grob

library(cowplot)

#convert the secondary graphic to a graphic object
degree_inset <- ggdraw() +
  draw_plot(log_log_combined) +  #Main plot
  draw_plot(degree_distr_x.p, 
            x = 0.4, y = 0.68, 
            width = 0.6,
            height = 0.3)  #Inset plot with positioning
degree_inset

#Save

# ggsave(filename = "degree_dis_inset.jpg",
#        plot = degree_inset,
#        width = 25,
#        height = 25,
#        units = "cm",
#        dpi = 300)

###########  k vs P(k) distribution ########

k_pk_AD <- ggplot(ADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
  labs(x = expression("<k>"), y = expression(p(k)), 
       title = "", 
       subtitle = "AD network \ndegree distribution") +
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  geom_smooth(method = NULL, se = FALSE, color = "black", linetype = "dashed", 
              size = 1.2, formula = y ~ x) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )

k_pk_noAD <-  ggplot(noADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
  labs(x = expression("<k>"), y = expression(p(k)), 
       title = "", 
       subtitle = "No AD network \ndegree distribution") +
  geom_smooth(method = NULL, se = FALSE, color = "black", linetype = "dashed", 
              size = 1.2, formula = y ~ x) +  
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )

k_pk_grid <- grid.arrange(k_pk_AD, k_pk_noAD, ncol =2)

###### k vs logp(k) #####

k_logpk_AD <- ggplot(ADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  
  scale_y_log10() +  # Escala logarítmica en el eje y
  labs(x = expression("<k>"), y = expression(logp(k)), 
       title = "", 
       subtitle = "AD network \ndegree distribution") +
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  geom_smooth(method = NULL, se = FALSE, color = "black", linetype = "dashed", 
              size = 1.2, formula = y ~ x) +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )

k_logpk_noAD <-  ggplot(noADdegree_freq, aes(x = degree, y = Prob)) +
  geom_point(color = "blue", size = 2, alpha = 0.8) +  # Puntos más grandes y semitransparentes
  labs(x = expression("<k>"), y = expression(logp(k)), 
       title = "", 
       subtitle = "No AD network \ndegree distribution") +
  scale_y_log10() +  # Escala logarítmica en el eje y
  geom_smooth(method = NULL, se = FALSE, color = "black", linetype = "dashed", 
              size = 1.2, formula = y ~ x) +  
  theme_minimal(base_size = 14) +  # Tamaño base de texto para mejor legibilidad
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),  # Centrar título y agrandar
    axis.title = element_text(size = 14),  # Tamaño de títulos de ejes
    axis.text = element_text(size = 12),   # Tamaño de etiquetas de ejes
    panel.grid = element_line(size = 0.5, color = "grey80")  # Líneas de cuadrícula más sutiles
  )

k_pk_grid <- grid.arrange(k_logpk_AD, k_logpk_noAD, ncol =2)

#Final panel

panel <- grid.arrange(log_log_combined, degree_distr_x.p, ncol = 2)

# ggsave(filename = "panel_loglog_accum_distr.jpg",
#        plot = panel,
#        width = 50,
#        height = 25,
#        units = "cm",
#        dpi = 300,
# )

#END