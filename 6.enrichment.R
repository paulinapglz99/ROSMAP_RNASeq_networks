#
#6.enrichment.R
#Script for a gene-list enrichment
#paulinapglz.99@gmail.com

pacman::p_load(clusterProfiler)

#Libraries

pacman::p_load(clusterProfiler, 
               tidyverse)

# Lectura de la tabla de genes diferencialemente expresados
degs = readRDS("data/degs.RDS")

# necesitamos el log2 fold change 
original_gene_list <- degs$logFC

# Nombramos el vector
names(original_gene_list) <- degs$ESGN

# eliminamos cualquier NA 
gene_list<-na.omit(original_gene_list)

# odernamos la lista en orden decreciente (requerido por clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)

# extraemos los genes significativos (p ajustada < 0.05)
sig_genes_df = subset(degs, adj.P.Val < 0.05)

# para los resultados significativos, queremos filtrar por log2fold change
genes <- sig_genes_df$logFC

# nombramos el vector
names(genes) <- sig_genes_df$ESGN

# omitimos posibles NAs
genes <- na.omit(genes)

# filtramos por mínimo log2fold change (log2FoldChange > 2)
genes <- names(genes)[abs(genes) > 2]