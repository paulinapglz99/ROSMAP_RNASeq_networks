#
#3.differential_expression_DEseq2.R
#this script processes RNA-seq data to identify differentially expressed genes. 

#Libraries --- ---

#BiocManager::install("DESeq2")

pacman::p_load("tidyverse", 
               "DESeq2", 
               "pheatmap", 
               "ggplot2", 
               "ggrepel", 
               "stringr",
               "biomaRt", 
               "RColorBrewer", 
               "gridExtra")

#Set seed --- --- 

set.seed(10)

#Define functions --- ---

#Function to translate gene names
ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Connecting to the Ensembl database through biomaRt

# Define function to convert from ENSMBL to SYMBOL
convert_ens_to_symbol <- function(ensembl_ids) {
 trad <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
        filters = "ensembl_gene_id",
        values = ensembl_ids,
        mart = ensembl)
  trad$external_gene_name <- ifelse(trad$external_gene_name == "", trad$ensembl_gene_id, trad$external_gene_name)
  return(trad)
}

#Delete duplicates
del_dupl <- function(counts) {
  # Verificar genes repetidos
  repeated_values <- counts %>%
    group_by(feature) %>%
    filter(n() > 1) %>%
    distinct(feature) %>%
    pull(feature)
  
  # Ver las filas duplicadas
  repeated_rows <- counts[counts$feature %in% repeated_values, ]
  
  # Ordenar y calcular la mediana de los valores duplicados
  repeated_rows <- repeated_rows[order(repeated_rows$feature),]
  repeated_rows <- repeated_rows %>%
    group_by(feature) %>%
    summarize(across(everything(), median, na.rm = TRUE))
  
  # Eliminar las filas duplicadas y agregar las filas con la mediana calculada
  counts <- counts %>% filter(!feature %in% repeated_rows$feature)
  counts <- bind_rows(counts, repeated_rows)
  
  # Convertir las columnas seleccionadas a enteros
  counts <- counts %>% mutate(across(-feature, as.integer))
  
  return(counts)
}

#Get data --- ---

counts <- readRDS(file = "~/redesROSMAP/RNAseq_QC_NOISeq/QC_nextflow/data/ROSMAP_RNAseq_rawcounts_DLPFC.rds") %>% as.data.frame()
dim(counts)
#[1]  60607  1142

# Get metadata
metadata <- vroom::vroom(file ="~/redesROSMAP/RNAseq_QC_NOISeq/QC_nextflow/metadata/RNA_seq_metadata_DLPFC.txt")
dim(metadata)
#[1] 1141   42

table(metadata$dicho_NIA_reagan, useNA = "ifany")

# 0    1 <NA> 
#   307  486  348 

#Filter to obtain only the ones with NIA-Reagan dicho

counts <- counts[-c(1:4),]

#Delete duplicates

counts <- counts %>% mutate(feature = str_remove(feature, "\\..*$"))

counts <- del_dupl(counts)

#Prepare feature column

rownames(counts)<- NULL
counts <- counts %>% column_to_rownames(var = "feature")

metadata <- metadata %>% filter(!is.na(dicho_NIA_reagan))

#Only samples with metadata
counts <- counts %>% dplyr::select(all_of(metadata$specimenID))
dim(counts)
#[1] 60558   880

#Differential expression --- ---

#Experimental design

coldata <- as.data.frame(colnames(counts))
colnames(coldata) <- "specimenID"

coldata <- coldata %>% 
  left_join(metadata, by = "specimenID")%>%
  dplyr::select("specimenID", "sequencingBatch","msex","cogdx", "ceradsc", "dicho_NIA_reagan")
coldata$dicho_NIA_reagan <- as.factor(coldata$dicho_NIA_reagan)  #Convert to factor

#DESeqData object
#DESeqData object con control para sequencingBatch
dds <- DESeqDataSetFromMatrix(countData = counts, 
                              colData = coldata,
                              design = ~ sequencingBatch + dicho_NIA_reagan) # Control sequencing batch and compare NIA-Reagan

#Differential expression analysis --- ---

#Specify conditions to compare
#Establishing 0 (not AD) as reference
dds$dicho_NIA_reagan <- relevel(dds$dicho_NIA_reagan, ref = "0")

dds <- DESeq(dds)  #Slow

#Results

res <- results(dds,
               contrast = c("dicho_NIA_reagan", "1", "0"), #contrast = c("condition", "problem", "control"))
               pAdjustMethod = 'BH',   #he method to use for adjusting p-values
               alpha = 0.05) #FDR

#How many adjusted p-values were less than 0.05?
sum(res$padj < 0.05, na.rm=TRUE)
#[1] 4423

res.df <- res %>% as.data.frame()
res.df <- res.df %>% filter(!is.na(padj))
dim(res.df)
#[1] 27683     6

# add a column of NAs
res.df$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
res.df$diffexpressed[res.df$log2FoldChange > 0.5 & res.df$padj < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
res.df$diffexpressed[res.df$log2FoldChange < -0.5 & res.df$padj < 0.05] <- "DOWN"

table(res.df$diffexpressed)

# DOWN    NO    UP 
#   18 34629   594 

#Add gene names in SYMBOL --- ---

#Create dictionary
symbol <- convert_ens_to_symbol(rownames(res))

# Make rownames columns
res.df <- res.df %>% mutate(ensembl_gene_id = rownames(.), .before = 1) 
res.df <- res.df %>% left_join(symbol, by = "ensembl_gene_id") #merge

# BaseMean by condition --- ---

#Extract counts normalized 

normcounts <- counts(dds, normalized = TRUE) %>% as.data.frame()

# Calcula los valores de baseMean por condición
baseMean_0 <- rowMeans(normcounts[, coldata$dicho_NIA_reagan == "0"]) %>% as.data.frame()
colnames(baseMean_0) <- "baseMean_NIA_R_0"

baseMean_1 <- rowMeans(normcounts[, coldata$dicho_NIA_reagan == "1"])%>% as.data.frame()
colnames(baseMean_1) <- "baseMean_NIA_R_1"

baseMean_genes <- cbind(baseMean_0, baseMean_1)
baseMean_genes <- baseMean_genes %>% rownames_to_column(var = "feature")

#Extract DEGs --- ---

DEGS <- res.df %>% filter(diffexpressed != "NO")
dim(DEGS)

baseMean_DEGS <- baseMean_genes %>% filter(feature %in% DEGS$ensembl_gene_id)
DEGS <-  DEGS %>% left_join(baseMean_DEGS, by = c("ensembl_gene_id" = "feature"))

#Save log2fold change full list --- ---
# 
vroom::vroom_write(res.df, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_DE_full_gene_list_dichoNIAReagan.txt")
# 
# #Save DEGs list --- ---
# 
vroom::vroom_write(DEGS, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_DEGS_dichoNIAReagan.txt")

#Vulcano plot --- ---

#Create labels for Vplot
res.df <- res.df %>% mutate(delabel = ifelse(diffexpressed != "NO", external_gene_name, NA ))

#
vplot <-  ggplot(data=res.df, aes(x=log2FoldChange, y= -log10(padj),  col = diffexpressed, label = delabel)) +
  geom_point() +
  # Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold 
  geom_vline(xintercept=c(0.5, -0.5), col="red") + #log2FoldChange threshold is 0.5
  geom_hline(yintercept=-log10(0.05), col="red") + #p-value threshold is 0.05
  scale_color_manual(values=c("#4E8098", "gray", "#A31621")) +
  theme_minimal() +
  geom_text_repel(max.overlaps = 25) +
  labs(title = "Differential expression", 
       subtitle = "Using dichotomic NIA-Reagan Criteria")

vplot
#Save Vulcano Plot ---

ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_dicho_NIAReagan_vulcano_plot.png", 
       plot = vplot, 
       width = 11,
       height = 10,
       units = "in",
       dpi = 300)

#Create matrix of DEG counts --- ---

DEG_mat <- normcounts %>% filter(rownames(normcounts) %in% DEGS$ensembl_gene_id)
dim(DEG_mat)
#[1] 28  880

htmap <- pheatmap(DEG_mat,
                  cluster_rows = F,   # Cluster the rows
                  cluster_cols = F,   # Do not cluster the columns
                  main = "Heatmap of Differentially Expressed Genes", 
                  #    annotation_col = metadata_DEG,
                  show_colnames = FALSE  # Hide the column names (sample names)
)

htmap

# Save heatmaps --- ---

# ggsave("/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_dicho_NIAReagan_heatmap.png", 
#        plot = htmap, 
#        width = 11,
#        height = 15,
#        units = "in",
#        dpi = 300)

#GSEA --- ---

#Libraries for GSEA

pacman::p_load(clusterProfiler,
  org.Hs.eg.db,
  fgsea, 
  enrichplot, 
  igraph)

#Set names of genes to search
DEGS.g <- setNames(DEGS$log2FoldChange, DEGS$ensembl_gene_id)
DEGS.g <- sort(DEGS.g, decreasing = TRUE)

#Biological process

gse_BP <- gseGO(geneList = DEGS.g, 
      ont ="BP", 
      keyType = "ENSEMBL", 
    #  nPerm = 10000, 
      minGSSize = 3, 
      maxGSSize = 800, 
      pvalueCutoff = 0.05, 
      verbose = TRUE, 
      OrgDb = "org.Hs.eg.db", 
      pAdjustMethod = "BH")

#dotplot
dotplot_BP <- dotplot(gse_BP, showCategory=10, split=".sign") + facet_grid(.~.sign) + 
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))

#emmap plot
simatrix_BP <- pairwise_termsim(gse_BP)
simatrix_BP.m <- as.matrix(simatrix_BP@termsim)
simatrix_BP.g <- graph_from_adjacency_matrix(simatrix_BP.m)
degree_BP <- degree(simatrix_BP.g) %>% as.data.frame()
degree_BP$. <- sort(degree_BP$., decreasing = TRUE)

ggraph(simatrix_BP.g, layout = 'kk') + 
  geom_edge_fan() + 
  geom_node_point()

emmap_BP <- emapplot(simatrix_BP, showCategory = 25, group_legend = T,
                     layout = "star", repel = T) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
emmap_BP

# 
vroom::vroom_write(simatrix_BP@result, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_GO_BP_results_dichoNIAReagan.txt")

#Molecular functin

gse_MF <- gseGO(geneList = DEGS.g, 
                ont ="MF", 
                keyType = "ENSEMBL", 
                #nPerm = 10000, 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = "org.Hs.eg.db", 
                pAdjustMethod = "BH")

#dotplot
dotplot_MF <- dotplot(gse_MF, showCategory=10, split=".sign") + facet_grid(.~.sign) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

#emmap plot
simatrix_MF <- pairwise_termsim(gse_MF)
simatrix_MF.df <- as.matrix(simatrix_MF@termsim)
simatrix_MF.g <- graph_from_adjacency_matrix(simatrix_MF.df)
degree_MF <- degree(simatrix_MF.g) %>% as.data.frame()
degree_MF$. <- sort(degree_MF$., decreasing = TRUE)

emmap_MF <- emapplot(simatrix_MF, showCategory = 19, group_legend = T,
                     layout = "circle", repel = T) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
emmap_MF

#Cellular process

gse_CC <- gseGO(geneList = DEGS.g, 
                ont ="CC", 
                keyType = "ENSEMBL", 
                #nPerm = 10000, 
                minGSSize = 3, 
                maxGSSize = 800, 
                pvalueCutoff = 0.05, 
                verbose = TRUE, 
                OrgDb = "org.Hs.eg.db", 
                pAdjustMethod = "BH")

#dotplot

dotplot_CC <- dotplot(gse_CC, showCategory=10, split=".sign") + facet_grid(.~.sign) + 
  theme(plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"))

#emmap plot
simatrix_CC <- pairwise_termsim(gse_CC)
simatrix_CC.df <- as.matrix(simatrix_CC@termsim)
simatrix_CC.g <- graph_from_adjacency_matrix(simatrix_CC.df)
degree_CC <- degree(simatrix_CC.g) %>% as.data.frame()
degree_CC$. <- sort(degree_CC$., decreasing = TRUE)

emmap_CC <- emapplot(simatrix_CC, showCategory = 19, group_legend = T,
                     layout = "circle", repel = T) +
  theme(plot.margin = unit(c(1,1,1,1), "cm"))
emmap_CC


#Save results

vroom::vroom_write(simatrix_CC@result, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/DEGs/ROSMAP_DLFPC_GO_CC_results_dichoNIAReagan.txt")

# Create titles

title_BP <- textGrob("Biological Process", rot = 90, gp = gpar(fontsize = 15))
title_CC <- textGrob("Cellular Component", rot = 90, gp = gpar(fontsize = 15))

# Crear los arreglos de gráficos con los títulos
grid_total <- grid.arrange(
  arrangeGrob(title_BP, dotplot_BP, emmap_BP, nrow = 1, widths = c(1, 10, 15)),
  arrangeGrob(title_CC, dotplot_CC, emmap_CC, nrow = 1, widths = c(1, 10, 15)),
  nrow = 2
)

#Save graphs

ggsave("gsea_DEGS_DLPFC_ROSMAP_circle.jpg", 
       plot = grid_total, 
       device = "jpg", 
       height = 15,
       width = 20,  
       dpi = 300)
