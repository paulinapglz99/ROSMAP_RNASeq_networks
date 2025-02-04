#differential_expression_limma.R
#This script performs differential expression with limma-voom 

#Libraries --- ---

# BiocManager::install("limma")
# BiocManager::install("Glimma")
# BiocManager::install("edgeR")
# BiocManager::install("EnhancedVolcano")
# BiocManager::install("sva")

pacman::p_load("tidyverse",
               "limma",
               "Glimma", 
               "edgeR", 
               "stringr", 
               "pheatmap",
               "edgeR",
               "EnhancedVolcano", 
               "sva")

#Get data --- ---

#Get counts

counts <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.txt") %>% as.data.frame()
counts <- counts[ -c(1:4),] #Delete alignment stats
counts$feature <- str_remove(counts$feature, "\\..*$")

# Get metadata
metadata <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
dim(metadata)
#[1] 1141   42

table(metadata$cogdx,  useNA = "ifany")
# AD  MCI noAD <NA> 
#   164  103  399  127 

# Combine and summarize duplicate features, then replace the original rows
# counts <- counts %>%
#   group_by(feature) %>% #group the data by the feature column.
#   summarize(across(everything(), #This combines the duplicate rows by calculating the median of their values.
#                    median, na.rm = TRUE),
#             .groups = 'drop')    #This ensures that the result does not maintain the grouping.
# dim(counts)
# #[1] 60558  1142

# Verify repeated genes in 'feature' column
repeated_values <- counts %>%
  group_by(feature) %>%
  filter(n() > 1) %>%
  distinct(feature) %>%
  pull(feature)
length(repeated_values)

# Ver las filas duplicadas
repeated_rows <- counts[counts$feature %in% repeated_values, ]
dim(repeated_rows)

repeated_rows <- repeated_rows[order(repeated_rows$feature),]

repeated_rows <- repeated_rows %>%
  group_by(feature) %>%
  summarize(across(everything(), median, na.rm = TRUE))
dim(repeated_rows)

#Delete rows in the matrix

counts <- counts %>% filter(!feature  %in% repeated_rows$feature)
counts <- bind_rows(counts, repeated_rows)
dim(counts)
#[1] 60558  1142

rownames(counts) <- counts$feature
#  Convert selected columns to integers
counts <- counts %>% mutate(across(-feature, as.integer))

#Filter to have only Samples with dicho_NIA_Reagan metadata

metadata <- metadata %>% filter(!is.na(cogdx))

counts <- counts %>% dplyr::select(one_of(metadata$specimenID))
dim(counts)
#[1] 60558   793

#QC --- ---

#Combat-seq to omit batch

counts <- ComBat_seq(counts = as.matrix(counts), batch = metadata$sequencingBatch)

# Create dge object --- ---

group <- factor(metadata$cogdx) # Ajusta esto según tu metadata
dge <- DGEList(counts, group = group)

#Filter data on lower count rate 

drop <- filterByExpr(y = dge, min.count = 10, min.prop = 0.8)

dge <- dge[drop,]
dim(dge)
#[1] 21544   793

# Library sizes

#Library sizes are a technical bias, to be corrected.
#This plot should look not that variant
barplot(dge$samples$lib.size)

#Normalize the data using the TMM scaling factor method.

dge <- calcNormFactors(dge, method = "TMM")

# MDS plot
mds <- plotMDS(dge, labels = group, col = as.numeric(group))

# Boxplot de valores de cuentas normalizadas
logCPM <- cpm(dge, log = TRUE)
boxplot(logCPM, las = 2, col = as.numeric(group))

# Model Matrix 

design <- model.matrix(~0 + group)
#colnames(design) <- c("dicho_NIA_Reagan_0", "dicho_NIA_Reagan_1")
rownames(design) <- metadata$specimenID
y <- estimateDisp(dge, design)

# Differential expression 

#voom for transformation and linear modeling
v <- voom(y, design, plot = TRUE)
fit <- lmFit(v, design)
#Use makeContrasts to define the contrasts you want to test
cont.matrix <- makeContrasts('group1 - group4', levels = design)
fit <- contrasts.fit(fit, cont.matrix )
fit <- eBayes(fit)
results <- topTable(fit, number = Inf, adjust.method = "BH")
results$feature <- rownames(results)

#Extract DEGs

# add a column of NAs
results$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and pvalue < 0.05, set as "UP" 
results$diffexpressed[results$logFC > 0.75 & results$adj.P.Val < 0.05] <- "UP"
# if log2Foldchange < -0.5 and pvalue < 0.05, set as "DOWN"
results$diffexpressed[results$logFC < -0.75 & results$adj.P.Val < 0.05] <- "DOWN"

#
DEGS <- results %>% filter(diffexpressed != "NO")
rownames(DEGS) <- DEGS$ensembl_gene_id
dim(DEGS)
#[1] 58  7

#Vulcano plot --- ---

volc <- EnhancedVolcano(results,
                        lab = rownames(results),
                        x = 'logFC',
                        y = 'adj.P.Val',
                        pCutoff = 0.05,
                        FCcutoff = 0.75,
                        title = 'Volcano plot',
                        subtitle = 'Differential Expression Analysis')

#Heatmap --- ---

#Extraer DEGS y hacer heatmap solo con eso 
normcounts_DEG <- as.data.frame(logCPM) %>% filter(rownames(logCPM) %in% DEGS$feature)

metadata_DEG <-  data.frame(
  cogdx = as.factor(metadata$cogdx),  # Replace with your actual metadata
  row.names = colnames(normcounts_DEG)  # Ensure row names match column names of counts_DEG
) %>% arrange(cogdx)

#Order

normcounts_DEG <- normcounts_DEG[, rownames(metadata_DEG)]

# Create the heatmap
pheatmap(normcounts_DEG,
         cluster_rows = T,   # Cluster the rows
         cluster_cols = T,   # Do not cluster the columns
         main = "Heatmap of Differentially Expressed Genes", 
         annotation_col = metadata_DEG,
         show_colnames = FALSE  # Hide the column names (sample names)
         )

#END