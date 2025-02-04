#differential_expression_limma.R
#This script performs differential expression with limma-voom 

#Libraries --- ---

# BiocManager::install("limma")
# BiocManager::install("Glimma")
# BiocManager::install("edgeR")
# BiocManager::install("EnhancedVolcano")

pacman::p_load("tidyverse",
               "limma",
               "Glimma", 
               "edgeR", 
               "stringr", 
               "edgeR",
               "EnhancedVolcano")

#Get data --- ---

counts <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.txt") %>% as.data.frame()
counts <- counts[ -c(1:4),] #Delete alignment stats
counts$feature <- str_remove(counts$feature, "\\..*$")

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
dim(counts)
counts <- bind_rows(counts, repeated_rows)
dim(counts)
#[1] 60558  1142

rownames(counts) <- counts$feature
#  Convert selected columns to integers
counts <- counts %>% mutate(across(-feature, as.integer))

# Get metadata
metadata <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
dim(metadata)
#[1] 1141   42

table(metadata$dicho_NIA_reagan, useNA = "ifany")
#  0    1 <NA> 
#  307  573  261  

#Filter to have only Samples with dicho_NIA_Reagan metadata

metadata <- metadata %>% filter(!is.na(dicho_NIA_reagan))

counts <- counts %>% dplyr::select(one_of(metadata$specimenID))
dim(counts)
#[1] 60558   880

# #Annotation 
# 
# mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
# 
# myannot <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "gene_biotype"),
#                  filters = "ensembl_gene_id",
#                  values = rownames(counts),
#                  mart = mart)

#Pre-processing --- ---

# Filtering out lowly expressed genes

mycpm <- cpm(counts)

filterByExpr(y = dge, min.count = 10, min.prop = 0.8)

plot(counts[,1],mycpm[,1],xlim=c(0,20),ylim=c(0,0.5))
abline(v=10,col=2)
abline(h=0.15,col=4)

thresh <- mycpm > 0.15
keep <- rowSums(thresh) >= 3
table(keep)

counts.keep <- counts[keep,]
dim(counts.keep)

## Convert to DGEList object
y <- DGEList(counts.keep)

#Sample names
snames <- colnames(counts)

## Convert to DGEList object
dge <- DGEList(counts)
d0 <- calcNormFactors(dge)

## Filtering out lowly expressed genes
# It must have at least 10 counts in at least 80% of the samples, which translates to 8 out of 10 samples.
counts.keep <- filterByExpr(y = dge, min.count = 10, min.prop = 0.8)
counts.keep <- dge[counts.keep, , keep.lib.sizes = FALSE]
dim(counts.keep)
#[1] 34291   880

snames <- colnames(counts.keep)


#Quality control --- ---

# Library sizes

#Library sizes are a technical bias, to be corrected.
#This plot should look not that variant
barplot(dge$samples$lib.size)

#Distribution of counts

# Get log2 counts per million
logcpm <- cpm(dge$counts,log=TRUE)
# Check distributions of samples using boxplots

metadata.dx <- metadata %>% dplyr::select(specimenID, dicho_NIA_reagan, cogdx, sequencingBatch)
metadata.dx$color_NIA_reagan <- ifelse(metadata.dx$dicho_NIA_reagan == 0, yes = "blue", no = "red")
metadata.dx$color_cogdx <- metadata.dx %>%
  mutate(color_cogdx = case_when(
    cogdx == 1 ~ "red",
    cogdx == 2 ~ "blue",
    cogdx == 3 ~ "green",
    cogdx == 4 ~ "yellow",
    cogdx == 5 ~ "purple",
    cogdx == 6 ~ "orange"
  ))

metadata.dx$color_cogdx <- metadata.dx %>%
  mutate(color_cogdx = case_when(
    cogdx == 1 ~ "red",
    cogdx == 2 ~ "blue",
    cogdx == 3 ~ "green",
    cogdx == 4 ~ "yellow",
    cogdx == 5 ~ "purple",
    cogdx == 6 ~ "orange"
  ))

# Reordena logcpm según el nuevo orden de metadata.dx
logcpm <- logcpm[, metadata.dx$specimenID]

boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,col=metadata.dx$color ,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")

# MDS plots
plotMDS(dge,col=metadata.dx$color_cogdx)
legend("topright", legend = levels(as.factor(metadata.dx$dicho_NIA_reagan)), fill=c("red","blue"))

prcomp(dge$counts)

# Normalization --- ---

## TMM normalization
dge <- calcNormFactors(dge)
dge$samples

# Differential expression analysis --- ---

#Design matrix

design <- data.frame(row.names = colnames(counts), dicho_NIA_reagan = metadata$dicho_NIA_reagan)

design <- model.matrix(~dicho_NIA_reagan, data = design)

# make the contrasts matrix
cont.matrix <- makeContrasts('', levels = design)

## Voom (variance modeling at the observational level) transform the data

v <- voom(dge, design, plot = TRUE)

#Fit linear models 
#Fit the linear model using the lmFit function and apply the contrasts.

fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

results <- decideTests(fit2)
summary(results)

topTable(fit2, coef = 1, sort.by = "p")

#Vulcano plot --- ---

tt<-topTable(fit2, coef = 1, adjust="fdr", number=Inf)

par(mfrow=c(1,1))
vulcano_plot <- EnhancedVolcano(tt, 
               lab = rownames(tt),
                x = "logFC",
                y = "adj.P.Val",
                pCutoff = 0.05,
                FCcutoff = 2,
                ylab = bquote(~-Log[10]~adjusted~italic(P)),
                title= "Differential Expression results")

#Heatmap --- ---

gene_list <- topTable(fit2, coef=1, number=dim(fit2)[1], sort.by="logFC")
selected <- gene_list[which(gene_list$adj.P.Val<0.05 & abs(gene_list$logFC)>=2),]
HM_matrix <- as.matrix(counts[which(rownames(counts) %in% rownames(selected)),])
dim(HM_matrix)

cases_colors <- rep(c("red", "blue"), each = 20)

par(mar=c(7,4,4,2)+0.1)

heatmap <- heatmap.2(HM_matrix, 
          col=greenred(75), 
          scale="row",
          key=TRUE,
          symkey=FALSE,
          density.info="none",
          trace="none",
          cexRow=0.005,
          ColSideColors = cases_colors,
          margins=c(12,8))

#Save results --- ---

# Write table with differentially expressed genes, logFC > 2, adj. P value < 0.05 
write.csv(selected, file = "DE_genes_results.csv", quote = FALSE, row.names = FALSE)
write.csv(gene_list, file = "full_DE_genes_results.csv", quote = FALSE, row.names = FALSE)
 
#END --- --- 


#Gente con cogdx == 4 y NIA Reagan 0

metadata.v <- metadata %>% filter(cogdx > 3 & dicho_NIA_reagan == 0 )


