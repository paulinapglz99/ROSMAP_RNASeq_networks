# 
#
#DE analysis and Vulcano plot

#Libraries --- ---

pacman::p_load('edgeR', 
               'dplyr')

#Get data --- --- 

#Counts data

counts <- vroom::vroom(file = '/datos/rosmap/FPKM_data/ROSMAP_QCed_count_matrixfiltered_090224.tsv')

#Metadata

metadata <- vroom::vroom(file = '/datos/rosmap/metadata/ROSMAP_QC_fitlered_annotation100224.tsv')

#Manage data --- ---

#Creat

counts_1 <- as.matrix(counts[2:ncol(counts)])
dim(counts_1)

d <- DGEList(counts = counts_1, genes = counts$ensembl_gene_id)

#Estimate dispersion

d <- estimateDisp(d)

#The fitGLM function fits a generalized linear model (GLM) to the data

group <- factor(c(1,1,2,2,3,3))
design <- model.matrix(~group)

d <- glmQLFit(d, design)

# Crear una matriz de diseño
design <- model.matrix(~ condición + factor(otra_variable), data = metadata)

# Ajustar por diseño
d <- fitGLM(d, design)


#Obtain p-values and log2 fold changes

results <- glmLRT(d)

ggplot(data=counts, aes(x=log2FoldChange, y=pvalue)) + geom_point()