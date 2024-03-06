# 
#
#DE analysis and Vulcano plot
#BIG NOTE: 

print("If FPKM is really all you have, then convert the values to a log2 scale (y = log2(FPKM+0.1) say) 
      and do an ordinary limma analysis as you would for microarray data, using eBayes() with trend=TRUE.
      Do not use voom, do not use edgeR, do not use DESeq. (Do not pass go and do not collect $200.) 
      This isn't 100% ideal, but is probably the best analysis available. You make this method somewhat better 
      by using arrayWeights() as well which, in this context, will attempt to estimate the library sizes the FPKMs 
      were originally computed from. Nevertheless, the mean-variance trend estimated by limma from the logFPKMs will
      never be as smooth or as informative as the trend that would have been estimated had you had the real counts.")

#Libraries --- ---

pacman::p_load('dplyr',
               'limma')
          
#Get data --- --- 

#Counts data

FPKM_counts <- vroom::vroom(file = '/datos/rosmap/FPKM_data/ROSMAP_QCed_count_matrixfiltered_090224.tsv')

#Metadata

metadata <- vroom::vroom(file = '/datos/rosmap/metadata/RNA_seq_metadata_080224.csv')

#Manage data --- ---

#Originalmente tenemos los datos en FPKMs,  pero los quiero  convertir a TPMs  con una  formula

log2data <- sapply(FPKM_counts[, 2:ncol(FPKM_counts)], function(x) log2(x + 0.1))

#Originalmente tenemos los datos en FPKMs,  pero los quiero  convertir a TPMs  con una  formula

fpkmToTpm <- function(fpkm)
{
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

#Convert FPKMs to TPMs

ex <- fpkmToTpm(counts[-1])

rownames(ex) <- counts$ensembl_gene_id

#Create a 





# log2 transform

ds <- DESeqDataSetFromMatrix(countData = ex,
                              colData = NULL,
                              design = ~1)

# Realiza el análisis de expresión diferencial
dds <- DESeq(dds)

# Extrae los resultados
results <- results(dds)

# Filtra los resultados significativos
results_sig <- subset(results, padj < 0.05)

# Muestra los resultados
print(results_sig[, c("gene_symbol", "log2FoldChange", "pvalue")])











