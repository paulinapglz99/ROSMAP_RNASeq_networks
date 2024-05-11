#
#3.differential_expression.R
#Script that makes differential expression of data 
#By paulinapglz.99@gmail.com
#For further details of how NOISeq works, go to https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf

#Libraries --- ---

pacman::p_load('dplyr', 
               'biomaRt',
               'NOISeq',
               'edgeR', 
               'EDASeq', 
               "ggplot2")

#Get data --- --- 

#Counts

counts <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/ROSMAP_RNAseq_filtered_counts_DLPFC.txt")
dim(counts)
#[1] 22094   445

#Metadata --- ---

metadata <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_filtered_DLPFC.txt")
dim(metadata)
#[1] 446  41

#Annotation --- ---

myannot <- vroom::vroom(myannot, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/RNA_seq_annotation_filtered_DLPFC.txt")

#readData object --- --- 

mydata <- 


mynoiseq = noiseq(counts[-1], k = 0.5, norm = "rpkm", factor = "Tissue", pnr = 0.2,
                  nss = 5, v = 0.02, lc = 1, replicates = "technical")
