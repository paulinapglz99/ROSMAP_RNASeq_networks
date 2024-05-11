#
#3.differential_expression.R
#Script that makes differential expression of data 
#By paulinapglz.99@gmail.com
#For further details of how NOISeq works, go to https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf

#Libraries --- ---

pacman::p_load('dplyr', 
               'NOISeq',
               'ggplot2')

#Get data --- --- 

#Counts

counts <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/ROSMAP_RNAseq_filteredQC_counts_DLPFC.txt") %>% as.data.frame()
dim(counts)
#[1] 22070   498
rownames(counts) <- counts$feature

#Metadata --- ---

metadata <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_filteredQC_DLPFC.txt")
dim(metadata)
#[1] 497  41

#Annotation --- ---

myannot <- vroom::vroom(file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_filteredQC_annotation_DLPFC.txt")
dim(myannot)
#[1] 22070     7

#readData object --- --- 

#Factors

#I do this to make sure the rowlength of factors match with the counts columns

factors <- data.frame("specimenID" = colnames(counts)[-1])

#Left join with metadata

factors <- factors %>% left_join(metadata, by = "specimenID")

#factors <- factors %>% dplyr::select(c(specimenID, libraryBatch, RIN , NIA_reagan_ADLikelihood, dicho_NIA_reagan, cogdx))
dim(factors)
#[1] 497  41

#Names of features characteristics

annot <- data.frame("feature" = counts$feature)
annot <- annot %>% left_join(myannot, by = "feature")

mylength <- setNames(annot$length, annot$feature)

mygc <- setNames(annot$percentage_gene_gc_content, annot$feature)

mybiotype <-setNames(annot$gene_biotype, annot$feature)

#Create NOISeq object

noiseqData <- NOISeq::readData(data = counts[-1],
                               factors = factors,           #variables indicating the experimental group for each sample
                               gc = mygc,                   #%GC in myannot
                               biotype = mybiotype,         #biotype
                               length =  mylength)          #gene length

#Differential expression --- ---
# NOISeq computes the following differentialexpression statistics for each feature:
#M (which is the log2-ratio of the two conditions) and D (the value of the difference between conditions). 

mynoiseq <- noiseq(noiseqData, k = 0.5, norm = "tmm", factor = "Tissue", pnr = 0.2,
                  nss = 5, v = 0.02, lc = 1, replicates = "technical")

mynoiseq = noiseq(mydata, k = 0.5, norm = "rpkm", factor = "Tissue", pnr = 0.2,
                  + nss = 5, v = 0.02, lc = 1, replicates = "technical")

#Hierarchical Clustering Heatmap --- --- 



