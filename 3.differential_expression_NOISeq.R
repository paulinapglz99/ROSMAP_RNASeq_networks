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

mynoiseq <- noiseqbio(noiseqData,
                   k = 0.5,     #Counts equal to 0 are replaced by k. By default, k = 0.5. 
                   norm = "n",  #. If the data have been previously normalized, norm parameter must be set to "n"
                   factor = "dicho_NIA_reagan",
                   pnr = 0.2,
                   nss = 5,
                   v = 0.02, 
                   lc = 1, 
                   replicates = "technical")

mynoiseqbio <- noiseqbio(noiseqData,
                         k = 0.5,           #Counts equal to 0 are replaced by k. By default, k = 0.5. 
                         norm = "n",        #. If the data have been previously normalized, norm parameter must be set to "n"
                         factor = "dicho_NIA_reagan", #name of factor whose levels are the conditions to be compared. 
                         lc = 0,            # Length correction is done by dividing expression by length^lc. By default, lc = 0. 
                         r = 20,            #Number of permutations to generate noise distribution by resampling. 
                         adj = 1.5,         #Smoothing parameter for the Kernel Density Estimation of noise distribution. Higher values produce smoother curves. 
                         plot = FALSE,      
                         a0per = 0.9, 
                         random.seed = 12345, #Random seed.
                         filter = 0)        #If filter=0, no filtering is performed.

#Select differential expression genes --- --- 
#Please remember that, when using NOISeq, the probability of dierential expression is not equivalent to 1 − pvalue. We recommend for q to use values around 0.8

#Differentially expressed features
mynoiseq.deg <- degenes(mynoiseqbio, q = 0.8, M = NULL)
#[1] "11841 differentially expressed features"

#Differentially expressed features that are more expressed in condition 1 than in condition 2 (M = "up"):
mynoiseq.deg_up = degenes(mynoiseqbio, q = 0.8, M = "up")
#[1] "9130 differentially expressed features (up in first condition)"

#Differentially expressed features that are under-expressed in condition 1 with regard to condition 2 (M = "down"):
mynoiseq.deg_down = degenes(mynoiseqbio, q = 0.8, M = "down")
#[1] "2711 differentially expressed features (down in first condition)"

#Plot --- ---

#Vulcano plot

ggplot(mynoiseq.deg, aes(x = log2FC, y = -log10(prob))) +
  geom_point(aes(color = ifelse(prob < 0.05 & abs(log2FC) > 1, "Significant", "Not Significant"))) +
  scale_color_manual(values = c("red", "black")) +  # Personalizar colores
  theme_minimal() +
  labs(title = "Volcano Plot",
       x = "log2 Fold Change",
       y = "-log10(p-value)")

##

DE.plot(mynoiseqbio, q = 0.9, graphic = "expr", log.scale = TRUE)
#[1] "1 differentially expressed features"

DE.plot(mynoiseqbio, chromosomes = NULL, log.scale = TRUE, join = FALSE,
 q = 0.8, graphic = "distr")
#[1] "11841 differentially expressed features"

#END
