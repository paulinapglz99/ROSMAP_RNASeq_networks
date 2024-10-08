#
# 1.1prepro-mRNA.R
#Script for annotation, bias detection, correction (QC) and normalization of RNA-seq data
#By paulinapglz.99@gmail.com, adapted from https://github.com/CSB-IG/SGCCA/blob/main/prepro-mRNA.R
#For further details of how NOISeq works, go to https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf

#Libraries --- ---

#Here we use the NOISeq, edgeR and EDAseq packages for bias detection and correction

pacman::p_load('tidyverse', 
               'biomaRt',
               'NOISeq',
               'edgeR', 
               'EDASeq', 
               "ggplot2")

#Seed --- ---

set.seed(10)

#Functions  --- ---

#Homogenize metadata for sex
homogenize_exclude <- function(metadata) {
  unique_values <- unique(metadata$exclude)
  
  if ("TRUE" %in% unique_values && any(is.na(unique_values))) {
    # Caso 1: TRUE para excluir, NA para mantener
    metadata <- metadata %>% mutate(exclude = ifelse(is.na(exclude), FALSE, TRUE))
  } else {
    # Caso 2: FALSE para mantener, TRUE para excluir
    metadata <- metadata %>% mutate(exclude = ifelse(exclude == FALSE, FALSE, TRUE))
  }
  
  return(metadata)
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

#Get the data --- ---

#Read counts data, this was already filtered by 1.QC_pre_analysis.R

counts <- readRDS(file = '/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_rawcounts_DLPFC.rds') %>%
  as.data.frame()
dim(counts)
#[1] 60558   1142

#Obtain factors from metadata --- --- 

metadata <- vroom::vroom(file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_DLPFC.txt")
dim(metadata)
#[1] 1141  41

# Acondicionar metadata

metadata <- homogenize_exclude(metadata)

#Replace msex with sex
names(metadata) <- ifelse(names(metadata) == "msex", "sex", names(metadata))

#Exclude specimens by metadata

metadata <- metadata %>% filter(exclude == FALSE)
metadata <- metadata %>% filter(if_all(c(ceradsc, sex), ~ !is.na(.)))
dim(metadata)

#Acondicionar conteos

counts <- counts[-c(1:4), ]
dim(counts)

#Delete duplicates

counts <- counts %>% mutate(feature = str_remove(feature, "\\..*$"))

counts <- del_dupl(counts)
dim(counts)

#Counts

counts <- counts %>% dplyr::select(1, all_of(metadata$specimenID))
dim(counts)

#I do this to make sure the rowlength of factors match with the counts columns, it is like metadata filtering

factors <- data.frame("specimenID" = colnames(counts)[-1])   

#Left join with metadata

factors <- factors %>% left_join(metadata, by = "specimenID") %>% dplyr::select(c(specimenID, ceradsc, sex))
dim(factors)

#Filter again counts

counts <- counts %>% dplyr::select(1,intersect(colnames(counts), factors$specimenID))
dim(counts)

#Set names to rows 
rownames(counts) <- NULL
counts <- counts %>% column_to_rownames(var = "feature")

#Annotation with ensembl --- ---

#Generate mart object

mart <- useEnsembl("ensembl", dataset="hsapiens_gene_ensembl")

#We create myannot, with GC content, biotype, info for length & names per transcript

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values =  rownames(counts),  #annotate the genes in the count matrix 
                 mart = mart)


myannot <- myannot %>%  rename(feature = ensembl_gene_id)
#Add length column

myannot$length <- abs(myannot$end_position-myannot$start_position)
dim(myannot)

#Create NOISeq object bias detect and bias correction --- ---

#For NOISeq, order of factors$specimenIDs and  colnames(expression_counts)[-1] must match
identical(colnames(counts), factors$specimenID)
#[1] TRUE

#Names of features characteristics

mylength <- setNames(myannot$length, myannot$feature)

mygc <- setNames(myannot$percentage_gene_gc_content, myannot$feature)

mybiotype <-setNames(myannot$gene_biotype, myannot$feature)

#Create NOISeq object

noiseqData <- NOISeq::readData(data = counts,
                               factors = factors,           #variables indicating the experimental group for each sample
                               gc = mygc,                   #%GC in myannot
                               biotype = mybiotype,         #biotype
                               length =  mylength)          #gene length

# 0) Diagnostic of data

#each sample "s" is compared to a reference "r" (which can be arbitrarily chosen).
#by computing M values=log2(counts = countsr). 

#Confidence intervals for the M median is computed by bootstrapping.
#If the median of M values for each comparison is not in the CI, the deviation
# of the sample is significant, therefore, normalization is needed 
#"cd" means "Cumulative Distribution."

mycd <- dat(noiseqData, type = "cd", norm = F) #slow

#[1] "Warning: 6321 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 01_120405"

#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#  438      7 

#FAILED PASSED 
#491      7 

#1)check expression bias for counts
#Use a same factor variable to all bias detection

mycountsbio <- dat(noiseqData, 
                   type =  "countsbio",  
                   norm = F,      #T when already normalized counts as input
                   factor = NULL) #When NULL, all factors are considered
#Plots

png("CountsOri.png")
explo.plot(mycountsbio,
           plottype = "boxplot", #type of plot
           samples = 1:20)  #only showing the first 15 samples
dev.off()

#2)check for low count genes

png("lowcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:15)
dev.off()

#3)check for transcript composition bias

#Plot for Mvalues

png("MvaluesOri.png")
explo.plot(mycd,samples=sample(1:100,10))
dev.off()

#4)check for length & GC bias

#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%), the expression depends on the feature

myGCcontent <- dat(noiseqData,
                   k = 0,            #A feature is considered to be detected if the corresponding number of read counts is > k. 
                   type = "GCbias", 
                   factor = NULL)

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."

png("GCbiasOri.png",width=1000)
explo.plot(myGCcontent,
           samples = 1:12,   
           toplot = "global")
dev.off()

#The GC-content of each gene does not change from sample to sample, so it can be expected to
#have little effect on differential expression analyses to a first approximation

mylengthbias <- dat(noiseqData, 
                    k = 0,
                    type = "lengthbias",
                    factor = NULL)

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."

#Plot length bias

png("lengthbiasOri.png", width=1000)
explo.plot(mylengthbias, 
           samples = 1:12, 
           toplot = "global")
dev.off()

#Display full QC report, don't run if not necessary

#QCreport(noiseqData, samples = NULL, factor = NULL, norm = FALSE)

#################SOLVE BIASES###################################

#1) filter low count genes.
#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 

#Hos do I know how to filter?

mycpm <- cpm(counts[-1])

#Plot cpm vs counts
#You will know where to cut according to the horizontal line, which must be orthogonal 
#to the horizontal line and coincide with the center point. 

plot(counts[[2]],mycpm[,1],xlim=c(0,20),ylim=c(0,0.5))
abline(v=10,col=2) #v hardcoded
abline(h=0.32,col=4)#h hardcoded

#This function normalizes and filters features to only have the ones with
#CPM>0.32
#Plotted cpm_counts vs counts and 0.32 is the intersection, this may change and it is hardcoded

countMatrixFiltered <- filtered.data(counts[,-1], 
                                     factor = factors$sequencingBatch,       #using all factors
                                     norm = F,            #counts are not normalized 
                                     depth = NULL,        #Sequencing depth of samples (column totals before normalizing the data). Depth only needs to be provided when method = 3 and norm = TRUE. 
                                     method = 1,          #Method 1 (CPM) removes those features that have an average expression per condition less than cpm value and a coefficient of variation per condition higher than cv.cutoff (in percentage) in all the conditions
                                     cpm = 1,             #Cutoff for the counts per million value to be used in methods 1 and 3. 
                                     p.adj = "fdr")       #Method for the multiple testing correction
dim(countMatrixFiltered)


#Filtering out low count features...
#22225 features are to be kept for differential expression analysis with filtering method 1

#Filtering out low count features...
#22191 features are to be kept for differential expression analysis with filtering method 1 <- 

#Create new annotation with filtered matrix

myannot <- myannot %>% filter(feature %in% rownames(countMatrixFiltered))
dim(myannot)
#[1] 22095     7
#[1] 22071     7? <- 
#[1] 25943     7

myannot <- myannot[!duplicated(myannot$feature), ]
dim(myannot)
#[1] 22094     7  1 duplicated feature was deleted
#[1] 22070     7? <- 

#Not all genes are annotated in the matrix

countMatrixFiltered <- countMatrixFiltered[rownames(countMatrixFiltered) %in% myannot$feature,]
dim(countMatrixFiltered)
#[1] 22094   446
#[1] 22070   499? <- 

#Feature data

featureData <-  data.frame("feature" = rownames(countMatrixFiltered)) %>% left_join(myannot, by = "feature")
rownames(featureData) <- featureData$feature
dim(featureData)
#[1] 22094     7
#[1] 22070     7? <- 
#[1] 25941   7

#Pheno data

phenoData <- data.frame("specimenID" = colnames(countMatrixFiltered)) %>% 
  left_join(factors, by = "specimenID")
rownames(phenoData) <- phenoData$specimenID

dim(phenoData)
#[1] 446   4
#[1] 499   6
#[1] 880   4 <-

##Create EDA object --- --

mydataEDA <- newSeqExpressionSet(
  counts = as.matrix(countMatrixFiltered),
  featureData = featureData[-1],
  phenoData = phenoData[-1])

#If you're re-normalizing, there's a warning, 
#Warning message:
#  In validityMethod(object) : 'counts' contains non-integer numbers

#order for less bias

#for gc content
gcFull <- withinLaneNormalization(mydataEDA, 
                                  "percentage_gene_gc_content",
                                  which = "full")#corrects GC bias 

#for length
lFull <- withinLaneNormalization(gcFull, 
                                 "length", 
                                 which = "full")#corrects length bias 

#cd Diagnostic test for length and gc correction

mycd_lessbias <- NOISeq::dat(lFull,
                             type = "cd",
                             norm = TRUE)
#[1] "Reference sample is: 01_120405"

#Table diagnostic

table(mycd_lessbias@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#   439      6   
#   493      5 ? <- 

#Solve batch effect --- --- 

#If your data is already batch solved, only renormalize if last table diagnostic tells you to do it

# Normalization --- ---
#TMM normalization adjusts library sizes based on the assumption that most genes are not differentially expressed.

norm_count <- NOISeq::tmm(normCounts(lFull),
                          long = 1000,  # If long == 1000, no length correction is applied (no matter the value of parameter lc). 
                          lc = 0, # If lc = 0, no length correction is applied
                          k = 0) # By default, k = 0. 

noiseqData_norm_count <- NOISeq::readData(data = norm_count, 
                                          factors = factors)

#cd has to preceed ARSyN or it won't work

mycd_norm <- NOISeq::dat(noiseqData_norm_count,
                         type="cd",
                         norm=TRUE)
#[1] "Warning: 94 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 594_120522"

#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd_norm@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#  76    369  
#    88    410? <-

#ARSyNseq for batch effect solution --- ---

#Only do this if your data is not batch solved

#When batch is identified with one of the factors described in the argument factor
#of the data object, ARSyNseq estimates this effect and removes it by estimating the
#main PCs of the ANOVA effects associated. 
#Selected PCs will be those that explain more than the variability proportion 
#specified in Variability. 

norm_ARSyn <- ARSyNseq(noiseqData_norm_count,      #Biobases eSet object
                       factor = "sequencingBatch",    #when NULL, all factors are considered
                       batch = T,                  #TRUE if factor argument is batch info
                       norm = "n",                 #type of normalization, "n" if already normalized
                       logtransf = F)              #If F, log-transformation will be applied before ARSyn

#############################FINAL QUALITY CHECK#######################################################

#New PCA for ARSyn data corrected data

myPCA_ARSyn <- dat(norm_ARSyn,
                   type = "PCA",
                   norm = T,
                   logtransf = T)

#Plot post-ARSyn

png("postArsynPCA.png")
explo.plot(myPCA_ARSyn, samples = c(1,2),
           plottype = "scores", 
           factor = "sequencingBatch")
dev.off()

#Save PCA file

#saveRDS(myPCA, "/datos/rosmap/PCAs/PCA_post_Arsyn.rsd")

#PCA with prcomp --- ---

norm_ARSyn <- exprs(norm_ARSyn)
dim(norm_ARSyn)
#[1] 22070   499

#Build matrix

pca_norm_ARSyn <- norm_ARSyn %>% 
  t()  # transpose the matrix so that rows = samples and columns = variables, this because dots in the PCA scatterplot will be the ones in the rows

# Look at the first 10 rows and first 5 columns of the matrix
pca_norm_ARSyn[1:5, 1:5]

# Perform the PCA

pca_norm_ARSyn <- prcomp(pca_norm_ARSyn, retx = TRUE, center = TRUE, scale. = FALSE)

#PCA to table

pca_norm_ARSyn_df <- pca_norm_ARSyn$x %>% as.data.frame() 

# Create a data frame with PC number and percentage of variance

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

variance_table_norm_ARSyn <- data.frame(
  PC = 1:length(pca_norm_ARSyn$sdev),
  Variance_Percentage = pca_norm_ARSyn$sdev^2 / sum(pca_norm_ARSyn$sdev^2) * 100,
  cumulative_percentage = cumsum(pca_norm_ARSyn$sdev^2 / sum(pca_norm_ARSyn$sdev^2) * 100))

#Plot by color

PC1_PC2_norm_ARSyn_batch <- pca_norm_ARSyn_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2, colour = as.factor(factors$sequencingBatch)) +
  geom_point() +
  geom_text(mapping = aes(label = factors$specimenID)) +
  labs(title = "PCA Scatterplot coloured by library batch",
       subtitle = "PC1 vs PC2", 
       x = paste("PC1 (", sprintf("%.2f", variance_table_norm_ARSyn$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table_norm_ARSyn$Variance_Percentage[3]), "%)")) +
  theme_minimal()

PC1_PC2_norm_ARSyn_NIA <- pca_norm_ARSyn_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2, colour = as.factor(factors$dicho_NIA_reagan)) +
  geom_point() +
  geom_text(mapping = aes(label = factors$specimenID)) +
  labs(title = "PCA Scatterplot coloured by library batch",
       subtitle = "PC1 vs PC2", 
       x = paste("PC1 (", sprintf("%.2f", variance_table_norm_ARSyn$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table_norm_ARSyn$Variance_Percentage[3]), "%)")) +
  theme_minimal()

#Create new noiseq object with normalized counts 

noiseqData_final <- NOISeq::readData(norm_ARSyn,
                                     factors = factors,
                                     gc = mygc,
                                     biotype = mybiotype,
                                     length = mylength)

#Check for bias with renormalized

mycountsbio_final <- dat(noiseqData_final, 
                         type = "countsbio", 
                         norm=T, 
                         factor = NULL)

#Plot final plots

png("CountsFinal.png")
explo.plot(mycountsbio_final,
           plottype = "boxplot",
           samples = 1:15)
dev.off()

#Low counts 

png("lowcountsFinal.png")
explo.plot(mycountsbio_final, plottype = "barplot", samples = 1:15)
dev.off()

#Plot for Mvalues

png("MvaluesFinal.png")
explo.plot(mycd_norm,samples=sample(1:ncol(counts),10))
dev.off()

#calculate final GC bias

myGCcontent_final <- dat(noiseqData_final,
                         k = 0, 
                         type = "GCbias", 
                         factor = NULL,
                         norm = T)

#Plot final GC bias

png("GCbiasFinal.png",width=1000)
explo.plot(myGCcontent_final, plottype = "boxplot", samples = 1:5)
dev.off()

#calculate final length bias

mylenBias <- dat(noiseqData_final, 
                 k = 0, 
                 type = "lengthbias", 
                 factor = "sequencingBatch",
                 norm=T)

#Plot final length bias

png("lengthbiasFinal.png",width=1000)
explo.plot(mylenBias, samples = 1:5)
dev.off()

#Final cuts --- ---

#In the DLPFC PCAs, there are 2 samples that are not cute, "792_180530", "15_120410"

norm_ARSyn_depl <- norm_ARSyn[, !(colnames(norm_ARSyn) %in% c("792_130530", "15_120410"))]
dim(norm_ARSyn)
# Pre - [1] 22094   446
# Post - [1] 22070   499 <-

factors_depleted <- factors %>% filter(specimenID != "792_130530" & specimenID != "15_120410")
dim(factors_depleted)
#[1] 444   4
#[1] 497   6 <-

#Build matrix

pca_norm_ARSyn_depl <- norm_ARSyn_depl %>% 
  t()  # transpose the matrix so that rows = samples and columns = variables, this because dots in the PCA scatterplot will be the ones in the rows

# Look at the first 10 rows and first 5 columns of the matrix
pca_norm_ARSyn_depl[1:5, 1:5]

# Perform the PCA

pca_norm_ARSyn_depl <- prcomp(pca_norm_ARSyn_depl, retx = TRUE, center = TRUE, scale. = FALSE)

#PCA to table

pca_norm_ARSyn_depl_df <- pca_norm_ARSyn_depl$x %>% as.data.frame() %>% rownames_to_column(var = 'specimenID')

# Create a data frame with PC number and percentage of variance

#Note: The percentage of variance is calculated as the squared singular value
#of each PC divided by the sum of squared singular values, multiplied by 100.

variance_table_depl <- data.frame(
  PC = 1:length(pca_norm_ARSyn_depl$sdev),
  Variance_Percentage = pca_norm_ARSyn_depl$sdev^2 / sum(pca_norm_ARSyn_depl$sdev^2) * 100,
  cumulative_percentage = cumsum(pca_norm_ARSyn_depl$sdev^2 / sum(pca_norm_ARSyn_depl$sdev^2) * 100))

#Plot

PC1_PC2_norm_ARSyn_depl_library_batch <- pca_norm_ARSyn_depl_df %>% 
  ggplot() +
  aes(x = PC1, y = PC2, colour = as.factor(factors_depleted$dicho_NIA_reagan)) +
  geom_point() +
  geom_text(mapping = aes(label = specimenID)) +
  labs(title = "PCA Scatterplot coloured by library batch",
       subtitle = "PC1 vs PC2", 
       x = paste("PC1 (", sprintf("%.2f", variance_table_depl$Variance_Percentage[1]), "%)"),
       y = paste("PC2 (",  sprintf("%.2f", variance_table_depl$Variance_Percentage[3]), "%)")) +
  theme_minimal()


#Final count matrix --- ---

final_counts <- norm_ARSyn_depl %>% as.data.frame() %>%  mutate(feature = rownames(norm_ARSyn), .before = 1) #add it as column so I can save it like tsv
dim(final_counts)
#[1] 22094   445  #This means 447 specimenIDs with 22094 features
#[1] 22070   498? <-

#Final metadata --- ---

final_metadata <- data.frame(
  "specimenID" = colnames(final_counts)[-1])   

#Don't run this more than once!
final_metadata <- final_metadata %>% 
  left_join(metadata, by = "specimenID")
dim(final_metadata)
#[1] 446  41
#[1] 497  41 <-

#Finally, save counts table --- ---

vroom::vroom_write(final_counts, file = "/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/DLFPC/full_counts/ROSMAP_RNAseq_filteredQC_counts_DLPFC.txt")

#Save filtered metadata --- ---

vroom::vroom_write(final_metadata, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_metadata_filteredQC_DLPFC.txt")

#Save annotation --- --- 

vroom::vroom_write(myannot, file ="/datos/rosmap/data_by_counts/ROSMAP_counts/counts_by_tissue/metadata/DLPFC/RNA_seq_filteredQC_annotation_DLPFC.txt")

#END