#
#Script 1.1prepro-mRNA.R
#Script for annotation, bias detection and correction (QC) 
#of already normalized gene expression data
#Previous script is 1.MatchCountsandClinicalMetadata.R
#By paulinapglz.99@gmail.com

####################### PACKAGES ############################## 

#Here we use the NOISeq, edgeR and EDAseq packages for bias detection and correction

pacman::p_load('dplyr', 
               'biomaRt',
               'NOISeq',
               'edgeR', 
               'EDASeq')

######################## A. Get the data #####################

#Read counts data
#This file was generated in 1.MatchFPKMandClinicalMetadata.R

counts <- vroom::vroom(file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_250124.csv') %>%   #counts for all cogdx but 6
  as.data.frame()
dim(counts)
#[1] 55889   625  #original expression counts have 55889 genes and 625 specimenIDs

colnames(counts)[1] <-"ensembl_gene_id" #change the name for further filtering

############################## B. Annotation ##############################

### Generate annotation with ensembl ------ ------
#First we generate mart object

mart <- useEnsembl("ensembl",                         
                   dataset="hsapiens_gene_ensembl")

#We create myannot, with GC content, biotype, info for length & names per transcript

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values =  counts$ensembl_gene_id,  #annotate the genes in the count matrix 
                 mart = mart)

#Add length column

myannot$length <- abs(myannot$end_position-myannot$start_position)
dim(myannot)
#[1] 49399     7   #The full annotation has 49410 genes, there's a difference between annotated genes and genes in counts

#Filtering   ------ ------

#left join to further filtering

counts <- myannot %>% left_join(counts, 
                                    by = "ensembl_gene_id")
dim(counts)
#[1] 49399   631

#Filter to obtain only protein coding 

counts <- counts %>% 
  filter(gene_biotype == "protein_coding" & hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
  distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column
dim(counts)
#[1] 18848   631

#Obtain new annotation after filtering

myannot <- counts %>% 
  dplyr::select(1:7)
dim(myannot)
#[1] 18848     7

#Obtain counts 

expression_counts <- counts %>% 
  dplyr::select(ensembl_gene_id, 8:ncol(counts))      
dim(expression_counts)
#[1] 18848   625

#Scale and center data --- --- 

#If your data is not already scaled, scale(scale = T) if you don't need it
#just center with scale(center = T)

scaled_expression_counts <- scale(expression_counts[-1], scale = T, center = T)  %>% #omited first column as it contains only names
  as.data.frame()

#This centering and scaling together has the effect of making the column standard 
#deviations equal to 1. 

apply(scaled_expression_counts, 2, sd)

#Obtain factors from metadata --- --- 

metadata <- vroom::vroom(file = "/datos/rosmap/metadata/RNA_seq_metadata_080224.csv")
dim(metadata)
#[1] 624   14

#I do this to make sure the rowlength of factors match with the counts columns

factors <- data.frame(
  "specimenID" = colnames(expression_counts)[-1])   

factors <- factors %>% 
  left_join(metadata, by = "specimenID")
dim(factors)
#[1] 624   14 # this means 624 specimen_IDs 

############################## C. NOISeq object ##############################

#Give format to table for NOIseq purposes ------ ------

gene_names <- myannot$ensembl_gene_id

rownames(expression_counts) <- gene_names    #for the counts
rownames(scaled_expression_counts) <- gene_names #for the NOISeq object to PCA

#Names of features characteristics

mylength <- setNames(myannot$length, myannot$ensembl_gene_id)

mygc <- setNames(myannot$percentage_gene_gc_content, myannot$ensembl_gene_id)

mybiotype <-setNames(myannot$gene_biotype, myannot$ensembl_gene_id)


#NOTE: I create here two NOISeq objects, one is only for building a PCA, the other one
#will be used for bias detection and solving, as dat() does not accept centered and scaled data

#Create NOISeq object FOR PCA --- ---

noiseqData_PCA <- NOISeq::readData(data = scaled_expression_counts,
                               factors = factors,           #variables indicating the experimental group for each sample
                               gc = mygc,                   #%GC in myannot
                               biotype = mybiotype,         #biotype
                               length =  mylength)          #gene length

#check for batch effect

myPCA <- dat(noiseqData_PCA, 
             type = "PCA",   #PCA type of plot
             norm = T,       #indicating data are already normalized
             logtransf = T)  #if T, it doesn't perform logtrasnf

#Plot PCA
#In this PCA, each point is a specimenID, and PCs are explained by gene variance

png("PCA_Ori.png")
explo.plot(myPCA,
           plottype = "scores",
           factor = NULL)   #use all factors
dev.off()

#Save PCA file

#saveRDS(myPCA, "/datos/rosmap/PCAs/PCA_ori_RNAseq_ROSMAP_proteincoding090224.rsd")

#Create NOISeq object bias detect and bias correction --- ---

#For NOISeq, order of factors$specimenIDs and  colnames(expression_counts)[-1] must match
#> identical(colnames(expression_counts)[-1], factors$specimenID)
#[1] TRUE

noiseqData <- NOISeq::readData(data = expression_counts[-1],
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

mycd <- dat(noiseqData, type = "cd", norm = T) #slooooow

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 594_120522"

#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
# 6    617 
#Good but not perfect

#1)check expression bias for counts
#Use a same factor variable to all bias detection

mycountsbio <- dat(noiseqData, 
                   type =  "countsbio",  
                   norm = T,      #T when already normalized counts as input
                   factor = NULL) #When NULL, all factors are considered
#Plots

png("CountsOri.png")
explo.plot(mycountsbio,
           plottype = "boxplot", #type of plot
           samples = 1:15)  #only showing the first 15 samples
dev.off()

#2)check for low count genes

png("lowcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:15)
dev.off()

#3)check for transcript composition bias

#Plot for Mvalues

png("MvaluesOri.png")
explo.plot(mycd,samples=sample(1:ncol(expression_counts),10))
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
           samples = NULL,   
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

#################SOLVE BIASES###################################

#1) filter low count genes.
#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 

#This function renormalizes and filters features to only have the ones with
#CPM>1

countMatrixFiltered <- filtered.data(expression_counts[-1], 
                                     factor = "batch",       #using all factors
                                     norm = T,           #already normalized counts
                                     depth = NULL,
                                     method = 1, 
                                     cpm = 0, 
                                     p.adj = "fdr")

#Filtering out low count features...
#14951 features are to be kept for differential expression analysis with filtering method 1

#Filter again myannot to have only genes after filtering

myannot <- myannot %>%
  filter(ensembl_gene_id %in% rownames(countMatrixFiltered))
dim(myannot)
#[1] 14951     7

###Save annotation for later --- ---

vroom::vroom_write(myannot, file = "/datos/rosmap/metadata/ROSMAP_QC_fitlered_annotation100224.tsv")

#Again regenerate names vector 
gene_names <- myannot$ensembl_gene_id

##Create EDA object

mydataEDA <- newSeqExpressionSet(
  counts = as.matrix(countMatrixFiltered),
  featureData = data.frame(myannot,
                           row.names = myannot$ensembl_gene_id),
  phenoData = data.frame(factors,
                         row.names=factors$specimenID))

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
#[1] "Reference sample is: 594_120522"

#Table diagnostic

table(mycd_lessbias@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#  448    175  

#############################SOLVE BATCH EFFECT#######################################################
#If your data is already batch solved, only renormalize if last table diagnostic tells you to do it
# Normalization --- ---

#Use Uqua (UpperQuartile) for renormalization

UquaNorm <-NOISeq::uqua(normCounts(lFull),
                        long = 1000, 
                        lc = 0,
                        k = 0)

noiseqData_Uqua <- NOISeq::readData(data = UquaNorm, 
                                    factors = factors)

#cd has to preceed ARSyN or won't work

mycd_Uqua <- NOISeq::dat(noiseqData_Uqua,
                         type="cd",
                         norm=TRUE)
#[1] "Warning: 94 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 594_120522"

#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd_Uqua@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#26    597

#ARSyNseq for batch effect solution --- ---
#Only do this if your data is not batch solved

#When batch is identified with one of the factors described in the argument factor
#of the data object, ARSyNseq estimates this effect and removes it by estimating the
#main PCs of the ANOVA effects associated. 
#Selected PCs will be those that explain more than the variability proportion 
#specified in Variability. 

norm_ARSyn <- ARSyNseq(noiseqData_Uqua,              #Biobases eSet object
                       factor = "batch",   #when NULL, all factors are considered
                       batch = T,      #TRUE if factor argument is batch info
                       norm = "n",            #type of normalization, "n" if already normalized
                       logtransf = F)      #If F, log-transformation will be applied before ARSyn

#New PCA for ARSyn data corrected data

myPCA_ARSyn <- dat(norm_ARSyn,
                   type = "PCA",
                   norm = T,
                   logtransf = T)

#Plot post-ARSyn

png("postArsynPCA.png")
explo.plot(myPCA_ARSyn, samples = c(1,2),
           plottype = "scores", 
           factor = "cogdx")
dev.off()

#Save PCA file

#saveRDS(myPCA, "/datos/rosmap/PCAs/PCA_post_Arsyn.rsd")

#############################FINAL QUALITY CHECK#######################################################

#Names of features characteristics to add to final data

mylength <- setNames(myannot$length, myannot$ensembl_gene_id)

mygc <- setNames(myannot$percentage_gene_gc_content, myannot$ensembl_gene_id)

mybiotype <-setNames(myannot$gene_biotype, myannot$ensembl_gene_id)


#Create new noiseq object with re-normalized counts 

noiseqData_final <- NOISeq::readData(exprs(noiseqData_Uqua),
                             gc = mygc,
                             biotype = mybiotype,
                             factors = factors,
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
explo.plot(mycountsbio, plottype = "barplot", samples = 1:15)
dev.off()

#calculate final GC bias

myGCcontent_final <- dat(noiseqData_final,
                         k = 0, 
                         type = "GCbias", 
                         factor = "cogdx",
                         norm = T)

#Plot final GC bias

png("GCbiasFinal.png",width=1000)
explo.plot(myGCcontent_final, plottype = "boxplot", samples = 1:5)
dev.off()

#calculate final length bias

mylenBias <- dat(noiseqData_final, 
                 k = 0, 
                 type = "lengthbias", 
                 factor = "cogdx",
                 norm=T)

#Plot final length bias

png("lengthbiasFinal.png",width=1000)
explo.plot(mylenBias, samples = 1:5)
dev.off()

#Save new count matriz --- ---

final_counts <- exprs(noiseqData_final) %>% as.data.frame() %>% 
  mutate(ensembl_gene_id = rownames(exprs(noiseqData_final)), .before = 1) #add it as column so I can save it like tsv
dim(final_counts)
#[1] 14951   624  #This means 624 specimenIDs with 14951 features

#Finally, save table
#vroom::vroom_write(final_counts, file = "/datos/rosmap/FPKM_data/ROSMAP_QCed_count_matrixfiltered_090224.tsv")
#END --- ---