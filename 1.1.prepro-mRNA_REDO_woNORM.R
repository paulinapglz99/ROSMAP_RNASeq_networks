#Script for annotation, bias detection and correction (QC) 
#of already normalized gene expression data

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

expression <- vroom::vroom(file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_250124.csv') %>%   #counts for all cogdx but 6
  as.data.frame()
dim(expression)
#[1] 55889   625  #original expression counts have 55889 genes and 625 specimenIDs

colnames(expression)[1] <-"ensembl_gene_id" #change the name for further filtering

############################## B. Annotation ##############################

#Generate annotation with ensembl.
#First we generate mart object

mart <- useEnsembl("ensembl",                         
                   dataset="hsapiens_gene_ensembl")

#We create myannot, with GC content, biotype, info for length & names per transcript

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values =  expression$ensembl_gene_id,  #annotate the genes in the count matrix 
                 mart = mart)

#Add length column

myannot$length <- abs(myannot$end_position-myannot$start_position)
dim(myannot)
#[1] 49399     7   #The full annotation has 49410 genes, there's a difference between annotated genes and genes in counts

#left join to further filtering

expression <- myannot %>% left_join(expression, 
                                    by = "ensembl_gene_id")
dim(expression)
#[1] 49399   631

#Filter to obtain only protein coding 

expression <- expression %>% 
  filter(gene_biotype == "protein_coding" & hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
  distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column

#Obtain new annotation after filtering

myannot <- expression %>% 
  dplyr::select(1:7)
dim(myannot)
#[1] 18848     7

#Obtain counts 

expression_counts <- expression %>% 
  dplyr::select(ensembl_gene_id, 8:ncol(expression))      
dim(expression_counts)
#[1] 18848   625

############################## C. NOISeq object ##############################

#Give format to table for NOIseq purposes

rownames(expression_counts) <- expression_counts$ensembl_gene_id

#Obtain factors from metadata

metadata <- vroom::vroom(file = "/datos/rosmap/metadata/RNA_seq_metadata_250124.csv") %>% 
  dplyr::select(specimenID, cogdx, ceradsc, braaksc)
dim(metadata)
#[1] 624   4

factors <- data.frame(
  "specimenID" = colnames(expression_counts)[-1])   

factors <- factors %>% 
  left_join(metadata, by = "specimenID")
dim(factors)
#[1] 624   4  # this means 624 specimen_IDs and only one factor

#Let's explore the metadata
library(ggplot2)

ggplot(factors, aes(x = factor(cogdx))) +
  geom_bar(fill = "#6495ed", color = "black") +
  labs(title = "Histograma de cogdx",
       x = "Cogdx",
       y = "Frecuencia") +
  theme_bw()
  
ggplot(factors, aes(x = factor(ceradsc))) +
  geom_bar(fill = "#6495ed", color = "black") +
  labs(title = "Histograma de ceradsc",
       x = "ceradsc",
       y = "Frecuencia") +
  theme_bw()


ggplot(factors, aes(x = factor(braaksc))) +
  geom_bar(fill = "#6495ed", color = "black") +
  labs(title = "Histograma de braaksc",
       x = "braaksc",
       y = "Frecuencia") +
  theme_bw()


#For NOISeq, order of factors$specimenIDs and  colnames(expression_counts)[-1] must match

#Names of features characteristics

mylength <- setNames(myannot$length, myannot$ensembl_gene_id)

mygc <- setNames(myannot$percentage_gene_gc_content, myannot$ensembl_gene_id)

mybiotype <-setNames(myannot$gene_biotype, myannot$ensembl_gene_id)

#For NOISeq, order of factors$specimenIDs and  colnames(expression_counts)[-1] must match
#> identical(colnames(expression_counts)[-1], factors$specimenID)
#[1] TRUE

noiseqData <- NOISeq::readData(data = expression_counts[-1],#not using 1st col 
                               factors = factors,           #variables indicating the experimental group for each sample
                               gc = mygc,                   #%GC in myannot
                               biotype = mybiotype,         #biotype
                               length =  mylength)          #gene length

#1)check expression bias per subtype

#Obtain the likely counts of genes, organized by subtype,
#from the noiseqData object

mycountsbio <- dat(noiseqData, 
                   type =  "countsbio",  
                   norm = T, 
                   factor = NULL)
#Plots

#patients with repeated measures
png("CountsOri.png")
explo.plot(mycountsbio, plottype = "boxplot", samples = 1:10)
dev.off()

#2)check for low count genes

png("lowcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:10)
dev.off()

#Histogram of row means

#esto no corre por alguna razon

png("lowCountThres.png")
hist(rowMeans(cpm(expression_counts,log=T)),
     ylab="genes",
     xlab="mean of log CPM",  #does this has any sense? it computes CPM values
     col="gray")
abline(v=0,col="red")
dev.off()

#3)check for transcript composition bias

#each sample "s" is compared to a reference "r" (which can be arbitrarily chosen).
#by computing M values=log2(counts = countsr). 

#Confidence intervals for the M median is computed by bootstrapping.
#If the median of M values for each comparison is not in the CI, the deviation
# of the sample is significant, therefore, normalization is needed 

mycd <- dat(noiseqData, type = "cd", norm = T) #slooooow

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 594_120522"

#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])

#When using real data

#FAILED PASSED 
# 6    617

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
                   factor = "cogdx")

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
                    factor = "cogdx")

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."

#Plot length bias

png("lengthbiasOri.png")
explo.plot(mylengthbias, 
           samples = NULL, 
           toplot = "global")
dev.off()

#5) check for batch effect

myPCA <- dat(noiseqData,
             type = "PCA", 
             norm = T,
             logtransf = F)

#Plot PCA

png("PCA_Ori.png")
explo.plot(myPCA, samples = c(1,2),
           plottype = "scores",
           factor = "cogdx")
dev.off()


#################SOLVE BIASES###################################

#1) filter low count genes.
#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 

#why? que me esta quitando? como calcula los CPM si ya estan normalizados 

countMatrixFiltered <- filtered.data(expression_counts[-1], 
                                     factor = "cogdx",
                                     norm = T, 
                                     depth = NULL,
                                     method = 1, 
                                     cpm = 0, 
                                     p.adj = "fdr")

#Filtering out low count features...
#14951 features are to be kept for differential expression analysis with filtering method 1

#Filter again myannot to have only genes after filtering

myannot <- myannot %>%
  filter(ensembl_gene_id %in% rownames(countMatrixFiltered))

##Create EDA object

mydataEDA <- newSeqExpressionSet(
  counts = as.matrix(countMatrixFiltered),
  featureData = data.frame(myannot,
                           row.names = myannot$ensembl_gene_id),
  phenoData = data.frame(factors,
                         row.names=factors$specimenID))

#order for less bias

#for gc content
gcFull <- withinLaneNormalization(mydataEDA, 
                                  "percentage_gene_gc_content",
                                  which = "full")#corrects GC bias 

#for length
lFull <- withinLaneNormalization(gcFull, 
                                 "length", 
                                 which = "full")#corrects length bias 

#cd has to preceed ARSyN or won't work

mycd_lessbias <- NOISeq::dat(lFull,
                         type = "cd",
                         norm = TRUE)

#[1] "Reference sample is: 594_120522"

table(mycd_lessbias@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#455    168 


#############################SOLVE BATCH EFFECT#######################################################

myPCAp_preARSyn <- dat(lFull, 
                       type = "PCA", 
                       norm = T, 
                       logtransf = F)

#plot preArsyn PCA

png("preArsyn.png")
explo.plot(myPCAp_preARSyn, samples = c(1,2),
           plottype = "scores",
           factor = "cogdx")
dev.off()

##################### TRY normalization ##################### 

#For Uqua

UquaNorm <-NOISeq::uqua(normCounts(lFull),    #OPCION A USAR TMM EN VEZ DE UQUA
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

table(mycd_Uqua@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#26    597

#ARSyNseq for batch effect solution

#When batch is identified with one of the factors described in the argument factor
#of the data object, ARSyNseq estimates this effect and removes it by estimating the
#main PCs of the ANOVA effects associated. 
#Selected PCs will be those that explain more than the variability proportion 
#specified in Variability. 

norm_ARSyn <- ARSyNseq(noiseqData_Uqua,              #Biobases eSet object
                       factor = "cogdx",   #when NULL, all factors are considered
                       batch = FALSE,      #TRUE if factor argument is batch info
                       norm = "n",            #type of normalization, "n" if already normalized
                       logtransf = F)      #If F, log-transformation will be applied before ARSyn

#New PCA for ARSyn data

myPCA_ARSyn <- dat(norm_ARSyn,
                   type = "PCA",
                   norm = T,
                   logtransf = T)

#Plot post-ARSyn

png("postArsyn.png")
explo.plot(myPCA_ARSyn, samples = c(1,2),
           plottype = "scores", 
           factor = "cogdx")
dev.off()


#############################FINAL QUALITY CHECK#######################################################

###Names of features characteristics to add to final data

mylength <- setNames(myannot$length, myannot$ensembl_gene_id)

mygc <- setNames(myannot$percentage_gene_gc_content, myannot$ensembl_gene_id)

mybiotype <-setNames(myannot$gene_biotype, myannot$ensembl_gene_id)

#Create new noiseq object with re-normalized counts

noiseqData_final <- readData(data = exprs(norm_ARSyn),
                             gc = mygc,
                             biotype = mybiotype,
                             factor = factors,
                             length = mylength)

mycountsbio_final <- dat(noiseqData_final, 
                         type = "countsbio", 
                         factor = "cogdx",
                         norm=T)

png("CountsFinal.png")
explo.plot(mycountsbio_final,
           plottype = "boxplot",
           samples=1:5)   #this doesnot run still
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

#Finally, save table
write.table(final,"filtered_FPKM_matrix_250124.tsv",sep='\t',quote=F)
#duplicates share everything except the plate
#Finally, save table
write.table(final,"RNAseqnormalized.tsv",sep='\t',quote=F)

#END
