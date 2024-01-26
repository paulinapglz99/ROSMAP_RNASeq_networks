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

factors <- vroom::vroom(file = "/datos/rosmap/metadata/RNA_seq_metadata_250124.csv") %>% 
  dplyr::select(specimenID, cogdx, ceradsc, braaksc)
dim(factors)
#[1] 624   4

#Redim factors to match with count data

factors <- factors %>% 
  filter(specimenID %in% colnames(expression_counts))
dim(factors)
#[1] 624   4  # this means 624 specimen_IDs and only one factor

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

png("lowCountThres.png")
hist(rowMeans(cpm(expression_counts,log=T)),  #esto no corre por alguna razon
     ylab="genes",
     xlab="mean of log CPM",
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
#[1] "Reference sample is: 525_120515"

#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])

#When only 1-type-factor 
#FAILED PASSED 
#10    611 

#When there's two categorical variants in factor

#FAILED PASSED 
#14    607 

#When using real data



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
                   factor = "group")

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."
#[1] "GC content bias detection is to be computed for:"
#[1] "1" "2"

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-9.5011 -1.5303 -0.0613  1.3300  6.4901 

#Residual standard error: 2.972 on 81 degrees of freedom
#Multiple R-squared:  0.8643,	Adjusted R-squared:  0.8459 
#F-statistic: 46.92 on 11 and 81 DF,  p-value: < 2.2e-16

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
                    factor = "group")

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Length bias detection information is to be computed for:"
#[1] "1"

#Residual standard error: 10.56 on 82 degrees of freedom
#Multiple R-squared:  0.2717,	Adjusted R-squared:  0.1829 
#F-statistic:  3.06 on 10 and 82 DF,  p-value: 0.002402

#Residuals:
# Min      1Q  Median      3Q     Max 
#-29.152  -2.099  -0.478   1.663  83.628 

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
           factor = "group")
dev.off()


#################SOLVE BIASES###################################

#1) filter low count genes.
#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 

countMatrixFiltered <- filtered.data(expression_counts[-1], 
                                     factor = "group",
                                     norm = T, 
                                     depth = NULL,
                                     method = 1, 
                                     cpm = 0, 
                                     p.adj = "fdr")
#why? que me esta quitando?

#Filtering out low count features...
#14952 features are to be kept for differential expression analysis with filtering method 1

#Filter again myannot to have only 

myannot <- myannot %>%
  filter(ensembl_gene_id %in% rownames(countMatrixFiltered))

##Create EDA object

#all names must match

mydataEDA <- newSeqExpressionSet(
  counts = as.matrix(countMatrixFiltered),
  featureData = data.frame(myannot,
                           row.names = myannot$ensembl_gene_id),
  phenoData = data.frame(factors,
                         row.names=factors$specimen_ID))

#order for less bias

#for gc content
gcFull <- withinLaneNormalization(mydataEDA, 
                                  "percentage_gene_gc_content",
                                  which = "full")#corrects GC bias 

#for length
lFull <- withinLaneNormalization(gcFull, 
                                 "length", 
                                 which = "full")#corrects length bias 
