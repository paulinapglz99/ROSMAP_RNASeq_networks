#Script for preprocessing (annot) and QC of gene expression data with NOISeq package
#paulinapglz.99@gmail.com
#Here we perform:
## Data preparation for NOISeq bias identification
## Quality Control & bias removal

####################### PACKAGES ############################## 

pacman::p_load('dplyr', 
               'biomaRt',
               'NOISeq')

######################## A. Get the data #####################

#Read counts data
#This file was generated in 1.MatchFPKMandClinicalMetadata.R

expression <- vroom::vroom(file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_new161223.csv') #counts for cogdx = 1, 2, 3, 4 and 5
dim(expression)
#[1] 55889   623  #original expression counts have 55889 genes and 623 specimenIDs

colnames(expression)[1] <-"ensembl_gene_id"

######################## B. Filter data ######################

#CPM=(counts/fragments sequenced)*one million.
#Filtering those genes with average CPM below 1, would be different
#to filtering by those with average counts below 1. 

countMatrixFiltered <- filtered.data(as.matrix(expression),
                                     norm = T, 
                                     depth = NULL,
                                     method = 1,   #method 1 (CPM)
                                     cpm = 1, 
                                     p.adj = "fdr")

############################## B. Annotation ##############################

#Generate annotation with ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

#annnotate GC content, length & biotype per transcript

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values =  expression$ensembl_gene_id,
                 mart = mart)

#Add length column

myannot$length <- abs(myannot$end_position-myannot$start_position)
dim(myannot)
#[1] 49399     7   #The full annotation has 49410 genes, there's a difference between annotated genes and genes in counts

#left join to further filtering

expression <- myannot %>% left_join(expression, 
                                    by = "ensembl_gene_id")
dim(expression)
#[1] 49399   629

#Filter to obtain only protein coding 

expression <- expression %>% 
  filter(gene_biotype == "protein_coding" & hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
      distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column

#Obtain new annotation after filtering

myannot <- expression %>% 
  dplyr::select(1:7)

#Obtain counts 

expression_counts <- expression %>% 
  dplyr::select(ensembl_gene_id, 8:ncol(expression))      
dim(expression_counts)
#[1] 18848   623

#Obtain factors      

factors <- data.frame(
  "specimen_ID" = colnames(expression_counts)[-1],
  "group" = 1)
dim(factors)
#[1] 622   2 # this means 621 specimen_IDs and only one factor

#for the factor format
#the order of the elements of the factor must coincide with the order of the samples (columns)
# in the expression data provided. Row number in factor must match with number of cols in data ("FPKM_exprots"). 

noiseqData <- readData(data = expression_counts[-1], 
                       gc = myannot$percentage_gene_gc_content,  #%GC in myannot
                       biotype = myannot$gene_biotype,          #biotype
                       factors = factors,                 #variables indicating the experimental group for each sample
                       length =  myannot$length)               #gene length

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
hist(rowMeans(cpm(FPKM_exprots,log=T)),
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

mycd_table <- table(mycd@dat)

table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#12    609 

#Plot for Mvalues

png("MvaluesOri.png")
explo.plot(mycd,samples=sample(1:ncol(FPKM_exprots),10))
dev.off()

#4)check for length & GC bias

#A cubic spline regression model is fitted. Both the model p-value and the coefficient
# of determination (R2) are shown. If the model p-value is significant and R2 value is
# high (more than 70%), the expression depends on the feature

myGCcontent <- dat(noiseqData,
                   k = 0, type = "GCbias",
                   factor = "ceradsc")

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."
#[1] "GC content bias detection is to be computed for:"
#[1] "1" "2" "3" "4"

#Residuals:
#  Min      1Q  Median      3Q     Max 
#-9.5011 -1.5303 -0.0613  1.3300  6.4901 

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
                    factor = "ceradsc")


#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Length bias detection information is to be computed for:"
#[1] "1" "2" "3" "4"

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
             norm = F,
             logtransf = F)

#Plot PCA

png("PCA_Ori.png")
explo.plot(myPCA, samples = c(1,2),
           plottype = "scores",
           factor = "ceradsc")
dev.off()
