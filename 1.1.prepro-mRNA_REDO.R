#Script for preprocessing, annot and QC of gene expression data
#Here we use the NOISeq, edgeR and EDAseq packages
#paulinapglz.99@gmail.com
#Here we perform:
## Data preparation for NOISeq bias identification
## Quality Control & bias removal

####################### PACKAGES ############################## 

pacman::p_load('dplyr', 
               'biomaRt',
               'NOISeq',
               'edgeR')

######################## A. Get the data #####################

#Read counts data
#This file was generated in 1.MatchFPKMandClinicalMetadata.R

expression <- vroom::vroom(file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_new161223.csv') %>%   #counts for cogdx = 1, 2, 3, 4 and 5
                           as.data.frame()
dim(expression)
#[1] 55889   623  #original expression counts have 55889 genes and 623 specimenIDs

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
#[1] 49399   629

#Filter to obtain only protein coding 

expression <- expression %>% 
  filter(gene_biotype == "protein_coding" & hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
      distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column

#Obtain new annotation after filtering

myannot <- expression %>% 
  dplyr::select(1:7)

############################## C. NOISeq object ##############################

#Obtain counts 

expression_counts <- expression %>% 
  dplyr::select(ensembl_gene_id, 8:ncol(expression))      
dim(expression_counts)
#[1] 18848   623

#Give format to table for NOIseq purposes

rownames(expression_counts) <- expression_counts$ensembl_gene_id

#Obtain factors      

factors <- data.frame(
  "specimen_ID" = colnames(expression_counts)[-1],
  "group" = 1)    #simulated factors
factors$group <- sample(c(1, 2), size = nrow(factors), replace = TRUE)  #Simulation to give random numbers 1 and 2 to the df de factors
dim(factors)
#[1] 622   2 # this means 621 specimen_IDs and only one factor

#Names of features characteristics

mylength <- setNames(myannot$length, myannot$ensembl_gene_id)

mygc <- setNames(myannot$percentage_gene_gc_content, myannot$ensembl_gene_id)

mybiotype <-setNames(myannot$gene_biotype, myannot$ensembl_gene_id)

#the order of the elements of the factor must coincide with the order of the samples (columns)
# in the expression data provided. Row number in factor must match with number of cols in data ("FPKM_exprots"). 

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

library(EDASeq)

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

#lFull ya es un objeto S3? con los datos corregidos por length y GC, podria volver a hacer el 
#diagnostico con NOISeq

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

#[1] "Warning: 95 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 525_120515"
#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd_Uqua@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#10    611 son los mismos?

#With FPKM normalization

RPKMNorm <- NOISeq::rpkm(normCounts(lFull),
                      long = 1000, 
                      lc = 0,
                      k = 0)


noiseqData_RPKM <- NOISeq::readData(data = RPKMNorm, 
                                    factors= factors)

#cd has to preceed ARSyN or won't work

mycd_RPKM <- NOISeq::dat(noiseqData_RPKM,
                         type="cd",
                         norm=TRUE)

#[1] "Warning: 95 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 525_120515"
#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd_RPKM@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
# 582     39 

#With TMM normalization

TMMNorm <-NOISeq::tmm(normCounts(lFull),    #OPCION A USAR TMM EN VEZ DE UQUA
                        long = 1000, 
                        lc = 0,
                        k = 0)

noiseqData_TMM <- NOISeq::readData(data = TMMNorm, 
                                    factors = factors)

#cd has to preceed ARSyN or won't work

mycd_TMM <- NOISeq::dat(noiseqData_TMM,
                         type="cd",
                         norm=TRUE)

#[1] "Warning: 95 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 525_120515"
#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd_TMM@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
#564     57 

#At this moment, any of the normalization makes data pass the Diagnostic Test

#############################SOLVE BATCH EFFECT#######################################################

myPCAp_preARSyn <- dat(noiseqData_Uqua, 
                      type = "PCA", 
                      norm = T, 
                      logtransf = F)

#plot preArsyn PCA

png("preArsyn.png")
explo.plot(myPCAp_preARSyn, samples = c(1,2),
           plottype = "scores",
           factor = "group")
dev.off()

#ARSyNseq for batch effect solution

#################ESTO AUN NO CORRE###################

#When batch is identified with one of the factors described in the argument factor
#of the data object, ARSyNseq estimates this effect and removes it by estimating the
#main PCs of the ANOVA effects associated. 
#Selected PCs will be those that explain more than the variability proportion 
#specified in Variability. 

norm_ARSyn <- ARSyNseq(noiseqData_Uqua,     #Biobases eSet object
                       factor = "group",  #when NULL, all factors are considered
                       batch = FALSE,      #TRUE if factor argument is batch info
                       norm = "n",     #type of normalization, "n" if already normalized
                       logtransf = F)  #If F, log-transformation will be applied before ARSyn

#ERROR, when factor = NULL
#Error in apply(sub, 2, mean) : dim(X) must have a positive length
#nor uqua, TMM or RPKM normalizations run

#I also tried changing the normalization type, and giving a different noiseqData (the normalized one)

#New PCA for ARSyn data

myPCA_ARSyn <- dat(norm_ARSyn,
             type = "PCA",
             norm = T,
             logtransf = T)

#Plot post-ARSyn

png("postArsyn.png")
explo.plot(myPCA_ARSyn, samples = c(1,2),
           plottype = "scores", 
           factor = "groups")
dev.off()

#############################FINAL QUALITY CHECK#######################################################

#Create new noiseq object with re-normalized counts
noiseqData_final <- readData(data = exprs(norm_ARSyn),
                      gc = myannot[,1:2],
                      biotype = myannot[,c(1,3)],
                      factor=designExp,
                      length=myannot[,c(1,8)])


mycountsbio_final <- dat(noiseqData, 
                   type = "countsbio", 
                   factor = "group",
                  norm=T)

png("CountsFinal.png")
explo.plot(mycountsbio_final,
           plottype = "boxplot",
           samples=1:5)
dev.off()

#calculate final GC bias

myGCcontent_final <- dat(noiseqData,
                         k = 0, 
                         type = "GCbias", 
                         factor = "group",
                         norm=T)

#Plot final GC bias

png("GCbiasFinal.png",width=1000)
par(mfrow=c(1,5))
sapply(1:5,function(x) explo.plot(myGCcontent_final, samples = x))
dev.off()

#calculate final length bias

mylenBias <- dat(noiseqData, k = 0, type = "lengthbias", 
                 factor = "group",
                 norm=T)

#Plot final length bias

png("lengthbiasFinal.png",width=1000)
par(mfrow=c(1,5))
sapply(1:5,function(x) explo.plot(mylenBias, samples = x))
dev.off()

############ ESTO NO LO EMPIEZO AUN ############ 
#############################RESOLVE DUPLICATES & SAVE##################################################
#get duplicates
i <- factors$specimen_ID %>% 
  filter(duplicated(.))

#get sample barcode per sample

i <- lapply(i,function(x) factors$specimen_ID [factors$specimen_ID == x])

#separate duplicates
final <- exprs(norm_ARSyn)
duplis <- final[,colnames(final)%in%unlist(i)]
prefi <- final[,!colnames(final)%in%unlist(i)]

#average duplicates
temp <- do.call(cbind,lapply(i,function(x) 
  rowMeans(duplis[,colnames(duplis)%in%x])))

#identify samples with barcode 
colnames(temp) <- factors$specimen_ID[duplicated(factors$specimen_ID)]
colnames(prefi) <-substr(colnames(prefi),1,19)

#joint matrices
final <- cbind(prefi,temp)
dim(final)
#[1] 17077   805

final <- final %>%
  select(order(match(colnames(.), subtype$samples)))

#Finally, save table
write.table(final,"RNAseqnormalized.tsv",sep='\t',quote=F)
#duplicates share everything except the plate
#Finally, save table
write.table(final,"RNAseqnormalized.tsv",sep='\t',quote=F)
#duplicates share everything except the plate
