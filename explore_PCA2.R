#Script que hace un PCA con NOISeq
#Experimental

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

### Generate annotation with ensembl ------ ------
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

#Filtering   ------ ------

#left join to further filtering

expression <- myannot %>% left_join(expression, 
                                    by = "ensembl_gene_id")
dim(expression)
#[1] 49399   631

#Filter to obtain only protein coding and delete repeated genes

expression <- expression %>% 
  #filter(gene_biotype == "protein_coding" & hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
  distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column
dim(expression)
#[1] 18848   631

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
#Obtain factors from metadata

metadata <- vroom::vroom(file = "/datos/rosmap/metadata/cli_bio_metadata.csv") %>% 
  dplyr::select(specimenID, msex, educ, apoe_genotype, race,cogdx, ceradsc, braaksc, exclude) %>% 
  filter(specimenID %in% colnames(expression_counts[-1]))
dim(metadata)
#[1] 624   4

#I do this to make sure the rowlength of factors match with the counts columns

factors <- data.frame(
  "specimenID" = colnames(expression_counts)[-1])   

factors <- factors %>%
  left_join(metadata, by = "specimenID") %>% 
  distinct()
dim(factors)
#[1] 624   4 # this means 624 specimen_IDs 

############################## C. NOISeq object ##############################

#Give format to table for NOIseq purposes ------ ------

rownames(expression_counts) <- expression_counts$ensembl_gene_id

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

#When using real data

#FAILED PASSED 
# 6    617 
#Good but not perfect

#1)check expression bias per subtype

#Obtain the likely counts of genes, organized by subtype,
#from the noiseqData object

#Use a same factor variable to all bias detection

mycountsbio <- dat(noiseqData, 
                   type =  "countsbio",  
                   norm = T,      #T when already normalized counts as input
                   factor = NULL) #When NULL, all factors
#Plots

png("CountsOri.png")
explo.plot(mycountsbio, plottype = "boxplot", samples = 1:10)
dev.off()

#2)check for low count genes

png("lowcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:10)
dev.off()

#Histogram of row means

#Warning, this histogram does not run properly
png("lowCountThres.png")
hist(rowMeans(cpm(as.matrix(expression_counts,log=T))), #renormalizes for CPM values
     ylab="genes",
     xlab="mean of log CPM",  
     col="gray")
abline(v=0,col="red")
dev.off()

#Trying w/ggplot

ggplot(data = data.frame(logCPM = rowMeans(cpm(as.matrix(expression_counts), log = TRUE))),
       aes(x = logCPM)) +
  geom_histogram(binwidth = 0.1, fill = "gray", color = "white") +
  labs(x = "Mean of log CPM", y = "Genes") +
  theme_minimal() +
  geom_vline(xintercept = 0, color = "red")

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

png("lengthbiasOri.png", width=1000)
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

#Save PCA file

#saveRDS(myPCA, "/datos/rosmap/PCAs/PCA_ori_preArsyn.rsd")