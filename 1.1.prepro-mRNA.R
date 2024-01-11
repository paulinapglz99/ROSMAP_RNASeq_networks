#Script for preprocessing (annot) and QC of gene expression data with NOISeq package

#paulinapglz.99@gmail.com

#Here we perform:

## Data preparation
## Quality Control & bias removal

pacman::p_load('dplyr', 
               'biomaRt')

############################## A. Get the data ############################## 

#Read normalized data
#This file was generated in 1.MatchFPKMandClinicalMetadata.R

FPKM <- vroom::vroom(file = '/datos/rosmap/FPKM_data/filtered_FPKM_matrix_new161223.csv') #counts for cogdx = 1, 2, 3, 4 and 5
dim(FPKM)
#[1] 55889   623  #original expression counts have 55889 genes and 623 specimenIDs

############################## B. Annotation ##############################

#Generate annotation with ensembl. Annotate gene_biotype, GC content

mart <- useEnsembl("ensembl",
                   dataset="hsapiens_gene_ensembl")

#annnotate GC content, length & biotype per transcript

myannot <- getBM(attributes = c("ensembl_gene_id", 
                                "percentage_gene_gc_content", "gene_biotype",
                                "start_position","end_position","hgnc_id","hgnc_symbol"),
                 filters = "ensembl_gene_id", 
                 values =  FPKM$gene_id,
                 mart = mart)

#Add lenght column

myannot$length <- abs(myannot$end_position-myannot$start_position)
dim(myannot)
#[1] 49410     8   #The full annotation has 49410 genes, there's a difference between annotated genes and genes in counts

#filter transcripts without annotation

myannot <- myannot %>% 
  dplyr::filter(gene_biotype == "protein_coding" & hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
  distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column

#Sort myannot's ensembl_gene_id column
myannot <- myannot[order(myannot$ensembl_gene_id), ]

dim(myannot)
#[[1] 18848     8  #18848 full annotated protein coding genes with 8 characteristics, difference of 30562 annotated genes

#Filter count matrix to have only values that match ensembl_gene_id in the myannot dataframe in a new object

FPKM_exprots <- FPKM %>%
  filter(FPKM$gene_id %in% myannot$ensembl_gene_id) %>% 
  as.data.frame() 

#Sorts the gene_id column of FPKM_exprots based on the sorted values of ensembl_gene_id:

FPKM_exprots <- FPKM_exprots[match(myannot$ensembl_gene_id, FPKM_exprots$gene_id), ] 

dim(FPKM_exprots)
#[1] 18848   623  #18848 protein coding genes from 622 specimenIDs, now annotation and gene counts match in order

#read metadata
RNA_seq_metadata <- vroom::vroom(file = "/datos/rosmap/metadata/ROSMAP_filtered_metadata_forRNAseq.csv")

#I will delete 2 non-necessary columns for now and make sure it converts to data.frame 

RNA_seq_metadata <- RNA_seq_metadata %>% 
  dplyr::select(-individualID, -msex) %>% 
  as.data.frame()

dim(RNA_seq_metadata)
#[1] 622   3  #622 individualIDs and specimenIDs and 3 characteristics


################## C. CHECK BIASES ########################################################

pacman::p_load("NOISeq", 
               "edgeR")

#The lengths of RNA_seq_metadata$specimenID and colnames(FPKM_exprots) must be the same for the noiseq::readData function
#They actually do at this point [622]

#for the factor format
#the order of the elements of the factor must coincide with the order of the samples (columns)
# in the expression data provided. Row number in factor must match with number of cols in data ("FPKM_exprots"). 

# Get indexes to reorder RNA_seq_metadata$specimenID
specimenID_positions <- match(colnames(FPKM_exprots)[-1], RNA_seq_metadata$specimenID)

# Use indexes to reorder RNA_seq_metadata$specimenID
RNA_seq_metadata$specimenID <- RNA_seq_metadata$specimenID[specimenID_positions]
names(RNA_seq_metadata$specimenID) <- colnames(FPKM_exprots)[-1] # sorted according to the order of colnames(FPKM_exprots)

#Assign gene_id rownames to rows in count data

rownames(FPKM_exprots) <- FPKM_exprots$gene_id

#And delete column

FPKM_exprots <- FPKM_exprots[,-1]

#NOISeq needs annotation in t() (don't really get why, but this works like this)

myannot <- t(myannot)
colnames(myannot) <- myannot[1,]

#Convert to a NOIseq object

noiseqData <- readData(data = FPKM_exprots, 
                       gc = myannot["percentage_gene_gc_content",],  #%GC in myannot
                       biotype = myannot["gene_biotype",],          #biotype
                       factors = RNA_seq_metadata,                 #variables indicating the experimental group for each sample
                       length =  myannot["length",])               #gene length

###########################HASTA AQUI VOY ###########################

#1)check expression bias per subtype

#Obtain the likely counts of genes, organized by subtype,
#from the noiseqData object


mycountsbio <- dat(noiseqData, 
                  type =  "countsbio",
                  norm = T, 
                  factor = NULL)
#Plot 

#patients with repeated measures
png("CountsOri.png")
explo.plot(mycountsbio, plottype = "boxplot", samples = 1:5)
dev.off()

#2)check for low count genes
png("lowcountsOri.png")
explo.plot(mycountsbio, plottype = "barplot", samples = 1:5)
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

#each sample s is compared to a reference r (which can be arbitrarily chosen).
#by computing M values=log2(counts = countsr). 

#Confidence intervals for the M median is computed by bootstrapping.
#If the median of M values for each comparison is not in the CI, the deviation
# of the sample is significant, therefore, normalization is needed 

mycd <- dat(noiseqData, type = "cd", norm = T) #slooooow

#[1] "Warning: 110 features with 0 counts in all samples are to be removed for this analysis."
#[1] "Reference sample is: 525_120515"

#[1] "Diagnostic test: FAILED. Normalization is required to correct this bias."

table(mycd@dat$DiagnosticTest[,  "Diagnostic Test"])

#FAILED PASSED 
# 11    610 

#Plot for Mvalues

png("MvaluesOri.png")
explo.plot(mycd,samples=sample(1:ncol(FPKM_exprots),10))
dev.off()


