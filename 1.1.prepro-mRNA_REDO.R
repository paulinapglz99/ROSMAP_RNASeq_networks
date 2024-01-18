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

#Obtain counts 

expression_counts <- expression %>% 
  dplyr::select(ensembl_gene_id, 9:ncol(expression))      
dim(expression_counts)
#[1] 18848   622

#Obtain new annotation after filtering



#Obtain factors      

factors <- data.frame(
  "specimen_ID" = colnames(expression_counts)[-1],
  "group" = 1)
dim(factors)

#[1] 621   2 # this means 621 specimen_IDs and only one factor


#for the factor format
#the order of the elements of the factor must coincide with the order of the samples (columns)
# in the expression data provided. Row number in factor must match with number of cols in data ("FPKM_exprots"). 

noiseqData <- readData(data = expression_counts, 
                       gc = myannot["percentage_gene_gc_content",],  #%GC in myannot
                       biotype = myannot["gene_biotype",],          #biotype
                       factors = RNA_seq_metadata,                 #variables indicating the experimental group for each sample
                       length =  myannot["length",])               #gene length

