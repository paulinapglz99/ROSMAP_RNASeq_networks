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
  dplyr::filter(gene_biotype == "protein_coding" | hgnc_symbol!="") %>% #only rows where gene_biotype is "protein_coding" and hgnc_symbol is not an empty string 
  distinct(ensembl_gene_id, .keep_all = TRUE) # Keeps only unique rows based on the ensembl_gene_id column

#Sort myannot's ensembl_gene_id column
myannot <- myannot[order(myannot$ensembl_gene_id), ]

dim(myannot)
#[1] 37963     8  #37963 full annotated protein coding genes with 8 characteristics, difference of 11447 annotated genes :c

#Filter count matrix to have only values that match ensembl_gene_id in the myannot dataframe.

FPKM_exprots <- FPKM %>%
  filter(FPKM$gene_id %in% myannot$ensembl_gene_id)

#Ordena la columna gene_id de FPKM_exprots en base a los valores ordenados de ensembl_gene_id:
FPKM_exprots <- FPKM_exprots[match(myannot$ensembl_gene_id, FPKM_exprots$gene_id), ] 

dim(FPKM_exprots)
#[1] 37963   623  #37963 protein coding genes from 622 specimenIDs, now annotation and gene counts match

#read metadata
RNA_seq_metadata <- vroom::vroom(file = "/datos/rosmap/metadata/ROSMAP_filtered_metadata_forRNAseq.csv")

#I will delete 2 non-necessary columns for now

RNA_seq_metadata$individualID <- NULL
RNA_seq_metadata$msex <- NULL
dim(RNA_seq_metadata)
#[1] 622   3  #622 individualIDs and specimenIDs and 3 characteristics


################## C. CHECK BIASES ########################################################

pacman::p_load("NOISeq", 
               "edgeR")

#assign design experiment format

#The lengths of RNA_seq_metadata$specimenID and colnames(FPKM_exprots) must be the same for the noiseq::readData function
#They actually do at this point [622]

#for the factor format
#the order of the elements of the factor must coincide with the order of the samples (columns)
# in the expression data provided.

#note: row number in factor must match with number of cols in data ("FPKM_exprots"). 

MAY DELETE 

# Get indexes to reorder RNA_seq_metadata$specimenID
#specimenID_positions <- match(colnames(FPKM_exprots)[-1], RNA_seq_metadata$specimenID)

# Use indexes to reorder RNA_seq_metadata$specimenID
RNA_seq_metadata$specimenID <- RNA_seq_metadata$specimenID[specimenID_positions]
names(RNA_seq_metadata$specimenID) <- colnames(FPKM_exprots)[-1] # sorted according to the order of colnames(FPKM_exprots)

#

FPKM_exprots <- as.data.frame(FPKM_exprots)

rownames(FPKM_exprots) <- FPKM_exprots$gene_id

FPKM_exprots <- FPKM_exprots[,-1]

RNA_seq_metadata <- as.data.frame(RNA_seq_metadata)

myannot <- t(myannot)
colnames(myannot) <- myannot[1,]

noiseqData <- readData(data = FPKM_exprots, 
                       gc = myannot["percentage_gene_gc_content",],  #%GC in myannot
                       biotype = myannot["gene_biotype",],          #biotype
                       factors = RNA_seq_metadata,                 #variables indicating the experimental group for each sample
                       length =  myannot["length",])               #gene length


#Meticulously verify that the sample names in FPKM_exprots (excluding the first column), RNA_seq_metadata$specimenID,
#and any other relevant data frames/objects are identical.

#1)check expression bias per subtype

mycountsbio = dat(noiseqData,
                  type = "countsbio",
                  factor = "subtype")

#patients with repeated measures

png("CountsOri.png")

explo.plot(mycountsbio, plottype = "boxplot",samples = 1:5)

##Next script is 2.dis_mat_coexpre.R
