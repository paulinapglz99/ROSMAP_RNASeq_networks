#Explore metadata
#paulinapglz.99@gmail.com

pacman::p_load(corrplot, 
              dplyr)

#Read metadata

metadata <- vroom::vroom(file = '/datos/rosmap/metadata/cli_bio_metadata.csv')  

#Filter to make a corplot

metadata <- metadata %>% 
  filter(assay == 'rnaSeq', 
         organ == 'brain') %>% 
  dplyr::select("individualID", "specimenID","msex", "educ",                
                "race", "spanish", "apoe_genotype",    
                 "braaksc", "ceradsc",
                "cogdx", "dcfdx_lv")

#Filter metadata just for corplot purposes  

cor_metadata <- metadata[3:11]

cor_metadata <- cor_metadata%>%
  filter(!is.na(cogdx) & !is.na(apoe_genotype))  #Here your variables with NAs

#Correlation table

C<- cor(cor_metadata)

#Corrplots

corrplot(C, method = 'number')

corrplot(C, method = 'square', order = 'FPC', type = 'lower', diag = FALSE)
