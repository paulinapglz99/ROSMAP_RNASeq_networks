#Explore metadata

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
                "cogdx", "dcfdx_lv") %>% 
  
cor_metadata <- metadata[3:11]

cor_metadata <- cor_metadata%>%
  filter(!is.na(cogdx) & !is.na(apoe_genotype))

C<- cor(cor_metadata)

corrplot(C)
