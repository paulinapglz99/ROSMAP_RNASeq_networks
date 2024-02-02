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
                "race", "spanish", "apoe_genotype", "age_at_visit_max",    
                 "age_death", "braaksc", "ceradsc",
                "cogdx", "dcfdx_lv")

cor_metadata<- metadata[3:13]

cor(cor_metadata)

