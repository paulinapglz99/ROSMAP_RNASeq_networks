#
#3.differential_expression.R
#Script that makes differential expression of data 
#By paulinapglz.99@gmail.com
#For further details of how NOISeq works, go to https://www.bioconductor.org/packages/release/bioc/vignettes/NOISeq/inst/doc/NOISeq.pdf

#Libraries --- ---

pacman::p_load('dplyr', 
               'biomaRt',
               'NOISeq',
               'edgeR', 
               'EDASeq', 
               "ggplot2")