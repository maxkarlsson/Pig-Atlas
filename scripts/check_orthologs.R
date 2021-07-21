

library(tidyverse)


gene_orthologs <- 
  read_csv("data/processed/gene_orthologs.csv")

all_gene_orthologs <- read_tsv("data/meta/Ensembl_orthologs_all.txt")

gene_orthologs %>% 
  filter(enssscg_id == "ENSSSCG00000040593") %>% View
all_gene_orthologs %>% 
  filter(`Gene stable ID` == "ENSSSCG00000040593")  %>% 
  View


gene_orthologs %>% 
  filter(enssscg_id == "ENSSSCG00000017296") %>%
  select(1, 2, 3, 4, 6)

all_gene_orthologs %>% 
  filter(`Gene stable ID` == "ENSSSCG00000017296")  %>% View
  select(1, 2, 3, 4, 5)
