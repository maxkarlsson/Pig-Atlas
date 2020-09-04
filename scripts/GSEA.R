get_GO_max <- 
  function(gene_id, 
           species,
           evidence.codes = c("EXP", "IDA", "IPI", "IMP", 
                              "IGI", "IEP", "TAS", "IC"), 
           ontologies = c("biological_process", 
                          "molecular_function", 
                          "cellular_component")) {
    require(biomaRt)
    
    ensembl <- 
      useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
              dataset = paste0(species, "_gene_ensembl"), 
              host = "http://apr2018.archive.ensembl.org")
    GOMap <- 
      getBM(mart = ensembl, 
            attributes = c("ensembl_gene_id", 
                           "go_id", 
                           "go_linkage_type", 
                           "namespace_1003", 
                           "name_1006"),
            filters = "ensembl_gene_id", 
            values = gene_id)
    
    subGOMap <- 
      GOMap[(GOMap$go_linkage_type %in% evidence.codes) & 
              (GOMap$namespace_1003 %in% ontologies), ]
    
    return(subGOMap)
    
  }

run_XGSA <- 
  function(species1, species2,
           genes1, genes2,
           universe1, universe2,
           ontologies = c("biological_process", 
                          "molecular_function", 
                          "cellular_component")) {
    
    XGSA_data1 <- 
      new_XGSA_dataset(species = species1, 
                       data = genes1, 
                       type = 'genesetlist', 
                       name = species1, 
                       universe = universe1)
    XGSA_data2 <- 
      new_XGSA_dataset(species = species2, 
                       data = genes2, 
                       type = 'genesetlist', 
                       name = species2, 
                       universe = universe2)
    
    GO_db <- get_GO_max(gene_id = universe2, 
                        species = species2, 
                        ontologies = ontologies)
    
    GO_dict <- 
      GO_db %>% 
      select(go_id, term = name_1006, ontology = namespace_1003) %>% 
      distinct()
    
    GO_list <- 
      GO_db %$%
      split(ensembl_gene_id, go_id) %>%
      collapseGO()
    
    GO_data <- 
      new_XGSA_dataset(species = species2, 
                       data = GO_list, 
                       type = 'genesetlist', 
                       name = "GO", 
                       universe = unique(unlist(GO_list)))
    
    XGSA_res1 <- run_XGSA_test(XGSA_data1, GO_data, max = 1e5)
    XGSA_res2 <- run_XGSA_test(XGSA_data2, GO_data, max = 1e5)
    
    bind_rows(lapply(XGSA_res1, function(X){ X[["pvals"]] }) %>%
                unlist() %>%
                enframe("GO_", "p_val") %>%
                mutate(species = species1),
              lapply(XGSA_res2, function(X){ X[["pvals"]] }) %>%
                unlist() %>%
                enframe("GO_", "p_val") %>%
                mutate(species = species2)) %>% 
      
      group_by(species, gene_set) %>%
      separate(GO_, into = c("gene_set", "GO"), sep = "\\.") %>%
      mutate(GO = gsub("genes.", "", GO),
             p_adj = p.adjust(p_val, method = "BH")) %>%
      ungroup() %>%
      left_join(GO_dict,
                by = c("GO" = "go_id")) %>%
      select(species, gene_set, GO, term, ontology, p_val, p_adj) %>%
      arrange(species, gene_set, p_adj) 
    
  }



library(xgsa)
library(tidyverse)
library(magrittr)
rename <- dplyr::rename
select <- dplyr::select

pig_gene_comparison_classification <- 
  read_delim("data/HPA/pig_compare_category_piggenome_92.tsv", delim = "\t") %>%
  rename(enssscg_id = ensg_id) %>%
  mutate(enhanced_tissues = gsub(", ", "_", enhanced_tissues),
         enhanced_tissues = gsub(",", ";", enhanced_tissues),
         enhanced_tissues = gsub("_", ", ", enhanced_tissues),
         enhanced_tissues = trimws(enhanced_tissues))

human_gene_comparison_classification <- 
  read_delim("data/HPA/tissue_compare_category_92.tsv", delim = "\t") %>%
  mutate(enhanced_tissues = gsub(", ", "_", enhanced_tissues),
         enhanced_tissues = gsub(",", ";", enhanced_tissues),
         enhanced_tissues = gsub("_", ", ", enhanced_tissues),
         enhanced_tissues = trimws(enhanced_tissues))


genlista1 <- 
  pig_gene_comparison_classification %>% 
  filter(!is.na(enhanced_tissues)) %>%
  separate_rows(enhanced_tissues, sep = ";") %$%
  split(enssscg_id, enhanced_tissues)


genlista2 <-
  human_gene_comparison_classification %>% 
  filter(!is.na(enhanced_tissues)) %>%
  separate_rows(enhanced_tissues, sep = ";") %$%
  split(ensg_id, enhanced_tissues)





XGSA_elevated_res <- 
  run_XGSA(species1 = "sscrofa", 
           species2 = "hsapiens",
           genes1 = genlista1, 
           genes2 = genlista2,
           universe1 = pig_gene_comparison_classification$enssscg_id, 
           universe2 = human_gene_comparison_classification$ensg_id,
           ontologies = c("biological_process", 
                          "molecular_function", 
                          "cellular_component"))






