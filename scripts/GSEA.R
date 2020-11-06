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
    
    GO_db1 <- get_GO_max(gene_id = universe1, 
                        species = species1, 
                        ontologies = ontologies)
    
    GO_dict1 <- 
      GO_db1 %>% 
      select(go_id, term = name_1006, ontology = namespace_1003) %>% 
      distinct()
    
    GO_list1 <- 
      GO_db1 %$%
      split(ensembl_gene_id, go_id) %>%
      collapseGO()
    
    GO_data1 <- 
      new_XGSA_dataset(species = species1, 
                       data = GO_list1, 
                       type = 'genesetlist', 
                       name = "GO", 
                       universe = unique(unlist(GO_list1)))
    
    ###
    
    
    XGSA_data2 <- 
      new_XGSA_dataset(species = species2, 
                       data = genes2, 
                       type = 'genesetlist', 
                       name = species2, 
                       universe = universe2)
    
    GO_db2 <- get_GO_max(gene_id = universe2, 
                        species = species2, 
                        ontologies = ontologies)
    
    GO_dict2 <- 
      GO_db2 %>% 
      select(go_id, term = name_1006, ontology = namespace_1003) %>% 
      distinct()
    
    GO_list2 <- 
      GO_db2 %$%
      split(ensembl_gene_id, go_id) %>%
      collapseGO()
    
    GO_data2 <- 
      new_XGSA_dataset(species = species2, 
                       data = GO_list2, 
                       type = 'genesetlist', 
                       name = "GO", 
                       universe = unique(unlist(GO_list2)))
    
    GO_dict <- 
      bind_rows(GO_dict1,
                GO_dict2) %>%
      distinct()
    ###
    
    XGSA_res1 <- run_XGSA_test(XGSA_data1, GO_data2, max = 1e5)
    XGSA_res2 <- run_XGSA_test(XGSA_data2, GO_data1, max = 1e5)
    
    bind_rows(lapply(XGSA_res1, function(X){ X[["pvals"]] }) %>%
                unlist() %>%
                enframe("GO_", "p_val") %>%
                mutate(species = species1),
              lapply(XGSA_res2, function(X){ X[["pvals"]] }) %>%
                unlist() %>%
                enframe("GO_", "p_val") %>%
                mutate(species = species2)) %>% 
      
      
      separate(GO_, into = c("gene_set", "GO"), sep = "\\.") %>%
      group_by(species, gene_set) %>%
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
library(enrichR)

library(rrvgo)
library(treemapify)
library(pbapply)

rename <- dplyr::rename
select <- dplyr::select
source("scripts/functions_utility.R")

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

load("data/processed/gene_overlap_comparison.Rdata")

gene_mapping <- 
  read_delim("data/meta/ensembl_pig_gene.tab", delim = "\t") %>% 
  filter(biotype == "protein_coding")
gene_chromosome <- read_delim("data/meta/ensembl92_gene_chromosome.txt", delim = "\t")

human_gene_mapping <- read_delim("data/meta/human/geninfo_92.tsv", delim = "\t")

gene_orthologs_all <- 
  read_delim("data/meta/Ensembl_orthologs_all.txt", delim = "\t") %>% 
  select(enssscg_id = 1,
         ensg_id = 3, 
         pig_gene_name = 2,
         human_gene_name = 4,
         ortholog_confidence = 10, 
         ortholog_type = 5,
         gene_order_conservation_score = 8, 
         whole_genome_alignment_coverage = 9, 
         sequence_identity = 6) %>%
  filter(!is.na(ortholog_confidence)) %>% 
  filter(ensg_id %in% human_gene_mapping$ensg_id) %>%
  filter(enssscg_id %in% gene_mapping$enssscg_id)

load("data/processed/pig_intraconsensus_anova2.Rdata")

intraconsens_gene_enrich_overlap_all <-
  read_csv("data/processed/intraconsensus results including elevated.csv")

# ------ OVerlap specificity GSEA -------

spec_overlap_enrichr_data <- 
  gene_overlap_data_all_orth %>%
  filter(comparison_id != "") %>%
  filter(spec_category_human %in% c("Group enriched", "Tissue enriched") |
           spec_category_pig %in% c("Group enriched", "Tissue enriched")) %>%
  # filter(comparison_id %in% c("kidney", "urinary bladder")) %>%
  select(comparison_id, enssscg_id, ensg_id, type) %>% 
  left_join(human_gene_mapping %>% 
              select(1:2))


if(file.exists("data/processed/GSEA_spec_overlap_enrichr_data.Rdata")) {
  load("data/processed/GSEA_spec_overlap_enrichr_data.Rdata")
} else {
  
  
  
  spec_overlap_enrichr_res <- 
    spec_overlap_enrichr_data %>% 
    group_by(comparison_id, type) %>%
    do({
      enrichr(unique(.$gene_name),
              databases = c("GO_Biological_Process_2018")) %>%
        bind_rows(.id = "database")
    }) %>%
    as_tibble() %>% 
    mutate(p_adj = p.adjust(P.value, method = "BH")) %>% 
    
    mutate(GO = str_extract(Term, "\\(.*\\)") %>%
             gsub("\\(|\\)", "", .)) %>%
    separate(Overlap, into = c("n_GO", "n_max"), sep = "/") 
  
  save(spec_overlap_enrichr_data,
       file = "data/processed/GSEA_spec_overlap_enrichr_data.Rdata")
  
}


  
  
spec_overlap_enrichr_res %>% 
  mutate(p_adj = p.adjust(P.value, method = "BH")) %>% 
  filter(p_adj < 0.05) %>%
  mutate(log_p = -log10(p_adj),
         log_p_dir = log_p * case_when(type == "Human" ~ 1,
                                       type == "Pig" ~ -1)) %>% 
  ggplot(aes(log_p, log10(Combined.Score))) +
  geom_hex(aes(fill = stat(log10(count))),
           bins = 100) +
  scale_fill_viridis_c()



plot_data_all <- 
  spec_overlap_enrichr_res %>%
  filter(p_adj < 0.05) %>%
  select(comparison_id, GO, n_GO, P = p_adj, type) %>% 
  left_join(spec_overlap_enrichr_data %>% 
              group_by(comparison_id, type) %>% 
              count()) %>%
  mutate(score = as.numeric(n_GO) / n)
  


plot_settings <- 
  plot_data_all %$% 
  expand.grid(type = c("Overlap", "Human", "Pig"),
              comparison_id = unique(comparison_id)) 

plots <- 
  pblapply(1:nrow(plot_settings),
           function(i) {
             cat(paste(i, "\n"))
             
             
             plot_data <- 
               plot_data_all %>% 
               filter(type == plot_settings$type[i]) %>%
               filter(comparison_id == plot_settings$comparison_id[i])
             
             if(nrow(plot_data) > 1) {
               plot_mat <- 
                 calculateSimMatrix(plot_data$GO, 
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel") 
               plot_mat2 <- 
                 reduceSimMatrix(plot_mat, 
                                 # setNames(-log10(plot_data$P), plot_data$GO),
                                 setNames(plot_data$score, plot_data$GO),
                                 threshold=0.7,
                                 orgdb="org.Hs.eg.db") 
               
               g <- 
                 plot_mat2 %>% 
                 as_tibble() %>%
                 group_by(parentTerm) %>%
                 mutate(n = row_number()) %>%
                 ggplot(aes(area = size, subgroup = parentTerm)) +
                 
                 geom_treemap(aes(fill = parentTerm,
                                  alpha = n), 
                              color = "black",
                              show.legend = F) +
                 geom_treemap_subgroup_border(color = "black") +
                 
                 geom_treemap_text(aes(label = term),
                                   colour = "black", 
                                   place = "centre",
                                   alpha = 0.4,
                                   grow = TRUE) +
                 
                 geom_treemap_subgroup_text(place = "centre", 
                                            grow = T, 
                                            reflow = T,
                                            alpha = 1, 
                                            colour = "white", 
                                            fontface = "bold",
                                            min.size = 0) +
                 scale_alpha_continuous(range = c(0.8, 1)) +
                 scale_fill_manual(values = rep(ggthemes::tableau_color_pal(palette = "Hue Circle", 
                                                                            type = c("regular"), 
                                                                            direction = 1)(19),
                                                10)) +
                 theme_void() +
                 ggtitle(paste(plot_settings$type[i], plot_settings$comparison_id[i]))
             } else {
               g <- 
                 ggplot() +
                 theme_void() +
                 ggtitle(paste(plot_settings$type[i], plot_settings$comparison_id[i]))
             }
           g
               
           })


pdf(savepath("Overlap GSEA treemaps.pdf"), width = 8, height = 5)
plots
# plots[rep(1:32, each = 2) + rep(c(0, 32), 32)]
dev.off()


spec_overlap_enrichr_res %>% 
  
  mutate(log_p = -log10(p_adj)) %>%
  select(1, 2, Term, GO, log_p) %>%
  
  spread(type, log_p, fill = 0) %>%
  filter(Human > -log10(0.05) |
           Pig > -log10(0.05)) %>%
  ggplot(aes(Human, Pig)) +
  geom_hline(yintercept = -log10(0.05)) +
  geom_vline(xintercept = -log10(0.05)) +
  geom_point() +
  # geom_hex(aes(fill = stat(log10(count))),
  #          bins = 50) +
  
  scale_fill_viridis_c() +
  facet_wrap(~comparison_id, scales = "free")


spec_overlap_enrichr_res %>%
  filter(p_adj < 0.05) %>%
  select(1, 2, Term, GO, p_adj, Genes) %>%
  write_csv(savepath("Overlap enrichr results.csv"))

# ------ Intraconsensus GSEA -------


intracon_overlap_enrichr_data <- 
  intraconsens_gene_enrich_overlap_all %>% 
  mutate(overlap_type = case_when(overlap == "not variable - elevated" ~ "Overlap",
                                  T ~ "Variable")) %>%
  select(1, 2, 4, overlap_type) %>% 
  inner_join(gene_orthologs %>% 
               select(1,2, 4)) %>%
  separate_rows(high_tissues, sep = ";")


if(file.exists("data/processed/GSEA_intracon_overlap_enrichr_res.Rdata")) {
  load("data/processed/GSEA_intracon_overlap_enrichr_res.Rdata")
} else {
  
  
  
  intracon_overlap_enrichr_res <- 
    intracon_overlap_enrichr_data %>% 
    group_by(consensus_tissue_name, overlap_type, high_tissues) %>%
    do({
      enrichr(unique(.$human_gene_name),
              databases = c("GO_Biological_Process_2018")) %>%
        bind_rows(.id = "database")
    }) %>%
    as_tibble() %>% 
    mutate(p_adj = p.adjust(P.value, method = "BH")) %>% 
    
    mutate(GO = str_extract(Term, "\\(.*\\)") %>%
             gsub("\\(|\\)", "", .)) %>%
    separate(Overlap, into = c("n_GO", "n_max"), sep = "/") 
  
  save(intracon_overlap_enrichr_res,
       file = "data/processed/GSEA_intracon_overlap_enrichr_res.Rdata")
  
}



plot_data_all <- 
  intracon_overlap_enrichr_res %>%
  filter(p_adj < 0.05) %>%
  select(consensus_tissue_name, type = overlap_type, high_tissues, GO, n_GO, P = p_adj) %>% 
  left_join(intracon_overlap_enrichr_data %>% 
              group_by(consensus_tissue_name, type = overlap_type) %>% 
              count()) %>%
  mutate(score = as.numeric(n_GO) / n,
         high_tissues = ifelse(is.na(high_tissues),
                               "",
                               high_tissues))



plot_settings <- 
  plot_data_all %>%
  select(1:3) %>% 
  distinct() 

plots <- 
  pblapply(1:nrow(plot_settings),
           function(i) {
             cat(paste(i, "\n"))
             
             
             plot_data <- 
               plot_data_all %>% 
               filter(type == plot_settings$type[i]) %>%
               filter(consensus_tissue_name == plot_settings$consensus_tissue_name[i]) %>%
               filter(high_tissues == plot_settings$high_tissues[i]) 
             
             if(nrow(plot_data) > 1) {
               plot_mat <- 
                 calculateSimMatrix(plot_data$GO, 
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method="Rel") 
               plot_mat2 <- 
                 reduceSimMatrix(plot_mat, 
                                 # setNames(-log10(plot_data$P), plot_data$GO),
                                 setNames(plot_data$score, plot_data$GO),
                                 threshold=0.7,
                                 orgdb="org.Hs.eg.db") 
               
               g <- 
                 plot_mat2 %>% 
                 as_tibble() %>%
                 group_by(parentTerm) %>%
                 mutate(n = row_number()) %>%
                 ggplot(aes(area = size, subgroup = parentTerm)) +
                 
                 geom_treemap(aes(fill = parentTerm,
                                  alpha = n), 
                              color = "black",
                              show.legend = F) +
                 geom_treemap_subgroup_border(color = "black") +
                 
                 geom_treemap_text(aes(label = term),
                                   colour = "black", 
                                   place = "centre",
                                   alpha = 0.4,
                                   grow = TRUE) +
                 
                 geom_treemap_subgroup_text(place = "centre", 
                                            grow = T, 
                                            reflow = T,
                                            alpha = 1, 
                                            colour = "white", 
                                            fontface = "bold",
                                            min.size = 0) +
                 scale_alpha_continuous(range = c(0.8, 1)) +
                 scale_fill_manual(values = rep(ggthemes::tableau_color_pal(palette = "Hue Circle", 
                                                                            type = c("regular"), 
                                                                            direction = 1)(19),
                                                10)) +
                 theme_void() +
                 ggtitle(paste(plot_settings$type[i], 
                               plot_settings$consensus_tissue_name[i],
                               plot_settings$high_tissues[i]))
             } else {
               g <- 
                 ggplot() +
                 theme_void() +
                 ggtitle(paste(plot_settings$type[i], 
                               plot_settings$consensus_tissue_name[i],
                               plot_settings$high_tissues))
             }
             g
             
           })


pdf(savepath("Intracon GSEA treemaps.pdf"), width = 8, height = 5)
plots
# plots[rep(1:32, each = 2) + rep(c(0, 32), 32)]
dev.off()

intracon_overlap_enrichr_res %>%
  filter(p_adj < 0.05) %>%
  select(1, 2, high_tissues, Term, GO, p_adj, Genes) %>%
  write_csv(savepath("Intraconsensus enrichr results.csv"))

######################


genlista1 <- 
  pig_gene_comparison_classification %>% 
  filter(!is.na(enhanced_tissues)) %>%
  separate_rows(enhanced_tissues, sep = ";") %>%
  filter(enhanced_tissues %in% c("kidney", "urinary bladder")) %$%
  split(enssscg_id, enhanced_tissues)


genlista2 <-
  human_gene_comparison_classification %>% 
  filter(!is.na(enhanced_tissues)) %>%
  separate_rows(enhanced_tissues, sep = ";") %>%
  filter(enhanced_tissues %in% c("kidney", "urinary bladder")) %$%
  split(ensg_id, enhanced_tissues)




t <- Sys.time()
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
t1 <- Sys.time() - t
t1

XGSA_elevated_res %>%
  filter(p_adj < 0.05) %>% View


