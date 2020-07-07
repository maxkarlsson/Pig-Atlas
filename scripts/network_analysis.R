
library(tidyverse)
library(magrittr)
library(ggraph)
library(igraph)
library(tidygraph)
library(patchwork)
library(pheatmap)
library(propr)
library(ggrepel)
source("scripts/theme.R")
source("scripts/functions_utility.R")
source("scripts/functions_graphics.R")
source("scripts/functions_network.R")

tissue_colors <- read_delim("data/meta/pig_tissue_colors.tsv", delim = "\t")


tissue_colors_palette_full <- 
  tissue_colors %$%
  tibble(color = c(tissue_color, region_tissue_color, consensus_tissue_color, organ_color),
         tissue = c(tissue_name, region_tissue_name, consensus_tissue_name, organ_name)) %>%
  unique() %$%
  c(setNames(color, tissue),
    setNames(color, str_to_sentence(tissue)))

xeno_ko_genes <- 
  c("GGTA1" = "ENSSSCG00000005518", 
    "CMAH" = "ENSSSCG00000001099", 
    "B4GALNT2" = "ENSSSCG00000040942",
    "VWF" = "ENSSSCG00000000712", 
    "ULBP1" = "ENSSSCG00000026042") %>% 
  enframe("gene_name", "enssscg_id") %>% 
  mutate(name_id = paste(gene_name, enssscg_id))

gene_mapping <- 
  read_delim("data/meta/ensembl_pig_gene.tab", delim = "\t") %>% 
  filter(biotype == "protein_coding")

pig_gene_classification <- 
  read_delim("data/HPA/pig_consensus_category_piggenome_92.tsv", delim = "\t") %>%
  rename(enssscg_id = ensg_id) 

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

tissue_mapping <- read_delim("data/meta/tissue_mapping.tsv", delim = "\t")

human_gene_class_blood <-
  read_delim("data/meta/human/bloodcells_hpa_category_92.tsv", delim = "\t") 


human_gene_class_bloodlineage <-
  read_delim("data/meta/human/bloodcells_hpa_regional_category_92.tsv", delim = "\t") 


human_gene_class <- 
  read_delim("data/human data/lims/consensus_all_category_92.tsv", delim = "\t") %>% 
  mutate(enhanced_tissues = gsub(", ", "_", enhanced_tissues),
         enhanced_tissues = gsub(",", ";", enhanced_tissues),
         enhanced_tissues = gsub("_", ", ", enhanced_tissues),
         enhanced_tissues = trimws(enhanced_tissues))

if(file.exists("data/processed/normalized_res.Rdata")) {
  load("data/processed/normalized_res.Rdata")
  
} else {
  source("scripts/normalization.R", local = T)
  load("data/processed/normalized_res.Rdata")
}

##Correlation

# load("data/processed/pig_gene_spearman_propr.Rdata")


  
###Proportionality



# load("data/processed/pig_gene_spearman_propr.Rdata")

if(file.exists("data/processed/pig_gene_propr_filtered.Rdata")) {
  load("data/processed/pig_gene_propr_filtered.Rdata")
  
} else {
  
  if(file.exists("data/processed/pig_gene_propr.Rdata")) {
    load("data/processed/pig_gene_propr_res.Rdata")
    
  } else {
    pig_gene_propr <- 
      pig_atlas_sample %>%
      select(1, 2, ptpm) %>% 
      spread(enssscg_id, ptpm) %>%
      column_to_rownames("sample_ID") %>% 
      {. + 1} %>%
      propr()
    
    pig_gene_propr_res <- 
      pig_gene_propr@matrix
    save(pig_gene_propr,
         file = "data/processed/pig_gene_propr.Rdata")
    save(pig_gene_propr_res,
         file = "data/processed/pig_gene_propr_res.Rdata")
    rm(pig_gene_propr)
    
  } 
  
  pig_gene_propr_res_05filtered <- 
    pig_gene_propr_res %>% 
    {temp <- .; temp[!upper.tri(temp)] <- NA; temp} %>%
    as_tibble(rownames = "enssscg_id1") %>%
    gather(enssscg_id2, prop, -1, na.rm = T) %>% 
    filter(abs(prop) > 0.5)
  
  save(pig_gene_propr_res_05filtered,
       file = "data/processed/pig_gene_propr_filtered.Rdata")
  
} 


# ------ Expression network 1 ------

if(file.exists("data/processed/pig_coexpression_network.Rdata")) {
  
  load("data/processed/pig_coexpression_network.Rdata")
} else {
  prop_filter <- 0.7
  
  test_filters <- 
    20:38/40
  
  net_n_per_filter <- 
    lapply(test_filters,
           function(filtr) {
             edges_ <- 
               pig_gene_propr_res_05filtered %>% 
               filter(prop > filtr)
             
             nodes_ <- 
               tibble(node = c(edges_$enssscg_id1, edges_$enssscg_id2)) %>% 
               distinct()
             
             tibble(edges = nrow(edges_),
                    nodes = nrow(nodes_))
           }) %>% 
    set_names(test_filters) %>%
    bind_rows(.id = "filter") %>% 
    mutate(ratio = edges/nodes)
  
  net_n_per_filter %>%
    gather(type, n, -1) %>%
    mutate(filter = as.numeric(filter)) %>%
    ggplot(aes(filter, n, color = type)) +
    geom_line() + 
    facet_wrap(~type, scales = "free_y") + 
    stripped_theme_facet
  ggsave(savepath("Cutoff metrics.pdf"), width = 6, height = 3)
  
  # xeno_ko_genes
  # pig_gene_propr_res_05filtered %>%
  #   filter(prop > prop_filter) %>%
  #   filter(enssscg_id1 %in% xeno_ko_genes$enssscg_id |
  #            enssscg_id2 %in% xeno_ko_genes$enssscg_id) -> A
  # 
  # A
  # 
  # A %>%
  #   ggplot(aes(prop)) +
  #   geom_density()
  # A %>%
  #   filter(enssscg_id1 %in% xeno_ko_genes$enssscg_id) %>%
  #   group_by(enssscg_id1) %>%
  #   summarise(max_prop = max(prop),
  #             mean_prop = mean(prop),
  #             median_prop = median(prop),
  #             n = n()) %>%
  #   arrange(max_prop)
  # 
  # A %>%
  #   filter(enssscg_id2 %in% xeno_ko_genes$enssscg_id) %>%
  #   group_by(enssscg_id2) %>%
  #   summarise(max_prop = max(prop),
  #             mean_prop = mean(prop),
  #             median_prop = median(prop),
  #             n = n()) %>%
  #   arrange(max_prop)
  
  zero_varying_genes <- 
    pig_atlas_sample_norm %>% 
    group_by(enssscg_id) %>% 
    summarise(sd = sd(ptpm), 
              mean = mean(ptpm)) %>% 
    arrange(sd) %>%
    filter(sd == 0)
  
  gene_coexpression_network <- 
    pig_gene_propr_res_05filtered %>% 
    filter(prop > prop_filter) %>%
    # filter(prop > 0.95) %>%
    filter(!(enssscg_id1 %in% zero_varying_genes$enssscg_id | 
               enssscg_id2 %in% zero_varying_genes$enssscg_id)) %>%
    
    #adjust prop to 1 - a few points are 1e-16 above or so:
    mutate(prop = ifelse(prop > 1, 1, prop), 
           weight = 1 - prop) %>%
    select(1, 2, weight, prop) %>%
    
    as_tbl_graph(directed = F, weighted = T) %>%
    
    
    activate(nodes) %>%
    left_join(gene_mapping %>% 
                select(1:2),
              by = c("name" = "enssscg_id"))
  
  


  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    add_node_feature("centrality_degree")
  # Centrality closeness takes a long time and tends to crash
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>% 
  #   add_node_feature("centrality_closeness")
  # Centrality betweenness takes a long time and tends to crash
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>% 
  #   add_node_feature("centrality_betweenness")
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    add_node_feature("centrality_eigen")
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    mutate(weighted_degree = centrality_degree / local_ave_degree())
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    add_node_feature("centrality_pagerank")
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>% 
  #   mutate(neighborhood_edges = map_local_dbl(.f = function(neighborhood, ...) {
  #     igraph::gsize(neighborhood)}))
  
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_fast_greedy")
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_louvain")
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_label_prop")
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_walktrap")
  
  # Group edge betweenness takes a long time!
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_edge_betweenness", directed = F)
  
  
  save(gene_coexpression_network,
       file = "data/processed/pig_coexpression_network.Rdata")
}

# ------ Expression network 2 ------

if(file.exists("data/processed/pig_coexpression_network2.Rdata")) {
  
  load("data/processed/pig_coexpression_network2.Rdata")
} else {
  
  if(file.exists("data/processed/pig_gene_pcor_res.Rdata")) {
    load("data/processed/pig_gene_pcor_res.Rdata")
    
  } else {
    
    if(file.exists("data/processed/pig_gene_pcor.Rdata")) {
      load("data/processed/pig_gene_pcor.Rdata")
      
    } else {
      
      
      pig_gene_pcor <- 
        pig_atlas_sample %>%
        select(1, 2, ptpm) %>% 
        spread(enssscg_id, ptpm) %>%
        column_to_rownames("sample_ID") %>% 
        {log10(. +1 )} %>%
        as.matrix() %>%
        GeneNet::ggm.estimate.pcor(method = "dynamic")
      
      
      
      save(pig_gene_pcor,
           file = "data/processed/pig_gene_pcor.Rdata")
      
    } 
    
    pig_gene_pcor_res <- 
      pig_gene_pcor[1:nrow(pig_gene_pcor), 1:nrow(pig_gene_pcor)] %>% 
      as.matrix() %>%
      {temp <- .; temp[!upper.tri(temp)] <- NA; temp} %>%
      as_tibble(rownames = "enssscg_id1") %>%
      gather(enssscg_id2, cor, -1, na.rm = T) 
    
    test_filters <- 
      seq(0.001, 0.01, length.out = 30)
    
    net_n_per_filter <- 
      lapply(test_filters,
             function(filtr) {
               edges_ <- 
                 pig_gene_pcor_res %>% 
                 filter(cor > filtr)
               
               nodes_ <- 
                 tibble(node = c(edges_$enssscg_id1, edges_$enssscg_id2)) %>% 
                 distinct()
               
               tibble(edges = nrow(edges_),
                      nodes = nrow(nodes_))
             }) %>% 
      set_names(test_filters) %>%
      bind_rows(.id = "filter") %>% 
      mutate(ratio = edges/nodes)
    
    net_n_per_filter %>%
      gather(type, n, -1) %>%
      mutate(filter = as.numeric(filter)) %>%
      ggplot(aes(filter, n, color = type)) +
      geom_line() + 
      facet_wrap(~type, ncol = 1, scales = "free_y") + 
      stripped_theme_facet
    ggsave(savepath("Cutoff metrics pcor.pdf"), width = 6, height = 10)
    
    net_n_per_filter %>%
      ggplot(aes(edges, nodes, label = round(edges / 1e6, 1))) + 
      geom_line() +
      geom_point() +
      geom_text() +
      scale_x_log10()
    
    
    
    save(pig_gene_pcor_res,
         file = "data/processed/pig_gene_pcor_res.Rdata")
    
  } 
  
  cor_filter <- 0.0025
  
  # pig_gene_pcor_res %>% 
  #   filter(cor > cor_filter)%>%
  #     filter(enssscg_id1 %in% xeno_ko_genes$enssscg_id |
  #              enssscg_id2 %in% xeno_ko_genes$enssscg_id) -> A
  # 
  # # xeno_ko_genes
  # 
  # A %>%
  #   filter(enssscg_id1 %in% xeno_ko_genes$enssscg_id) %>%
  #   group_by(enssscg_id1) %>%
  #   summarise(max_prop = max(cor),
  #             mean_prop = mean(cor),
  #             median_prop = median(cor),
  #             n = n()) %>%
  #   arrange(max_prop)
  # 
  # A %>%
  #   filter(enssscg_id2 %in% xeno_ko_genes$enssscg_id) %>%
  #   group_by(enssscg_id2) %>%
  #   summarise(max_prop = max(prop),
  #             mean_prop = mean(prop),
  #             median_prop = median(prop),
  #             n = n()) %>%
  #   arrange(max_prop)
  # 
  zero_varying_genes <- 
    pig_atlas_sample_norm %>% 
    group_by(enssscg_id) %>% 
    summarise(sd = sd(ptpm), 
              mean = mean(ptpm)) %>% 
    arrange(sd) %>%
    filter(sd == 0)
  
  gene_coexpression_network <- 
    pig_gene_pcor_res %>% 
    filter(cor > cor_filter) %>%
    
    filter(!(enssscg_id1 %in% zero_varying_genes$enssscg_id | 
               enssscg_id2 %in% zero_varying_genes$enssscg_id)) %>%
    
    select(1, 2, cor) %>%
    
    as_tbl_graph(directed = F, weighted = F) %>%
    
    
    activate(nodes) %>%
    left_join(gene_mapping %>% 
                select(1:2),
              by = c("name" = "enssscg_id"))
  
  
  
  
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    add_node_feature("centrality_degree")
  # Centrality closeness takes a long time and tends to crash
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>% 
  #   add_node_feature("centrality_closeness")
  # Centrality betweenness takes a long time and tends to crash
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>% 
  #   add_node_feature("centrality_betweenness")
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    add_node_feature("centrality_eigen")
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    mutate(weighted_degree = centrality_degree / local_ave_degree())
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    add_node_feature("centrality_pagerank")
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>% 
  #   mutate(neighborhood_edges = map_local_dbl(.f = function(neighborhood, ...) {
  #     igraph::gsize(neighborhood)}))
  
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_fast_greedy")
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_louvain")
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_label_prop")
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_walktrap")
  
  # Group edge betweenness takes a long time!
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>%
  #   add_node_feature("group_edge_betweenness", directed = F)
  
  
  network_hclust <- 
    gene_coexpression_network %>% 
    activate(edges) %>% 
    as_tibble() %>%
    select(from, to) %>% 
    bind_rows(tibble(from = .$to,
                     to = .$from, 
                     weight = 0)) %>% 
    spread(to, weight, fill = 1) %>%
    column_to_rownames("from") %>%
    as.dist() %>%
    hclust(method = "ward.D2")
  
  plot_data <- 
    tibble(h = seq(1, 6, length.out = 30)) %>% 
    group_by(h) %>% 
    do({
      cutree(network_hclust, h = .$h) %>% 
        enframe()
    }) 
  
  plot_data %>% 
    group_by(h) %>%
    summarise(n = n_distinct(value)) %>%
    ggplot(aes(as.factor(round(h, 2)), n, label = n)) +
    geom_col() + 
    geom_text(hjust = 0, 
              vjust = 0.35,
              angle = 90) +
    ylab("Number of clusters") +
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    
    
    plot_data %>% 
    group_by(h, value) %>%
    summarise(n = n()) %>%
    ggplot(aes(as.factor(round(h, 2)), n)) +
    geom_boxplot() +
    ylab("Number of genes per cluster") +
    xlab("Cutoff") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    scale_y_log10() + 
    plot_layout(ncol = 1)
  
  ggsave(savepath("network hclust clustering.pdf"), height = 6, width = 6)
  
  
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    activate(nodes) %>%
    mutate(group_hclust2p5 = network_hclust %>%
             cutree(h = 2.5) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust3 = network_hclust %>%
             cutree(h = 3) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust3p5 = network_hclust %>%
             cutree(h = 3.5) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust4 = network_hclust %>%
             cutree(h = 4) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust4p5 = network_hclust %>%
             cutree(h = 4.5) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust5 = network_hclust %>%
             cutree(h = 5) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2))
  
  
    
  
  save(gene_coexpression_network,
       file = "data/processed/pig_coexpression_network2.Rdata")
}


# ------ Expression network 3 ------

if(file.exists("data/processed/pig_coexpression_network3.Rdata")) {
  
  load("data/processed/pig_coexpression_network3.Rdata")
} else {
  
  if(file.exists("data/processed/pig_gene_spearman_res.Rdata")) {
    load("data/processed/pig_gene_spearman_res.Rdata")
    
  } else {
    
    if(file.exists("data/processed/pig_gene_spearman.Rdata")) {
      load("data/processed/pig_gene_spearman.Rdata")
      
    } else {
      
      
      pig_gene_spearman <- 
        pig_atlas_sample %>%
        select(1, 2, ptpm) %>% 
        spread(enssscg_id, ptpm) %>%
        column_to_rownames("sample_ID") %>% 
        as.matrix() %>%
        cor(method = "spearman")
        
      
      save(pig_gene_spearman,
           file = "data/processed/pig_gene_spearman.Rdata")
      
    } 
    
    pig_gene_spearman_res <- 
      pig_gene_spearman %>%
      {temp <- .; temp[!upper.tri(temp)] <- NA; temp} %>%
      as_tibble(rownames = "enssscg_id1") %>%
      gather(enssscg_id2, cor, -1, na.rm = T) 
    
    test_filters <- 
      seq(0.3, 0.9, length.out = 30)
    
    net_n_per_filter <- 
      lapply(test_filters,
             function(filtr) {
               edges_ <- 
                 pig_gene_spearman_res %>% 
                 filter(cor > filtr)
               
               nodes_ <- 
                 tibble(node = c(edges_$enssscg_id1, edges_$enssscg_id2)) %>% 
                 distinct()
               
               tibble(edges = nrow(edges_),
                      nodes = nrow(nodes_))
             }) %>% 
      set_names(test_filters) %>%
      bind_rows(.id = "filter") %>% 
      mutate(ratio = edges/nodes)
    
    net_n_per_filter %>%
      gather(type, n, -1) %>%
      mutate(filter = as.numeric(filter)) %>%
      ggplot(aes(filter, n, color = type)) +
      geom_line() + 
      facet_wrap(~type, ncol = 1, scales = "free_y") + 
      stripped_theme_facet
    ggsave(savepath("Cutoff metrics spearman.pdf"), width = 6, height = 10)
    
    net_n_per_filter %>%
      ggplot(aes(edges, nodes, label = round(edges / 1e6, 1))) + 
      geom_line() +
      geom_point() +
      geom_text() +
      scale_x_log10()
    
    
    
    save(pig_gene_spearman_res,
         file = "data/processed/pig_gene_spearman_res.Rdata")
    
  } 
  
  cor_filter <- 0.5
  # 
  # pig_gene_spearman_res %>%
  #   filter(cor > cor_filter)%>%
  #     filter(enssscg_id1 %in% xeno_ko_genes$enssscg_id |
  #              enssscg_id2 %in% xeno_ko_genes$enssscg_id) -> A
  # 
  # # xeno_ko_genes
  # 
  # A %>%
  #   filter(enssscg_id1 %in% xeno_ko_genes$enssscg_id) %>%
  #   group_by(enssscg_id1) %>%
  #   summarise(max_prop = max(cor),
  #             mean_prop = mean(cor),
  #             median_prop = median(cor),
  #             n = n()) %>%
  #   arrange(max_prop)
  # 
  # A %>%
  #   filter(enssscg_id2 %in% xeno_ko_genes$enssscg_id) %>%
  #   group_by(enssscg_id2) %>%
  #   summarise(max_prop = max(cor),
  #             mean_prop = mean(cor),
  #             median_prop = median(cor),
  #             n = n()) %>%
  #   arrange(max_prop)

  zero_varying_genes <- 
    pig_atlas_sample_norm %>% 
    group_by(enssscg_id) %>% 
    summarise(sd = sd(ptpm), 
              mean = mean(ptpm)) %>% 
    arrange(sd) %>%
    filter(sd == 0)
  
  gene_coexpression_network <- 
    pig_gene_spearman_res %>% 
    filter(cor > cor_filter) %>%
    
    filter(!(enssscg_id1 %in% zero_varying_genes$enssscg_id | 
               enssscg_id2 %in% zero_varying_genes$enssscg_id)) %>%
    
    select(1, 2, cor) %>%
    
    as_tbl_graph(directed = F, weighted = F) %>%
    
    
    activate(nodes) %>%
    left_join(gene_mapping %>% 
                select(1:2),
              by = c("name" = "enssscg_id"))
  
  
  
  
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    add_node_feature("centrality_degree")
  # Centrality closeness takes a long time and tends to crash
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>% 
  #   add_node_feature("centrality_closeness")
  # Centrality betweenness takes a long time and tends to crash
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>% 
  #   add_node_feature("centrality_betweenness")
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    add_node_feature("centrality_eigen")
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    mutate(weighted_degree = centrality_degree / local_ave_degree())
  gene_coexpression_network <- 
    gene_coexpression_network %>% 
    add_node_feature("centrality_pagerank")
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>% 
  #   mutate(neighborhood_edges = map_local_dbl(.f = function(neighborhood, ...) {
  #     igraph::gsize(neighborhood)}))
  
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_fast_greedy")
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_louvain")
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_label_prop")
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    add_node_feature("group_walktrap")
  
  # Group edge betweenness takes a long time!
  # gene_coexpression_network <- 
  #   gene_coexpression_network %>%
  #   add_node_feature("group_edge_betweenness", directed = F)
  
  
  network_hclust <- 
    gene_coexpression_network %>% 
    activate(edges) %>% 
    as_tibble() %>%
    select(from, to) %>% 
    bind_rows(tibble(from = .$to,
                     to = .$from, 
                     weight = 0)) %>% 
    spread(to, weight, fill = 1) %>%
    column_to_rownames("from") %>%
    as.dist() %>%
    hclust(method = "ward.D2")
  
  plot_data <- 
    tibble(h = seq(1, 6, length.out = 30)) %>% 
    group_by(h) %>% 
    do({
      cutree(network_hclust, h = .$h) %>% 
        enframe()
    }) 
  
  plot_data %>% 
    group_by(h) %>%
    summarise(n = n_distinct(value)) %>%
    ggplot(aes(as.factor(round(h, 2)), n, label = n)) +
    geom_col() + 
    geom_text(hjust = 0, 
              vjust = 0.35,
              angle = 90) +
    ylab("Number of clusters") +
    theme(axis.text.x = element_blank(), 
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()) + 
    
    
    plot_data %>% 
    group_by(h, value) %>%
    summarise(n = n()) %>%
    ggplot(aes(as.factor(round(h, 2)), n)) +
    geom_boxplot() +
    ylab("Number of genes per cluster") +
    xlab("Cutoff") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
    scale_y_log10() + 
    plot_layout(ncol = 1)
  
  ggsave(savepath("network hclust clustering.pdf"), height = 6, width = 6)
  
  
  gene_coexpression_network <- 
    gene_coexpression_network %>%
    activate(nodes) %>%
    mutate(group_hclust2p5 = network_hclust %>%
             cutree(h = 2.5) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust3 = network_hclust %>%
             cutree(h = 3) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust3p5 = network_hclust %>%
             cutree(h = 3.5) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust4 = network_hclust %>%
             cutree(h = 4) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust4p5 = network_hclust %>%
             cutree(h = 4.5) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2),
           group_hclust5 = network_hclust %>%
             cutree(h = 5) %>% 
             enframe("gene_i", "cluster") %>%
             arrange(gene_i) %>% 
             pull(2))
  
  
  
  
  save(gene_coexpression_network,
       file = "data/processed/pig_coexpression_network2.Rdata")
}
# ------ small network ------

if(file.exists("data/processed/small_net.Rdata")) {
  
  load("data/processed/small_net.Rdata")
} else {
  
  
  small_net <- 
    gene_coexpression_network %>% 
    filter(group_fast_greedy == 3)
  
  
  set.seed(42)
  small_net <- 
    small_net %>%
    mutate(group_louvain = group_louvain(weights = weight))
  small_net <- 
    small_net %>%
    mutate(group_fast_greedy = group_fast_greedy(weights = weight))
  small_net <- 
    small_net %>%
    mutate(group_infomap = group_infomap(weights = weight))
  small_net <- 
    small_net %>%
    mutate(group_label_prop = group_label_prop(weights = weight))
  # Spinglass takes a long time:
  # small_net <- 
  #   small_net %>%
  #   mutate(group_spinglass = group_spinglass(weights = prop))
  # small_net <- 
  #   small_net %>%
  #   mutate(group_edge_betweenness = group_edge_betweenness(weights = weight))
  small_net <- 
    small_net %>%
    mutate(group_leading_eigen = group_leading_eigen(weights = weight))
  # Optimal takes a long time:
  # small_net <- 
  #   small_net %>%
  #   mutate(group_optimal = group_optimal(weights = prop))
  small_net <- 
    small_net %>%
    mutate(group_walktrap = group_walktrap(weights = weight))
  
  # #####
  # 
  # set.seed(42)
  # small_net2 <- 
  #   small_net %>%
  #   mutate(group_louvain = group_louvain(weights = prop))
  # small_net2 <- 
  #   small_net2 %>%
  #   mutate(group_fast_greedy = group_fast_greedy(weights = prop))
  # small_net2 <- 
  #   small_net2 %>%
  #   mutate(group_infomap = group_infomap(weights = prop))
  # small_net2 <- 
  #   small_net2 %>%
  #   mutate(group_label_prop = group_label_prop(weights = prop))
  # # Spinglass takes a long time:
  # # small_net <- 
  # #   small_net %>%
  # #   mutate(group_spinglass = group_spinglass(weights = prop))
  # # small_net2 <- 
  # #   small_net2 %>%
  # #   mutate(group_edge_betweenness = group_edge_betweenness(weights = prop))
  # small_net2 <- 
  #   small_net2 %>%
  #   mutate(group_leading_eigen = group_leading_eigen(weights = prop))
  # # Optimal takes a long time:
  # # small_net <- 
  # #   small_net %>%
  # #   mutate(group_optimal = group_optimal(weights = prop))
  # small_net2 <- 
  #   small_net2 %>%
  #   mutate(group_walktrap = group_walktrap(weights = prop))
  # 
  # 
  # 
  
  #####
  
  
  save(small_net,
       file = "data/processed/small_net.Rdata")
}

small_net_specs <- 
  small_net %>% 
  as_tibble() %>% 
  select(name, contains("group")) %>% 
  gather(group_type, cluster, -1) %>% 
  group_by(group_type, cluster) %>% 
  summarise(n = n()) %>% 
  mutate(n_cluster = n_distinct(cluster)) %>% 
  ungroup()
  
small_net_specs %>%
  ggplot(aes(cluster, n)) +
  geom_col() + 
  facet_wrap(~group_type)

small_net %>% 
  activate(nodes) %>%
  rename(group = group_col) %>%
  morph(to_contracted, group, 
        edge.attr.comb= "mean",
        simplify = F) %>%
  crystallise() %>%
  pull(graph) %>%
  {.[[1]]} %>%
  mutate(n = sapply(.tidygraph_node_index,
                    length)) %>%
  activate(edges) %>% 
  mutate(weight = sapply(.orig_data, 
                         function(x) {mean(x$weight)}), 
         prop = sapply(.orig_data, 
                       function(x) {mean(x$prop)})) 

bind_rows(get_edge_groups(small_net, 
                cols_ <- c('group_fast_greedy', 
                           'group_louvain', 
                           'group_infomap', 
                           'group_label_prop', 
                           # 'group_edge_betweenness', 
                           # 'group_leading_eigen', 
                           'group_walktrap')) %>%
            mutate(type = "weight"),
          get_edge_groups(small_net2, 
                          cols_ <- c('group_fast_greedy', 
                                     'group_louvain', 
                                     'group_infomap', 
                                     'group_label_prop', 
                                     # 'group_edge_betweenness', 
                                     # 'group_leading_eigen', 
                                     'group_walktrap')) %>%
            mutate(type = "prop")) %>%
  mutate(within = new_group_from == new_group_to) %>%
  filter(within) %>%
  filter(new_group_from %in% 1:5) %>%
  ggplot(aes(as.factor(new_group_from), prop, fill = type)) +
  geom_boxplot() + 
  
  facet_grid(name ~ .)

plots <- 
  lapply(grep("group", colnames(as_tibble(small_net)), value = T),
         function(x) {plot_network_groups(small_net, x)})

pdf(savepath("Small network groupings.pdf"), width = 6, height = 6)
plots
dev.off()



# ----- Analysis of non-participating genes -----
pig_gene_classification %>%
  filter(!enssscg_id %in% as_tibble(gene_coexpression_network)$name) -> A



pig_atlas_sample_norm %>%
  filter(enssscg_id %in% A$enssscg_id) %>% 
  group_by(enssscg_id) %>% 
  summarise(max_ptpm = max(ptpm)) %>% 
  ggplot(aes(log10(max_ptpm + 1))) + 
  geom_density()

A %>% 
  group_by(specificity_category) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(specificity_category, n)) +
  geom_col()

#  
# pig_gene_spearman_propr %>%
#   filter(enssscg_id1 %in% A$enssscg_id |
#            enssscg_id2 %in% A$enssscg_id) -> B

B %>% 
  gather(type, cor, prop, cor) %>%
  ggplot(aes(cor, color = type)) +
  geom_density()

pig_gene_spearman_propr %>% 
  ggplot(aes(cor, prop)) +
  geom_hex(aes(fill = stat(log10(count))),
           bins = 100)

B %>%
  filter(prop < 0.55 & prop > 0 & cor > 0.9) %>%
  arrange(-cor) -> B_x

plots <- 
  lapply(1:10,
         function(i) {
           
           genes_ <- c(B_x$enssscg_id1[i], B_x$enssscg_id2[i])
           
           pig_atlas_sample %>%
             filter(enssscg_id %in% genes_) %>%
             select(-ptpm) %>%
             spread(enssscg_id, tpm) %>% 
             rename(gen1 = genes_[1], gen2 = genes_[2]) %>%
             ggplot(aes(log10(gen1 + 1), log10(gen2 + 1))) +
             geom_point()+
             geom_text(aes(label = sample_ID)) +
             xlab(genes_[1]) + 
             ylab(genes_[2])
         })

pdf(savepath("high spearman low propr examples.pdf"), width = 6, height = 6)
plots
dev.off()


pig_gene_spearman_propr %>%
  filter(prop > 0.9 & cor < 0.1) %>%
  arrange(-cor) -> B_x2

plots <- 
  lapply(1:10,
         function(i) {
           
           genes_ <- c(B_x2$enssscg_id1[i], B_x2$enssscg_id2[i])
           
           pig_atlas_sample %>%
             filter(enssscg_id %in% genes_) %>%
             select(-ptpm) %>%
             spread(enssscg_id, tpm) %>% 
             rename(gen1 = genes_[1], gen2 = genes_[2]) %>%
             ggplot(aes(log10(gen1 + 1), log10(gen2 + 1))) +
             geom_point()+
             geom_text(aes(label = sample_ID)) +
             xlab(genes_[1]) + 
             ylab(genes_[2])
         })

pdf(savepath("low spearman high propr examples.pdf"), width = 6, height = 6)
plots
dev.off()


pig_gene_spearman_propr %>%
  filter(prop > 0.9 & cor > 0.9) %>%
  arrange(-cor) -> B_x3

plots <- 
  lapply(1:10,
         function(i) {
           
           genes_ <- c(B_x3$enssscg_id1[i], B_x3$enssscg_id2[i])
           
           pig_atlas_sample %>%
             filter(enssscg_id %in% genes_) %>%
             select(-ptpm) %>%
             spread(enssscg_id, tpm) %>% 
             rename(gen1 = genes_[1], gen2 = genes_[2]) %>%
             ggplot(aes(log10(gen1 + 1), log10(gen2 + 1))) +
             geom_point()+
             geom_text(aes(label = sample_ID)) +
             xlab(genes_[1]) + 
             ylab(genes_[2])
         })

pdf(savepath("high spearman high propr examples.pdf"), width = 6, height = 6)
plots
dev.off()

pig_atlas_sample %>%
  filter(enssscg_id %in% c("ENSSSCG00000040769", "ENSSSCG00000037479", "ENSSSCG00000033978", "ENSSSCG00000032291")) %>%
  select(-ptpm) %>%
  spread(enssscg_id, tpm) 

gene_orthologs_all %>% 
  filter(enssscg_id %in% c("ENSSSCG00000040769", "ENSSSCG00000037479", "ENSSSCG00000033978", "ENSSSCG00000032291")) 
  
  

  pig_atlas_sample %>%
  filter(enssscg_id %in% c("ENSSSCG00000009687", "ENSSSCG00000040042")) %>%
  select(-ptpm) %>%
  spread(enssscg_id, tpm) %>% 
  ggplot(aes(log10(ENSSSCG00000009687 + 1), log10(ENSSSCG00000040042 + 1))) +
  geom_point()+
  geom_text(aes(label = sample_ID))

# ----- Network visualization -----

plots <- 
  lapply(grep("group", colnames(as_tibble(gene_coexpression_network)), value = T),
         function(x) {plot_network_groups(gene_coexpression_network, x, cor_col = "cor", min_nodes = 10)})

pdf(savepath("Network groupings.pdf"), width = 6, height = 6)
plots
dev.off()




# ----- Network metrics ------

gene_coexpression_network %>%
  as_tibble() %>% 
  select(-contains("group")) %>% 
  gather(metric, value, -1, -2) %>% 
  ggplot(aes(value)) +
  geom_histogram() +
  facet_wrap(~metric, scales = "free") +
  stripped_theme_facet

# ------ Hypergeometric test ------

network_hyper_test_data <-
  gene_coexpression_network %>%
  as_tibble() %>%
  select(enssscg_id = name,
         contains("group")) %>%
  gather(group_type, cluster, -1) %>%
  left_join(pig_gene_classification %>%
              mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues), "not enriched", enhanced_tissues)) %>%
              separate_rows(enhanced_tissues, sep = ",")) %>% 
   group_by(group_type) %>%
  
  do({
    # q is the number of successes
    q <- 
      group_by(., enhanced_tissues, cluster, .drop = F) %>%
      summarise(q = n_distinct(enssscg_id))
    # k is the number of tries - i.e. the cluster size
    k <- 
      group_by(., cluster, .drop = F) %>%
      summarise(k = n_distinct(enssscg_id))
    
    # m is the number of possible successes
    m <- 
      group_by(., enhanced_tissues, .drop = F) %>%
      summarise(m = n_distinct(enssscg_id))
    
    # n is the population size - i.e. the number of genes
    n <- n_distinct(.$enssscg_id)
    
    q %>% 
      left_join(k) %>% 
      left_join(m) %>%
      # n is the population size - i.e. the number of genes
      mutate(n = n - m) 
    
  })

network_enrich_hyper <- 
  network_hyper_test_data %>%
  group_by(group_type, enhanced_tissues, cluster, q, k, .drop = F) %>%
  # Testing the chance of getting the observed number or higher per cluster and tissue type
  summarise(p_value = phyper(q - 1, m, n, k, lower.tail = F)) %>%
  ungroup() %>%
  mutate(p_value = ifelse(p_value == 0, 1e-300, p_value),
         adj_pval = p.adjust(p_value, method = "BH"),
         sign_level = case_when(adj_pval < 0.00000001 ~ 8,
                                adj_pval < 0.0000001 ~ 7,
                                adj_pval < 0.000001 ~ 6,
                                adj_pval < 0.00001 ~ 5,
                                adj_pval < 0.0001 ~ 4,
                                adj_pval < 0.001 ~ 3,
                                adj_pval < 0.01 ~ 2,
                                adj_pval < 0.05 ~ 1,
                                T ~ 0))



for(group_type_ in unique(network_enrich_hyper$group_type)) {
  # print(group_type_)
  plot_data <- 
    network_enrich_hyper %>% 
    filter(q >= 10) %>%
    filter(group_type == group_type_) %>%
    select(enhanced_tissues, cluster, sign_level) %>%
    spread(cluster, sign_level, fill = 0) %>% 
    column_to_rownames("enhanced_tissues") 
  n_u <- 
    unlist(plot_data) %>% 
    unique() %>% 
    length()
  
  if(n_u > 1) {
    plot_data %>% 
      pheatmap(clustering_method = "ward.D2",
               color = colorRampPalette(c("white", "orange", "orangered"))(9),
               filename = savepath(paste0("network overlap heatmap ", group_type_, ".pdf")),
               width = 16, 
               height = 6)
  }
  
}

dev.off()
# ----- Xeno-KO genes -----

pig_atlas_tissue %>% 
  inner_join(xeno_ko_genes) %>%
  left_join(tissue_mapping,
            by = c("tissue_ID" = "tissue")) %>%
  select(gene_name, tissue_name, ptpm) %>% 
  spread(tissue_name, ptpm) %>% 
  column_to_rownames("gene_name") %>% 
  {log10(. + 1)} %>%
  pheatmap(clustering_method = "ward.D2",
           color = heatmap_palette,
           border_color = NA,
           filename = savepath("xeno-ko genes expression heatmap.pdf"),
           width = 12, 
           height = 4)

xeno_subgraph_settings <- 
  expand.grid(enssscg_id = xeno_ko_genes$enssscg_id,
              cor_lim = seq(0.5, 0.95, 0.05))

if(file.exists("data/processed/pig_xeno_subgraphs.Rdata")) {
  
  load("data/processed/pig_xeno_subgraphs.Rdata")
} else {
  xeno_subgraphs_full <- 
    lapply(unique(xeno_subgraph_settings$enssscg_id),
           function(gene_id) {
             gene_coexpression_network %>%
               morph(to_local_neighborhood, node = which(.N()$name == gene_id), order = 1) %>% 
               crystallize() %>%
               pull(2) %>%
               {.[[1]]}
           }) %>%
    set_names(unique(xeno_subgraph_settings$enssscg_id))
  
  xeno_subgraphs <- 
    lapply(1:nrow(xeno_subgraph_settings),
           function(i) {
             gene_id <- xeno_subgraph_settings$enssscg_id[i]
             cor_lim <- xeno_subgraph_settings$cor_lim[i]
             
             xeno_subgraphs_full[[gene_id]] %>%
               activate(edges) %>% 
               filter(cor > cor_lim) %>%
               activate(nodes) %>%
               morph(to_local_neighborhood, node = which(.N()$name == gene_id), order = 1) %>% 
               crystallize() %>%
               pull(2) %>%
               {.[[1]]}
           })
  
  save(xeno_subgraphs,
       file = "data/processed/pig_xeno_subgraphs.Rdata")
}



xeno_subgraph_genes <- 
  xeno_subgraphs %>% 
  map(as_tibble) %>%
  set_names(as.character(1:length(.))) %>%
  bind_rows(.id = "settings") %>% 
  left_join(xeno_subgraph_settings %>% 
              mutate(settings = as.character(row_number()))) %>% 
  select(xeno_id = enssscg_id, cor_lim, enssscg_id = name, everything()) %>% 
  select(-settings)

xeno_subgraph_genes %>% 
  group_by(xeno_id, cor_lim) %>% 
  summarise(n = n()) %>% 
  ggplot(aes(n, xeno_id, fill = as_factor(cor_lim), 
             label = n)) + 
  geom_col(position = "dodge") + 
  geom_text(position = position_dodge(width = 1),
            hjust = 0)

xeno_hyper_test_data <-
  xeno_subgraph_genes %>% 
  left_join(pig_gene_classification %>%
              mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues), "not enriched", enhanced_tissues)) %>%
              separate_rows(enhanced_tissues, sep = ",") %>% 
              mutate(enhanced_tissues = factor(enhanced_tissues, unique(enhanced_tissues)))) %>% 
  mutate(cor_lim = factor(cor_lim, unique(cor_lim))) %>%
  # group_by() %>%
  
  do({
    # q is the number of successes
    q <- 
      group_by(., enhanced_tissues, xeno_id, cor_lim, .drop = F) %>%
      summarise(q = n_distinct(enssscg_id)) 
    # k is the number of tries - i.e. the cluster size
    k <- 
      group_by(., xeno_id, cor_lim, .drop = F) %>%
      summarise(k = n_distinct(enssscg_id))
    
    # m is the number of possible successes
    
    total_net <- 
      gene_coexpression_network %>% 
      as_tibble()
    
    m <- 
      pig_gene_classification %>%
      filter(enssscg_id %in% total_net$name) %>%
      mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues), "not enriched", enhanced_tissues)) %>%
      separate_rows(enhanced_tissues, sep = ",") %>%
      group_by(., enhanced_tissues, .drop = F) %>%
      summarise(m = n_distinct(enssscg_id))
    
    # n is the population size - i.e. the number of genes
    n <- 
      total_net %>% 
      nrow()
    
    q %>% 
      left_join(k) %>% 
      left_join(m) %>%
      # n is the population size - i.e. the number of genes
      mutate(n = n - m) 
    
  })

xeno_enrich_hyper <- 
  xeno_hyper_test_data %>%
  group_by(xeno_id, cor_lim, enhanced_tissues, q, k, .drop = F) %>%
  # Testing the chance of getting the observed number or higher per cluster and tissue type
  summarise(p_value = phyper(q - 1, m, n, k, lower.tail = F)) %>%
  ungroup() %>%
  mutate(p_value = ifelse(p_value == 0, 1e-300, p_value),
         adj_pval = p.adjust(p_value, method = "BH"),
         sign_level = case_when(adj_pval < 0.00000001 ~ 8,
                                adj_pval < 0.0000001 ~ 7,
                                adj_pval < 0.000001 ~ 6,
                                adj_pval < 0.00001 ~ 5,
                                adj_pval < 0.0001 ~ 4,
                                adj_pval < 0.001 ~ 3,
                                adj_pval < 0.01 ~ 2,
                                adj_pval < 0.05 ~ 1,
                                T ~ 0)) 

xeno_hyper_cluster_data <-
  xeno_enrich_hyper %>% 
  left_join(xeno_ko_genes,
            by = c("xeno_id" = "enssscg_id")) %>%
  filter(cor_lim == 0.6) %>%
  select(gene_name, enhanced_tissues, sign_level) %>% 
  spread(gene_name, sign_level) %>%
  column_to_rownames("enhanced_tissues")

xeno_hyper_tis_cluster <- 
  xeno_hyper_cluster_data %>% 
  dist() %>%
  hclust(method = "ward.D2")

xeno_hyper_gene_cluster <- 
  xeno_hyper_cluster_data %>% 
  t() %>%
  dist() %>%
  hclust(method = "ward.D2")

xeno_enrich_hyper %>%
  
  group_by(xeno_id, enhanced_tissues) %>% 
  mutate(any_sign = any(sign_level != 0)) %>% 
  ungroup() %>%
  filter(any_sign) %>%
  ggplot(aes(cor_lim, -log10(p_value), group = enhanced_tissues, color = enhanced_tissues)) +
  geom_line(show.legend = F) + 
  geom_text_repel(data = . %>% 
                    group_by(enhanced_tissues, xeno_id) %>%
                    summarise(p_value = min(p_value)),
                  aes(label = enhanced_tissues), 
                  size = 1.5,
                  show.legend = F) +
  geom_line(show.legend = F) + 
  geom_hline(yintercept = -log10(0.05)) +
  facet_wrap(~xeno_id, scales = "free_y") + 
  theme_minimal()+ 
  scale_color_manual(values = c(tissue_colors_palette_full, "not enriched" = "black"))
ggsave(savepath("xeno subnetwork tis p-val per corlim.pdf"), width = 10, height = 6)

xeno_enrich_hyper %>% 
  left_join(xeno_ko_genes,
            by = c("xeno_id" = "enssscg_id")) %>%
  mutate(enhanced_tissues = factor(enhanced_tissues, ggdendro::dendro_data(xeno_hyper_tis_cluster)$labels$label),
         gene_name = factor(gene_name, ggdendro::dendro_data(xeno_hyper_gene_cluster)$labels$label)) %>%
  filter(sign_level > 0) %>%
  ggplot(aes(enhanced_tissues, gene_name, fill = sign_level)) +
  geom_tile() +
  facet_wrap(~cor_lim) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "orange", "orangered"))(9)) +
  stripped_theme_facet +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major.x = element_line(color = "lightgray"))

ggsave(savepath("Xeno-KO genes network neighborhood enrichment 1.pdf"), width = 10, height = 6)

xeno_enrich_hyper %>% 
  filter(cor_lim == 0.6) %>%
  left_join(xeno_ko_genes,
            by = c("xeno_id" = "enssscg_id")) %>%
  arrange(adj_pval) %>%
  mutate(enhanced_tissues = factor(enhanced_tissues, ggdendro::dendro_data(xeno_hyper_tis_cluster)$labels$label),
         gene_name = factor(gene_name, ggdendro::dendro_data(xeno_hyper_gene_cluster)$labels$label),
         adj_pval = ifelse(adj_pval < 1e-20, 1e-20, adj_pval)) %>%
  filter(sign_level != 0) %>%
  ggplot(aes(enhanced_tissues, gene_name, 
             fill = -log10(adj_pval), 
             size = -log10(adj_pval))) +
  geom_point(shape = 21,
             color = "black") +
  geom_text(aes(enhanced_tissues, gene_name, label = q), 
            inherit.aes = F,
            vjust = -1) +
  # facet_wrap(~order) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "orange", "orangered"))(21), 
                       values = (0:20)/20) +
  stripped_theme_facet +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major.x = element_line(color = "lightgray"))

ggsave(savepath("Xeno-KO genes network neighborhood enrichment 2.pdf"), width = 6, height = 4)

# ----- Go analysis -----

library(enrichR)


xeno_subgraphs_enrichr_res <- 
  lapply(1:nrow(xeno_subgraph_settings),
         function(i) {
           
           gene_id <- xeno_subgraph_settings$enssscg_id[i]
           cor_lim_ <- xeno_subgraph_settings$cor_lim[i]
           
           enrichR_genes <- 
             xeno_subgraph_genes %>% 
             select(1:4) %>%
             left_join(gene_orthologs_all) %>%
             filter(xeno_id == gene_id & cor_lim == cor_lim_) %>%
             filter(!is.na(human_gene_name)) %>%
             pull(human_gene_name) %>% 
             unique()
           
           enrichr(enrichR_genes, 
                   databases = c("GO_Biological_Process_2018")) %>%
             bind_rows(.id = "database") 
         }) %>%
  map(as_tibble) %>%
  set_names(as.character(1:length(.))) %>%
  bind_rows(.id = "settings") %>% 
  left_join(xeno_subgraph_settings %>% 
              mutate(settings = as.character(row_number()))) %>% 
  select(xeno_id = enssscg_id, cor_lim, everything()) %>% 
  select(-settings)


xeno_subgraphs_enrichr_res %>%
  
  group_by(xeno_id, Term) %>% 
  mutate(any_sign = any(Adjusted.P.value < 0.01)) %>% 
  ungroup() %>%
  filter(any_sign) %>%
  ggplot(aes(cor_lim, -log10(Adjusted.P.value), group = Term)) +
  geom_line(show.legend = F) + 
  geom_hline(yintercept = -log10(0.01)) +
  facet_wrap(~xeno_id, scales = "free_y") + 
  theme_minimal()
ggsave(savepath("xeno subnetwork term p-val per corlim 1.pdf"), width = 10, height = 6)


xeno_subgraphs_enrichr_res %>%
  
  group_by(xeno_id, Term) %>% 
  mutate(any_sign = any(Adjusted.P.value < 0.01)) %>% 
  ungroup() %>%
  filter(any_sign) %>%
  ggplot(aes(cor_lim, -log10(Adjusted.P.value), group = cor_lim)) +
  geom_boxplot() +
  geom_hline(yintercept = -log10(0.01)) +
  geom_vline(data = . %>% 
               group_by(xeno_id, cor_lim) %>% 
               summarise(median_p = median(Adjusted.P.value)) %>% 
               group_by(xeno_id) %>%
               filter(median_p == min(median_p)), 
             aes(xintercept = cor_lim)) +
  facet_wrap(~xeno_id, scales = "free_y") + 
  theme_minimal()
ggsave(savepath("xeno subnetwork term p-val per corlim 2.pdf"), width = 10, height = 6)


xeno_subgraphs_enrichr_res %>%
  
  group_by(xeno_id, Term) %>% 
  mutate(any_sign = any(Adjusted.P.value < 0.01)) %>% 
  ungroup() %>%
  filter(any_sign) %>%
  ggplot(aes(cor_lim, -log10(Adjusted.P.value), group = cor_lim)) +
  geom_boxplot() +
  geom_hline(yintercept = -log10(0.01)) +
  geom_vline(data = . %>% 
               group_by(cor_lim) %>% 
               summarise(median_p = median(Adjusted.P.value)) %>% 
               filter(median_p == min(median_p)), 
             aes(xintercept = cor_lim)) +
  theme_minimal()
ggsave(savepath("xeno subnetwork term p-val per corlim 3.pdf"), width = 3, height = 4)


xeno_subgraphs_enrichr_res_filtered <-
  xeno_subgraphs_enrichr_res %>% 
  filter(Adjusted.P.value < 0.01) %>%
  filter(cor_lim == 0.6) 

xeno_subgraphs_enrichr_cluster_data <-
  xeno_subgraphs_enrichr_res_filtered %>%
  select(xeno_id, Term, Adjusted.P.value) %>% 
  spread(xeno_id, Adjusted.P.value, fill = 1) %>%
  column_to_rownames("Term") %>% 
  {-log10(.)}

xeno_subgraphs_enrichr_term_cluster <- 
  xeno_subgraphs_enrichr_cluster_data %>% 
  dist() %>%
  hclust(method = "ward.D2")

xeno_subgraphs_enrichr_gene_cluster <- 
  xeno_subgraphs_enrichr_cluster_data %>% 
  t() %>%
  dist() %>%
  hclust(method = "ward.D2")

xeno_subgraphs_enrichr_res_filtered %>%
  arrange(Adjusted.P.value) %>% 
  left_join(xeno_ko_genes,
            by = c("xeno_id" = "enssscg_id")) %>%
  arrange(-Combined.Score) %>%
  mutate(Term = factor(Term, unique(Term)),
         xeno_id = factor(xeno_id, ggdendro::dendro_data(xeno_subgraphs_enrichr_gene_cluster)$labels$label)) %>%
  group_by(xeno_id) %>%
  top_n(20, Combined.Score) %>%
  ggplot(aes(Term, gene_name, fill = -log10(Adjusted.P.value))) +
  geom_tile() +
  # facet_wrap(~order) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "orange", "orangered"))(9)) +
  stripped_theme_facet +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_line(color = "lightgray"))

ggsave(savepath("Xeno neighborhood GO result heatmap.pdf"), width = 10, height = 8)

# Text mining results

library("tm")
library("ggwordcloud")
# docs <- 

xeno_subgraphs_enrichr_res_filtered %>%
  select(xeno_id, Term, Adjusted.P.value, Combined.Score) %>%
  mutate(term = gsub(" \\(.*\\)$", "", Term)) %>% 
  mutate(term = gsub("B cell", "B-cell", term)) %>% 
  mutate(term = gsub("T cell", "T-cell", term)) %>% 
  mutate(term = tolower(term)) %>% 
  separate_rows(term, sep = " ") %>% 
  filter(grepl("t", term)) %>% 
  filter(!term %in% c(stopwords("english"), 
                      c("regulation", "positive", "cell", "response", 
                        "negative", "signaling", "cellular", "pathway", "process"))) %>% 
  
  group_by(xeno_id, term) %>% 
  summarise(n = n()) %>% 
  
  arrange(-n) %>% 
  top_n(20, n) %>%
  left_join(xeno_ko_genes, 
            by = c("xeno_id" = "enssscg_id")) %>%
  # mutate(freq = n / sum(n)) %>%
  
  ggplot(aes(label = term, size = n)) +
  geom_text_wordcloud(eccentricity = 0.8) +
  theme_minimal() +
  facet_wrap(~gene_name)
ggsave(savepath("Xeno neighborhood GO wordclouds.pdf"), width = 10, height = 8)

# ----- Blood atlas overlap analysis -----

human_gene_class_blood_lymph_enr <- 
  human_gene_class_blood %>% 
  mutate(enhanced_tissues = ifelse(ensg_id %in% filter(human_gene_class, 
                                                       grepl("blood|lymph", enhanced_tissues))$ensg_id,
                                   enhanced_tissues, 
                                   NA))


xeno_blood_hyper_test_data <-
  xeno_subgraph_genes %>% 
  inner_join(gene_orthologs_all) %>%
  left_join(human_gene_class_blood_lymph_enr %>%
              mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues), "not enriched", enhanced_tissues)) %>%
              separate_rows(enhanced_tissues, sep = ",") %>% 
              mutate(enhanced_tissues = factor(enhanced_tissues, unique(enhanced_tissues)))) %>% 
  mutate(cor_lim = factor(cor_lim, unique(cor_lim))) %>%
  # filter(xeno_id == "ENSSSCG00000005518" & cor_lim == 0.55) -> a
  # group_by() %>%
  
  do({
    # q is the number of successes
    q <- 
      group_by(., enhanced_tissues, xeno_id, cor_lim, .drop = F) %>%
      summarise(q = n_distinct(ensg_id)) 
    # k is the number of tries - i.e. the cluster size
    k <- 
      group_by(., xeno_id, cor_lim, .drop = F) %>%
      summarise(k = n_distinct(ensg_id))
    
    # m is the number of possible successes
    
    total_net <- 
      gene_coexpression_network %>% 
      as_tibble()
    
    m <- 
      human_gene_class_blood_lymph_enr %>%
      mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues), "not enriched", enhanced_tissues)) %>%
      separate_rows(enhanced_tissues, sep = ",") %>%
      group_by(., enhanced_tissues, .drop = F) %>%
      summarise(m = n_distinct(ensg_id))
    
    # n is the population size - i.e. the number of genes
    n <- 
      human_gene_class_blood_lymph_enr %>% 
      nrow()
    
    q %>% 
      left_join(k) %>% 
      left_join(m) %>%
      # n is the population size - i.e. the number of genes
      mutate(n = n - m) 
    
  })

xeno_blood_hyper <- 
  xeno_blood_hyper_test_data %>%
  group_by(xeno_id, cor_lim, enhanced_tissues, q, k, .drop = F) %>%
  # Testing the chance of getting the observed number or higher per cluster and tissue type
  summarise(p_value = phyper(q - 1, m, n, k, lower.tail = F)) %>%
  ungroup() %>%
  mutate(p_value = ifelse(p_value == 0, 1e-300, p_value),
         adj_pval = p.adjust(p_value, method = "BH"),
         sign_level = case_when(adj_pval < 0.00000001 ~ 8,
                                adj_pval < 0.0000001 ~ 7,
                                adj_pval < 0.000001 ~ 6,
                                adj_pval < 0.00001 ~ 5,
                                adj_pval < 0.0001 ~ 4,
                                adj_pval < 0.001 ~ 3,
                                adj_pval < 0.01 ~ 2,
                                adj_pval < 0.05 ~ 1,
                                T ~ 0)) 

xeno_blood_hyper_cluster_data <-
  xeno_blood_hyper %>% 
  left_join(xeno_ko_genes,
            by = c("xeno_id" = "enssscg_id")) %>%
  filter(cor_lim == 0.6) %>%
  select(gene_name, enhanced_tissues, sign_level) %>% 
  spread(gene_name, sign_level) %>%
  column_to_rownames("enhanced_tissues")

xeno_hyper_tis_cluster <- 
  xeno_blood_hyper_cluster_data %>% 
  dist() %>%
  hclust(method = "ward.D2")

xeno_hyper_gene_cluster <- 
  xeno_blood_hyper_cluster_data %>% 
  t() %>%
  dist() %>%
  hclust(method = "ward.D2")

xeno_blood_hyper %>%
  
  group_by(xeno_id, enhanced_tissues) %>% 
  mutate(any_sign = any(sign_level != 0)) %>% 
  ungroup() %>%
  filter(any_sign) %>%
  ggplot(aes(cor_lim, -log10(p_value), group = enhanced_tissues)) +
  geom_line(show.legend = F) + 
  geom_text_repel(data = . %>% 
                    group_by(enhanced_tissues, xeno_id) %>%
                    summarise(p_value = min(p_value)),
                  aes(label = enhanced_tissues), 
                  size = 1.5,
                  show.legend = F) +
  geom_line(show.legend = F) + 
  geom_hline(yintercept = -log10(0.05)) +
  facet_wrap(~xeno_id, scales = "free_y") + 
  theme_minimal()+ 
  scale_color_manual(values = c(tissue_colors_palette_full, "not enriched" = "black"))
ggsave(savepath("xeno subnetwork blood p-val per corlim.pdf"), width = 10, height = 6)

xeno_blood_hyper %>% 
  left_join(xeno_ko_genes,
            by = c("xeno_id" = "enssscg_id")) %>%
  mutate(enhanced_tissues = factor(enhanced_tissues, ggdendro::dendro_data(xeno_hyper_tis_cluster)$labels$label),
         gene_name = factor(gene_name, ggdendro::dendro_data(xeno_hyper_gene_cluster)$labels$label)) %>%
  filter(sign_level > 0) %>%
  ggplot(aes(enhanced_tissues, gene_name, fill = sign_level)) +
  geom_tile() +
  facet_wrap(~cor_lim) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "orange", "orangered"))(9)) +
  stripped_theme_facet +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major.x = element_line(color = "lightgray"))

ggsave(savepath("Xeno-KO genes network neighborhood blood 1.pdf"), width = 10, height = 6)

xeno_blood_hyper %>% 
  filter(cor_lim == 0.6) %>%
  left_join(xeno_ko_genes,
            by = c("xeno_id" = "enssscg_id")) %>%
  arrange(adj_pval) %>%
  mutate(enhanced_tissues = factor(enhanced_tissues, ggdendro::dendro_data(xeno_hyper_tis_cluster)$labels$label),
         gene_name = factor(gene_name, ggdendro::dendro_data(xeno_hyper_gene_cluster)$labels$label),
         adj_pval = ifelse(adj_pval < 1e-20, 1e-20, adj_pval)) %>%
  filter(sign_level != 0) %>%
  ggplot(aes(enhanced_tissues, gene_name, 
             fill = -log10(adj_pval), 
             size = -log10(adj_pval))) +
  geom_point(shape = 21,
             color = "black") +
  geom_text(aes(enhanced_tissues, gene_name, label = q), 
            inherit.aes = F,
            vjust = -1) +
  # facet_wrap(~order) +
  scale_fill_gradientn(colours = colorRampPalette(c("white", "orange", "orangered"))(21), 
                       values = (0:20)/20) +
  stripped_theme_facet +
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        panel.grid.major.x = element_line(color = "lightgray"))

ggsave(savepath("Xeno-KO genes network neighborhood blood 2.pdf"), width = 6, height = 4)

# ----- Visualization Xeno neighborhoods -----


xeno_subgraphs_n <- 
  lapply(1:nrow(xeno_subgraph_settings),
         function(i) {
           net_ <-
             xeno_subgraphs[[i]]
           
           edges_ <- 
             net_ %>% 
             activate(edges) %>%
             as_tibble()
           
           nodes_ <- 
             net_ %>% 
             activate(nodes) %>%
             as_tibble()
           
           tibble(edges = nrow(edges_),
                  nodes = nrow(nodes_))
         }) %>% 
  bind_rows() %>% 
  mutate(ratio = edges/nodes) %>%
  bind_cols(xeno_subgraph_settings %>% 
              left_join(xeno_ko_genes))


xeno_subgraphs_n %>% 
  gather(metric, number, 1:3) %>%
  ggplot(aes(cor_lim, number, color = gene_name)) + 
  geom_line() + 
  facet_wrap(~metric, scales = "free_y", ncol = 1) +
  stripped_theme_facet
ggsave(savepath("metric per subgraph.pdf"), width = 5, height = 7)

xeno_subgraphs_n %>% 
  ggplot(aes(edges, nodes, color = gene_name)) + 
  geom_line() + 
  # facet_wrap(~gene_name, scales = "free", nrow = 1) +
  stripped_theme_facet
ggsave(savepath("metric per subgraph2.pdf"), width = 4, height = 3)

# Bundled

plots <-
  lapply(1:nrow(xeno_subgraph_settings),
         function(i) {
           print(i)
           
           enssscg_id_ <- xeno_subgraph_settings$enssscg_id[i]
           cor_lim_ <- xeno_subgraph_settings$cor_lim[i]
           gene_name_ <- 
             xeno_ko_genes %>% 
             filter(enssscg_id == enssscg_id_) %>% 
             pull(gene_name)
           
           
           filtered_edges <- 
             xeno_subgraphs[[i]] %>%
             activate(edges) %>%
             as_tibble() %>% 
             {bind_rows(mutate(., 
                              type = "original"),
                       rename(., 
                              from = to, 
                              to = from) %>%
                         mutate(type = "flipped"))} %>% 
             group_by(from) %>%
             mutate(rank = rank(1 - cor)) %>% 
             ungroup() %>%
             filter(rank <= 10) %>%
             mutate(from_ = ifelse(type == "flipped", 
                                   to, from),
                    to_ = ifelse(type == "flipped", 
                                 from, to)) %>%
             select(from = from_,
                    to = to_) %>% 
             distinct()
           
           graph_data <-
             xeno_subgraphs[[i]]  %>%
             
             activate(edges) %>%
             inner_join(filtered_edges) %>% 
             activate(nodes) %>%
             
             mutate(xeno_gene = name %in% xeno_ko_genes$enssscg_id,
                    display_name = ifelse(display_name == "NULL", name, display_name)) 
           
           if(length(pull(graph_data, name)) == 1) {
             
             g_plot <- 
               ggplot() + 
               ggtitle(gene_name_, 
                       paste("thresh =", cor_lim_))
             
           } else {
             
             graph_clust <- 
               pig_atlas_sample_wide_tmm[pull(graph_data, name),] %>%
               t() %>%
               cor(method = "spearman") %>%
               {1 - .} %>%
               as.dist() %>% 
               
               hclust(method = "ward.D2")
             
             dendrograph <-
               graph_clust %>%
               as.dendrogram() %>%
               as_tbl_graph(directed = T) %>% 
               left_join(cutree(graph_clust, 
                                min(c(ceiling(length(graph_clust$labels) / 10), 10))) %>% 
                           enframe("label", "cluster"))
             
             
             g1 <- 
               dendrograph %>% 
               ggraph(layout = 'dendrogram', circular = T) 
             
             
             
             matched_graph <- 
               graph_data %>% 
               # as.igraph()
               activate(edges) %>%
               mutate(from_name = pull(activate(graph_data, nodes), name)[from],
                      to_name = pull(activate(graph_data, nodes), name)[to],
                      new_from = g1$data$.ggraph.orig_index[match(from_name, pull(g1$data, label))],
                      new_to = g1$data$.ggraph.orig_index[match(to_name, pull(g1$data, label))])
             
             
             from <- 
               matched_graph %>%
               pull(new_from)
             
             to <- 
               matched_graph %>%
               pull(new_to)
             
             #### Color non-leaf nodes as connected leaves:
             reps <- 
               get_nodes()(g1$data) %>%
               as_tibble() %>%
               mutate(radius = x^2 + y^2) %>%
               pull(radius) %>% 
               unique()
             
             
             
             temp <- 
               get_edges()(g1$data) %>%
               as_tibble() %>% 
               select(from = node1..ggraph.orig_index, 
                      to = node2..ggraph.orig_index, 
                      node1.cluster,
                      node2.cluster) 
             
             for(i in reps) {
               temp <- 
                 temp %>%
                 mutate(node1.cluster = node2.cluster)
               
               temp_clusters <- 
                 temp %$%
                 tibble(node = c(from, to),
                        cluster = c(node1.cluster, node2.cluster)) %>% 
                 group_by(node) %>% 
                 summarise(cluster = cluster[1])
               
               temp <- 
                 temp %>%
                 select(-node2.cluster) %>% 
                 left_join(temp_clusters,
                           by = c("to" = "node")) %>%
                 rename(node2.cluster = cluster)
             }
             
             #### Make new graph with assigned clusters on non-leaf nodes
             g <- 
               dendrograph %>% 
               left_join(gene_mapping, 
                         by = c("label" = "enssscg_id")) %>%
               mutate(cluster = temp_clusters$cluster,
                      display_name2 = ifelse(display_name == "NULL", label, display_name),
                      display_name = ifelse(display_name == "NULL", NA, display_name)) %>%
               ggraph(layout = 'dendrogram', circular = T) 
             # ggraph::
             
             node_data <- 
               get_nodes()(g$data)
             
             g_plot <- 
               g +
               
               geom_conn_bundle2(data = get_con(from = from, to = to), 
                                 aes(color = cluster),
                                 width=0.1, alpha=0.2, 
                                 tension = 0.95, 
                                 show.legend = F) +
               geom_node_point(data = node_data %>% 
                                 filter(label == enssscg_id_),
                               show.legend = F, 
                               color = "red",
                               size = 4) +
               geom_node_point(data = node_data %>% 
                                 filter(leaf),
                               aes(color = cluster), 
                               show.legend = F) +
               node_data %>%
               filter(label != "") %>%
               mutate(degree = case_when(x >= 0 ~ asin(y) * 180 / pi,
                                         x < 0 ~ 360 - asin(y) * 180 / pi)) %>%
                         {geom_node_text(data = .,
                                         aes(label = display_name),
                                         angle = .$degree,
                                         hjust = ifelse(.$x < 0, 
                                                        1, 
                                                        0),
                                         vjust = 0.5,
                                         size = 1.5)} +
               theme_void() + 
               scale_edge_colour_distiller(palette = "RdPu") +
               scale_colour_distiller(palette = "RdPu") +
               
               ggtitle(gene_name_, 
                       paste("thresh =", cor_lim_))
           }
           
           
           g_plot
           
           # ggsave(savepath("sa.pdf"), plot = g_plot, width = 10, height = 10)
           
         })

pdf(savepath("xeno-ko neighborhood bundled networks.pdf"), width = 10, height = 10)
plots
dev.off()

# Simple

plots <-
  lapply(1:nrow(xeno_subgraph_settings),
         function(i) {
           print(i)
           
           enssscg_id_ <- xeno_subgraph_settings$enssscg_id[i]
           cor_lim_ <- xeno_subgraph_settings$cor_lim[i]
           gene_name_ <- 
             xeno_ko_genes %>% 
             filter(enssscg_id == enssscg_id_) %>% 
             pull(gene_name)

           filtered_edges <-
             xeno_subgraphs[[i]] %>%
             activate(edges) %>%
             as_tibble() %>%
             {bind_rows(mutate(.,
                               type = "original"),
                        rename(.,
                               from = to,
                               to = from) %>%
                          mutate(type = "flipped"))} %>%
             group_by(from) %>%
             mutate(rank = rank(1 - cor)) %>%
             ungroup() %>%
             filter(rank <= 50) %>%
             mutate(from_ = ifelse(type == "flipped",
                                   to, from),
                    to_ = ifelse(type == "flipped",
                                 from, to)) %>%
             select(from = from_,
                    to = to_) %>%
             distinct()

           graph_data <-
             xeno_subgraphs[[i]]  %>%


             activate(edges) %>%
             inner_join(filtered_edges) %>%
             activate(nodes) %>%
             mutate(xeno_gene = name %in% xeno_ko_genes$enssscg_id,
                    display_name = ifelse(display_name == "NULL", name, display_name))
           
           if(length(pull(graph_data, name)) == 1) {
             
             g_plot <- 
               ggplot() + 
               ggtitle(gene_name_, 
                       paste("thresh =", cor_lim_))
             
           } else {
             
             layout_ <- 
               pig_atlas_sample_wide_tmm[pull(graph_data, name),] %>%
               
               {
                 dt <- log10(. + 1)
                 uwot::umap(dt, n_neighbors = max(round(sqrt(dim(dt)[1])), 2)) %>%
                   as_tibble() %>% 
                   mutate(name = rownames(dt))
               }
             
             
             
             
             graph_ <-
               ggraph(graph_data,
                      layout = "manual",
                      x = layout_$V1,
                      y = layout_$V2)
             
             
             
             G <-
               graph_ +
               geom_edge_arc(color = "gray",
                             width = 0.1,strength = 0.1,
                             alpha = 0.1,
                             show.legend = F)  +
               
               geom_node_point(data = graph_$data %>%
                                 filter(name == enssscg_id_),
                               size = 10,
                               color = "red",
                               show.legend = F) +
               geom_node_point(aes(size = centrality_eigen,
                                   color = weighted_degree),
                               stroke = 3,
                               show.legend = F) +
               geom_node_text(data = graph_$data %>%
                                filter(name == enssscg_id_),
                              aes(label = display_name),
                              size = 4,
                              color = "black",
                              show.legend = F) +
               # geom_node_text(aes(label = display_name),
               #                show.legend = F,
               #                color = "black") +
               theme_graph(base_family="sans") +
               scale_colour_distiller(palette = "RdPu") +
               # scale_color_gradientn(colours = ggsci::pal_material("purple")(10)) +
               scale_size_continuous(range = c(0.2, 5)) +
               # scale_color_manual(values = c("TRUE" = "red", "FALSE" = NA))+
               
               ggtitle(gene_name_, 
                       paste("thresh =", cor_lim_))
           }
           

           # ggsave(savepath("ex.pdf"), plot = G, width = 10, height = 10)
           
           G


         })

pdf(savepath("xeno-ko neighborhood networks.pdf"), width = 10, height = 10)
plots
dev.off()


####


read_csv("doc/network_res_summary.csv") %>%
  mutate_at(-1, .funs = function(x) ifelse(is.na(x), 0, 1)) %>% 
  column_to_rownames("gene_name") %>% 
  pheatmap(clustering_method = "ward.D2")

library(GeneNet)

lapply((.packages()),
       function(x) {
         capture.output(utils:::print.bibentry(citation(x), style = "Bibtex"),
                        file = savepath(paste0("endnote_", x, ".bib")))
       })




