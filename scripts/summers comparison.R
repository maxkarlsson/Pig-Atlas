

library(tidyverse)
library(clusterProfiler)
library(aricode)

summer_res_all <- 
  list(t2_all = readxl::read_excel("data/temp/Table_2_Functional Annotation of the Transcriptome of the Pig, Sus scrofa, Based Upon Network Analysis of an RNAseq Transcriptional Atlas.xlsx",
                                   sheet = 2)%>% 
         slice(-(1:5)) %>% 
         select(gene_name = 1, 
                enssscg_id = 2, 
                cluster = 6),
       t3_CNS = readxl::read_excel("data/temp/Table_3_Functional Annotation of the Transcriptome of the Pig, Sus scrofa, Based Upon Network Analysis of an RNAseq Transcriptional Atlas.xlsx",
                                   sheet = 2) %>% 
         slice(-(1:5)) %>% 
         select(gene_name = 1, 
                enssscg_id = 3, 
                cluster = 7),
       t4_liver = readxl::read_excel("data/temp/Table_4_Functional Annotation of the Transcriptome of the Pig, Sus scrofa, Based Upon Network Analysis of an RNAseq Transcriptional Atlas.xlsx",
                                     sheet = 2)%>% 
         slice(-(1:5)) %>% 
         select(gene_name = 1, 
                enssscg_id = 3, 
                cluster = 7)) %>% 
  bind_rows(.id = "table") %>% 
  
  mutate(cluster = gsub("Cluster0", "", cluster)) %>% 
  mutate(term = cluster, 
         term_id = cluster)

summer_res <- 
  summer_res_all %>% 
  filter(enssscg_id %in% pig_gene_umap$enssscg_id)
  

summer_anno <- 
  readxl::read_excel("data/temp/Table_2_Functional Annotation of the Transcriptome of the Pig, Sus scrofa, Based Upon Network Analysis of an RNAseq Transcriptional Atlas.xlsx",
                     sheet = 3)%>% 
  slice(-(1:5)) %>% 
  select(cluster = 1, 
         n_genes = 2, 
         cluster_name = 9)

UMAP_annotations <-
  readxl::read_xlsx("data/processed/UMAP cluster annotation 20201126.xlsx", sheet = 2) %>%
  rename(merged_cluster = `Cluster nr`, 
         n_genes = `n genes`, 
         cluster_name = `Cluster name`) %>%
  mutate(merged_cluster = as.factor(merged_cluster))
# readxl::read_excel("data/temp/Table_2_Functional Annotation of the Transcriptome of the Pig, Sus scrofa, Based Upon Network Analysis of an RNAseq Transcriptional Atlas.xlsx",
#                    sheet = 1)
# 
summer_res %>%
  filter(table == "t2_all") %>% 
  filter(cluster %in% summer_anno$cluster) %>% 
  left_join(pig_gene_umap %>%
              filter(enssscg_id %in% summer_res$enssscg_id) %>% 
              select(enssscg_id, merged_cluster)) %>%
  group_by(table) %>%
  summarise(ari = ARI(cluster, merged_cluster))
  

summer_enrichment_res <- 
  summer_res %>% 
  filter(cluster %in% summer_anno$cluster) %>% 
  filter(table == "t2_all") %>% 
  group_by(table) %>% 
  do({
    # For each database
    database <- .
    
    term2gene <- 
      database %>% 
      select(term_id, enssscg_id)
    
    # For 
    pig_gene_umap %>% 
      filter(enssscg_id %in% summer_res$enssscg_id) %>% 
      group_by(merged_cluster) %>% 
      do({
        g_data <- .
        
        pull(g_data, enssscg_id) %>% 
          enricher(maxGSSize = Inf, 
                   universe = unique(pig_gene_umap$enssscg_id),
                   TERM2GENE = term2gene) %>% 
          as_tibble()
      }) %>% 
      ungroup()
    
    
  }) %>% 
  ungroup()

summer_enrichment_res %>% 
  arrange(p.adjust) 




# summer_res %>% 
#   filter(table == "t2_all") %>% 
#   group_by(cluster) %>% 
#   count() %>% 
#   ungroup() %>% 
#   mutate(n_cum = cumsum(n)) %>% View
#   n_distinct()


plot_data <- 
  summer_enrichment_res %>% 
  select(-1, -Description) %>% 
  rename(summer_cluster = ID) %>% 
  left_join(UMAP_annotations %>% 
              select(merged_cluster, cluster_name)) %>% 
  left_join(summer_anno %>% 
              select(summer_cluster = cluster,
                     summer_cluster_name = cluster_name))

plot_data %>% 
  pull(2) %>% 
  n_distinct

plot_data_wide <- 
  plot_data %>% 
  select(merged_cluster, summer_cluster ) %>% 
  mutate(i = 1) %>% 
  spread(summer_cluster , i, fill = 0) %>% 
  column_to_rownames("merged_cluster") 
  
plot_clusts <- 
  plot_data_wide %>% 
  list(merged_cluster = ., 
       ID = t(.)) %>% 
  map(. %>% 
        dist(method = "binary") %>% 
        hclust(method = "ward.D2") %>% 
        with(labels[order]))

plot_data %>% 
  mutate(summer_cluster = factor(summer_cluster, plot_clusts$ID),
         merged_cluster = factor(merged_cluster, plot_clusts$merged_cluster))%>% 
    ggplot(aes(merged_cluster, summer_cluster, size = -log10(p.adjust))) +
    geom_point()
  


plot_data %>% 
  select(cluster_name, summer_cluster_name) %>% 
  distinct() %>% 
  arrange(cluster_name) %>% 
  View

plot_data %>% 
  separate(GeneRatio, into = c("overlap", "pot_overlap"), sep = "/", convert = T, remove = F) %>% 
  separate(BgRatio, into = c("bg", "tot_bg"), sep = "/", convert = T, remove = F) %>% 
  mutate(overlap_frac = overlap / pot_overlap,
         odds_ratio = overlap_frac / (bg / tot_bg)) %>% 
  filter(overlap > 10) %>% 
  select(1,2, GeneRatio, BgRatio, qvalue, cluster_name, 
         summer_cluster_name,
         overlap_frac, odds_ratio) %>% 
  arrange(-odds_ratio) %>% View



library(tidygraph)
library(ggraph)

graph_data <- 
  plot_data %>% 
  separate(GeneRatio, into = c("overlap", "pot_overlap"), sep = "/", convert = T, remove = F) %>% 
  separate(BgRatio, into = c("bg", "tot_bg"), sep = "/", convert = T, remove = F) %>% 
  mutate(overlap_frac = overlap / pot_overlap,
         odds_ratio = overlap_frac / (bg / tot_bg)) %>% 
  filter(overlap > 10) %>% 
  select(1,2, GeneRatio, BgRatio, qvalue, 
         cluster_name, 
         summer_cluster_name,
         overlap_frac, odds_ratio) %>% 
  select(cluster_name, 
         summer_cluster_name,
         everything()) 

graph_data %>%
  arrange(-odds_ratio)


graph_data %>% 
  as_tbl_graph() %>% 
  left_join(graph_data %>% 
              select(1,2) %>% 
              gather(type, name) %>% 
              mutate(type = ifelse(type == "cluster_name",
                                   "Ours",
                                   "Summers")) %>% 
              distinct()) %>% 
  ggraph(layout = "kk") +
  geom_edge_link(aes(alpha = log2(odds_ratio))) +
  geom_node_point(aes(color = type)) +
  geom_node_text(aes(label = name),
                 size = 1.5)+
  scale_edge_width(range = c(1, 2)) +
  scale_edge_alpha(range = c(0.1, 1)) +
  theme_void()
ggsave(savepath("Summers network.pdf"),
       width = 10, height = 7)



graph_data %>% 
  write_tsv(savepath("Summers network data.tsv"))

graph_data %>% 
  arrange(-overlap_frac) %>% 
  slice(7) %>% 
  head


summer_res_all %>%
  filter(table == "t2_all") %>% 
  group_by(cluster) %>% 
  count %>% 
  arrange(n)

pig_gene_umap %>% 
  filter(merged_cluster == 11)

summer_res %>%
  filter(table == "t2_all") %>%
  filter(cluster == "003")
  
plot_data$merged_cluster %>% n_distinct()

plot_data %>% 
  unite(our_cluster, merged_cluster, cluster_name) %>% 
  unite(summers_cluster, summer_cluster, summer_cluster_name) %>% 
  select(-Count) %>% 
  arrange(qvalue) %>% 
  write_tsv(savepath("Summers hypergeom.tsv"))
         
