

sscrofa92 <- 
  read_tsv("data/meta/sscrofa_92.txt") %>% 
  select(gene_id = 1,
         gene_id_version = 2,
         gene_name = `Gene name`,
         gene_type = `Gene type`)

sscrofa103 <- 
  read_tsv("data/meta/sscrofa_103.txt") %>% 
  select(gene_id = 1,
         gene_id_version = 2,
         gene_name = `Gene name`,
         gene_type = `Gene type`)



sscrofa_comb <- 
  list(s92 = sscrofa92,
       s103 = sscrofa103) %>% 
  map(. %>% 
        filter(gene_type == "protein_coding"))


sscrofa_comb$s92 %>% 
  filter(gene_id %in% sscrofa_comb$s103$gene_id)
plot_data <- 
  list(gene_id = sscrofa_comb$s92 %>% 
         select(gene_id, gene_type) %>% 
         full_join(sscrofa103 %>% 
                     select(gene_id, gene_type),
                   by = c("gene_id"),
                   suffix = c("_92", "_103")) %>% 
         group_by(gene_type_92, gene_type_103) %>% 
         count(),
       gene_id_version = sscrofa_comb$s92 %>% 
         select(gene_id_version, gene_type) %>% 
         full_join(sscrofa103 %>% 
                     select(gene_id_version, gene_type),
                   by = c("gene_id_version"),
                   suffix = c("_92", "_103")) %>% 
         group_by(gene_type_92, gene_type_103) %>% 
         count(),
       gene_name = sscrofa_comb$s92 %>% 
         select(gene_name, gene_type) %>% 
         distinct() %>% 
         full_join(sscrofa103 %>% 
                     select(gene_name, gene_type) %>% 
                     distinct(),
                   by = c("gene_name"),
                   suffix = c("_92", "_103")) %>% 
         group_by(gene_type_92, gene_type_103) %>% 
         count()) %>% 
  bind_rows(.id = "type") %>% 
  mutate(overlap = case_when(is.na(gene_type_103) ~ "92 only",
                             is.na(gene_type_92) ~ "103 only",
                             T ~ "overlap"),
         gene_type_switch = case_when(gene_type_103 == "protein_coding" ~ "protein coding",
                                      is.na(gene_type_103) ~ "92 only",
                                      gene_type_92 == "protein_coding" & 
                                        gene_type_103 != "protein_coding" ~ "not protein coding in 103",
                                      
                                      T ~ "not protein coding")) %>% 
  filter(!(is.na(gene_type_92) & gene_type_103 != "protein_coding"))

plot_data %>%
  filter(type == "gene_id") %>% 
  group_by(type, overlap, gene_type_switch) %>% 
  summarise(n = sum(n)) %>% 
  ggplot(aes(gene_type_switch, n, fill = overlap, label = n)) +
  geom_col() +
  # facet_wrap(~gene_type_switch) +
  geom_text(position = position_stack(vjust = 0.5))
ggsave(savepath("92 vs 103 overlap.pdf"),
       width = 5, height = 4)

E92_103_comp <- 
  sscrofa_comb$s92 %>% 
  select(gene_id, gene_type) %>% 
  full_join(sscrofa103 %>% 
              select(gene_id, gene_type),
            by = c("gene_id"),
            suffix = c("_92", "_103")) 
  
protein_coding_E92_only <- 
  E92_103_comp %>% 
  filter(is.na(gene_type_103))

protein_coding_reclassed <- 
  E92_103_comp %>% 
  filter(gene_type_92 == "protein_coding" & 
           gene_type_103 != "protein_coding")

protein_coding_reclassed %>% 
  group_by(gene_type_103) %>% 
  count

plot_data <- 
  protein_coding_E92_only %>% 
  full_join(pig_gene_classification,
            by = c("gene_id" = "enssscg_id")) 



plot_data %>% 
  filter(!is.na(gene_type_92)) %>% 
  group_by(specificity_category) %>% 
  count %T>%
  print() %>% 
  mutate(specificity_category = factor(specificity_category,
                                       spec_category_levels)) %>% 
  ggplot(aes(1, n, fill = specificity_category, label = n)) +
  geom_col(show.legend = F) +
  geom_text(position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = gene_category_pal)
ggsave(savepath("genes spec removed from 103.pdf"),
       width = 3, height = 4)

plot_data %>% 
  separate_rows(enhanced_tissues, sep = ";") %>% 
  group_by(gene_type_92, enhanced_tissues) %>% 
  filter(!is.na(enhanced_tissues)) %>% 
  count %>%  
  ggplot(aes(n, enhanced_tissues, fill = gene_type_92)) +
  geom_col(show.legend = F, 
           position = "fill") 




plot_data <- 
  protein_coding_reclassed %>% 
  full_join(pig_gene_classification,
            by = c("gene_id" = "enssscg_id")) 



plot_data %>% 
  filter(!is.na(gene_type_92)) %>% 
  group_by(specificity_category) %>% 
  count %T>%
  print() %>% 
  mutate(specificity_category = factor(specificity_category,
                                       spec_category_levels)) %>% 
  ggplot(aes(1, n, fill = specificity_category, label = n)) +
  geom_col(show.legend = F) +
  geom_text(position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = gene_category_pal)
ggsave(savepath("genes spec reclassed.pdf"),
       width = 3, height = 4)

plot_data %>% 
  separate_rows(enhanced_tissues, sep = ";") %>% 
  group_by(gene_type_92, enhanced_tissues) %>% 
  filter(!is.na(enhanced_tissues)) %>% 
  count %>%  
  ggplot(aes(n, enhanced_tissues, fill = gene_type_92)) +
  geom_col(show.legend = F, 
           position = "fill") 




list(reclassed = protein_coding_reclassed,
     E92_only = protein_coding_E92_only) %>% 
  bind_rows(.id = "type") %>% 
  left_join(gene_mapping,
            by = c("gene_id" = "enssscg_id")) %>% 
  filter(display_name != "NULL") %>% 
  arrange(display_name) %>% 
  mutate(name_in_103 = display_name %in% sscrofa103$gene_name) %>% 
  filter(name_in_103)


sscrofa92 %>% 
  filter(gene_type == "protein_coding") %>% 
  left_join(protein_coding_reclassed %>% 
              select(1, 3)) %>% 
  mutate(comment = case_when(gene_id %in% protein_coding_E92_only$gene_id ~ "Not in Ensembl 103",
                             !is.na(gene_type_103) ~ paste("Changed to", gene_type_103, "in Ensembl 103"),
                             T ~ "")) %>% 
  select(gene_id, gene_name, gene_type, comment) %>% 
  write_tsv(savepath("E92 E103 comparison comments.tsv"))

#################################################
# umap & uttrycksnivå & chromosom
pig_ge
gene_mapping


gene_chromosome %>% 
  select(-2) %>% 
  filter(`Gene stable ID` %in% protein_coding_E92_only$gene_id) %>%
  distinct() %>% 
  mutate(chromosome_match = `Chromosome/scaffold name` %in% c(1:40, "MT", "X", "Y")) %>% 
  group_by(chromosome_match) %>% 
  count()


pig_atlas_tissue %>% 
  filter(enssscg_id %in% protein_coding_E92_only$gene_id) %>%
  group_by(enssscg_id) %>% 
  summarise(tmm = max(tmm)) %>%
  left_join(pig_gene_classification) %>% 
  mutate(specificity_category = factor(specificity_category, spec_category_levels)) %>% 
  ggplot(aes(tmm, specificity_category, fill = specificity_category)) +
  geom_violin(scale = "width", 
              show.legend = F,
              draw_quantiles = 0.5) +
  scale_x_log10() +
  scale_fill_manual(values = gene_category_pal)
ggsave(savepath("Not E103 genes expression.pdf"),
       width = 4, height = 4)  
  

UMAP_annotations <-
  readxl::read_xlsx("data/processed/UMAP cluster annotation 20201126.xlsx", sheet = 2) %>%
  rename(merged_cluster = `Cluster nr`, 
         n_genes = `n genes`, 
         cluster_name = `Cluster name`) %>%
  mutate(merged_cluster = as.factor(merged_cluster))

pig_gene_umap %>%
  select(enssscg_id, merged_cluster) %>% 
  mutate(E103_removed = ifelse(enssscg_id %in% protein_coding_E92_only$gene_id,
                               "removed", 
                               "both")) %>% 
  group_by(merged_cluster, E103_removed) %>% 
  count %>% 
  ungroup() %>% 
  spread(E103_removed, n) %>% 
  mutate(frac = removed / (removed + both)) %>% 
  arrange(-frac) %>% 
  left_join(UMAP_annotations) %>% 
  select(Fraction_lost = frac, 
         cluster_name) %>% 
  head(12)


  

pig_gene_umap %>%
  select(enssscg_id, merged_cluster) %>% 
  mutate(E103_removed = ifelse(enssscg_id %in% protein_coding_E92_only$gene_id,
                               "removed", 
                               "both")) %>% 
  group_by(merged_cluster, E103_removed) %>% 
  count %>% 
  ungroup() %>% 
  spread(E103_removed, n) %>% 
  mutate(frac = removed / (removed + both)) %>% 
  arrange(-frac) %>% 
  left_join(UMAP_annotations) %>% 
  mutate(cluster_type = ifelse(grepl("uncharacterized|Not detected", cluster_name),
                               "not annotated",
                               "annotated")) %>% 
  group_by(cluster_type) %>% 
  summarise(n = sum(removed))


#########################


protein_coding_E92_only %>% 
  left_join(gene_mapping,
            by = c("gene_id" = "enssscg_id")) %>% 
  filter(display_name != "NULL") %T>% 
  print() %>% 
  filter(display_name %in% sscrofa103$gene_name)

