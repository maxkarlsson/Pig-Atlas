

# # Load data
pig_atlas_sample_temp <-
  read_delim("data/processed/pig_filtered_sample_data.tab", delim = "\t")

pig_atlas_tissue <-
  read_delim("data/HPA/consensus_pig_all_tissues_92.tsv", delim = "\t") %>%
  rename(enssscg_id = ensg_id,
         tissue_name = tissue) %>%
  mutate(tissue_name = gsub(" 1$", "", tissue_name) %>%
           gsub("facial adipose tissue", "orbital adipose tissue", .) %>%
           trimws())

## Tissue data
human_atlas_tissue <-
  bind_rows(read_delim("data/human data/lims/consensus_hpa_92.tsv", delim = "\t") %>%
              mutate(source = "hpa"),
            read_delim("data/human data/lims/consensus_gtex_92.tsv", delim = "\t") %>%
              mutate(source = "gtex"),
            read_delim("data/human data/lims/consensus_fantom_92.tsv", delim = "\t") %>%
              mutate(source = "fantom")) %>% 
  rename(tissue_name = celltype) %>%
  group_by(tissue_name, ensg_id) %>% 
  summarise(tpm = mean(tpm, na.rm = T),
            ptpm = mean(ptpm, na.rm = T),
            tmm = mean(tmm, na.rm = T),
            nx = mean(nx, na.rm = T)) %>%
  ungroup()

##Normalization

pig_atlas_sample <- 
  pig_atlas_sample_temp %>% 
  group_by(sample_ID) %>% 
  group_by(sample_ID) %>% 
  mutate(sum_tpm = sum(tpm), 
         ptpm = tpm * 1e6 / sum_tpm) %>% 
  select(-sum_tpm) %>% 
  ungroup()

# Create wide format sample data
pig_atlas_sample_wide_ptpm <- 
  pig_atlas_sample %>%
  select(1, 2, ptpm) %>%
  spread(key = sample_ID, value = ptpm) %>%
  column_to_rownames("enssscg_id")

# TMM normalize sample data
pig_atlas_sample_wide_tmm <- 
  pig_atlas_sample_wide_ptpm %>%
  tmm_norm_median_ref()

# Create long format normalized sample data
pig_atlas_sample_norm <- 
  pig_atlas_sample %>% 
  left_join(pig_atlas_sample_wide_tmm %>% 
              as_tibble(rownames = "enssscg_id") %>%
              gather("sample_ID", "tmm", -1),
            by = c("enssscg_id", "sample_ID")) %>%
  group_by(enssscg_id) %>%
  mutate(nx = case_when(tpm == 0 ~ 0, 
                        T ~ tmm/sqrt(sd(tmm)))) %>% 
  ungroup()


# Aggregate to tissue data


pig_atlas_comparison <- 
  pig_atlas_tissue %>% 
  inner_join(tissue_mapping %>% 
               select(tissue_name, comparison_tissue), 
             by = "tissue_name") %>%
  group_by(enssscg_id, comparison_tissue) %>% 
  summarise(tpm = max(tpm), 
            ptpm = max(ptpm), 
            tmm = max(tmm),
            nx = max(nx)) %>% 
  ungroup() %>% 
  filter(!is.na(comparison_tissue))

human_atlas_comparison <- 
  human_atlas_tissue %>% 
  inner_join(human_tissue_mapping, 
             by = "tissue_name") %>% 
  group_by(comparison_tissue, ensg_id) %>% 
  summarise(tpm = max(tpm, na.rm = T),
            ptpm = max(ptpm, na.rm = T),
            tmm = max(tmm, na.rm = T),
            nx = max(nx, na.rm = T)) %>% 
  ungroup() %>% 
  filter(!is.na(comparison_tissue)) %>%
  filter(comparison_tissue %in% pig_atlas_comparison$comparison_tissue)



# ----- Joined human pig tables

joined_atlas_comparison_temp <- 
  pig_atlas_comparison %>%
  inner_join(gene_orthologs_all) %>% 
  unite(mutual_id, enssscg_id, ensg_id, sep = "_") %>% 
  mutate(species = "pig") %>%
  select(mutual_id, comparison_tissue, species, ptpm) %>%
  
  
  bind_rows(human_atlas_comparison %>%
              inner_join(gene_orthologs_all) %>% 
              unite(mutual_id, enssscg_id, ensg_id, sep = "_") %>% 
              mutate(species = "human") %>%
              select(mutual_id, comparison_tissue, species, ptpm))

joined_atlas_comparison_temp2 <- 
  joined_atlas_comparison_temp %>%
  unite(id, comparison_tissue, species) %>% 
  spread(id, ptpm) %>% 
  gather(id, ptpm, -1) %>% 
  separate(id, into = c("comparison_tissue", "species"), sep = "_", remove = F) %>%
  group_by(mutual_id, species) %>%
  mutate(mean_ptpm = mean(ptpm, na.rm = T),
         imputed = is.na(ptpm),
         imp_ptpm = ifelse(imputed, mean_ptpm, ptpm)) %>%
  ungroup()


joined_atlas_comparison_temp3_tmm <-
  joined_atlas_comparison_temp2 %>% 
  select(mutual_id, id, imp_ptpm) %>% 
  spread(id, imp_ptpm) %>% 
  column_to_rownames("mutual_id") %>%
  tmm_norm_median_ref(doWeighting = F)

joined_atlas_comparison_temp3_limma <-
  joined_atlas_comparison_temp3_tmm %>% 
  {log10(. + 1)} %>% 
  limma::removeBatchEffect(batch = colnames(joined_atlas_comparison_temp3_tmm) %>% 
                             str_extract("_.*$")) %>%
  {10^. - 1}

joined_atlas_comparison <- 
  joined_atlas_comparison_temp2 %>%
  left_join(joined_atlas_comparison_temp3_tmm %>% 
              as_tibble(rownames = "mutual_id") %>% 
              gather(id, mutual_tmm, -1)) %>% 
  left_join(joined_atlas_comparison_temp3_limma %>% 
              as_tibble(rownames = "mutual_id") %>% 
              gather(id, mutual_limma, -1)) %>% 
  mutate(mutual_tmm = ifelse(imputed, NA, mutual_tmm),
         mutual_limma = ifelse(imputed, NA, mutual_limma)) %>% 
  select(1:4, ptpm, mutual_tmm, mutual_limma)





# Save
save(list = c("human_atlas_tissue",
              "pig_atlas_sample",
              "pig_atlas_sample_wide_ptpm",
              "pig_atlas_sample_wide_tmm",
              "pig_atlas_sample_norm",
              "pig_atlas_tissue",
              "pig_atlas_comparison",
              "human_atlas_comparison",
              "joined_atlas_comparison"),
     file = "data/processed/normalized_res.Rdata")
