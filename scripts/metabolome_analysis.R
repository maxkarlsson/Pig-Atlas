


library(tidyverse)
library(magrittr)
library(broom)
library(pathview)
library(pheatmap)
library(patchwork)
select <- dplyr::select
rename <- dplyr::rename
source("scripts/theme.R")
source("scripts/functions_utility.R")
load("data/processed/gene_overlap_comparison.Rdata")

if(file.exists("data/processed/normalized_res.Rdata")) {
  load("data/processed/normalized_res.Rdata")

} else {
  source("scripts/normalization.R", local = T)
  load("data/processed/normalized_res.Rdata")
}



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

######



metabol_pws <- 
  read_delim("data/HPA/metabolic atlas pw genes.txt", delim = "\t") %>% 
  rename(ensg_id = ensembl_gene_id, 
         pathway = location_name)





######

metabol_pws %>% 
  left_join(gene_orthologs_all) %>% 
  select(ensg_id, enssscg_id, ortholog_type) %>% 
  distinct() %>%
  group_by(ensg_id) %>% 
  summarise(n = n_distinct(enssscg_id[which(!is.na(enssscg_id))])) %>% 
  group_by(n) %>% 
  summarise(n_genes = n()) %>% 
  ggplot(aes(n_genes, as.factor(n), label = n_genes)) +
  geom_col() +
  geom_text(hjust = 0) +
  xlab("Number of human genes") +
  ylab("Number of pig orthologs") +
  ggtitle("Overlap of genes in metabolic pathways") +
  stripped_theme +
  scale_x_continuous(expand = expansion(c(0.01, 0.1)))

ggsave(savepath("Number gene overlap metabolism.pdf"), width = 5, height = 5)  

metabol_pws %>% 
  left_join(gene_orthologs_all) %>% 
  group_by(pathway, ensg_id) %>% 
  summarise(n_orths = n_distinct(enssscg_id[which(!is.na(enssscg_id))])) %>%
  mutate(has_orth = ifelse(n_orths > 0, "orth", "no_orth")) %>% 
  group_by(pathway, has_orth, .drop = F) %>% 
  summarise(n = n()) %>%
  mutate(n_tot = sum(n),
         fract = n/n_tot) %>% 
  ungroup() %>% 
  mutate(pathway = factor(pathway, filter(., has_orth == "orth") %$%
                            pathway[order(fract)])) %>% 
  ggplot(aes(fract, pathway, fill = has_orth)) +
  geom_col() +
  stripped_theme +
  theme(legend.position = "bottom") +
  ggtitle("Fraction of genes with pig ortholog per pathway")

ggsave(savepath("fract of genes having ortholog.pdf"), width = 6, height = 15)  

#------ Correlation analysis ------

cor_data <- 
  joined_atlas_comparison %>% 
  select(mutual_id, comparison_tissue, species, mutual_tmm) %>% 
  separate(mutual_id, into = c("enssscg_id", "ensg_id")) %>%
  inner_join(metabol_pws) %>% 
  left_join(gene_orthologs_all %>% 
              select(1:6)) %>%
  mutate(mutual_tmm = log10(mutual_tmm + 1)) %>%
  spread(species, mutual_tmm)


####------- Gene correlation per pw -----


metabol_pws_gene_cor <-
  cor_data %>% 
  group_by(pathway, enssscg_id, ensg_id, pig_gene_name, human_gene_name, ortholog_type, ortholog_confidence) %>% 
  do({
    cor.test(.$human, .$pig, method = "pearson") %>% 
      tidy()
  }) %>%
  ungroup



metabol_pws_gene_cor %>% 
  group_by(pathway) %>%
  mutate(median = median(estimate, na.rm = T)) %>% 
  ungroup() %>% 
  arrange(median) %>%
  mutate(pathway = factor(pathway, unique(pathway))) %>%
  ggplot(aes(pathway, estimate)) +
  geom_boxplot(fill = "lightgray") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
 
ggsave(savepath("Pearson per pw.pdf"), width = 16, height = 10)

####------- Pathway correlation per tissue -----

metabol_pws_tis_cor <-
  cor_data %>% 
  group_by(pathway) %>% 
  mutate(n = n_distinct(ensg_id)) %>% 
  ungroup() %>% 
  filter(n >= 5) %>%
  filter(pathway != "Isolated") %>%
  group_by(pathway, comparison_tissue) %>% 
  do({
    cor.test(.$human, .$pig, method = "pearson") %>% 
      tidy()
  }) %>%
  ungroup %>% 
  mutate(adj_p = p.adjust(p.value, method = "BH"))

metabol_pws_tis_cor %>% 
  mutate(estimate = ifelse(adj_p < 0.05, estimate, 0)) %>%
  select(pathway, comparison_tissue, estimate) %>% 
  
  spread(comparison_tissue, estimate) %>% 
  column_to_rownames("pathway") %>%
  pheatmap(clustering_method = "ward.D2",
           border_color = NA,
           main = "Pearson correlation of pathways in tissues",
           filename = savepath("Metabolic spearman per pw and tissue.pdf"),
           width = 9, 
           height = 16
           )




#########

metabol_pws_tis_prop <-
  cor_data %>% 
  group_by(pathway) %>% 
  mutate(n = n_distinct(ensg_id)) %>% 
  ungroup() %>% 
  filter(n >= 5) %>%
  filter(pathway != "Isolated") %>%
  group_by(pathway, comparison_tissue) %>% 
  summarise(prop = (2*cov(human, pig)) / (var(human) + var(pig))) %>%
  ungroup()

metabol_pws_tis_prop %>% 
  
  select(pathway, comparison_tissue, prop) %>% 
  
  spread(comparison_tissue, prop) %>% 
  column_to_rownames("pathway") %>%
  pheatmap(clustering_method = "ward.D2",
           border_color = NA,
           main = "",
           filename = savepath("Metabolic prop per pw and tissue.pdf"),
           width = 9, 
           height = 16)

# ------ Pathway comparison heatmap ------


gene_overlap_data_summarized <- 
  gene_overlap_data %>%
  left_join(gene_orthologs_all) %>%
  group_by(enssscg_id,
           ensg_id, 
           gene_name = human_gene_name) %>%
  summarise(type = case_when(any(type == "Overlap") ~ "Overlap",
                             any(type == "Pig") & 
                               any(type == "Human") &
                               (!any(spec_category_human %in% c("Low tissue specificity", 
                                                                "Not detected")) &
                                  !any(spec_category_pig %in% c("Low tissue specificity", 
                                                                "Not detected"))) ~ "Different tissues",
                             any(spec_category_human %in% c("Low tissue specificity", 
                                                            "Not detected")) ~ "Pig",
                             any(spec_category_pig %in% c("Low tissue specificity", 
                                                          "Not detected")) ~ "Human")) %>%
  ungroup()

plot_data <-
  joined_atlas_comparison %>%
  select(mutual_id, comparison_tissue, species, mutual_tmm) %>%
  separate(mutual_id, into = c("enssscg_id", "ensg_id")) %>%
  left_join(gene_orthologs_all) %>%
  filter(ortholog_type == "ortholog_one2one" & ortholog_confidence == 1) %>%
  
  # mutate(mutual_tmm = ifelse(mutual_tmm < 1, 0, mutual_tmm)) %>%
  
  group_by(enssscg_id, species) %>% 
  mutate(max_tmm = max(mutual_tmm, na.rm = T),
         mutual_tmm = ifelse(max_tmm == 0, 0, mutual_tmm / max_tmm)) %>%
  left_join(human_gene_mapping) 

plot_data2 <-
  gene_overlap_data_summarized 

plot_data3 <-
  human_pig_gene_class_comparison %>% 
  inner_join(plot_data %>% 
               select(enssscg_id, ensg_id) %>% 
               distinct()) %>%
  select(enssscg_id, ensg_id, spec_category_pig, spec_category_human) %>%
  left_join(human_gene_mapping) %>%
  gather(species, spec_category, spec_category_pig, spec_category_human) %>%
  mutate(species = str_to_sentence(gsub("spec_category_", "", species)))


#####

plots <- 
  lapply(sort(unique(metabol_pws$pathway)),
         function(pw) {
           cat(paste(pw, "\n"))
           filter_ensg <- 
             metabol_pws %>% 
             filter(pathway == pw) %>% 
             pull(ensg_id)
           
           cluster_data <-
             plot_data %>%
             ungroup() %>%
             filter(ensg_id %in% filter_ensg) %>%
             select(comparison_tissue, species, mutual_tmm, gene_name) %>%
             group_by(gene_name, comparison_tissue) %>%
             summarise(mutual_tmm = mean(mutual_tmm)) %>%
             spread(gene_name, mutual_tmm) %>%
             column_to_rownames("comparison_tissue")
           
           # cluster_data <-
           #   plot_data %>% 
           #   ungroup() %>%
           #   filter(ensg_id %in% filter_ensg) %>%
           #   select(comparison_tissue, species, mutual_tmm, gene_name) %>%
           #     spread(species, mutual_tmm) %>%
           #   mutate(diff = human - pig) %>% 
           #     select(-human, -pig) %>% 
           #     spread(gene_name, diff) %>%
           #   column_to_rownames("comparison_tissue") 
           
           tis_cluster <- 
             cluster_data %>%
             dist() %>% 
             hclust(method = "ward.D2")
           
           tis_order <- 
             tis_cluster$labels[tis_cluster$order]
           
           if(dim(cluster_data)[2] >= 2) {
             # gene_cluster <- 
             #   cluster_data %>%
             #   t() %>%
             #   dist() %>% 
             #   hclust(method = "ward.D2")
             
             gene_order <- 
               plot_data %>% 
                 ungroup() %>%
                 filter(ensg_id %in% filter_ensg) %>%
                 select(comparison_tissue, species, mutual_tmm, gene_name) %>%
                   spread(species, mutual_tmm) %>%
               group_by(gene_name) %>% 
               summarise(dist = sqrt(sum((human - pig) ^ 2))) %>% 
               arrange(dist) %>% 
               pull(gene_name)
                 
             
           } else {
             gene_cluster <- 
               tibble(labels = colnames(cluster_data),
                      order = 1)
             
             gene_order <- 
               gene_cluster$labels[gene_cluster$order]
             
           }
           
           
           
           # Plot:
           
           plot_data %>% 
             ungroup() %>%
             filter(ensg_id %in% filter_ensg) %>%
             select(gene_name, species) %>% 
             distinct() %>%
             mutate(species = toupper(substr(species, 1, 1))) %>%
             mutate(gene_name = factor(gene_name, gene_order)) %>%
             
             
             ggplot(aes(1, gene_name, group = species, fill = species, label = species)) +
             geom_tile(show.legend = F, 
                       position = "dodge", 
                       height = 0.9) +
             geom_text(position = position_dodge(width = 1), 
                       size = 2, 
                       fontface = "bold",
                       color = "white") +
             scale_fill_manual(values = c("H" = "#DB1F48",
                                          "P" = "#01949A")) +
             scale_x_discrete(limits = "Species") +
             theme_void() +
             theme(plot.margin = unit(c(0,0,0,0), "mm"),
                   axis.text.x = element_text(angle = 90, hjust = 0, size = 6, vjust = 0.2)) +
             
             # Specificity
             plot_data3 %>%
             filter(ensg_id %in% filter_ensg) %>%
             mutate(gene_name = factor(gene_name, gene_order)) %>%
             
             ggplot(aes(1, gene_name, group = species, fill = spec_category, label = species)) +
             geom_tile(show.legend = F,
                       position = "dodge", 
                       height = 0.9, width = 1) +
             # geom_text(position = position_dodge(width = 1), 
             #           size = 2, 
             #           fontface = "bold",
             #           color = "white") +
             scale_fill_manual(values = gene_category_pal) +
             scale_x_discrete(limits = "Specificity") +
             theme_void() +
             theme(plot.margin = unit(c(0,0,0,0), "mm"),
                   axis.text.x = element_text(angle = 90, hjust = 0, size = 6, vjust = 0.2)) +

             
             # Overlap
             plot_data2 %>%
             filter(ensg_id %in% filter_ensg) %>%
             mutate(gene_name = factor(gene_name, gene_order)) %>%
             
             ggplot(aes(1, gene_name, fill = type)) +
             geom_tile(show.legend = F,
                       height = 0.9) +
             scale_fill_manual(values = overlap_type_pal) +
             scale_x_discrete(limits = "Overlap") +
             theme_void() +
             theme(axis.text.x = element_text(angle = 90, hjust = 0, size = 6, vjust = 0.2)) + 
             
             # Heatmap
             plot_data %>% 
             ungroup() %>%
             filter(ensg_id %in% filter_ensg) %>%
             mutate(gene_name = factor(gene_name, gene_order),
                    comparison_tissue = factor(comparison_tissue, tis_order)) %>%
             
             ggplot(aes(gene_name, comparison_tissue, group = species, fill = mutual_tmm)) +
             geom_tile(position= "dodge", width = 0.9) +
             coord_flip() +
             scale_fill_viridis_c(option = "B", direction = -1, name = "Expression") +
             # facet_grid(ensg_id ~ .) + 
             # facet_wrap(~gene_name, ncol = 1, strip.position = "left") + 
             stripped_theme_facet +
             theme(axis.text.x = element_text(angle = 60, hjust = 1),
                   strip.placement = "outside",
                   strip.text.y.left = element_text(angle = 0),
                   panel.spacing = unit(0, "mm"),
                   legend.position = "right",
                   panel.border = element_rect(fill = NA, color = "white"),
                   axis.title = element_blank(), 
                   axis.line = element_blank()) +
             theme(plot.margin = unit(c(0,0,0,0), "mm")) +
             ggtitle(pw) +
             
             
             plot_layout(widths = c(0.025, 0.025, 0.025, 1), nrow = 1)
         })
  
plot_heights <- 
  plot_data %>% 
  ungroup() %>% 
  select(ensg_id) %>% 
  distinct() %>% 
  inner_join(metabol_pws)  %>%
  group_by(pathway) %>%
  summarise(n = n()) %>%
  arrange(pathway) %>% 
  mutate(height = 2 + 0.3*n / 2)  
  

pdf_multiple(plots, 
             filepath = savepath("metabolism pathways expression heatmaps.pdf"), 
             widths = rep(6, nrow(plot_heights)), 
             heights = plot_heights$height)




######



# ----- Visualize KEGG pathways ------
  
mapIds(org.Hs.eg.db, 
       keys = metabol_pws %>% 
         filter(pathway == "Steroid metabolism") %>% 
         pull(ensg_id),
       column=c("ENTREZID"), 
       keytype ="ENSEMBL")

# keytypes(org.Hs.eg.db)
# head(keys(org.Hs.eg.db, keytype = "GENENAME"))
data("gse16873.d")

pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id = demo.paths$sel.paths[i],
                   species = "hsa", out.suffix = "gse16873", kegg.native = F,
                   sign.pos = demo.paths$spos[i])


