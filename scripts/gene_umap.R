
library(uwot)
library(multidplyr)
library(FNN)

umap_data_temp <- 
  pig_atlas_sample_norm %>% 
  mutate(log_tmm = log10(tmm + 1)) %>% 
  group_by(enssscg_id) %>% 
  mutate(mean = mean(log_tmm), 
         sd = sd(log_tmm)) %>% 
  ungroup() %>% 
  mutate(z_score = (log_tmm - mean) / sd) %>%
  select(enssscg_id, sample_ID, z_score) %>% 
  mutate(z_score = ifelse(is.na(z_score), 0, z_score)) %>%
  spread(sample_ID, z_score) %>% 
  column_to_rownames("enssscg_id")




pig_gene_pca <-
  umap_data_temp %>%
  pca_calc(npcs = 200)

pig_gene_pca$scores %>%
  as_tibble() %>%
  ggplot(aes(PC1, PC3)) +
  geom_hex(aes(fill = stat(log10(count))),
           bins = 100) +
  scale_fill_viridis_c() +
  theme_minimal()


#####

# cluster <- new_cluster(6)
# 
# cluster_copy(cluster, names = c("%>%", "umap_data_temp", "pig_gene_pca", "set_rownames") )
# 


pig_gene_umap_all <- 
  expand.grid(n_neighbors = c(15), #3*(2:15),
              scale = c("none"),
              type = c("original"),
              n_epochs = 1000,
              n_components = c(2, 3, 10)) %>%
  group_by(n_neighbors, scale, type, n_epochs, n_components) %>%
  # partition(cluster) %>%
  do({
    n_neigh = .$n_neighbors
    scale_ = .$scale
    n_epochs_ = .$n_epochs
    type_ = .$type
    n_comps_ = .$n_components
    
    if(type_ == "PCA") {
      dat_ <- pig_gene_pca$scores[,1:50]#pig_gene_pca$scores[,1:pig_gene_pca$pc_95]
    } else if(type_ == "original") {
      dat_ <- umap_data_temp
    }
    
    
    set.seed(42)
    
    dat_ %>%
      uwot::umap(n_neighbors = n_neigh,
                 scale = scale_,
                 n_epochs = n_epochs_,
                 n_components = n_comps_) %>%
      set_rownames(rownames(dat_)) %>%
      tibble::as_tibble(rownames = "enssscg_id")
  })  %>%
  # collect() %>%
  ungroup()



# 
# rm(cluster)


####

pig_gene_umap_temp <- 
  pig_gene_umap_all %>% 
  filter(n_components == 2)

pig_gene_umap_cl <- 
  pig_gene_umap_all %>% 
  filter(n_components == 2) %>%
  select(enssscg_id, V1, V2) %>%
  # select(enssscg_id, matches("^V\\d")) %>%
  column_to_rownames("enssscg_id")


# k = 10
# set.seed(42)
# lc.norm <- 
#   get.knn(pig_gene_umap_cl, k = k) %$%
#   data.frame(from = rep(1:nrow(nn.index), 
#                         k), 
#              to = as.vector(nn.index), 
#              weight = 1/(1 + as.vector(nn.dist))) %>%
#   graph_from_data_frame(directed = FALSE) %>%
#   igraph::simplify() %>%
#   cluster_louvain()
# 
# pig_gene_louvain <- 
#   tibble(enssscg_id = rownames(pig_gene_umap_cl),
#          louvain = as.factor(membership(lc.norm)))




pig_gene_cluster_density <- 
  pig_gene_umap_cl %>% 
  # head(1000) %>%
  fpc::dbscan(0.1) %$%
  cluster %>%
  enframe("enssscg_id", "cluster") %>% 
  mutate(cluster = factor(cluster),
         enssscg_id = rownames(pig_gene_umap_cl))

set.seed(42)
km_coordinates <- 
  pig_gene_umap_cl %>% 
  as_tibble(rownames = "enssscg_id") %>%
  left_join(pig_gene_cluster_density) %>%
  group_by(cluster) %>% 
  mutate(range_V1 = max(V1) - min(V1),
         range_V2 = max(V2) - min(V2),
         approx_size = range_V1 * range_V2,
         mean_V1 = mean(V1),
         mean_V2 = mean(V2)) %>%
  ungroup() %>%
  filter(approx_size != max(approx_size)) 



km_centers_1 <-
  km_coordinates %>% 
  select(cluster, mean_V1, mean_V2) %>%
  distinct()
set.seed(42)
km_centers_2 <- 
  km_coordinates %>% 
  filter(approx_size > 10) %>% 
  
  group_by(cluster) %>%
  
  mutate(row_n = row_number()) %>% 
  filter(row_n %in% sample(1:max(row_n), 3 * unique(round(approx_size/10)))) %>%
  select(cluster, 
         mean_V1 = V1, 
         mean_V2 = V2) %>% 
  ungroup() %>%
  mutate(cluster = factor(row_number() + max(as.numeric(km_centers_1$cluster))))

km_centers <-
  bind_rows(km_centers_1,
            km_centers_2) %>%
  column_to_rownames("cluster")

set.seed(42)
pig_gene_cluster_km <-
  pig_gene_umap_cl %>% 
  # head(1000) %>%
  kmeans(centers = km_centers, iter.max = 50) %$%
    cluster %>%
  enframe("enssscg_id", "cluster") %>% 
  mutate(cluster = factor(cluster))


pig_gene_cluster <- pig_gene_cluster_km
                                                     
# pig_gene_cluster <- 
#   pig_gene_cluster_km %>% 
#   left_join(pig_gene_cluster_density,
#             by = "enssscg_id",
#             suffix = c("_km", "_density")) %>% 
#   mutate(cluster = factor(unclass(factor(paste(cluster_km, cluster_density)))))

pig_gene_umap <- 
  pig_gene_umap_temp %>%
  left_join(pig_gene_cluster) %>% 
  group_by(cluster_old = cluster) %>%
  mutate(V1_mean = mean(V1),
         V2_mean = mean(V2)) %>% 
  ungroup() %>%
  arrange(-V2_mean,
          V1_mean) %>%
  mutate(cluster = factor(unclass(factor(cluster_old, unique(cluster_old)))))
  

# plot_data <-
#   pig_gene_umap %>%
#   left_join(pig_gene_classification)
# # 
# # ####
# plot_data %>%
#   group_by(cluster) %>%
#   mutate(mean_V1 = mean(V1),
#          mean_V2 = mean(V2),
#          dist = sqrt((V1 - mean_V1)^2 + (V2 - mean_V2)^2),
#          n = n()) %>%
#   # select(V1, V2, dist, louvain, n) %>%
#   arrange(cluster, dist) %>%
#   mutate(cum_member = cumsum(1/n)) %>%
#   # filter(cum_member < 0.9) %>%
#   mutate(hull = row_number() %in% chull(V1, V2)) %>%
#   filter(hull) %>%
#   mutate(hull = order(chull(V1, V2))) %>%
#   arrange(cluster, hull) %>%
#   ggplot(aes(V1, V2)) +
#   geom_hex(aes(fill = 1),
#            data = plot_data,
#            fill = "gray90",
#            bins = plot_bins) +
#   geom_encircle(aes(fill = cluster,
#                     color = cluster),
#                 alpha = 0.1,
#                 s_shape=1,
#                 expand=0,
#                 # color = "white",
#                 show.legend = F) +
#   geom_encircle(aes(color = cluster),
#                 fill = NA,
#                 alpha = 0.7,
#                 s_shape=1,
#                 expand=0,
#                 # color = "white",
#                 show.legend = F) +
#   geom_text(data = . %>%
#               select(cluster, mean_V1, mean_V2, n) %>%
#               distinct(),
#             aes(mean_V1, mean_V2, label = cluster),
#             show.legend = F) +
#   geom_point(data = km_centers,
#              aes(mean_V1, mean_V2),
#              inherit.aes = F) +
#   stripped_theme

# ----- Hypergeometric test ------

pig_gene_cluster_enrich_hyper_test_data <- 
  pig_gene_umap %>%
  # head(1000) %>%
  select(enssscg_id, cluster) %>%
  left_join(pig_gene_classification) %>% 
  mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues), "not enriched", enhanced_tissues)) %>%
  separate_rows(enhanced_tissues, sep = ";") %>%
  
  {
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
    
  } 

pig_gene_cluster_enrich_hyper <- 
  pig_gene_cluster_enrich_hyper_test_data %>%
  group_by(enhanced_tissues, cluster) %>%
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
                                T ~ 0),
         log_p = -log10(adj_pval))

# 
# a <- 
#   pig_gene_cluster_enrich_hyper %>% 
#   select(1:2, log_p) %>%
#   spread(enhanced_tissues, log_p) %>% 
#   column_to_rownames("cluster") %>% 
#   umap_calc(npcs = 2)
# 
# a %>%
#   as_tibble(rownames = "cluster") %>% 
#   ggplot(aes(V1, V2, label = cluster)) +
#   geom_text()
# 

# ----- Cluster merging -----

pig_gene_umap_cluster_props <- 
  pig_gene_umap %>% 
  select(cluster, enssscg_id, V1, V2) %>%
  group_by(cluster) %>% 
  summarise(radius = sqrt((max(V1) - min(V1))^2 +
                            (max(V2) - min(V2))^2) / 2,
            V1 = mean(V1),
            V2 = mean(V2)) 

max(pig_gene_umap_cluster_props$radius)

pig_gene_umap_cluster_props %>%
  select(cluster, V1, V2) %>%
  column_to_rownames("cluster") %>% 
  dist() %>%
  as.matrix() %>% 
  as_tibble(rownames = "cluster1") %>% 
  gather(cluster2, dist, -1) %>%
  ggplot(aes(dist)) +
  geom_density()

pig_gene_umap_cluster_close <- 
  pig_gene_umap_cluster_props %>%
  select(cluster, V1, V2) %>%
  column_to_rownames("cluster") %>% 
  dist() %>%
  as.matrix() %>% 
  as_tibble(rownames = "cluster1") %>% 
  gather(cluster2, dist, -1) %>%
  filter(cluster1 < cluster2 & dist < 2 & cluster1 != cluster2)
  
pig_gene_umap_hull <- 
  pig_gene_umap %>% 
  group_by(cluster) %>% 
  mutate(row_n = row_number()) %>%
  filter(row_n %in% chull(V1, V2)) %>% 
  select(enssscg_id, V1, V2, cluster)

pig_gene_umap_cluster_adjacent <- 
  pig_gene_umap_cluster_close %>% 
  select(cluster1, cluster2) %>% 
  left_join(pig_gene_umap_hull ,
            by = c("cluster1" = "cluster")) %>% 
  left_join(pig_gene_umap_hull,
            by = c("cluster2" = "cluster"), 
            suffix = c("1", "2")) %>%
  mutate(dist = sqrt((V11 - V12)^2 + (V21 - V22)^2)) %>%
    filter(dist < 0.1) %>%
    select(cluster1, cluster2) %>% 
    distinct() %>% 
    bind_rows(tibble(cluster1 = .$cluster2,
                     cluster2 = .$cluster1)) %>%
    rename(cluster = cluster1, 
           adjacent_cluster = cluster2)

plot_data <- 
  pig_gene_umap %>% 
  left_join(pig_gene_umap_cluster_adjacent)

plots <- 
  lapply(unique(plot_data$cluster),
         function(clus_) {
           # clus_ <- 1
           
           adjacent_clusters <- 
             pig_gene_umap_cluster_adjacent %>% 
             filter(cluster == clus_)
           
           plot_data %>% 
             mutate(cluster_type = case_when(cluster == clus_ ~ "selected",
                                             cluster %in% adjacent_clusters$adjacent_cluster ~ "adjacent",
                                             T ~ "not adjacent")) %>%
             group_by(cluster) %>% 
             mutate(mean_V1 = mean(V1),
                    mean_V2 = mean(V2),
                    dist = sqrt((V1 - mean_V1)^2 + (V2 - mean_V2)^2),
                    n = n()) %>% 
             # select(V1, V2, dist, cluster, n) %>% 
             arrange(cluster, dist) %>% 
             # mutate(cum_member = cumsum(1/n)) %>% 
             # filter(cum_member < 0.9) %>%
             # mutate(hull = row_number() %in% chull(V1, V2)) %>% 
             # filter(hull) %>%
             # mutate(hull = order(chull(V1, V2))) %>%
             # arrange(cluster, hull) %>%
             ggplot(aes(V1, V2)) +
             geom_hex(aes(fill = 1),
                      data = plot_data,
                      fill = "gray90",
                      bins = plot_bins) +
             geom_encircle(aes(group = cluster, 
                               fill = cluster_type, 
                               color = cluster_type),
                           alpha = 0.1, 
                           s_shape=1, 
                           expand=0,
                           # color = "white",
                           show.legend = F) +
             geom_encircle(aes(group = cluster, 
                               color = cluster_type),
                           fill = NA,
                           alpha = 0.7, 
                           s_shape=1, 
                           expand=0,
                           # color = "white",
                           show.legend = F) +
             geom_text(data = . %>% 
                         select(cluster, mean_V1, mean_V2, n) %>% 
                         distinct(),
                       aes(mean_V1, mean_V2, label = cluster),
                       size = 2,
                       alpha = 0.6,
                       show.legend = F) +
             stripped_theme +
             scale_fill_manual(values = c("selected" = "green3", "adjacent" = "red", "not adjacent" = "gray"))
         })

pdf(savepath("UMAP adjacent clusters.pdf"), width = 6, height = 6)
plots
dev.off()



pig_gene_cluster_enrich_hyper %>% 
  select(enhanced_tissues, cluster, log_p) %>%
  group_by(enhanced_tissues, cluster) %>%
  mutate(log_p = min(log_p, 20)) %>%
  ungroup() %>% 
  mutate(log_p = ifelse(log_p < 3, 0, log_p)) %>%
  spread(cluster, log_p) %>% 
  column_to_rownames("enhanced_tissues") %>%
  dist() %>% 
  hclust(method = "ward.D2")


# Save
save(list = c("pig_gene_umap",
              "pig_gene_cluster_enrich_hyper_test_data",
              "pig_gene_cluster_enrich_hyper"),
     file = "data/processed/gene_umap_res.Rdata")

