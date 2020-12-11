
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
  mutate(z_score = (tmm)) %>%
  select(enssscg_id, sample_ID, z_score) %>% 
  mutate(z_score = ifelse(is.na(z_score), 0, z_score)) %>%
  spread(sample_ID, z_score) %>% 
  column_to_rownames("enssscg_id")




pig_gene_pca <-
  umap_data_temp %>%
  t() %>%
  pca_calc(npcs = 200)

pig_gene_pca$loadings

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
              type = c("PCA"),
              n_epochs = 3000,
              n_components = c(2, 3)) %>%
  group_by(n_neighbors, scale, type, n_epochs, n_components) %>%
  # partition(cluster) %>%
  do({
    n_neigh = .$n_neighbors
    scale_ = .$scale
    n_epochs_ = .$n_epochs
    type_ = .$type
    n_comps_ = .$n_components
    
    if(type_ == "PCA") {
      dat_ <- pig_gene_pca$loadings[,1:20]#pig_gene_pca$scores[,1:pig_gene_pca$pc_95]
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

pig_gene_umap_all %>%
  ggplot(aes(V1, V2)) +
  geom_hex(aes(fill = stat(log10(count))),
           bins = 120) +
  scale_fill_viridis_c() +
  stripped_theme +
  facet_grid(n_neighbors ~ n_components)

pig_gene_umap_all %>%
  filter(n_components == 3 & n_neighbors == 15) %$%
  car::scatter3d(V1, V2, V3,
                 surface=FALSE, point.col = "blue")


x_range <- 
  pig_gene_umap_all$V1 %>% 
  range %>%
  {.[2] - .[1]}

y_range <- 
  pig_gene_umap_all$V2 %>% 
  range %>%
  {.[2] - .[1]}

plot_density <- 
  pig_gene_umap_all %>%
  filter(n_components == 2 & n_neighbors == 15) %$% 
  MASS::kde2d(V1, V2, n = 1000, h = 0.5, 
              lims = c(min(V1) - 0.05*x_range, max(V1) + 0.05*x_range,
                       min(V2) - 0.05*y_range, max(V2) + 0.05*y_range)) %>% 
  {.$z} %>%
  as_tibble(rownames = "x") %>%
  gather(y, density, -1) 

x_s_range <- 
  pig_gene_umap_all %>%
  filter(n_components == 2 & n_neighbors == 15) %>%
  pull(V1) %>%
  range

y_s_range <- 
  pig_gene_umap_all %>%
  filter(n_components == 2 & n_neighbors == 15) %>%
  pull(V2) %>%
  range

plot_data <-
  plot_density %>% 
  mutate(y = gsub("V", "", y),
         x = as.numeric(x),
         y = as.numeric(y),
         group = as.numeric(cut(log10(density + 1), breaks = 10)),
         group_factor = factor(group, 1:10),
         color = colorRampPalette(c('red','blue'))(10)[group]) %>%
  arrange(x,y) %>%
  mutate(x_s = scales::rescale(x, x_s_range + c(-0.05*x_range, 0.05*x_range)),
         y_s = scales::rescale(y, y_s_range + c(-0.05*y_range, 0.05*y_range)))


library(rgl)

plot_data%>%
  
  # mutate(density = round(density, 4)) %>%
  # filter(density > 0.005) %>%
  # head(1000) %>%
  select(x, y, density) %>%
  spread(y, density, fill = NA) %>% 
  column_to_rownames("x") %>%
  as.matrix() %>%
  persp3d(col = "green", aspect = c(1, 1, 0.1),
          front = "lines",
          back = "lines")


plot_x <- unique(plot_data$x_s)

plot_y <- unique(plot_data$y_s)

plot_z <-
  plot_data %>%
  select(x, y, density) %>%
  spread(y, density, fill = NA) %>% 
  column_to_rownames("x") %>%
  as.matrix()

plot_col <-
  ggthemes::tableau_color_pal(palette = "Hue Circle", 
                              type = c("regular"),
                              direction = 1)(19) %>%
  prismatic::color() %>%
  {.[cut(plot_z, 19)]}
# zcol2 = as.numeric(apply(z,2, mycut, breaks=nbcol))
persp3d(x = plot_x, 
        y= plot_y, 
        z = plot_z,
        xlim = c(-10, max(plot_x)),
        ylim = c(-10, max(plot_y)),
        zlim = c(0.0025, max(plot_z)),background = "black",
        color = plot_col,
        aspect = c(1, 1, 0.05),
        
        front = "lines",
        back = "lines")

pig_gene_umap_all %>%
  filter(n_components == 2 & n_neighbors == 15) %$%
  plot3d(x = V1,
         y = V2, 
         z = 0, add = T,
         xlim = c(-10, max(plot_x)),
         ylim = c(-10, max(plot_y)),
         size = 2)
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

pig_gene_umap_temp2 <- 
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
  pig_gene_umap_temp2 %>%
  # head(1000) %>%
  select(enssscg_id, cluster) %>%
  left_join(pig_gene_classification) %>% 
  mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues), 
                                   "not enriched", enhanced_tissues)) %>%
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
  mutate(p_value = ifelse(p_value == 0, .Machine$double.xmin, p_value),
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
  pig_gene_umap_temp2 %>% 
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
    filter(dist < 0.2) %>%
    select(cluster1, cluster2) %>% 
    distinct() %>% 
    bind_rows(tibble(cluster1 = .$cluster2,
                     cluster2 = .$cluster1)) %>%
    rename(cluster = cluster1, 
           adjacent_cluster = cluster2)

plot_data <- 
  pig_gene_umap_temp2 %>% 
  left_join(pig_gene_umap_cluster_adjacent) %>%
  group_by(cluster) %>% 
  mutate(mean_V1 = mean(V1),
         mean_V2 = mean(V2),
         dist = sqrt((V1 - mean_V1)^2 + (V2 - mean_V2)^2),
         n = n()) %>% 
  arrange(cluster, dist) 

plot <- 
  plot_data %>%
  ggplot(aes(V1, V2)) +
  geom_hex(aes(fill = 1),
           data = plot_data,
           fill = "gray90",
           bins = plot_bins) +
  geom_text(data = . %>% 
              select(cluster, mean_V1, mean_V2, n) %>% 
              distinct(),
            aes(mean_V1, mean_V2, label = cluster),
            size = 2,
            alpha = 0.6,
            show.legend = F) +
  stripped_theme +
  scale_fill_manual(values = c("selected" = "green3", "adjacent" = "red", "not adjacent" = "gray")) +
  scale_color_manual(values = c("selected" = "green3", "adjacent" = "red", "not adjacent" = "gray"))
  
plots <- 
  lapply(unique(plot_data$cluster),
         function(clus_) {
           # clus_ <- 1
           
           adjacent_clusters <- 
             pig_gene_umap_cluster_adjacent %>% 
             filter(cluster == clus_)
           
           plt_dat <- 
             plot_data %>% 
             mutate(cluster_type = case_when(cluster == clus_ ~ "selected",
                                             cluster %in% adjacent_clusters$adjacent_cluster ~ "adjacent",
                                             T ~ "not adjacent")) 
             
           plot +
             geom_encircle(data = plt_dat,
                           aes(group = cluster, 
                               fill = cluster_type, 
                               color = cluster_type),
                           alpha = 0.1, 
                           s_shape=1, 
                           expand=0,
                           # color = "white",
                           show.legend = F) +
             geom_encircle(data = plt_dat,
                           aes(group = cluster, 
                               color = cluster_type),
                           fill = NA,
                           alpha = 0.7, 
                           s_shape=1, 
                           expand=0,
                           # color = "white",
                           show.legend = F)
             
         })

pdf(savepath("UMAP adjacent clusters.pdf"), width = 6, height = 6)
plots
dev.off()



pig_gene_cluster_enrich_hyper_clust <- 
  pig_gene_cluster_enrich_hyper %>% 
  select(enhanced_tissues, cluster, adj_pval, log_p) %>%
  mutate(log_p = ifelse(adj_pval < 0.05, 1, 0)) %>%
  select(-adj_pval) %>%
  spread(cluster, log_p, fill = 0) %>% 
  column_to_rownames("enhanced_tissues") %>%
  t() %>%
  dist() %>% 
  hclust(method = "ward.D2")



similar_clusters <- 
  cutree(pig_gene_cluster_enrich_hyper_clust, h = 0.1) %>% 
  enframe("cluster", "similar_cluster") %>% 
  mutate(similar_cluster = as.character(similar_cluster))
  

pig_gene_cluster_enrich_hyper %>% 
  select(enhanced_tissues, cluster, adj_pval, log_p) %>%
  group_by(enhanced_tissues, cluster) %>%
  mutate(log_p = min(log_p, 20)) %>%
  ungroup() %>% 
  mutate(log_p = ifelse(adj_pval < 0.05, 1, 0)) %>%
  select(-adj_pval) %>%
  spread(cluster, log_p, fill = 0) %>% 
  column_to_rownames("enhanced_tissues") %>%
  pheatmap(clustering_method = "ward.D2",
           cutree_cols = n_distinct(similar_clusters$similar_cluster),
           annotation_col = similar_clusters %>% 
             column_to_rownames("cluster"), 
           filename = savepath("UMAP similar clusters heatmap.pdf"),
           width = 15, 
           height = 6)



pig_gene_umap_merged_clusters <-
  pig_gene_umap_cluster_adjacent %>%
  left_join(similar_clusters) %>% 
  left_join(similar_clusters,
            by = c("adjacent_cluster" = "cluster"),
            suffix = c("", "_adj")) %>%
  left_join(pig_gene_umap_temp2 %>% 
              group_by(cluster) %>% 
              count()) %>%
  
  mutate(similar_adj = similar_cluster == similar_cluster_adj) %>%
  
  filter(similar_adj | n < 25) %>% 
  select(1:2) %>% 
  graph_from_data_frame(directed = F) %>%
  components() %$%
  membership %>%
  enframe("cluster", "merged_cluster") %>%
    right_join(pig_gene_umap_temp2 %>% 
                 select(cluster) %>% 
                 distinct()) %>% 
    mutate(merged_cluster = ifelse(is.na(merged_cluster),
                                   cluster,
                                   paste("mg", merged_cluster)))


# # pig_gene_umap_merged_clusters
# pig_gene_umap_cluster_adjacent %>%
#   left_join(pig_gene_umap %>% 
#               group_by(cluster) %>% 
#               count()) %>%
#   filter(n < 20) %>% 
#   graph_from_data_frame(directed = F) %>%
#   components() %$%
#   membership %>%
#   enframe("cluster", "merged_cluster") %>%
#   arrange(merged_cluster)


pig_gene_umap <- 
  pig_gene_umap_temp2 %>% 
  left_join(pig_gene_umap_merged_clusters) %>%
  group_by(merged_cluster_old = merged_cluster) %>%
  mutate(V1_mean = mean(V1),
         V2_mean = mean(V2)) %>% 
  ungroup() %>%
  arrange(-V2_mean,
          V1_mean) %>%
  mutate(merged_cluster = factor(unclass(factor(merged_cluster_old, unique(merged_cluster_old)))))

plot_data <-
  pig_gene_umap %>%
  left_join(pig_gene_classification)
#
# ####
# plot_data %>%
#   group_by(merged_cluster) %>%
#   mutate(mean_V1 = mean(V1),
#          mean_V2 = mean(V2),
#          dist = sqrt((V1 - mean_V1)^2 + (V2 - mean_V2)^2),
#          n = n()) %>%
#   # select(V1, V2, dist, louvain, n) %>%
#   arrange(merged_cluster, dist) %>%
#   mutate(cum_member = cumsum(1/n)) %>%
#   # filter(cum_member < 0.9) %>%
#   mutate(hull = row_number() %in% chull(V1, V2)) %>%
#   filter(hull) %>%
#   mutate(hull = order(chull(V1, V2))) %>%
#   arrange(merged_cluster, hull) %>%
#   ggplot(aes(V1, V2)) +
#   geom_hex(aes(fill = 1),
#            data = plot_data,
#            fill = "gray90",
#            bins = plot_bins) +
#   geom_encircle(aes(fill = merged_cluster,
#                     color = merged_cluster),
#                 alpha = 0.1,
#                 s_shape=1,
#                 expand=0,
#                 # color = "white",
#                 show.legend = F) +
#   geom_encircle(aes(color = merged_cluster),
#                 fill = NA,
#                 alpha = 0.7,
#                 s_shape=1,
#                 expand=0,
#                 # color = "white",
#                 show.legend = F) +
#   geom_text(data = . %>%
#               select(merged_cluster, mean_V1, mean_V2, n) %>%
#               distinct(),
#             aes(mean_V1, mean_V2, label = merged_cluster),
#             show.legend = F) +
#   geom_point(data = km_centers,
#              aes(mean_V1, mean_V2),
#              inherit.aes = F) +
#   stripped_theme


pig_gene_umap %>%
  group_by(merged_cluster) %>%
  mutate(mean_V1 = mean(V1),
         mean_V2 = mean(V2),
         dist = sqrt((V1 - mean_V1)^2 + (V2 - mean_V2)^2),
         n = n()) %>%
  ggplot(aes(V1, V2)) +
  geom_hex(aes(fill = 1),
           fill = "gray90",
           bins = plot_bins) +
  geom_text(data = . %>% 
              select(merged_cluster, mean_V1, mean_V2, n) %>% 
              distinct(),
            aes(mean_V1, mean_V2, label = merged_cluster),
            size = 2,
            alpha = 0.6,
            show.legend = F) +
  geom_encircle(aes(group = cluster),
                alpha = 0.1, 
                s_shape=1, 
                expand=0,
                fill = "gray",
                color = "black",
                show.legend = F) +
  geom_encircle(data = . %>%
                  filter(grepl("^mg", merged_cluster_old)),
                aes(group = merged_cluster, 
                    color = merged_cluster,
                    fill = merged_cluster),
                
                alpha = 0.5, 
                s_shape=1, 
                expand=0,
                # color = "white",
                show.legend = F) +
  stripped_theme 
ggsave(savepath("UMAP merged clusters.pdf"), width = 6, height = 6)



# ----- Hypergeometric test ------

pig_gene_cluster_enrich_hyper_test_data <- 
  pig_gene_umap %>%
  # head(1000) %>%
  select(enssscg_id, merged_cluster) %>%
  left_join(pig_gene_classification) %>% 
  mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues), 
                                   "not enriched", enhanced_tissues)) %>%
  separate_rows(enhanced_tissues, sep = ";") %>%
  
  {
    # q is the number of successes
    q <- 
      group_by(., enhanced_tissues, merged_cluster, .drop = F) %>%
      summarise(q = n_distinct(enssscg_id))
    # k is the number of tries - i.e. the cluster size
    k <- 
      group_by(., merged_cluster, .drop = F) %>%
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
  group_by(enhanced_tissues, merged_cluster) %>%
  # Testing the chance of getting the observed number or higher per cluster and tissue type
  summarise(p_value = phyper(q - 1, m, n, k, lower.tail = F)) %>%
  ungroup() %>%
  mutate(p_value = ifelse(p_value == 0, .Machine$double.xmin, p_value),
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

# pig_gene_cluster_enrich_hyper %>% 
#   select(enhanced_tissues, merged_cluster, log_p) %>%
#   group_by(enhanced_tissues, merged_cluster) %>%
#   mutate(log_p = min(log_p, 20)) %>%
#   ungroup() %>% 
#   mutate(log_p = ifelse(log_p < 3, 0, 1)) %>%
#   spread(merged_cluster, log_p) %>% 
#   column_to_rownames("enhanced_tissues") %>%
#   pheatmap(clustering_method = "ward.D2")

# 
# 
# pig_gene_cluster_only_enrich_hyper_test_data <- 
#   pig_gene_umap %>%
#   # head(1000) %>%
#   select(enssscg_id, merged_cluster) %>%
#   left_join(pig_gene_classification) %>% 
#   mutate(enhanced_tissues = ifelse(is.na(enhanced_tissues) | specificity_category == "Tissue enhanced", 
#                                    "not enriched", enhanced_tissues)) %>%
#   separate_rows(enhanced_tissues, sep = ";") %>%
#   
#   {
#     # q is the number of successes
#     q <- 
#       group_by(., enhanced_tissues, merged_cluster, .drop = F) %>%
#       summarise(q = n_distinct(enssscg_id))
#     # k is the number of tries - i.e. the cluster size
#     k <- 
#       group_by(., merged_cluster, .drop = F) %>%
#       summarise(k = n_distinct(enssscg_id))
#     
#     # m is the number of possible successes
#     m <- 
#       group_by(., enhanced_tissues, .drop = F) %>%
#       summarise(m = n_distinct(enssscg_id))
#     
#     # n is the population size - i.e. the number of genes
#     n <- n_distinct(.$enssscg_id)
#     
#     q %>% 
#       left_join(k) %>% 
#       left_join(m) %>%
#       # n is the population size - i.e. the number of genes
#       mutate(n = n - m) 
#     
#   } 
# 
# pig_gene_cluster_only_enrich_hyper <- 
#   pig_gene_cluster_only_enrich_hyper_test_data %>%
#   group_by(enhanced_tissues, merged_cluster) %>%
#   # Testing the chance of getting the observed number or higher per cluster and tissue type
#   summarise(p_value = phyper(q - 1, m, n, k, lower.tail = F)) %>%
#   ungroup() %>%
#   mutate(p_value = ifelse(p_value == 0, 1e-300, p_value),
#          adj_pval = p.adjust(p_value, method = "BH"),
#          sign_level = case_when(adj_pval < 0.00000001 ~ 8,
#                                 adj_pval < 0.0000001 ~ 7,
#                                 adj_pval < 0.000001 ~ 6,
#                                 adj_pval < 0.00001 ~ 5,
#                                 adj_pval < 0.0001 ~ 4,
#                                 adj_pval < 0.001 ~ 3,
#                                 adj_pval < 0.01 ~ 2,
#                                 adj_pval < 0.05 ~ 1,
#                                 T ~ 0),
#          log_p = -log10(adj_pval))
#####


# Save
save(list = c("pig_gene_umap",
              "pig_gene_cluster_enrich_hyper_test_data",
              "pig_gene_cluster_enrich_hyper"),
     file = "data/processed/gene_umap_res.Rdata")

