# ----- Classification plots -----

class_spec_dist_piechart <- function(class_table) {
  spec_plot <-
    class_table %>%
    group_by(spec_category) %>%
    summarise(n_genes = n()) %>%
    mutate(spec_category = factor(spec_category, levels = spec_category_levels),
           perc = paste0("(", round(100 * n_genes / sum(n_genes), digits = 1), "%)")) %>%

    ggplot(aes("", n_genes, fill = spec_category)) +
    geom_bar(stat = "identity",
             show.legend = F,
             color = "white",
             size = 1)+
    geom_text(aes(y = n_genes, label = paste(n_genes, perc, sep = "\n")),
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4)+
    geom_text(aes(x = 1.5, y = n_genes, label = spec_category),
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4)+

    coord_polar("y", start = 0)+
    scale_fill_manual(values = gene_category_pal) +
    theme_void()


  ### Distribution
  dist_plot <-
    class_table %>%
    group_by(dist_category) %>%
    summarise(n_genes = n()) %>%
    mutate(dist_category = factor(dist_category, levels = dist_category_levels),
           perc = paste0("(", round(100 * n_genes / sum(n_genes), digits = 1), "%)")) %>%

    ggplot(aes("", n_genes, fill = dist_category)) +
    geom_bar(stat = "identity",
             show.legend = F,
             color = "white",
             size = 1)+
    geom_text(aes(y = n_genes, label = paste(n_genes, perc, sep = "\n")),
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4)+
    geom_text(aes(x = 1.5, y = n_genes, label = dist_category),
              position = position_stack(vjust = 0.5),
              color = "black",
              size = 4)+

    coord_polar("y", start = 0)+
    scale_fill_manual(values = gene_category_pal) +
    theme_void()

  list(spec_plot, dist_plot)
}

class_elevated_bar_plot <- function(class_table){

  class_table %>%
  {names <- rownames(.); as.tibble(.) %>% mutate(tissue = names)} %>%
    gather(key = "Classification", value = "Number of genes", -tissue) %>%
    mutate(tissue = factor(tissue, levels = rev(unique(tissue[order(mapply(tissue, FUN = function(x) sum(`Number of genes`[tissue == x & Classification %in% c("Tissue enriched","Celltype enriched",
                                                                                                                                                               "Group enriched","Tissue enhanced","Celltype enhanced")])))]))),
           Classification = factor(Classification, levels = c("Not detected in any tissues","Not detected in any celltypes","Not detected in this tissue","Not detected in this celltype","Mixed in this tissue", "Mixed in this celltype","Expressed in all tissues","Expressed in all celltypes","Tissue enhanced", "Celltype enhanced","Group enriched","Tissue enriched", "Celltype enriched")),
           Classification = gsub(pattern = translate_categories, replacement = names(translate_categories), Classification)) %>%
    filter(Classification %in% c("Tissue enriched","Celltype enriched",
                                 "Group enriched","Tissue enhanced","Celltype enhanced")) %>%
    ggplot(aes(tissue, `Number of genes`, fill = Classification))+
    geom_bar(stat = "identity") +
    scale_fill_manual(name = "",values = cat2.cols)+
    simple_theme+
    xlab("")+
    ylab("Number of genes")+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = c(0.8, 0.8))

}


class_spec_dist_barplot <- function(class_table) {
  cat <-
    class_table %>%
    group_by(spec_category, dist_category) %>%
    summarise(n_genes=n()) %>%
    ungroup() %>%
    mutate(spec_category = factor(spec_category, levels = spec_category_levels),
           dist_category = factor(dist_category, levels = dist_category_levels))


  list(ggplot(cat, aes(x=dist_category, y=n_genes, fill=spec_category)) +
         geom_bar(stat = "identity")+
         ggtitle("Distribution category")+
         scale_fill_manual(values = gene_category_pal, name = "Specificity category")+
         ylab("Number of genes")+
         xlab("")+
         simple_theme+
         theme(axis.text.x = element_text(angle=60, hjust=1)),

       ggplot(cat, aes(x=spec_category, y=n_genes, fill=dist_category)) +
         geom_bar(stat = "identity")+
         ggtitle("Specificity category")+
         scale_fill_manual(values = gene_category_pal, name = "Distribution category")+
         ylab("Number of genes")+
         xlab("")+
         simple_theme+
         theme(axis.text.x = element_text(angle=60, hjust=1)))
}

chord_classification <- function(from, to, sizes, grid.col, groups, plot.order, size_labels = F){
  require(circlize)

  factors.from <- unique(from)
  factors.to <- unique(to)
  factors <- c(factors.from, factors.to)


  tb <-
    tibble(from, to, sizes)

  #groups <- groups[plot.order]
  gap.after.par <- c()
  for(i in 1:(length(groups)-1)) {
    if(groups[i] == groups[i+1]) {
      gap.after.par <- c(gap.after.par, 2)
    } else {
      gap.after.par <- c(gap.after.par, 15)
    }
  }

  if(groups[length(groups)] == groups[1]) {
    gap.after.par <- c(gap.after.par, 2)
  } else {
    gap.after.par <- c(gap.after.par, 15)
  }

  circos.par(gap.after = gap.after.par)

  chord <-
    tb %>%
    chordDiagram(grid.col = grid.col,
                 directional = 0,
                 annotationTrack="grid",
                 annotationTrackHeight = 0.05,
                 preAllocateTracks = 1,
                 order = plot.order)

  if(size_labels) {
    for(i in 1:nrow(chord)) {
      value <- chord$value[i]
      if(is.null(value)) value <- chord$value1[i]
      x1 <- chord$x1[i] - value / 2
      x2 <- chord$x2[i] - value / 2
      to_ <- chord$cn[i]
      from_ <- chord$rn[i]
      circos.text(x = x1, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = from_, niceFacing = T)
      circos.text(x = x2, y = -1, track.index = 2, labels = value, cex = 0.7, sector.index = to_, niceFacing = T)
    }
  }



  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector.name <- get.cell.meta.data("sector.index")
    sector.index <- get.cell.meta.data("sector.numeric.index")

    adjustment <- ifelse(sector.index %% 2 == 1, 0.3, -0.2)

    circos.segments(x0 = mean(xlim), x1 = mean(xlim),
                    y0 = min(ylim), y1 = mean(ylim)-0.2 + adjustment,
                    sector.name)


    circos.text(mean(xlim), mean(ylim) + adjustment, sector.name, niceFacing = TRUE, facing = "bending")
  }, bg.border = NA)

  circos.clear()
}

class_chord_plot <- function(class_table) {

  class_table %>%
    group_by(spec_category, dist_category) %>%
    summarise(n_genes=n()) %>%
    ungroup() %>%
    mutate(spec_category = case_when(spec_category == "not detected" ~ "not detected ",
                                     T ~ spec_category)) %$%
    chord_classification(from = spec_category,
                         to = dist_category,
                         sizes = n_genes,
                         grid.col = gene_category_pal,
                         groups = c(rep(1, 5), rep(2, 5)),
                         plot.order = c(c(spec_category_levels[-5], "not detected "),
                                        dist_category_levels),
                         size_labels = T)

}

class_tissue_n_enriched_barplot <- function(class_table) {

  class_table_temp <-
    class_table %>%
    select(gene, spec_category, enriched_tissues) %>%
    separate_rows(enriched_tissues, sep = ", ")


  tissues_not_in_plot <- with(tissue_mapping, consensus_tissue_name[which(!consensus_tissue_name %in% class_table_temp$enriched_tissues)])

  if(length(tissues_not_in_plot) != 0) warning(paste0("These tissues are not in the plot: ",
                                                      paste(tissues_not_in_plot, collapse = ", ")))

  class_table_temp %>%
    filter(!is.na(enriched_tissues)) %>%
    group_by(enriched_tissues, spec_category) %>%
    summarise(n_genes = n()) %>%
    ungroup() %>%
    mutate(enriched_tissues = factor(enriched_tissues, levels = {group_by(., enriched_tissues) %>%
        summarise(n_genes = sum(n_genes)) %$%
        enriched_tissues[order(n_genes)]})) %>%
    ggplot(aes(enriched_tissues, n_genes, fill = spec_category)) +
    geom_col(width = 0.9, color = "white", size = 0.1) +
    simple_theme +
    scale_fill_manual(values = gene_category_pal, name = "Specificity") +
    coord_flip() +
    ggtitle("Number of genes enriched per tissue") +
    xlab("Tissue") +
    ylab("Number of genes") +
    scale_y_continuous(position = "bottom")
}

class_tissue_n_expressed_barplot <- function(class_table) {

  class_table_temp <-
    class_table %>%
    select(gene, dist_category, tissues_detected) %>%
    separate_rows(tissues_detected, sep = ", ")


  tissues_not_in_plot <- with(tissue_mapping, consensus_tissue_name[which(!consensus_tissue_name %in% class_table_temp$tissues_detected)])

  if(length(tissues_not_in_plot) != 0) warning(paste0("These tissues are not in the plot: ",
                                                      paste(tissues_not_in_plot, collapse = ", ")))

  class_table_temp %>%
    filter(!is.na(tissues_detected)) %>%
    group_by(tissues_detected, dist_category) %>%
    summarise(n_genes = n()) %>%
    ungroup() %>%
    mutate(tissues_detected = factor(tissues_detected, levels = {group_by(., tissues_detected) %>%
        summarise(n_genes = sum(n_genes)) %$%
        tissues_detected[order(n_genes)]})) %>%
    ggplot(aes(tissues_detected, n_genes, fill = dist_category)) +
    geom_col(width = 0.9, color = "white", size = 0.1) +
    simple_theme +
    scale_fill_manual(values = gene_category_pal, name = "Distribution") +
    coord_flip() +
    ggtitle("Number of genes detected per tissue") +
    xlab("Tissue") +
    ylab("Number of genes") +
    scale_y_continuous(position = "bottom")
}

## Retinagram

map_colors_to_edge <- function(edge_data, color_mapping, label_col, color_col, mean_color = F, default_color = "gray80") {
  mapped_colors <-
    edge_data %>%
    as_tibble() %>%
    left_join(color_mapping %>%
                select(label = label_col,
                       color = color_col),
              by = c("node2.label" = "label")) %>%
    mutate(radius = xend^2 + yend^2,
           edge.id = as.character(edge.id)) %>%
    arrange(-radius)

  mapped_colors %>%
    filter(!near(radius, 1)) %$%
    sapply(edge.id,
           FUN = function(edge_id_) {
             xend_ = mapped_colors$xend[which(mapped_colors$edge.id == edge_id_)]
             yend_ = mapped_colors$yend[which(mapped_colors$edge.id == edge_id_)]

             new_color <-
               mapped_colors %>%
               filter(near(x, xend_) & near(y, yend_)) %$%
               ifelse(mean_color,
                      colorRampPalette(color)(3)[2],
                      ifelse(length(unique(color)) == 1,
                             unique(color), default_color))


             new_color_column <- mapped_colors$color
             new_color_column[which(mapped_colors$edge.id == edge_id_)] <- new_color

             mapped_colors$color <<- new_color_column
             NULL


           })


  mapped_colors
}

circular_dendrogram_retinastyle <-
  function(clust, color_mapping, label_col, color_col, mean_color = F) {
    require(ggraph)
    require(igraph)
    require(viridis)
    require(tidyverse)
    require(magrittr)

    dendrogram <-
      clust %>%
      as.dendrogram()



    g <-
      ggraph(dendrogram, layout = 'dendrogram', circular = T)

    edge_data <- get_edges()(g$data)

    edge_data_colors <- map_colors_to_edge(edge_data, color_mapping, label_col, color_col, mean_color = mean_color)


    g +
      scale_edge_width(range = c(1.5, 6))+
      geom_edge_diagonal(aes(edge_color = as.character(edge_data$edge.id),
                             edge_width = 1 - sqrt(xend^2 + yend^2)),
                         strength = 0.8,
                         show.legend = F) +
      scale_edge_color_manual(values = edge_data_colors %$%
                                set_names(c(color, "gray80"), c(edge.id, "")))  +
      g$data %>%
      filter(label != "") %>%
      mutate(degree = case_when(x >= 0 ~ asin(y) * 180 / pi,
                                x < 0 ~ 360 - asin(y) * 180 / pi)) %>%
      left_join(color_mapping %>%
                  select(label = label_col,
                         color = color_col),
                by = "label") %>%
                {geom_node_text(data = .,
                                aes(label = label),
                                angle = .$degree,
                                hjust = ifelse(.$x < 0, 1, 0),
                                vjust = 0.5,
                                size = 3)}  +
      scale_x_continuous(expand = expand_scale(c(0.25, 0.25))) +
      scale_y_continuous(expand = expand_scale(c(0.25, 0.25))) +

      coord_fixed() +
      theme_void()
  }


# ----- tissue clustering plots -----

evolutionary_tree_plot <- function(data,
                                   dist_fun = function(x) cor(x, method = 'spearman', use="pairwise.complete.obs") %>%
                                     {as.dist(1-.)},
                                   tree_fun = function(x) nj(),
                                   color_mapping,
                                   mapping_col,
                                   color_col) {

  # calculate distance
  if(is.function(dist_fun) & is.function(tree_fun)) {
    data_dist <- dist_fun(data)
  } else {
    data_dist <- dist_fun
  }

  # neighbor-joining tree estimation
  if(is.function(tree_fun)) {
    nj_tree <- tree_fun(data_dist)
  } else {
    nj_tree <- tree_fun
  }



  plot.phylo(as.phylo(nj_tree),
             type="u",
             lab4ut = "axial",
             font = 1,
             cex = 0.8,
             tip.color = color_mapping %>%
               rename(mapping_col = mapping_col,
                      color_col = color_col) %>%
               {.[match(nj_tree$tip.label, .$mapping_col),]} %$%
               color_col,
             edge.col="black",
             edge.width=2,
             show.node.label=TRUE, no.margin=TRUE,
             use.edge.length = T)

}

pca_calc <- function(data, npcs) {

  pca_res <-
    data %>%
    pca(nPcs = npcs)

  pca_stats <-
    tibble(PC = 1:npcs,
           R2cum = R2cum(pca_res))

  informative_pcs <- pca_stats$PC[which(pca_stats$R2cum > 0.95)[1]]

  pca_stats_plot <-
    pca_stats %>%
    select(PC, R2cum) %>%
    ggplot(aes(PC, R2cum)) +
    geom_point() +
    geom_line() +
    simple_theme +
    geom_vline(xintercept = informative_pcs, linetype = "dashed") +
    annotate("text",
             x = informative_pcs,
             y = 0.55,
             label = paste0("PC ", informative_pcs,
                            "\nR2 = ", round(pca_stats[informative_pcs, ]$R2cum, 3)),
             hjust = 1,
             vjust = 0)

  list(pca = pca_res,
       scores = scores(pca_res),
       loadings = loadings(pca_res),
       stats = pca_stats,
       pc_95 = informative_pcs,
       plot = pca_stats_plot)
}

pca_score_plot <- function(pca_scores, group_mapping, mapping_col, group_col, pal, xpc = 1, ypc = 2) {

  pca_scores_mapped <-
    pca_scores %>%
    as_tibble(rownames = "mapping_col") %>%
    left_join(group_mapping,
              by = c("mapping_col" = mapping_col)) %>%
    rename(xpc = xpc + 1,
           ypc = ypc + 1)

  if(group_col == mapping_col) {
    pca_scores_mapped <-
      pca_scores_mapped %>%
      mutate(group_col = mapping_col)
  } else {
    pca_scores_mapped <-
      pca_scores_mapped %>%
      rename(group_col = group_col)
  }


  pca_scores_mapped %>%
    group_by(group_col) %>%
    mutate(mean_x = mean(xpc),
           mean_y = mean(ypc)) %>%
    ungroup() %>%

    {ggplot(., aes(xpc, ypc, color = group_col)) +
        geom_point(show.legend = F) +
        scale_color_manual(values = pal) +
        simple_theme +
        geom_segment(aes(xend = mean_x, yend = mean_y),
                     show.legend = F) +
        geom_text(data = select(., group_col, mean_x, mean_y) %>%
                    unique(),
                  aes(x = mean_x, y = mean_y, label = group_col),
                  show.legend = F) +
        xlab(paste0("PC", xpc)) +
        ylab(paste0("PC", ypc))}


}

# ----- Not currently in use -----





tissue_distributions_plot <- function(data, expression_col, content_col, det_lim, pal, do.tissues = "all") {

  data %>%
    filter(consensus_tissue_name %in% sample(unique(consensus_tissue_name), 5)) %>%
    rename(content = content_col,
           expression = expression_col) %>%
        {ggplot(., aes(content, expression, fill = content, color = content))+
            stat_summary(aes(content, expression, group = 1),
                         fun.y = "median",
                         fun.args = c("na.rm" = T),
                         geom = "line",
                         size = 2,
                         alpha = 0.5)+
            stat_summary(aes(content, expression, group = 1),
                         fun.y = "min",
                         geom = "line",
                         inherit.aes = F,
                         size = 1,
                         alpha = 0.5)+
            stat_summary(aes(content, expression, group = 1),
                         fun.y = "max",
                         geom = "line",
                         inherit.aes = F,
                         size = 1,
                         alpha = 0.5)+
            # stat_summary(fun.y = "min", geom = "line", aes(group = method), size = 0.5, alpha = 0.5)+
            # stat_summary(fun.y = "max", geom = "line", aes(group = method), size = 0.5, alpha = 0.5)+


            geom_violin(draw_quantiles = 0.5, alpha = 0.5, position = "identity")+
            geom_hline(yintercept = det_lim, linetype = "dashed")+
            simple_theme+
            scale_y_log10()+
            theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
            theme(axis.title = element_blank())+
            scale_fill_manual(values = pal)+
            scale_color_manual(values = pal)}

}


