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


chord_with_title <- function(from, to, sizes, grid.col, groups, plot.order, titles, 
                             from_labels, to_labels, size_labels = F, label_style = "horisontal", gap_adj = 1){
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
      gap.after.par <- c(gap.after.par, 2 * gap_adj)
    } else {
      gap.after.par <- c(gap.after.par, 15 * gap_adj)
    }
  }
  
  if(groups[length(groups)] == groups[1]) {
    gap.after.par <- c(gap.after.par, 2 * gap_adj)
  } else {
    gap.after.par <- c(gap.after.par, 15 * gap_adj)
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
  
  
  
  group_coordinates <- 
    sapply(plot.order, FUN = function(si) get.cell.meta.data("xlim",sector.index = si)) %>% 
    t() %>% 
    as_tibble(rownames = "sector_name") %>% 
    mutate(group = groups) %>% 
    group_by(group) %>% 
    mutate(x = cumsum(max.data)) %>%
    
    mutate(group_n = row_number(), 
           group_first = group_n == 1) 
  
  
  if(label_style == "horisontal") {
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      sector.name <- get.cell.meta.data("sector.index")
      sector_label <- c(from_labels, to_labels)[match(sector.name, c(from, to))]
      sector.index <- get.cell.meta.data("sector.numeric.index")
      
      adjustment <- ifelse(sector.index %% 2 == 1, 0.3, -0.2)
      
      circos.segments(x0 = mean(xlim), x1 = mean(xlim),
                      y0 = min(ylim), y1 = mean(ylim)-0.2 + adjustment,
                      sector.name)
      
      
      circos.text(mean(xlim), mean(ylim) + adjustment, sector_label, niceFacing = TRUE, facing = "bending")
      
      sector_i <- match(sector.name, group_coordinates$sector_name)
      
      if(group_coordinates$group_first[sector_i]) {
        
        group_coordinates_lim <- 
          group_coordinates %>% 
          filter(group == group_coordinates$group[sector_i]) %>% 
          filter(group_n == max(group_n) | group_n == min(group_n) )
        
        sector_indices <- which(gr_ == groups)
        
        group_n <- unique(group_coordinates_lim$group)
        
        max_gr <- max(group_coordinates_lim$x)
        
        titl_x <-
          max_gr / 2
        
        titl_lab <- titles[match(group_n, unique(group_coordinates$group))]
        
        circos.text(titl_x, 1.2, titl_lab, niceFacing = TRUE, facing = "bending", cex = 2)
      }
    }, bg.border = NA)
  } else if(label_style == "vertical") {
    circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      sector.name <- get.cell.meta.data("sector.index")
      sector_label <- c(from_labels, to_labels)[match(sector.name, c(from, to))]
      sector.index <- get.cell.meta.data("sector.numeric.index")
      
      
      
      circos.text(mean(xlim), min(ylim), sector_label, niceFacing = TRUE, facing = "clockwise", adj = c(0, 0))
      
      sector_i <- match(sector.name, group_coordinates$sector_name)
      
      if(group_coordinates$group_first[sector_i]) {
        
        group_coordinates_lim <- 
          group_coordinates %>% 
          filter(group == group_coordinates$group[sector_i]) %>% 
          filter(group_n == max(group_n) | group_n == min(group_n) )
        
        sector_indices <- which(gr_ == groups)
        
        group_n <- unique(group_coordinates_lim$group)
        
        max_gr <- max(group_coordinates_lim$x)
        
        titl_x <-
          max_gr / 2
        
        titl_lab <- titles[match(group_n, unique(group_coordinates$group))]
        
        circos.text(titl_x, 1.2, titl_lab, niceFacing = TRUE, facing = "bending", cex = 2)
      }
    }, bg.border = NA)
  }
  
  
  
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

pca_R2_plot <- function(pca_stats, mark_pc) {
  pca_stats %>% 
    ggplot(aes(PC, R2cum)) +
    geom_point() + 
    geom_line() + 
    simple_theme +
    geom_vline(xintercept = mark_pc, linetype = "dashed") + 
    annotate("text", 
             x = mark_pc,
             y = pca_stats$R2cum[mark_pc], 
             label = paste0("PC ", mark_pc, ": R2 = ", round(pca_stats$R2cum[mark_pc], 3)),
             hjust = -0.1, 
             vjust = -0.5, 
             angle = -90) 
}

plot_dendrogram <- 
  function(clust, 
           color_mapping, label_col, color_col, pal, do_label = F) {
    
    dendr <- dendro_data(clust)
    
    
    dendro_plot_data <- 
      left_join(dendr$segments, 
                dendr$labels, 
                by = c("x" = "x", "yend" = "y")) 
    
    dendro_plot <- 
      dendro_plot_data %>%
      ggplot() +
      geom_segment(aes(x=y, y=x, xend=yend, yend=xend, group = label))+
      
      scale_color_manual(values = pal)+
      scale_fill_manual(values = pal)+
      scale_x_reverse(expand = expand_scale(mult = 0.25))+
      theme(axis.text.y = element_blank(), 
            axis.title = element_blank(), 
            axis.ticks.y = element_blank(),
            plot.margin = unit(c(1,1,1,1), units = "mm"), 
            panel.background = element_blank()) 
    
    if(do_label) {
      dendro_plot + 
        geom_label(data = ggdendro::label(dendr) %>%
                    as_tibble(),
                  aes(label=label,
                      x=0,
                      y=x,
                      fill = label),
                  size = 3,
                  hjust = 0,
                  label.r = unit(0, units = "cm"),
                  nudge_x = 0.005,
                  show.legend = F)
    } else {
      dendro_plot + 
        geom_text(data = ggdendro::label(dendr) %>%
                    as_tibble(),
                  aes(label=label,
                      x=0,
                      y=x,
                      color = label),
                  size = 3,
                  hjust = 0,
                  nudge_x = 0.005,
                  show.legend = F)
    }
    
  }

umap_calc <- function(data, npcs = 2, n_epochs = 400, n_neighbors = round(sqrt(dim(data)[1]))) {
  
  data %>%
    umap(n_components = npcs,
         n_epochs = n_epochs,
         n_neighbors = n_neighbors)
  
}

umap_score_plot <- function(umap_scores, group_mapping, mapping_col, group_col, pal, xpc = 1, ypc = 2) {
  
  umap_scores_mapped <-
    umap_scores %>%
    as_tibble(rownames = "mapping_col") %>%
    left_join(group_mapping,
              by = c("mapping_col" = mapping_col)) %>%
    rename(xpc = xpc + 1,
           ypc = ypc + 1)
  
  if(group_col == mapping_col) {
    umap_scores_mapped <-
      umap_scores_mapped %>%
      mutate(group_col = mapping_col)
  } else {
    umap_scores_mapped <-
      umap_scores_mapped %>%
      rename(group_col = group_col)
  }
  
  
  umap_scores_mapped %>%
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
        xlab(paste0("V", xpc)) +
        ylab(paste0("V", ypc))}
  
  
}

# ----- correlation plots -----

basic_network_plot <- function(data, 
                         cor_col, 
                         var1_col, var2_col, 
                         var1_label_col, var2_label_col, 
                         var1_fill_col, var2_fill_col, 
                         fill_pal = NA, 
                         var1_color = "black",
                         var2_color = "black",
                         scale_edge = F, 
                         edge_color = "black") {
  
  if(is.na(fill_pal)) fill_pal <- rep("white", 1000)
  
  cor_links <- 
    data %>%
    do(tibble(var1 = .[,var1_col][[1]],
              var2 = .[,var2_col][[1]], 
              cor = .[,cor_col][[1]], 
              var1_label = .[,var1_label_col][[1]], 
              var2_label = .[,var2_label_col][[1]], 
              var1_fill = .[,var1_fill_col][[1]], 
              var2_fill = .[,var2_fill_col][[1]]))
  
  link_mapping <- 
    cor_links %>% 
    do(bind_rows(select(., 
                        var = var1, 
                        label = var1_label, 
                        fill = var1_fill) %>% 
                   mutate(var_n = "var1"), 
                 select(., 
                        var = var2, 
                        label = var2_label, 
                        fill = var2_fill) %>%
                   mutate(var_n = "var2"))) %>% 
    unique()

  g <-
    cor_links %>%
    graph_from_data_frame(directed = FALSE) %>%
    ggraph(layout = "kk") 
  
  edge_data <- get_edges()(g$data)
  node_data <- 
    get_nodes()(g$data) %>% 
    left_join(link_mapping, 
              by = c("name" = "var"))
  
  if (scale_edge) {
    g <-
      g + 
      geom_edge_fan(aes(width = cor, 
                        alpha = cor), 
                    color = edge_color) + 
      scale_edge_alpha_continuous(range = c(0.3, 1)) +
      scale_edge_width_continuous(range = c(1, 3))
    
  } else {
    g <- 
      g + 
      geom_edge_fan(color = edge_color)
  }
  
  
  g + 
    geom_node_point(data = node_data,
                    aes(fill = fill, 
                        color = var_n), 
                    stroke = 2,
                    size = 10,
                    shape = 21,
                    show.legend = F)+
    geom_node_text(data = node_data,
                   aes(label = label),
                   size = 4) +
    scale_fill_manual(values = fill_pal) +
    scale_color_manual(values = c("var1" = var1_color, "var2" = var2_color)) +
      # scale_edge_colour_gradientn(limits = c(-1, 1), colors = c("dodgerblue2", "white", "firebrick2")) +
      
      
    theme_void()
  
  
  
    
}

# ----- tissue grouping plots -----

tissue_connection_plot <- 
  function(mapping_table, col1, col2, col1_name, col2_name, title, pal = rep("black", 1000), na.rm = F) {
    tissue_overlap <- 
      mapping_table %>% 
      rename(col1 = col1,
             col2 = col2)
    
    if(na.rm) {
      tissue_overlap <- 
        tissue_overlap %>%
        filter(!(is.na(col1) | is.na(col2)))
    }
    tissue_overlap <- 
      tissue_overlap %>% 
      select(col1, col2) %>% 
      distinct() %>% 
      mutate(col1_i = unclass(factor(col1, levels = unique(col1))),
             col2_i = unclass(factor(col2, levels = unique(col2)))) %>% 
      group_by(col2) %>%
      mutate(col2_y = case_when(!is.na(col2) ~ mean(col1_i))) %>% 
      ungroup()
    
    
    
    tissue_overlap  %>%
      ggplot(aes(x = 1, xend = 2, 
                 y = col1_i, yend = col2_y)) + 
      annotate("text",
               x = 1:2,
               y = max(tissue_overlap$col2_y) + 1,
               label = c(col1_name, col2_name), 
               vjust = 0, 
               hjust = 1:0,
               fontface = "bold", 
               size = 4) +
      geom_segment(aes(color = col1), 
                   size = 2, 
                   alpha = 1, 
                   show.legend = F) + 
      geom_text(aes(x = 1, y = col1_i, label = col1, color = col1), 
                inherit.aes = F,
                hjust = 1, 
                show.legend = F) + 
      geom_text(aes(x = 2, y = col2_y, label = col2, color = col1), 
                inherit.aes = F, 
                hjust = 0, 
                show.legend = F) + 
      theme_void() + 
      scale_x_continuous(expand = expand_scale(3)) + 
      scale_color_manual(values = pal) + 
      ggtitle(title) + 
      theme(plot.title = element_text(hjust = 0.5))
  }

# ----- ggplot2 utility -----

label_unique <-  function (labels) {
  apply(labels, MARGIN = 1, 
        FUN = function(row_) unique(row_) %>% paste(collapse = " - ")) %>%
    list()
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


