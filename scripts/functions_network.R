

plot_network_groups <- 
  function(net_, group_col, cor_col, min_nodes = 0) {
    
    net_nodes <- 
      net_ %>% 
      activate(nodes) %>%
      rename(group = group_col) %>% 
      mutate(row_n = row_number()) %>%
      as_tibble() %>% 
      select(row_n, group)
    
    intranode_ <- 
      net_ %>% 
      activate(edges) %>% 
      left_join(net_nodes, 
                by = c("from" = "row_n")) %>%
      left_join(net_nodes, 
                by = c("to" = "row_n"),
                suffix = c("_from", "_to")) %>%
    filter(group_from == group_to) %>% 
      as_tibble() %>% 
      rename(cor = cor_col) %>%
      group_by(group = group_from) %>%
      summarise(cor = mean(cor))
    
    
    net_ %>% 
      activate(nodes) %>%
      rename(group = group_col) %>%
      morph(to_contracted, group, 
            edge.attr.comb= "mean",
            simplify = T) %>%
      crystallise() %>%
      pull(graph) %>%
      {.[[1]]} %>%
      mutate(n = sapply(.tidygraph_node_index,
                        length)) %>%
      filter(n > min_nodes) %>% 
      left_join(intranode_, 
                by = "group") %>%
      activate(edges) %>% 
      mutate(cor = sapply(.orig_data, 
                           function(x) {mean(x$cor)})) %>%
      ggraph(layout = "kk") +
      geom_edge_fan(aes(width = cor),
                    color = "gray",
                    show.legend = F)  +
      geom_node_point(aes(size = sqrt(n),
                          color = cor),
                      show.legend = F) +
      geom_node_text(aes(label = paste(n, round(cor, 2), sep = "\n")),
                     show.legend = F,
                     color = "black") +
      theme_graph(base_family="sans") +
      scale_size_continuous(range = c(6, 20))+
      scale_color_gradientn(colours = ggsci::pal_material("red")(10)) + 
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(group_col)
    
    # nodes_ <- 
    #   net_ %>% 
    #   activate(nodes) %>%
    #   as_tibble() %>% 
    #   select(group = group_col) %>%
    #   mutate(row_number = row_number())
    # 
    # edges_ <- 
    #   net_ %>% 
    #   activate(edges) %>%
    #   as_tibble() %>%
    #   select(1:2)
    # 
    # 
    # group_members <- 
    #   nodes_ %>%
    #   group_by(group) %>%
    #   summarise(n_nodes = n())
    # 
    # group_connections <- 
    #   edges_ %>%
    #   
    #   left_join(nodes_, 
    #             by = c("from" = "row_number")) %>%
    #   left_join(nodes_, 
    #             by = c("to" = "row_number"),
    #             suffix= c("_from", "_to")) %>%
    #   filter(group_from != group_to) %>%
    #   
    #   mutate(swap = group_from > group_to,
    #          new_from = ifelse(swap, to, from),
    #          new_to = ifelse(swap, from, to),
    #          new_group_from = ifelse(swap, group_to, group_from),
    #          new_group_to = ifelse(swap, group_from, group_to)) %>%
    #   
    #   select(from = new_from, 
    #          to = new_to, 
    #          group_from = new_group_from, 
    #          group_to = new_group_to) %>%
    #   
    #   group_by(group_from, group_to) %>% 
    #   summarise(n_connections = n()) %>% 
    #   ungroup() 
    # 
    # 
    # 
    # group_connections %>% 
    #   as_tbl_graph(directed = F) %>%
    #   left_join(group_members %>% 
    #               mutate(group = as.character(group)),
    #             by = c("name" = "group")) %>%
    #   ggraph(layout = "kk") +
    #   geom_edge_fan(aes(width = sqrt(n_connections),
    #                     alpha = sqrt(n_connections)),
    #                 color = "gray",
    #                 show.legend = F) + 
    #   geom_node_point(aes(size = sqrt(n_nodes),
    #                       color = log10(n_nodes)),
    #                   show.legend = F) +
    #   geom_node_text(aes(label = n_nodes),
    #                  show.legend = F) +
    #   theme_graph(base_family="sans") +
    #   ggtitle(group_col) +
    #   # scale_color_binned(type = "viridis") +
    #   scale_color_viridis_c() + 
    #   theme(plot.title = element_text(hjust = 0.5))
             
  }

get_edge_groups <- 
  function(net_, cols_) {
    nodes_ <- 
      net_ %>% 
      activate(nodes) %>%
      as_tibble() %>% 
      select(all_of(cols_)) %>%
      mutate(row_number = row_number())
    
    net_ %E>% 
      left_join(nodes_, 
                by = c("from" = "row_number")) %>%
      left_join(nodes_, 
                by = c("to" = "row_number"),
                suffix= c("_from", "_to")) %>%
      as_tibble() %>%
      pivot_longer(cols = contains("group"), values_to = "cluster") %>%
      mutate(origin = paste0("group_", str_extract(name, "from$|to$")),
             name = gsub("_from$|_to$", "", name)) %>% 
      spread(origin, cluster) %>%
      mutate(swap = group_from > group_to,
             new_from = ifelse(swap, to, from),
             new_to = ifelse(swap, from, to),
             new_group_from = ifelse(swap, group_to, group_from),
             new_group_to = ifelse(swap, group_from, group_to))
    
  }

add_node_feature <- 
  function(net_, node_feat_type, ...) {
    
    txt_ <- 
      paste0("net_ %<>% mutate(", node_feat_type, " = ", node_feat_type, "(...))")
    
    eval(parse(text = txt_))
    
    net_
  }


# add_centrality_to_network <- 
#   function(net_, centrality_type) {
#     
#     txt_ <- 
#       paste0("net_ %<>% mutate(", centrality_type, " = ", centrality_type, "())")
#     
#     eval(parse(text = txt_))
#     
#     net_
#   }

